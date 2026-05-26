use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::time::Duration;

use serde_json::{Value, json};
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader, Lines};
use tokio::process::{Child, ChildStdin, ChildStdout, Command};

// ---------- shared harness ----------

struct McpHarness {
    child: Child,
    stdin: ChildStdin,
    reader: Lines<BufReader<ChildStdout>>,
    next_id: i64,
}

impl McpHarness {
    async fn start(cfg_path: &Path) -> Self {
        let bin = env!("CARGO_BIN_EXE_vcf-mcp");
        let mut child = Command::new(bin)
            .arg("serve")
            .arg("--config")
            .arg(cfg_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .kill_on_drop(true)
            .spawn()
            .expect("spawn server");
        let mut stdin = child.stdin.take().expect("server stdin");
        let stdout = child.stdout.take().expect("server stdout");
        let mut reader = BufReader::new(stdout).lines();

        write_to(
            &mut stdin,
            &json!({
                "jsonrpc": "2.0",
                "id": 0,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": { "name": "vcf-mcp-it", "version": "0.0.0" }
                }
            }),
        )
        .await;
        let _ = recv_line(&mut reader).await;

        write_to(
            &mut stdin,
            &json!({"jsonrpc": "2.0", "method": "notifications/initialized"}),
        )
        .await;

        Self {
            child,
            stdin,
            reader,
            next_id: 1,
        }
    }

    async fn call_tool(&mut self, name: &str, arguments: Value) -> Value {
        let id = self.next_id;
        self.next_id += 1;
        let req = json!({
            "jsonrpc": "2.0",
            "id": id,
            "method": "tools/call",
            "params": { "name": name, "arguments": arguments }
        });
        write_to(&mut self.stdin, &req).await;
        let resp = recv_line(&mut self.reader).await;
        assert_eq!(resp["id"], id, "mismatched response id: {resp}");
        resp
    }

    async fn shutdown(mut self) {
        // Closing stdin causes the server to exit cleanly.
        drop(self.stdin);
        let _ = tokio::time::timeout(Duration::from_secs(5), self.child.wait()).await;
    }
}

async fn write_to(stdin: &mut ChildStdin, msg: &Value) {
    let line = format!("{msg}\n");
    stdin.write_all(line.as_bytes()).await.expect("write stdin");
    stdin.flush().await.expect("flush stdin");
}

async fn recv_line(reader: &mut Lines<BufReader<ChildStdout>>) -> Value {
    let line = tokio::time::timeout(Duration::from_secs(30), reader.next_line())
        .await
        .expect("timed out reading server response")
        .expect("server stdout error")
        .expect("server closed stdout");
    serde_json::from_str(&line)
        .unwrap_or_else(|e| panic!("server emitted non-JSON line: {line:?} ({e})"))
}

fn fixture_config_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/test_config.toml")
}

fn extract_text(resp: &Value) -> &str {
    resp["result"]["content"][0]["text"]
        .as_str()
        .unwrap_or_else(|| panic!("missing tool text content: {resp}"))
}

fn is_error_response(resp: &Value) -> bool {
    resp["error"].is_object()
}

// ---------- tests ----------

#[tokio::test]
async fn list_samples_returns_configured_samples() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h.call_tool("list_samples", json!({})).await;
    let text = extract_text(&resp);
    let samples: Vec<Value> = serde_json::from_str(text).unwrap();
    assert_eq!(samples.len(), 1);
    assert_eq!(samples[0]["name"], "NA12878");
    assert_eq!(samples[0]["build"], "GRCh38");
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_returns_variants_in_range() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chr17", "start": 50180000, "end": 50210000}),
        )
        .await;
    let text = extract_text(&resp);
    let body: Value = serde_json::from_str(text).unwrap();

    assert_eq!(body["sample"], "NA12878");
    assert_eq!(body["chrom"], "chr17");
    assert_eq!(body["truncated"], false);
    let count = body["count"].as_u64().unwrap();
    assert!(
        count > 0,
        "expected >0 variants in COL1A1 region, got {count}"
    );
    let variants = body["variants"].as_array().unwrap();
    assert_eq!(variants.len() as u64, count);

    let mut saw_genotype = false;
    let mut saw_alleles = false;
    for v in variants {
        assert_eq!(v["chrom"], "chr17");
        let pos = v["pos"].as_u64().unwrap();
        assert!(
            (50180000..=50210000).contains(&pos),
            "pos {pos} out of range"
        );
        assert!(v["ref"].is_string());
        assert!(v["alt"].is_string());
        assert!(v["filter"].is_string());
        if let Some(gt) = v["genotype"].as_str() {
            // Expect canonical "<idx>(/|<idx>)*" form, e.g. "0/1", "1|1", "./."
            assert!(
                gt.chars()
                    .all(|c| c.is_ascii_digit() || matches!(c, '/' | '|' | '.')),
                "unexpected genotype char in {gt}"
            );
            saw_genotype = true;
        }
        if v["genotype_alleles"].is_string() {
            saw_alleles = true;
        }
    }
    assert!(
        saw_genotype,
        "GIAB benchmark should yield at least one parsed genotype"
    );
    assert!(saw_alleles, "and at least one allele rendering");
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_accepts_unprefixed_chrom() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "17", "start": 50180000, "end": 50210000}),
        )
        .await;
    let text = extract_text(&resp);
    let body: Value = serde_json::from_str(text).unwrap();
    // Server normalizes to file convention; should hit the same data as chr17.
    assert_eq!(body["chrom"], "chr17");
    assert!(body["count"].as_u64().unwrap() > 0);
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_empty_when_no_variants() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    // Position outside the slice (slice covers ~50.1-50.3 Mbp, this is at the very edge).
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chr17", "start": 50299999, "end": 50299999}),
        )
        .await;
    let text = extract_text(&resp);
    let body: Value = serde_json::from_str(text).unwrap();
    assert_eq!(body["count"], 0);
    assert_eq!(body["variants"].as_array().unwrap().len(), 0);
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_single_position() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    // A 1-bp window — verify it works and produces a well-formed (possibly empty) response.
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chr17", "start": 50190000, "end": 50190000}),
        )
        .await;
    let text = extract_text(&resp);
    let body: Value = serde_json::from_str(text).unwrap();
    assert_eq!(body["start"], 50190000);
    assert_eq!(body["end"], 50190000);
    assert!(body["count"].is_number());
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_unknown_sample_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NOPE", "chrom": "chr17", "start": 50180000, "end": 50210000}),
        )
        .await;
    assert!(
        is_error_response(&resp),
        "expected error response, got {resp}"
    );
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains("NOPE"),
        "error should mention the bad sample name: {msg}"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_invalid_range_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chr17", "start": 50200000, "end": 50100000}),
        )
        .await;
    assert!(is_error_response(&resp), "expected error response");
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_too_large_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chr17", "start": 1, "end": 50_000_000}),
        )
        .await;
    assert!(is_error_response(&resp), "expected error response");
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(msg.contains("max") || msg.contains("region"), "msg: {msg}");
    h.shutdown().await;
}

#[tokio::test]
async fn query_region_unknown_chrom_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "query_region",
            json!({"sample": "NA12878", "chrom": "chrZZ", "start": 1, "end": 1000}),
        )
        .await;
    assert!(is_error_response(&resp), "expected error response");
    h.shutdown().await;
}
