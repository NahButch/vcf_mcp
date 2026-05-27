use std::fs::File;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::time::Duration;

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use serde_json::{Value, json};
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader, Lines};
use tokio::process::{Child, ChildStdin, ChildStdout, Command};

// ---------- fixture generators ----------
//
// `examples/annotate_rsids.rs` is the documented standalone path, but Windows
// Smart App Control on some boxes silently blocks freshly-built unsigned
// example exes. The test binary itself is reliably trusted, so we expose the
// same annotation logic as an #[ignore]d test. Run it manually with:
//   cargo test --release --test integration -- --ignored generate_rsid_fixture
// The resulting `tests/data/na12878_chr17_slice_rsids.vcf.gz` (+ .tbi) is
// what the lookup_rsids tests consume; commit the regenerated pair if the
// underlying slice ever changes.

fn annotate_rsids_to_file(src: &Path, dst: &Path) -> anyhow::Result<u32> {
    use anyhow::Context as _;

    let idx_dst = PathBuf::from(format!("{}.tbi", dst.display()));

    let mut reader = bgzf::io::Reader::new(File::open(src)?);
    let mut writer = bgzf::io::Writer::new(File::create(dst)?);

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut buf = String::new();
    let mut record_idx: u32 = 0;
    let mut line_start = writer.virtual_position();

    loop {
        buf.clear();
        if reader.read_line(&mut buf)? == 0 {
            break;
        }
        if buf.starts_with('#') {
            writer.write_all(buf.as_bytes())?;
            line_start = writer.virtual_position();
            continue;
        }

        let trimmed = buf.trim_end_matches('\n');
        let mut parts = trimmed.splitn(5, '\t');
        let chrom = parts.next().context("missing CHROM")?;
        let pos_str = parts.next().context("missing POS")?;
        let _orig_id = parts.next().context("missing ID")?;
        let ref_bases = parts.next().context("missing REF")?;
        let rest = parts.next().unwrap_or("");

        record_idx += 1;
        let new_id = format!("rs{}", 1_000_000 + record_idx);
        let line_out = format!("{chrom}\t{pos_str}\t{new_id}\t{ref_bases}\t{rest}\n");
        writer.write_all(line_out.as_bytes())?;

        let line_end = writer.virtual_position();
        let pos: usize = pos_str.trim().parse()?;
        let start = Position::try_from(pos)?;
        let end_pos = pos + ref_bases.len().saturating_sub(1);
        let end = Position::try_from(end_pos.max(pos))?;
        indexer.add_record(chrom, start, end, Chunk::new(line_start, line_end))?;
        line_start = line_end;
    }

    writer.finish()?;
    let mut idx_writer = tabix::io::Writer::new(File::create(&idx_dst)?);
    idx_writer.write_index(&indexer.build())?;
    Ok(record_idx)
}

// Preprocess an Ensembl GTF into the compact gene-coordinate TSV that gets
// embedded into the binary via `include_str!`. Filters to feature=gene
// + gene_biotype=protein_coding, drops gene symbols that map to more than
// one Ensembl ID (paralog ambiguity), and sorts by chrom then position.
fn preprocess_gtf(src: &Path, dst: &Path) -> anyhow::Result<()> {
    use std::collections::{HashMap, HashSet};
    use std::io::{BufRead as _, BufReader};

    #[derive(Clone)]
    struct Row {
        chrom: String,
        start: u32,
        end: u32,
        strand: char,
        gene_id: String,
    }

    let f = File::open(src)?;
    let gz = flate2::read::GzDecoder::new(f);
    let reader = BufReader::new(gz);

    let mut by_name: HashMap<String, Vec<Row>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }
        if parts[2] != "gene" {
            continue;
        }
        let attributes = parts[8];
        if !attributes.contains("gene_biotype \"protein_coding\"") {
            continue;
        }
        let Some(gene_name) = extract_attr(attributes, "gene_name") else {
            continue;
        };
        let Some(gene_id) = extract_attr(attributes, "gene_id") else {
            continue;
        };
        let chrom = parts[0].to_string();
        let Ok(start) = parts[3].parse::<u32>() else {
            continue;
        };
        let Ok(end) = parts[4].parse::<u32>() else {
            continue;
        };
        let strand = parts[6].chars().next().unwrap_or('.');

        by_name.entry(gene_name).or_default().push(Row {
            chrom,
            start,
            end,
            strand,
            gene_id,
        });
    }

    let mut out: Vec<(String, Row)> = Vec::with_capacity(by_name.len());
    for (name, rows) in by_name {
        let unique_ids: HashSet<&str> = rows.iter().map(|r| r.gene_id.as_str()).collect();
        if unique_ids.len() != 1 {
            continue; // ambiguous symbol — multiple Ensembl IDs share this name
        }
        out.push((name, rows.into_iter().next().unwrap()));
    }
    out.sort_by(|a, b| {
        a.1.chrom
            .cmp(&b.1.chrom)
            .then(a.1.start.cmp(&b.1.start))
            .then(a.0.cmp(&b.0))
    });

    let mut w = std::io::BufWriter::new(File::create(dst)?);
    writeln!(w, "gene_symbol\tchrom\tstart\tend\tstrand\tensembl_id")?;
    for (name, r) in &out {
        writeln!(
            w,
            "{name}\t{}\t{}\t{}\t{}\t{}",
            r.chrom, r.start, r.end, r.strand, r.gene_id
        )?;
    }
    w.flush()?;
    eprintln!("Wrote {} genes to {}", out.len(), dst.display());
    Ok(())
}

fn extract_attr(attributes: &str, key: &str) -> Option<String> {
    let pat = format!("{key} \"");
    let i = attributes.find(&pat)?;
    let start = i + pat.len();
    let rel_end = attributes[start..].find('"')?;
    Some(attributes[start..start + rel_end].to_string())
}

#[test]
#[ignore = "run explicitly: cargo test --release -- --ignored generate_gene_tables"]
fn generate_gene_tables() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    preprocess_gtf(
        &root.join("data/ensembl/GRCh38.115.gtf.gz"),
        &root.join("data/genes_grch38.tsv"),
    )
    .expect("GRCh38 preprocess");
    preprocess_gtf(
        &root.join("data/ensembl/GRCh37.87.gtf.gz"),
        &root.join("data/genes_grch37.tsv"),
    )
    .expect("GRCh37 preprocess");
}

#[test]
#[ignore = "run explicitly: cargo test --release -- --ignored generate_rsid_fixture"]
fn generate_rsid_fixture() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = root.join("tests/data/na12878_chr17_slice.vcf.gz");
    let dst = root.join("tests/data/na12878_chr17_slice_rsids.vcf.gz");
    let n = annotate_rsids_to_file(&src, &dst).expect("annotate");
    eprintln!(
        "Wrote {n} records to {} (rs1000001..rs{})",
        dst.display(),
        1_000_000 + n
    );
}

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
            // --ephemeral so each test starts with an empty registry and
            // doesn't share / clobber the user's real state.json.
            .arg("--ephemeral")
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
    assert_eq!(samples.len(), 2);
    let names: Vec<&str> = samples
        .iter()
        .map(|s| s["name"].as_str().unwrap())
        .collect();
    assert!(names.contains(&"NA12878"));
    assert!(names.contains(&"NA12878_copy"));
    assert!(samples.iter().all(|s| s["build"] == "GRCh38"));
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

// ---------- lookup_rsids ----------
//
// The committed NA12878 GIAB slice has ID="." on every record, so the
// integration tests here exercise the not-found path, empty-cache build,
// input-validation errors, and order preservation. The found-path is
// exercised by `cargo test --release vcf::tests::*` (unit-level cache
// tests) and by manual verification against a real annotated VCF.

#[tokio::test]
async fn lookup_rsids_all_not_found() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "lookup_rsids",
            json!({"sample": "NA12878", "rsids": ["rs429358", "rs7412", "rs1800012"]}),
        )
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["sample"], "NA12878");
    assert_eq!(body["count"], 3);
    assert_eq!(body["found_count"], 0);
    let results = body["results"].as_array().unwrap();
    assert_eq!(results.len(), 3);
    for r in results {
        assert_eq!(r["found"], false);
        assert!(r["variant"].is_null() || !r.as_object().unwrap().contains_key("variant"));
    }
    // Order preserved
    assert_eq!(results[0]["rsid"], "rs429358");
    assert_eq!(results[1]["rsid"], "rs7412");
    assert_eq!(results[2]["rsid"], "rs1800012");
    h.shutdown().await;
}

#[tokio::test]
async fn lookup_rsids_cache_reused_across_calls() {
    // Two back-to-back calls. Second should be fast (cache hit). We can't
    // assert wall-clock from inside the protocol, but a second call
    // succeeding with identical shape demonstrates the cache survives.
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let first = h
        .call_tool(
            "lookup_rsids",
            json!({"sample": "NA12878", "rsids": ["rs429358"]}),
        )
        .await;
    let second = h
        .call_tool(
            "lookup_rsids",
            json!({"sample": "NA12878", "rsids": ["rs429358"]}),
        )
        .await;
    let b1: Value = serde_json::from_str(extract_text(&first)).unwrap();
    let b2: Value = serde_json::from_str(extract_text(&second)).unwrap();
    assert_eq!(b1["count"], b2["count"]);
    assert_eq!(b1["found_count"], b2["found_count"]);
    h.shutdown().await;
}

#[tokio::test]
async fn lookup_rsids_empty_list_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("lookup_rsids", json!({"sample": "NA12878", "rsids": []}))
        .await;
    assert!(is_error_response(&resp));
    h.shutdown().await;
}

#[tokio::test]
async fn lookup_rsids_too_many_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let rsids: Vec<String> = (1..=101).map(|i| format!("rs{i}")).collect();
    let resp = h
        .call_tool("lookup_rsids", json!({"sample": "NA12878", "rsids": rsids}))
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(msg.contains("max") || msg.contains("100"), "msg: {msg}");
    h.shutdown().await;
}

// ---------- query_gene ----------

#[tokio::test]
async fn query_gene_col1a1_against_slice() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("query_gene", json!({"sample": "NA12878", "gene": "COL1A1"}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["gene"], "COL1A1");
    assert_eq!(body["build"], "GRCh38");
    assert_eq!(body["ensembl_id"], "ENSG00000108821");
    assert_eq!(body["flank_bp"], 0);
    // Gene coords come from the embedded table — chrom in the table is bare "17"
    // but the response.chrom is normalized to whatever the file uses ("chr17").
    assert_eq!(body["chrom"], "chr17");
    let start = body["start"].as_u64().unwrap();
    let end = body["end"].as_u64().unwrap();
    assert!(
        start >= 50_180_000 && end <= 50_205_000,
        "coords: {start}-{end}"
    );
    assert!(
        body["count"].as_u64().unwrap() > 0,
        "expected variants in COL1A1"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn query_gene_with_flank_expands_window() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let r0 = h
        .call_tool(
            "query_gene",
            json!({"sample": "NA12878", "gene": "COL1A1", "flank_bp": 0}),
        )
        .await;
    let r1 = h
        .call_tool(
            "query_gene",
            json!({"sample": "NA12878", "gene": "COL1A1", "flank_bp": 5000}),
        )
        .await;
    let b0: Value = serde_json::from_str(extract_text(&r0)).unwrap();
    let b1: Value = serde_json::from_str(extract_text(&r1)).unwrap();
    // Gene coords reported should be identical (they're the gene's, not the window's).
    assert_eq!(b0["start"], b1["start"]);
    assert_eq!(b0["end"], b1["end"]);
    // But the variant count with flank should be >= without.
    assert!(
        b1["count"].as_u64().unwrap() >= b0["count"].as_u64().unwrap(),
        "flanked count should not decrease"
    );
    assert_eq!(b1["flank_bp"], 5000);
    h.shutdown().await;
}

#[tokio::test]
async fn query_gene_unknown_gene_errors_with_suggestion() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("query_gene", json!({"sample": "NA12878", "gene": "BRCA22"}))
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains("BRCA22"),
        "msg should echo the bad symbol: {msg}"
    );
    // Distance-1: "BRCA22" → "BRCA2" should be in suggestions
    assert!(
        msg.contains("BRCA2") || msg.contains("Did you mean"),
        "expected a suggestion: {msg}"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn query_gene_case_insensitive() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("query_gene", json!({"sample": "NA12878", "gene": "col1a1"}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    // Canonical case from the table
    assert_eq!(body["gene"], "COL1A1");
    h.shutdown().await;
}

#[tokio::test]
async fn query_gene_gene_outside_slice_returns_empty() {
    // BRCA1 is on chr17 too but well outside the 50.1–50.3 Mbp slice.
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("query_gene", json!({"sample": "NA12878", "gene": "BRCA1"}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["gene"], "BRCA1");
    assert_eq!(body["count"], 0);
    assert_eq!(body["variants"].as_array().unwrap().len(), 0);
    h.shutdown().await;
}

#[tokio::test]
async fn lookup_rsids_unknown_sample_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool("lookup_rsids", json!({"sample": "NOPE", "rsids": ["rs1"]}))
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains("NOPE"),
        "error should mention bad sample: {msg}"
    );
    h.shutdown().await;
}

// ---------- compare_samples ----------
//
// Both fixture samples (NA12878 and NA12878_copy) point at the same VCF, so
// any query should return identical per-sample genotypes — that's what these
// tests assert.

#[tokio::test]
async fn compare_samples_region_identical_for_duplicate() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "compare_samples",
            json!({
                "samples": ["NA12878", "NA12878_copy"],
                "query": {"chrom": "chr17", "start": 50180000, "end": 50210000}
            }),
        )
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["query_type"], "region");
    assert_eq!(body["samples"][0], "NA12878");
    assert_eq!(body["samples"][1], "NA12878_copy");
    assert_eq!(body["truncated"], false);
    let entries = body["results"].as_array().unwrap();
    assert!(!entries.is_empty(), "expected variants in COL1A1 region");

    for v in entries {
        let results = v["results"].as_array().unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0]["sample"], "NA12878");
        assert_eq!(results[1]["sample"], "NA12878_copy");
        assert_eq!(results[0]["found"], true);
        assert_eq!(results[1]["found"], true);
        assert_eq!(results[0]["genotype"], results[1]["genotype"]);
        assert_eq!(
            results[0]["genotype_alleles"],
            results[1]["genotype_alleles"]
        );
        assert_eq!(results[0]["depth"], results[1]["depth"]);
    }
    let positions: Vec<u64> = entries.iter().map(|e| e["pos"].as_u64().unwrap()).collect();
    let mut sorted = positions.clone();
    sorted.sort();
    assert_eq!(positions, sorted, "entries should be position-sorted");
    h.shutdown().await;
}

#[tokio::test]
async fn compare_samples_rsids_preserves_input_order() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "compare_samples",
            json!({
                "samples": ["NA12878", "NA12878_copy"],
                "query": {"rsids": ["rs1", "rs2", "rs3"]}
            }),
        )
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["query_type"], "rsids");
    let entries = body["results"].as_array().unwrap();
    assert_eq!(entries.len(), 3);
    assert_eq!(entries[0]["rsid"], "rs1");
    assert_eq!(entries[1]["rsid"], "rs2");
    assert_eq!(entries[2]["rsid"], "rs3");
    // GIAB slice has ID="." so all samples report not-found.
    for e in entries {
        let per_sample = e["results"].as_array().unwrap();
        assert_eq!(per_sample.len(), 2);
        assert_eq!(per_sample[0]["found"], false);
        assert_eq!(per_sample[1]["found"], false);
    }
    h.shutdown().await;
}

#[tokio::test]
async fn compare_samples_too_few_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "compare_samples",
            json!({
                "samples": ["NA12878"],
                "query": {"chrom": "chr17", "start": 50180000, "end": 50210000}
            }),
        )
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(msg.contains('2'), "msg: {msg}");
    h.shutdown().await;
}

// ---------- add_sample / add_samples_from_folder / remove_sample ----------

fn fixture_vcf_path() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/na12878_chr17_slice.vcf.gz")
        .to_string_lossy()
        .into_owned()
}

#[tokio::test]
async fn add_sample_registers_with_explicit_args() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    // Existing registry has NA12878 and NA12878_copy from the TOML import.
    let path = fixture_vcf_path();
    let resp = h
        .call_tool(
            "add_sample",
            json!({"path": path, "name": "freshly_added", "build": "GRCh38",
                   "description": "added via tool call"}),
        )
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    assert_eq!(body["name"], "freshly_added");
    assert_eq!(body["build"], "GRCh38");
    assert_eq!(body["description"], "added via tool call");

    // list_samples should now show three entries.
    let resp2 = h.call_tool("list_samples", json!({})).await;
    let samples: Vec<Value> = serde_json::from_str(extract_text(&resp2)).unwrap();
    let names: Vec<&str> = samples
        .iter()
        .map(|s| s["name"].as_str().unwrap())
        .collect();
    assert!(names.contains(&"freshly_added"), "got {names:?}");
    h.shutdown().await;
}

#[tokio::test]
async fn add_sample_auto_derives_name_and_detects_build() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let path = fixture_vcf_path();
    let resp = h
        .call_tool("add_sample", json!({"path": path, "build": "GRCh38"}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    // The slice file is "na12878_chr17_slice.vcf.gz". The TOML bootstrap
    // already registered it as "NA12878"; calling add_sample WITHOUT an
    // explicit name + same canonical path triggers idempotency, returning
    // the existing entry. So we may see either the derived name or the
    // pre-registered one — both are correct outcomes.
    let name = body["name"].as_str().unwrap();
    assert!(
        name == "NA12878"
            || name == "na12878_chr17_slice"
            || name.starts_with("na12878_chr17_slice_"),
        "got {name}"
    );
    assert_eq!(body["build"], "GRCh38");
    h.shutdown().await;
}

#[tokio::test]
async fn add_sample_idempotent_without_explicit_name() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let path = fixture_vcf_path();
    // First call without name registers (or returns existing under canonical
    // path — TOML import already registered it as "NA12878").
    let r1 = h
        .call_tool("add_sample", json!({"path": &path, "build": "GRCh38"}))
        .await;
    let r2 = h
        .call_tool("add_sample", json!({"path": &path, "build": "GRCh38"}))
        .await;
    let b1: Value = serde_json::from_str(extract_text(&r1)).unwrap();
    let b2: Value = serde_json::from_str(extract_text(&r2)).unwrap();
    assert_eq!(
        b1["name"], b2["name"],
        "calls without explicit name should be idempotent"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn add_sample_rejects_bad_path() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "add_sample",
            json!({"path": "Z:\\does\\not\\exist.vcf.gz", "build": "GRCh38"}),
        )
        .await;
    assert!(is_error_response(&resp));
    h.shutdown().await;
}

#[tokio::test]
async fn add_sample_rejects_non_vcf_extension() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    // Point at a real file that isn't .vcf.gz (use this very test source).
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("Cargo.toml")
        .to_string_lossy()
        .into_owned();
    let resp = h
        .call_tool("add_sample", json!({"path": path, "build": "GRCh38"}))
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains(".vcf.gz") || msg.contains("filename"),
        "msg: {msg}"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn remove_sample_works() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    // Pre-add a removable entry.
    let path = fixture_vcf_path();
    h.call_tool(
        "add_sample",
        json!({"path": &path, "name": "to_be_removed", "build": "GRCh38"}),
    )
    .await;
    let r = h
        .call_tool("remove_sample", json!({"name": "to_be_removed"}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&r)).unwrap();
    assert_eq!(body["name"], "to_be_removed");

    // Removing again returns the not-removed sentinel.
    let r2 = h
        .call_tool("remove_sample", json!({"name": "to_be_removed"}))
        .await;
    let body2: Value = serde_json::from_str(extract_text(&r2)).unwrap();
    assert_eq!(body2["removed"], false);
    h.shutdown().await;
}

#[tokio::test]
async fn add_samples_from_folder_finds_the_slice() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let folder = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .to_string_lossy()
        .into_owned();
    let resp = h
        .call_tool("add_samples_from_folder", json!({"folder": folder}))
        .await;
    let body: Value = serde_json::from_str(extract_text(&resp)).unwrap();
    // tests/data contains the slice + .tbi + config.toml (the .tbi and
    // config.toml don't end in .vcf.gz so they're skipped at collection time).
    // The slice was already registered via the TOML import on harness start,
    // so the auto-derive idempotency path applies — it returns the existing
    // entry (NA12878) without erroring.
    let scanned = body["scanned"].as_u64().unwrap();
    assert!(scanned >= 1, "expected at least one .vcf.gz");
    let registered = body["registered"].as_array().unwrap();
    let skipped = body["skipped"].as_array().unwrap();
    assert_eq!(
        registered.len() + skipped.len(),
        scanned as usize,
        "registered + skipped should equal scanned"
    );
    h.shutdown().await;
}

#[tokio::test]
async fn compare_samples_unknown_sample_errors() {
    let mut h = McpHarness::start(&fixture_config_path()).await;
    let resp = h
        .call_tool(
            "compare_samples",
            json!({
                "samples": ["NA12878", "GHOST"],
                "query": {"chrom": "chr17", "start": 50180000, "end": 50210000}
            }),
        )
        .await;
    assert!(is_error_response(&resp));
    let msg = resp["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains("GHOST"),
        "error should mention bad sample: {msg}"
    );
    h.shutdown().await;
}
