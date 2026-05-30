//! Append-only error log for `Category::Unexpected` errors.
//!
//! Writes a deliberately minimal, scrubbed JSON Lines record to a
//! platform-stable path (the same `directories::ProjectDirs` scheme as the
//! state file — e.g. `%LOCALAPPDATA%\vcf-mcp\error.log` on Windows). On an
//! Unexpected error the assistant tells the user this path and offers a
//! `report_url`; the user reviews before attaching/submitting.
//!
//! Records contain ONLY: timestamp, vcf-mcp version, error kind, category, and
//! the curated static `hint` text. No paths, sample names, regions, arg
//! values, or error `Display` strings — those can carry PII (file paths
//! especially) and are intentionally excluded.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::error::Error;

const FULL_VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    "+build.",
    env!("VCF_MCP_BUILD"),
    ".",
    env!("VCF_MCP_COMMIT"),
);

/// Platform-stable error log location, mirroring `registry::default_state_path`.
pub fn default_path() -> Option<PathBuf> {
    directories::ProjectDirs::from("", "", "vcf-mcp").map(|p| p.data_dir().join("error.log"))
}

fn record_json(error: &Error) -> serde_json::Value {
    let ts_secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    serde_json::json!({
        "timestamp_unix": ts_secs,
        "version": FULL_VERSION,
        "kind": error.kind(),
        "category": error.category().as_str(),
        "hint": error.hint(),
    })
}

/// Append a scrubbed record for the given error to `path`. Best-effort I/O
/// error; the caller decides whether to log/discard.
pub fn append_record_to(path: &Path, error: &Error) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let line = record_json(error).to_string();
    let mut f = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(f, "{line}")
}

/// Append a scrubbed record to the default platform path. Failures are
/// warn-logged and never propagated — the underlying error is what the
/// caller is already handling; we never want logging to mask it.
pub fn append_record(error: &Error) {
    let Some(path) = default_path() else {
        return;
    };
    if let Err(e) = append_record_to(&path, error) {
        tracing::warn!(
            error_log_path = %path.display(),
            error = %e,
            "failed to append to error log"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn temp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "vcf-mcp-error-test-{}-{}.log",
            std::process::id(),
            name
        ))
    }

    #[test]
    fn default_path_ends_in_error_log() {
        let p = default_path().expect("project dirs should resolve");
        assert!(p.ends_with("error.log"), "got {p:?}");
    }

    #[test]
    fn append_record_writes_jsonl_with_only_safe_fields() {
        let path = temp_path("safe-fields");
        let _ = fs::remove_file(&path);

        append_record_to(&path, &Error::NoSamples).unwrap();
        append_record_to(&path, &Error::QueryTimeout { secs: 30 }).unwrap();

        let content = fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2);

        for line in &lines {
            let v: serde_json::Value = serde_json::from_str(line).unwrap();
            // Required scrubbed fields
            assert!(v["timestamp_unix"].is_u64());
            assert!(v["version"].as_str().unwrap().starts_with("0.1.0+build."));
            assert!(v["kind"].is_string());
            assert!(v["category"].is_string());
            // PII bright lines — these keys must never appear
            for forbidden in ["message", "path", "vcf_path", "sample", "args"] {
                assert!(
                    v.get(forbidden).is_none(),
                    "forbidden key {forbidden:?} in record: {v}"
                );
            }
        }

        // Second line is the QueryTimeout — confirm kind matches.
        let v2: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
        assert_eq!(v2["kind"], "QueryTimeout");
        assert_eq!(v2["category"], "unexpected");

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn append_record_appends_does_not_truncate() {
        let path = temp_path("append");
        let _ = fs::remove_file(&path);

        for _ in 0..3 {
            append_record_to(&path, &Error::NoSamples).unwrap();
        }
        let content = fs::read_to_string(&path).unwrap();
        assert_eq!(content.lines().count(), 3);

        let _ = fs::remove_file(&path);
    }
}
