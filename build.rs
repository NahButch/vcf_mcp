//! Build script: stamp the binary with the current git commit count and
//! short SHA so `server_info` (and `--version`) can report exactly which
//! build is running. Re-runs whenever `.git/HEAD` changes (i.e., on every
//! new commit or checkout), but not on every file save.
//!
//! If git isn't available or this isn't a git checkout, falls back to
//! `build=0` / `commit=unknown` so the build still succeeds.

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");

    let commit_count = run_git(&["rev-list", "--count", "HEAD"]).unwrap_or_else(|| "0".to_string());

    let short_sha =
        run_git(&["rev-parse", "--short", "HEAD"]).unwrap_or_else(|| "unknown".to_string());

    println!("cargo:rustc-env=VCF_MCP_BUILD={commit_count}");
    println!("cargo:rustc-env=VCF_MCP_COMMIT={short_sha}");
}

fn run_git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    if !output.status.success() {
        return None;
    }
    let s = String::from_utf8(output.stdout).ok()?.trim().to_string();
    if s.is_empty() { None } else { Some(s) }
}
