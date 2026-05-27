# vcf-mcp

A [Model Context Protocol](https://modelcontextprotocol.io) server that lets an
AI assistant query bgzipped, tabix-indexed VCF files — your personal genome, a
cohort sample, anything in VCF format — through natural-language prompts.
Built in Rust on top of [`noodles`](https://crates.io/crates/noodles) for VCF
I/O and [`rmcp`](https://crates.io/crates/rmcp) for the MCP protocol.

## Disclaimer

vcf-mcp is research / engineering tooling. It is **not a medical device** and
must not be used for clinical diagnosis, treatment decisions, or any other
purpose involving patient care. Output may be incomplete, incorrect, or
out-of-date relative to the underlying VCF data. Always verify results against
the source files and authoritative annotations before drawing biological
conclusions. Use of this software is at your own risk; see the [LICENSE](LICENSE)
for the full warranty disclaimer.

## What it does

Exposes five tools to an MCP client:

| Tool | Description |
|---|---|
| `list_samples` | List configured VCF samples (name, build, description). |
| `query_region` | Return variant calls overlapping a chromosomal region. |
| `lookup_rsids` | Look up variants by dbSNP rsid; lazy per-sample cache. |
| `query_gene` | Return variants in a gene's coordinates using the embedded Ensembl 115/87 gene table. |
| `compare_samples` | Run the same query across multiple samples and merge per-sample results. |

The tools return structured JSON. An MCP client like Claude Desktop calls them
on demand as the model reasons about your prompts.

## Quick start

```bash
git clone <repo> vcf-mcp && cd vcf-mcp
cargo build --release
# Generate a tabix index if your VCF doesn't have one:
cargo run --release --example index_vcf -- /path/to/your.vcf.gz
# Edit config.example.toml to point at your VCF, save as config.toml,
# then wire the binary into your MCP client (see "Claude Desktop" below).
```

## Prerequisites

- A **bgzipped, tabix-indexed VCF** (or several). If you don't have a `.tbi`,
  vcf-mcp ships a helper to build one (see below).
- **Rust ≥ 1.85** (transitive deps use Rust 2024 edition features). Install
  from [rustup.rs](https://rustup.rs/).
- An **MCP-capable client**. Tested with Claude Desktop on Windows.

On Windows you'll need either the MSVC build tools or the GNU toolchain
(`rustup target add x86_64-pc-windows-msvc` is the default).

## Build

```bash
cargo build --release
```

The binary lands at `target/release/vcf-mcp` (or `.exe` on Windows).

> If you have `CARGO_TARGET_DIR` set globally, the project's
> [`.cargo/config.toml`](.cargo/config.toml) tries to override it back to
> `./target/`. Cargo's env var still wins over the config file though, so in a
> shell that has the env set, either clear it
> (`unset CARGO_TARGET_DIR` / `$env:CARGO_TARGET_DIR=''`) or pass
> `--target-dir target` on the command line.

## Configure

vcf-mcp reads a TOML config file. Default path:

| OS | Path |
|---|---|
| Linux | `~/.config/vcf-mcp/config.toml` |
| macOS | `~/Library/Application Support/vcf-mcp/config.toml` |
| Windows | `%APPDATA%\vcf-mcp\config.toml` |

Or pass `--config <path>` on the command line. Example (also see
[`config.example.toml`](config.example.toml)):

```toml
[[samples]]
name = "me"                                    # short identifier you'll use in prompts
vcf_path = "patient_001.snp-indel.vcf.gz"      # absolute, or relative to this config file
build = "GRCh38"                               # or "GRCh37"
description = "Whole-genome sequencing, 2025-06-24"

[[samples]]
name = "spouse"
vcf_path = "spouse_001.vcf.gz"
build = "GRCh38"
description = "..."
```

Each `vcf_path` must point at a bgzipped VCF with a matching `.tbi` alongside.

To validate config without starting the server:

```bash
vcf-mcp serve --config /path/to/config.toml --check
```

### Building a `.tbi` index

If you don't have `tabix` or `bcftools` installed:

```bash
cargo run --release --example index_vcf -- /path/to/your.vcf.gz
```

This streams the file and writes `<your.vcf.gz>.tbi`. Throughput is roughly
3 million records/second on a typical desktop; a 30× WGS file with ~12M
variants indexes in under 5 seconds.

## Wire into Claude Desktop

Locate `claude_desktop_config.json`. The cleanest way is **Settings → Developer
→ Edit Config** from inside Claude Desktop — that button opens the file at
whatever path your install actually uses (it differs between the standalone
download and the Microsoft Store / MSIX build on Windows).

Add the `vcf-mcp` entry under `mcpServers` (merge with any existing entries —
don't replace the whole file):

```json
{
  "mcpServers": {
    "vcf-mcp": {
      "command": "/absolute/path/to/vcf-mcp",
      "args": ["serve", "--config", "/absolute/path/to/config.toml"]
    }
  }
}
```

On Windows, escape backslashes (`"D:\\code\\vcf-mcp\\target\\release\\vcf-mcp.exe"`)
or use forward slashes.

**Restart Claude Desktop fully** — right-click the tray icon → Quit. Just
closing the window leaves the process running, and the running process won't
pick up the config change. On the Windows Store build you may need to manually
kill any lingering `claude.exe` processes:

```powershell
Get-Process -Name claude -ErrorAction SilentlyContinue |
  Where-Object { $_.Path -like '*WindowsApps*' } |
  Stop-Process -Force
```

Then open a new chat and try:

```
What samples are configured?
What variants do I have at chr19:44905000-44910000?
Look up rs429358 and rs7412 for me.
What variants do I have in the APOE gene?
What's my ApoE genotype?
```

The last prompt is interesting — it requires composing tool output with
biological knowledge (the absence of rs429358 / rs7412 in the VCF implies
homozygous reference at both, which is the ε3/ε3 genotype).

## Tool reference

All tools take a `sample` argument (the `name` from your config) and return
JSON in the `text` field of an MCP `CallToolResult`. Coordinate conventions are
**1-based inclusive** throughout, matching the VCF spec and command-line tools
like `tabix`.

### `list_samples`

No arguments.

```json
[
  {"name": "me", "build": "GRCh38", "description": "..."},
  {"name": "spouse", "build": "GRCh38", "description": "..."}
]
```

### `query_region`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `sample` | string | ✓ | Configured sample name |
| `chrom` | string | ✓ | With or without `chr` prefix; server normalizes |
| `start` | integer | ✓ | 1-based inclusive |
| `end` | integer | ✓ | 1-based inclusive |

Maximum 10 Mb region. Result is capped at 500 records; when truncated, the
`truncated` flag is `true`.

```json
{
  "sample": "me",
  "chrom": "chr19",
  "start": 44905000,
  "end": 44910000,
  "count": 10,
  "truncated": false,
  "variants": [
    {
      "chrom": "19",
      "pos": 44905579,
      "ref": "T",
      "alt": "G",
      "rsid": "rs405509",
      "genotype": "1/1",
      "genotype_alleles": "GG",
      "depth": 33,
      "gq": 73,
      "filter": "PASS"
    }
  ]
}
```

- `alt` is comma-separated for multi-allelic records.
- gVCF reference-call records have `alt: "."` and `genotype: "0/0"`.
- `filter` is `"PASS"` for variants that pass all filters, `"."` for unfiltered
  records, or a semicolon-joined filter name list.

### `lookup_rsids`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `sample` | string | ✓ | |
| `rsids` | array of string | ✓ | 1–100 entries |

First call against a sample triggers a one-time scan of the VCF to build an
in-memory `rsid → position` cache (a few seconds on a 30× WGS file).
Subsequent calls reuse the cache for O(1) lookup. The cache lives for the
lifetime of the server process.

```json
{
  "sample": "me",
  "count": 3,
  "found_count": 2,
  "results": [
    {"rsid": "rs405509", "found": true, "variant": {/* same shape as query_region */}},
    {"rsid": "rs9999999", "found": false},
    {"rsid": "rs440446", "found": true, "variant": {...}}
  ]
}
```

Results are in input order. Entries not found yield `{rsid, found: false}`.

**Note on "not found"**: VCFs typically carry records only at sites where the
sample differs from reference. An rsid not present usually means the sample
is homozygous reference at that position — *not* that the lookup is broken.
Common variants like rs429358 / rs7412 are absent for everyone with the
common haplotype.

### `query_gene`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `sample` | string | ✓ | |
| `gene` | string | ✓ | HGNC symbol, case-insensitive |
| `flank_bp` | integer | | Bases added to each side of the gene's coords. Default 0. |

Gene coordinates come from the embedded Ensembl release pinned to the
sample's `build`: **release 115** for GRCh38, **release 87** (the GRCh37
archive freeze) for GRCh37. Protein-coding genes only. Symbols that resolve
to more than one Ensembl ID are dropped to avoid ambiguity.

```json
{
  "gene": "APOE",
  "ensembl_id": "ENSG00000130203",
  "chrom": "19",
  "start": 44903787,
  "end": 44909396,
  "build": "GRCh38",
  "flank_bp": 0,
  "count": 9,
  "truncated": false,
  "variants": [/* same shape as query_region */]
}
```

Unknown symbols return an error with Levenshtein-1 distance suggestions:

```
gene "BRCA22" not found in GRCh38. Did you mean: BRCA2?
```

### `compare_samples`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `samples` | array of string | ✓ | Minimum 2 configured sample names |
| `query` | object | ✓ | Either `{"rsids": [...]}` or `{"chrom": ..., "start": ..., "end": ...}` |

Composes `lookup_rsids` or `query_region` across multiple samples and merges
the per-sample results into one position-keyed table.

```json
{
  "samples": ["me", "spouse"],
  "query_type": "rsids",
  "count": 1,
  "truncated": false,
  "results": [
    {
      "rsid": "rs405509",
      "chrom": "19",
      "pos": 44905579,
      "ref": "T",
      "alt": "G",
      "results": [
        {"sample": "me",     "found": true, "genotype": "1/1", "genotype_alleles": "GG", "depth": 33, "gq": 73, "filter": "PASS"},
        {"sample": "spouse", "found": false}
      ]
    }
  ]
}
```

Each merged entry's `results` array contains one item per input sample, in
input order. Samples without a matching `(chrom, pos, ref, alt)` get
`{found: false}`.

## Architecture

- [`src/main.rs`](src/main.rs): clap CLI; the `serve` subcommand bootstraps
  tracing (stderr only — stdout is reserved for MCP JSON-RPC) and runs the
  MCP server over stdio.
- [`src/config.rs`](src/config.rs): TOML config loader; validates VCF + index
  existence at startup; resolves relative `vcf_path` against the config file's
  directory.
- [`src/vcf.rs`](src/vcf.rs): query implementations. All blocking VCF I/O
  happens inside `tokio::task::spawn_blocking` with a 30s timeout. The
  `RsidCache` is `Arc<RwLock<HashMap<sample, Arc<RsidIndex>>>>` for lock-free
  reads after the per-sample cache is warm.
- [`src/genes.rs`](src/genes.rs): embedded gene tables via `include_str!`,
  parsed lazily through `OnceLock` on first lookup per build.
- [`src/tools.rs`](src/tools.rs): rmcp tool definitions. Each tool is an
  `async` method on `VcfServer` annotated with `#[tool]`; the `#[tool_router]`
  and `#[tool_handler]` macros wire them into the `ServerHandler` impl.

Built on:

- [`rmcp`](https://crates.io/crates/rmcp) — official Rust MCP SDK.
- [`noodles`](https://crates.io/crates/noodles) — bioinformatics I/O (VCF,
  bgzf, tabix, csi).
- [`tokio`](https://crates.io/crates/tokio) — async runtime.

## Development

```bash
cargo test --release
cargo clippy --release --all-targets --all-features -- -D warnings
cargo fmt --all -- --check
```

> **Smart App Control on Windows 11**: SAC can silently block freshly-compiled
> unsigned exes. Release binaries placed in `target/release/` are typically
> allowed; debug-mode build scripts in `target/debug/build/.../build-script-build`
> may be blocked, which breaks `cargo test` without `--release`. Either run
> tests in release mode (the CI does this implicitly) or turn SAC off — note
> that disabling SAC is a one-way change without an OS reinstall.

### Examples

- [`examples/index_vcf.rs`](examples/index_vcf.rs) — build a tabix `.tbi`
  index for an existing bgzipped VCF.
- [`examples/slice_vcf.rs`](examples/slice_vcf.rs) — extract a subregion of a
  bgzipped VCF (+ its index) using `noodles`.
- [`examples/annotate_rsids.rs`](examples/annotate_rsids.rs) — rewrite a
  VCF's ID column with synthetic deterministic rsids (used to build a
  `lookup_rsids` test fixture on hosts where SAC permits it).

### Regenerating the gene tables

The committed `data/genes_grch3{7,8}.tsv` files are pre-built from Ensembl
release 115 / 87. To regenerate from a different release:

1. Download fresh GTFs into `data/ensembl/` (gitignored). The URLs encoded in
   the source comments are stable.
2. Run the preprocessor:
   ```bash
   cargo test --release -- --ignored generate_gene_tables
   ```

This re-emits `data/genes_grch3{7,8}.tsv`, which is then embedded into the
binary via `include_str!` at the next build.

## Data sources

- **Gene coordinates**: [Ensembl release 115](https://ftp.ensembl.org/pub/release-115/)
  (GRCh38) and the [GRCh37 archive release 87](https://ftp.ensembl.org/pub/grch37/current/).
  Protein-coding genes only.
- **Test fixture**: [GIAB NA12878 / HG001 v4.2.1 benchmark VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/),
  sliced to chr17:50100000-50300000 (COL1A1 region).

## License

[MIT](LICENSE).
