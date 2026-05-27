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

Exposes eight tools to an MCP client:

| Tool | Description |
|---|---|
| `add_sample` | Register a single VCF file at runtime. Validates BGZF magic, header structure, and runs a tabix probe before accepting. Auto-detects genome build from the VCF header; auto-derives sample name from the filename. |
| `add_samples_from_folder` | Scan a folder for `*.vcf.gz` and register each that passes validation. Optional recursive walk. Default cap 50 files, hard cap 200. |
| `remove_sample` | Unregister a sample by name (doesn't delete the file). |
| `list_samples` | List currently registered samples. |
| `query_region` | Return variant calls overlapping a chromosomal region. |
| `lookup_rsids` | Look up variants by dbSNP rsid; lazy per-sample cache. |
| `query_gene` | Return variants in a gene's coordinates using the embedded Ensembl 115/87 gene table. |
| `compare_samples` | Run the same query across multiple samples and merge per-sample results. |

The tools return structured JSON. An MCP client like Claude Desktop calls them
on demand as the model reasons about your prompts. **Zero-config startup is
supported** — point Claude Desktop at the binary, then in a chat paste a VCF
path (or a folder of VCFs) and ask Claude to register them. The server
auto-detects genome build and persists registrations to a state file so they
survive restarts.

## Quick start

```bash
git clone <repo> vcf-mcp && cd vcf-mcp
cargo build --release
# (Optional) generate a tabix index if your VCF doesn't have one:
cargo run --release --example index_vcf -- /path/to/your.vcf.gz
# Wire the binary into your MCP client (see "Claude Desktop" below).
# No config file needed — just chat with Claude and paste your VCF path
# when prompted; the server registers it via add_sample.
```

## Prerequisites

- A **bgzipped VCF** (`.vcf.gz`), or several. No `.tbi` index file is
  required — vcf-mcp builds the index in memory at registration time. (If a
  `.tbi` already sits alongside the file it's used as a fast-load shortcut,
  but it's optional.)
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

## Configure (optional)

There are three ways to register samples; pick whichever fits your workflow:

**1. In chat (recommended)** — say "add this VCF: D:\path\file.vcf.gz" and
Claude calls `add_sample` for you. The server validates the file and remembers
it. Registrations persist across server restarts via a small auto-managed
state file:

| OS | Default state file |
|---|---|
| Linux | `~/.local/share/vcf-mcp/state.json` |
| macOS | `~/Library/Application Support/vcf-mcp/state.json` |
| Windows | `%LOCALAPPDATA%\vcf-mcp\data\state.json` |

You never edit this file directly; the server writes it. Override the location
with `--state-file <path>` or disable persistence entirely with `--ephemeral`.

**2. TOML bootstrap** — for an existing setup, or to seed several samples at
once. Pass `--config /path/to/config.toml`; the file is imported into the
registry on startup. After import, normal `add_sample` / `remove_sample`
operations still work and persist to the state file. Example (also see
[`config.example.toml`](config.example.toml)):

```toml
[[samples]]
name = "me"                                    # short identifier you'll use in prompts
vcf_path = "person_genome_001.snp-indel.vcf.gz"  # absolute, or relative to this config file
build = "GRCh38"                               # or "GRCh37"
description = "Whole-genome sequencing, 2025-06-24"
```

**3. Folder scan** — say "register everything in D:\cohort\" and Claude calls
`add_samples_from_folder`. Default cap 50 files, hard cap 200; see the tool
reference below.

Whichever path you choose, each VCF must be **bgzipped** (`.vcf.gz`). No
`.tbi` is required — the server builds the index in memory when the sample is
registered.

To validate startup without serving:

```bash
vcf-mcp serve --check
```

### Building a `.tbi` index (optional)

vcf-mcp builds the tabix index in memory on first registration when a `.tbi`
isn't present, so you usually don't need to pre-index anything. If you'd
prefer to write a `.tbi` to disk anyway (e.g. for use with `tabix` or
`bcftools`):

```bash
cargo run --release --example index_vcf -- /path/to/your.vcf.gz
```

This streams the file and writes `<your.vcf.gz>.tbi`. Throughput is roughly
3 million records/second on a typical desktop; a 30× WGS file with ~12M
variants indexes in under 5 seconds. vcf-mcp will use that file directly
on next `add_sample` rather than rebuilding.

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
      "args": ["serve"]
    }
  }
}
```

The minimal `serve` form uses state-file persistence — registrations made via
chat survive restarts. To bootstrap from a TOML config, add
`"--config", "/path/to/config.toml"` to the args. To run entirely without
persistence (each session blank), add `"--ephemeral"`.

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

### `add_sample`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `path` | string | ✓ | Absolute path to a bgzipped, tabix-indexed VCF |
| `name` | string | | Auto-derived from filename if omitted |
| `build` | string | | "GRCh37" or "GRCh38"; auto-detected from VCF header if omitted |
| `description` | string | | Human-readable note |

Validation chain (short-circuits on the first failure, cheapest first):

1. Path canonicalizes (resolves symlinks)
2. Path exists, is a regular file, readable
3. Filename ends `.vcf.gz` (case-insensitive)
4. First 4 bytes match BGZF magic `1F 8B 08 04` (rejects plain gzip)
5. Tabix index is acquired: `.tbi` is loaded from disk if present, **otherwise built in memory** — no `.tbi` is ever written next to the VCF
6. `noodles_vcf::io::IndexedReader` opens the file with the acquired index
7. Header has `#CHROM` line, ≥1 sample column, ≥1 `##contig`
8. Tabix probe on the first indexed contig proves data ↔ index consistency

The in-memory tabix index is cached on the server (per sample) for the
lifetime of the process, so subsequent queries don't rebuild. Indexes are
NOT persisted across restarts — they're rebuilt on first query against
each sample, which is fast (~3–4 s for a 30× WGS file).

Build detection cascade: `##reference=` substring → `##contig=<...assembly=...>`
field → chr1 length heuristic (GRCh38: 248,956,422; GRCh37: 249,250,621).

Re-adding the same canonical path **without** an explicit name returns the
existing entry (idempotent). Re-adding **with** an explicit name registers a
new entry — useful for viewing the same file under multiple labels.

### `add_samples_from_folder`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `folder` | string | ✓ | Absolute path to a directory |
| `recursive` | bool | | Walk subdirectories. Default false |
| `max_files` | integer | | Max files to register; default 50, hard cap 200 |

If the folder contains more `.vcf.gz` files than `max_files`, the call errors
with a hint telling the caller exactly which `max_files` value would let the
scan proceed. Per-file validation failures **don't** abort the scan — they
collect in a `skipped` array with the per-file reason:

```json
{
  "folder": "D:\\cohort",
  "scanned": 12,
  "registered": [
    {"name": "person_genome_001", "build": "GRCh38", "description": "", "vcf_path": "..."}
  ],
  "skipped": [
    {"path": "D:\\cohort\\malformed.vcf.gz", "reason": "BGZF magic mismatch ..."}
  ]
}
```

### `remove_sample`

| Arg | Type | Required | Notes |
|---|---|---|---|
| `name` | string | ✓ | Registered sample name |

Drops the sample from the registry (and from the state file). Returns the
removed sample's summary, or `{"removed": false, "name": "..."}` if no such
sample was registered.

### `list_samples`

No arguments.

```json
[
  {"name": "me", "build": "GRCh38", "description": "...", "vcf_path": "..."},
  {"name": "spouse", "build": "GRCh38", "description": "...", "vcf_path": "..."}
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
