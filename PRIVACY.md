# Privacy Policy

_Last updated: 2026-05-28_

vcf-mcp is a **local** Model Context Protocol server. It runs as a process on
the local machine and communicates with the MCP client (e.g. Claude Desktop)
over local stdio. This document describes exactly what it does and does not do
with user data.

## Summary

- **No network transmission.** vcf-mcp makes no outbound network connections.
  It does not phone home, send telemetry, or upload any data anywhere.
- **No data collection by the author.** The author/maintainer of vcf-mcp
  receives nothing. There is no analytics, no crash reporting, no usage
  tracking.
- **Genomic data never leaves the local machine.** VCF file contents are read
  from local disk, queried in memory, and returned only to the local MCP
  client.

## What data vcf-mcp handles

| Data | Where it comes from | Where it goes |
|---|---|---|
| VCF file **contents** (variants, genotypes, etc.) | Local files registered | Read into memory to answer queries; results returned to the local MCP client only |
| Tabix / rsID **indexes** | Built in memory from the registered VCF files | Held in RAM for the life of the server process; never written to the source folder |
| **Registry metadata**: sample name, file path, genome build, and an optional description | Provided via `add_sample` / config | Persisted to a local state file (see below) |

vcf-mcp does **not** read, infer, or store any data beyond what is needed to
answer the queries made.

## What is stored on disk, and where

The only thing vcf-mcp writes to disk is a small **state file** recording the
registered samples, so they survive a restart. It contains, per
sample: the **name**, the **absolute file path**, the **genome build**, and
an optional **description**. It does **not** contain genomic variant data.

Default location:

| OS | Path |
|---|---|
| Linux | `~/.local/share/vcf-mcp/state.json` |
| macOS | `~/Library/Application Support/vcf-mcp/state.json` |
| Windows | `%LOCALAPPDATA%\vcf-mcp\data\state.json` |

- Override the location with `--state-file <path>`.
- Disable persistence entirely with `--ephemeral` (nothing is written to disk;
  each session starts empty).

vcf-mcp **never** writes a `.tbi` index, or any other file, next to the source
VCFs. Index construction happens in memory.

## Data sharing

vcf-mcp shares data with **no one**. The only recipient of any output is the
local MCP client it is connected to, on the same machine. What that client
(e.g. Claude Desktop) subsequently does with the query results is governed by
**that client's** own privacy policy, not this one. In particular, an AI
assistant may transmit tool results to its own model provider as part of
generating a response — review the client's privacy terms to understand that
boundary.

## Data retention and deletion

- **Indexes/caches** live only in RAM and vanish when the server process exits.
- **Registry metadata** persists in the state file until removed:
  - `remove_sample` drops a single sample,
  - `reset_samples` (with `confirm: true`) drops all samples and rewrites the
    state file empty,
  - or delete the state file directly.
- vcf-mcp never deletes or modifies the underlying VCF files.

## Permissions and access scope

When started with one or more `--allowed-root` paths, vcf-mcp will only register
files located under those directories; any path outside them is refused. Without
that flag, it can register any readable file path provided. It only ever
reads files explicitly registered or found by an
`add_samples_from_folder` scan.

## Medical disclaimer

vcf-mcp is research / engineering tooling and is **not a medical device**. It
must not be used for clinical diagnosis, treatment decisions, or any purpose
involving patient care. See the [README](README.md#disclaimer) and
[LICENSE](LICENSE) for the full disclaimer.

## Contact

Questions about this policy or vcf-mcp's data handling:
**`NahButch@users.noreply.github.com`**

## Changes

Material changes to this policy will be reflected in this file with an updated
"Last updated" date.
