# Known Issues

## Cold-cache first-query latency (re-measured — original "4 min" was a misdiagnosis)

**Status:** mostly resolved / re-scoped
**Severity:** low (was filed as medium based on a wrong build-time estimate)
**Filed:** 2026-05-27 — first `query_region` on a cold sample appeared to take
~4 min for a 200bp window.
**Re-measured:** 2026-05-28 — index construction is **seconds, not minutes**.
The original diagnosis (below) was wrong about the cause.

### What the perf logs actually show

Measured on two real 30× WGS `snp-indel.genome.vcf.gz` files (~4M variants
each), from `mcp-server-vcf-mcp.log`:

| Phase | Time | Notes |
|---|---|---|
| In-memory tabix index build | **~3.3–3.6 s** | what `add_sample` does synchronously |
| rsID cache build | **~3.3–4.1 s** | 4.0–4.3M entries; ~1.0–1.3M entries/sec |
| Per-sample warmup (tabix + rsID) | **~7 s** | |
| Full startup warmup (2 samples) | **~14 s** | sequential |

So a cold first query waits **single-digit seconds** for index construction,
well under any client timeout — not the 4 minutes originally filed.

### What likely caused the original 4-minute observation

Index build was never the bottleneck. The two plausible culprits, both visible
in the environment:

1. **OneDrive on-demand hydration.** One sample lives under a OneDrive folder.
   If the VCF is a "files on-demand" placeholder, the *first* read forces a
   cloud download of a multi-GB file — minutes on a slow link — while every
   subsequent build is ~3.5 s once the file is local. Fix is environmental:
   mark the file **"Always keep on this device."**
2. **Server respawn churn.** The logs show the server starting several times
   within seconds during a session (here, coinciding with our own
   rebuild/sync work). Each respawn re-warms from cold, multiplying the
   apparent wait.

### Residual structural note (low priority)

The single-flight gap is real but now low-impact: a query arriving mid-build
sees `get_tabix_index() == None` and starts its *own* build rather than
awaiting the in-flight warmup build (`set_tabix_index` only caches *after* the
build completes). With ~3.5 s builds this wastes a few seconds of duplicate
work, not minutes. Worth fixing only if cold-start latency ever becomes a
measured problem again; not worth the complexity today.

Also note `spawn_blocking` tasks are not cancellable, so a `query_region`
that times out at `query_timeout_secs` leaves its build running to completion
in the background (which then warms the cache for the next call anyway).

### Mitigations already in place

- **Background warmup at startup** (`warmup_samples_in_background`) builds both
  caches for all registered samples before the first user query.
- **Background rsID warm after `add_sample`** (`warm_rsids_in_background`): the
  tabix index is built synchronously at registration; the rsID cache is now
  warmed in the background for the new sample and any siblings, so the first
  `lookup_rsids` on a freshly added sample doesn't pay the build cost.

### If it ever regresses

- Confirm the VCF is not a OneDrive/cloud placeholder (hydrate it locally).
- Check the logs for repeated `phase="warmup" event="all_start"` within
  seconds — that's respawn churn, an MCP-client/process-lifecycle problem, not
  a build-speed problem.
- Only then consider single-flighting the build or persisting indexes to a
  cache dir (never the source folder).

---

## Building on Windows 11 with Smart App Control (SAC)

**Status:** environmental, not a vcf-mcp bug
**Affects:** developers building from source on Windows 11 with SAC enforced

Windows 11's Smart App Control nondeterministically blocks freshly-compiled,
unsigned executables with `os error 4551` ("An Application Control policy has
blocked this file"). It hits two distinct surfaces in this project's dev loop:

1. **The built `vcf-mcp.exe`** — a clean release build can fail to launch
   immediately after `cargo build --release`, even though the same file ran
   minutes earlier. SAC's verdict is per file-hash, time-dependent, and not
   user-overridable via a SmartScreen-style "Run anyway" button.
2. **Crate build scripts** — every dependency's `build-script-build.exe` under
   `target/release/build/<crate>-<hash>/` is also a fresh unsigned exe. After a
   `rustup update`, all of them get recompiled and any one of them
   (we've observed `proc-macro2`) can hit SAC and abort the build mid-compile.

### Workarounds (no perfect fix while SAC is enforced)

- **After a Rust toolchain bump, do a clean build:**
  ```powershell
  cargo clean
  cargo build --release
  ```
  Fresh compiles produce new file hashes that SAC re-evaluates from scratch.
  Retrying without `cargo clean` keeps hitting the same blocked hashes.
- **After each rebuild, verify the binary actually launches** — don't trust the
  "Finished" message alone:
  ```powershell
  & .\target\release\vcf-mcp.exe --version
  ```
  If it errors with "Application Control policy has blocked this file," another
  `cargo clean && cargo build --release` (or a short wait) usually clears it.
- **Kill child processes before rebuilding** — Windows locks running `.exe`s,
  so MCP clients holding the server process must be told to exit:
  ```powershell
  Get-Process -Name vcf-mcp -ErrorAction SilentlyContinue | Stop-Process -Force
  ```

### What does *not* reliably help

- **Self-signing** the binary with a `New-SelfSignedCertificate` cert plus
  Trusted-Publishers install: in our testing the self-signed signature ended up
  `UnknownError` (no machine-wide trust chain), and SAC's verdict was still
  driven by hash reputation rather than the embedded signature. Some builds ran
  signed; others ran unsigned; the signing wasn't the deciding factor.
- **Retrying without `cargo clean`** after the first block — same hash, same
  verdict.

### The reliable cures (not recommended for a personal dev box)

- **EV (Extended Validation) code-signing certificate** — chains to a CA SAC's
  intelligence recognizes; deterministic allow. Significant cost; appropriate
  for distributing signed releases, overkill for self-built dev binaries.
- **Turn SAC off** — eliminates all of the above. **This is irreversible
  without a Windows reset/reinstall**; never do it casually.
