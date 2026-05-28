# Known Issues

## Cold-cache: first query against an unindexed sample is gated by full index build, not by `query_timeout_secs`

**Status:** open
**Severity:** medium (correctness of the timeout contract / first-query UX)
**Found:** 2026-05-27, via stress testing — first `query_region` on a cold sample took ~4 min for a 200bp window despite `query_timeout_secs = 30`.

### Symptom

The first `query_region` (or any indexed query) against a sample whose tabix
index is not yet cached takes as long as the full in-memory index build
(~4 min for a whole-genome VCF) before any result is returned — even for a
tiny window like 200bp. This exceeds the stated 30s `query_timeout_secs`.

### Root cause

The 30s timeout is *not* bypassed. `query_region` wraps
`spawn_blocking(query_region_blocking)` in `tokio::time::timeout(30s, …)`
(`src/vcf.rs` ~L95–L108), and the in-memory tabix build
(`build_tabix_index_in_memory`, which streams the entire bgzipped VCF,
`src/vcf.rs` ~L1875) runs *inside* that blocking task. At 30s the client does
receive `QueryTimeout`. But three structural facts make that protection
ineffective on a cold index:

1. **`spawn_blocking` tasks are not cancellable.** When the timeout future
   elapses and we return `QueryTimeout`, the blocking closure keeps running to
   completion (~4 min). The work is not abandoned — it's just no longer awaited.

2. **The index is cached only *after* the build completes**
   (`set_tabix_index`, `src/vcf.rs` ~L1940). Any retry issued during the build
   window finds `get_tabix_index() == None` and starts *another* independent
   full-file build. There is no single-flight / dedup, so retries also contend
   with the background warmup build (`warmup_one_sample`) for the same file.

3. **Net effect:** no query against that sample can succeed until *some* build
   finishes, regardless of how small the requested window is. A 200bp window
   pays the whole-genome index-construction cost. Time-to-first-successful
   -result is gated by index build time, not by `query_timeout_secs`.

The background warmup (`warmup_samples_in_background`) does build and cache the
index, but it runs sequentially across samples (~4 min each) and a user query
can easily race ahead of it, triggering a redundant build.

### Why the two natural hypotheses are/aren't right

- "The in-memory tabix build bypasses the timeout." — **No.** The build is
  inside the timeout-wrapped task; the timeout fires at 30s.
- "The timeout applies only post-build." — **Effectively yes, in user-impact
  terms.** The 30s bound governs the RPC response, but because the build is
  uncancellable, uncached-until-complete, and not single-flighted, a successful
  query result cannot be obtained until the full build completes.

### Fix options (not yet implemented)

- **Single-flight the index build.** Track a per-sample "build in progress"
  state so concurrent queries await the one in-flight build instead of each
  starting their own. This alone removes the redundant-build contention.
- **Distinguish "index building" from "query timed out."** Return an
  actionable, distinct error during a cold build (e.g. "index for `{sample}`
  is building, retry shortly") rather than a generic `QueryTimeout`.
- **Separate `index_build_timeout` from `query_timeout`.** The 30s query budget
  conceptually applies to querying a *ready* index; first-time construction is a
  different operation and arguably deserves its own (larger) budget or no
  client-facing timeout at all (build proceeds in background, query returns a
  "building" status).
- **Persist the built index to a cache directory** (NOT the source folder — see
  the in-memory/no-source-write constraint) so cold starts after a Claude
  Desktop respawn don't rebuild from scratch.
