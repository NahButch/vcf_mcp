# vcf-mcp — MCP Submission Roadmap

Status as of 2026-05-28. This tracks getting vcf-mcp into discoverable MCP
catalogs. Check items off as they land.

## Reality check (read first)

vcf-mcp is a **local stdio server** that reads **personal genomic VCF files**
off the local disk. Two facts shape everything below:

1. **It is not a remote/hosted connector.** It speaks JSON-RPC over stdio and
   runs as a local process. So the Anthropic path is the **Desktop Extension
   (MCPB)** submission — *not* the remote-MCP / OAuth directory form. It is also
   independently eligible for the **open MCP Registry**.
2. **It touches health/genetic data.** Anthropic's submission form has an
   explicit *health-data access* compliance question, and reviewers test with
   real prompts. Expect scrutiny. Our story is strong but must be stated
   plainly: **no network egress, no telemetry, no writes to the source folder,
   state limited to a local paths-and-names file, explicit non-medical-device
   disclaimer.** This is a privacy-by-architecture argument — make it loud.

The two biggest rejection causes (from Anthropic's published guidance) are
**missing tool annotations** (~30% of rejections) and **missing/incomplete
privacy policy**. Both are fully in our control.

## Two target catalogs (not mutually exclusive)

| | Anthropic Desktop Extensions Directory | Open MCP Registry |
|---|---|---|
| Format | MCPB bundle (`manifest.json` + binary) | `server.json` metadata |
| Audience | Claude Desktop "Extensions" surface | Any MCP client / aggregators |
| Gatekeeping | Human review (~2 weeks) | Automated namespace-ownership check |
| Submit via | Desktop extension submission form (`clau.de/desktop-extention-submission`) | `mcp-publisher` CLI |
| Hard requirement | Privacy policy in **both** README **and** manifest.json; tool annotations; production-ready | Publicly available install method + `mcpName` validation metadata |
| Namespace | n/a | `io.github.<you>/vcf-mcp` (GitHub-auth'd) |

Recommended order: **Registry first** (fast, automated, low-risk), then the
**Desktop Extensions directory** (slower human review, where the health-data
question lives).

---

## Phase 0 — Make it public (prerequisite for both)

- [ ] Create a **public GitHub repo** and push (`git remote add origin …; git push`).
      Distributing the tool is fine — it transmits no user data; each user runs
      it on their own machine and their own files.
- [ ] Decide the namespace owner: `io.github.<github-username>/vcf-mcp`.
- [ ] Confirm `LICENSE` (MIT already present) and that no personal genome path
      or PII is committed anywhere in history.
- [ ] Tag a release (`v0.1.0`) and attach built binaries (at minimum Windows;
      ideally also macOS + Linux — see cross-platform note in Phase 2).

## Phase 1 — Make every tool submission-grade

- [ ] **Add tool annotations to all 10 tools.** Each needs a human-readable
      `title` plus the correct hint. Proposed map:

  | Tool | `title` | Hints |
  |---|---|---|
  | `server_info` | "Get server info" | `readOnlyHint: true` |
  | `list_samples` | "List samples" | `readOnlyHint: true` |
  | `query_region` | "Query region" | `readOnlyHint: true` |
  | `query_gene` | "Query gene" | `readOnlyHint: true` |
  | `lookup_rsids` | "Look up rsIDs" | `readOnlyHint: true` |
  | `compare_samples` | "Compare samples" | `readOnlyHint: true` |
  | `add_sample` | "Register sample" | `readOnlyHint: false`, `idempotentHint: true`, `openWorldHint: false` |
  | `add_samples_from_folder` | "Register folder" | `readOnlyHint: false`, `openWorldHint: false` |
  | `remove_sample` | "Remove sample" | `readOnlyHint: false`, `destructiveHint: true` |
  | `reset_samples` | "Reset registry" | `readOnlyHint: false`, `destructiveHint: true` |

  (Verify the exact rmcp 1.7 annotation syntax — the `#[tool]` macro takes an
  `annotations(...)` argument.)

- [ ] **Write `PRIVACY.md`** and mirror a privacy section in the README. Must
      cover: data **collected** (none beyond local file paths/names you supply),
      **storage** (local `state.json` only; index/rsid caches in memory),
      **sharing** (none — no network calls), **retention** (until you
      `remove_sample`/`reset_samples` or delete the state file), and a
      **contact**. This is also a manifest.json field for MCPB.
- [ ] **Expand README usage docs** beyond install: concrete example prompts and
      what a successful response looks like (reviewers test with real prompts).
      Largely done — audit against the directory checklist.
- [ ] Resolve / acknowledge the cold-cache caveat in
      [KNOWN_ISSUES.md](../KNOWN_ISSUES.md) so "production-ready" holds up.

## Phase 2 — Package

### MCP Registry (`server.json`)
- [ ] `mcp-publisher init` to scaffold `server.json`; set namespace
      `io.github.<you>/vcf-mcp`, `$schema`, version matching the release.
- [ ] Add the required ownership-validation metadata (`mcpName` in a README
      mention / package label).
- [ ] `mcp-publisher login github` then `mcp-publisher publish --dry-run`.

### MCPB bundle (`manifest.json`)
- [ ] Author `manifest.json`: server name, version, description, the 10 tools
      with annotations, the privacy section, and the run command.
- [ ] Decide **cross-platform support.** Today's binary is Windows-only. MCPB
      bundles a binary/runtime; either ship per-OS bundles or document the
      supported platform. Set up CI release builds for the OSes you commit to.
- [ ] Validate the bundle loads in Claude Desktop's Extensions UI locally.

## Phase 3 — Submit

- [ ] **Registry:** `mcp-publisher publish` (automated; near-instant).
- [ ] **Anthropic directory:** complete the Desktop extension submission form.
      Have ready: tagline + description + use cases, tool inventory with
      annotations confirmation, **health-data handling answer**, privacy-policy
      link, support contact, example prompts for reviewers, logo/icon, GA date,
      and the "production-ready" confirmation.

## Phase 4 — Review

- [ ] Registry: live immediately once validation passes.
- [ ] Anthropic: ~2-week human review. Watch for follow-ups on the health-data
      question and tool annotations specifically.

---

## Document inventory

| Artifact | Needed for | Status |
|---|---|---|
| Public GitHub repo + release | both | ☐ todo (no remote today) |
| Tool annotations in `src/tools.rs` | both | ☐ todo |
| `PRIVACY.md` + README privacy section | Anthropic (req), Registry (nice) | ☐ todo |
| `server.json` | Registry | ☐ todo |
| `manifest.json` (MCPB) | Anthropic | ☐ todo |
| README example prompts | Anthropic | ◑ partial (have prompts; audit) |
| Logo / 256px icon | Anthropic | ☐ todo |
| Cross-platform release binaries | MCPB (depends on scope) | ☐ todo |
| `SECURITY.md` (vuln contact) | recommended | ☐ todo |

## Inputs I need from you

1. **GitHub username/org** for the repo URL and `io.github.*` namespace.
2. **Support/contact channel** for the privacy policy + form (email? GitHub
   issues?).
3. **Which OSes** to support in the MCPB bundle (Windows only, or +macOS/Linux).
4. **Logo/icon** (or I can generate a placeholder spec).
5. Confirmation you're comfortable **listing a genomics tool publicly**, given
   the health-data review — the architecture supports it, but it's your call.
