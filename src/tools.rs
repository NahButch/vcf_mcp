use std::path::PathBuf;
use std::sync::Arc;

use rmcp::{
    ErrorData as McpError, ServerHandler,
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{
        CallToolResult, Content, Implementation, ProtocolVersion, ServerCapabilities, ServerInfo,
    },
    schemars, tool, tool_handler, tool_router,
};
use serde::{Deserialize, Serialize};

use crate::config::Sample;
use crate::error::{Category, Error};
use crate::registry::SampleRegistry;
use crate::vcf::{
    self, AddSampleArgs, AddSamplesFromFolderArgs, CompareQuery, CompareSamplesArgs,
    LookupRsidsArgs, QueryGeneArgs, QueryRegionArgs, RsidCache,
};

#[derive(Clone)]
pub struct VcfServer {
    registry: Arc<SampleRegistry>,
    rsid_cache: Arc<RsidCache>,
    allowed_roots: Arc<Vec<PathBuf>>,
    // Populated and read by the #[tool_router] / #[tool_handler] proc-macros.
    #[allow(dead_code)]
    tool_router: ToolRouter<VcfServer>,
}

#[derive(Serialize)]
struct SampleSummary<'a> {
    name: &'a str,
    build: &'a str,
    description: &'a str,
    vcf_path: String,
}

impl<'a> From<&'a Sample> for SampleSummary<'a> {
    fn from(s: &'a Sample) -> Self {
        SampleSummary {
            name: &s.name,
            build: &s.build,
            description: &s.description,
            vcf_path: s.vcf_path.to_string_lossy().into_owned(),
        }
    }
}

#[tool_router]
impl VcfServer {
    pub fn new(registry: Arc<SampleRegistry>, allowed_roots: Vec<PathBuf>) -> Self {
        Self {
            registry,
            rsid_cache: Arc::new(RsidCache::default()),
            allowed_roots: Arc::new(allowed_roots),
            tool_router: Self::tool_router(),
        }
    }

    #[tool(
        description = "Return the running server's identity: vcf-mcp version, embedded Ensembl gene-table release per build, and the current count of registered samples. Use this when the user asks what version they're talking to or what reference data is in play."
    )]
    async fn server_info(&self) -> Result<CallToolResult, McpError> {
        tool_ok(&vcf::server_info(&self.registry))
    }

    #[tool(
        description = "List the VCF samples currently registered on this server. Returns name, genome build, description, and absolute vcf_path for each. Returns an empty list if no samples are registered yet — use add_sample to register one."
    )]
    async fn list_samples(&self) -> Result<CallToolResult, McpError> {
        let samples = self.registry.list();
        let summaries: Vec<SampleSummary<'_>> = samples.iter().map(SampleSummary::from).collect();
        tool_ok(&summaries)
    }

    #[tool(
        description = "Register a single bgzipped, tabix-indexed VCF file as a queryable sample. Pass an absolute file path (must end in .vcf.gz and have a matching .tbi alongside). The server validates BGZF magic, opens the file with noodles, checks header structure, and runs a tabix probe before accepting. Genome build is auto-detected from the VCF header when possible; pass `build` explicitly if detection fails. Name is auto-derived from the filename; pass `name` to override. Re-adding the same path is idempotent — the existing entry is returned, no duplicate. On name collision with a different file, a random suffix is appended."
    )]
    async fn add_sample(
        &self,
        Parameters(args): Parameters<AddSampleParams>,
    ) -> Result<CallToolResult, McpError> {
        let result = vcf::add_sample(
            self.registry.clone(),
            self.allowed_roots.clone(),
            AddSampleArgs {
                path: args.path,
                name: args.name,
                build: args.build,
                description: args.description,
            },
        )
        .await;
        match result {
            Ok(s) => tool_ok(&SampleSummary::from(&s)),
            Err(e) => Ok(domain_error_to_result(e)),
        }
    }

    #[tool(
        description = "Scan a folder for .vcf.gz files and register each that passes validation. Set recursive=true to walk subdirectories. Default max_files is 50; server hard cap is 200. If the folder contains more VCFs than max_files, the call errors with a hint telling the caller exactly which value to re-call with. Per-file validation failures don't abort the scan — they collect in a `skipped` array with the reason per file. Use this when the user pastes a folder path instead of file paths."
    )]
    async fn add_samples_from_folder(
        &self,
        Parameters(args): Parameters<AddSamplesFromFolderParams>,
    ) -> Result<CallToolResult, McpError> {
        into_tool_result(
            vcf::add_samples_from_folder(
                self.registry.clone(),
                self.allowed_roots.clone(),
                AddSamplesFromFolderArgs {
                    folder: args.folder,
                    recursive: args.recursive,
                    max_files: args.max_files,
                },
            )
            .await,
        )
    }

    #[tool(
        description = "Unregister a sample by name. Does not delete the underlying VCF file. Returns the removed sample's details if it existed."
    )]
    async fn remove_sample(
        &self,
        Parameters(args): Parameters<RemoveSampleParams>,
    ) -> Result<CallToolResult, McpError> {
        match vcf::remove_sample(self.registry.clone(), args.name.clone()).await {
            Ok(Some(s)) => tool_ok(&SampleSummary::from(&s)),
            Ok(None) => tool_ok(&serde_json::json!({"removed": false, "name": args.name})),
            Err(e) => Ok(domain_error_to_result(e)),
        }
    }

    #[tool(
        description = "Reset the sample registry — drops every registered sample, the per-sample tabix and rsid caches, and rewrites the state file as empty. Destructive: requires `confirm: true`; refuses otherwise. Useful as a workflow / cowork step (clean slate before re-registering a new cohort, end-of-session cleanup, scheduled refresh). Returns the list of what was removed so the caller can mirror it elsewhere or re-register if needed."
    )]
    async fn reset_samples(
        &self,
        Parameters(args): Parameters<ResetSamplesParams>,
    ) -> Result<CallToolResult, McpError> {
        if !args.confirm {
            return Ok(domain_error_to_result(Error::ResetNotConfirmed));
        }
        into_tool_result(vcf::reset_samples(self.registry.clone(), self.rsid_cache.clone()).await)
    }

    #[tool(
        description = "Return variant calls overlapping a chromosome region. Coordinates are 1-based inclusive. Chromosome may be given with or without the 'chr' prefix; the server normalizes to the file's convention. Multi-allelic ALTs are comma-separated. The result is capped at 500 records; when more would match, `truncated` is true."
    )]
    async fn query_region(
        &self,
        Parameters(args): Parameters<QueryRegionParams>,
    ) -> Result<CallToolResult, McpError> {
        into_tool_result(
            vcf::query_region(
                self.registry.clone(),
                QueryRegionArgs {
                    sample: args.sample,
                    chrom: args.chrom,
                    start: args.start,
                    end: args.end,
                },
            )
            .await,
        )
    }

    #[tool(
        description = "Return variants in a gene's coordinates using the embedded Ensembl 115 (GRCh38) / 87 (GRCh37) gene table. Gene symbol is case-insensitive HGNC (e.g. COL1A1, brca1). Optional `flank_bp` extends the queried window on each side. If the gene is unknown, the error suggests close matches."
    )]
    async fn query_gene(
        &self,
        Parameters(args): Parameters<QueryGeneParams>,
    ) -> Result<CallToolResult, McpError> {
        into_tool_result(
            vcf::query_gene(
                self.registry.clone(),
                QueryGeneArgs {
                    sample: args.sample,
                    gene: args.gene,
                    flank_bp: args.flank_bp,
                },
            )
            .await,
        )
    }

    #[tool(
        description = "Run the same query across multiple samples and merge the per-sample results into one position-keyed table. The `query` field is either {\"rsids\": [...]} (uses lookup_rsids; max 100) or {\"chrom\": ..., \"start\": ..., \"end\": ...} (uses query_region; max 10 Mb, 500 records per sample). For each variant any sample has, every sample gets an entry — `found: false` when that sample has no matching call. Minimum 2 samples."
    )]
    async fn compare_samples(
        &self,
        Parameters(args): Parameters<CompareSamplesParams>,
    ) -> Result<CallToolResult, McpError> {
        let inner_query = match args.query {
            CompareQueryParams::Rsids { rsids } => CompareQuery::Rsids(rsids),
            CompareQueryParams::Region { chrom, start, end } => {
                CompareQuery::Region { chrom, start, end }
            }
        };
        into_tool_result(
            vcf::compare_samples(
                self.registry.clone(),
                self.rsid_cache.clone(),
                CompareSamplesArgs {
                    samples: args.samples,
                    query: inner_query,
                },
            )
            .await,
        )
    }

    #[tool(
        description = "Look up variants by dbSNP rsid. First call against a sample triggers a one-time cache build (streams the full VCF, takes seconds-to-minutes depending on file size). Returns one entry per input rsid in input order; entries not found yield {rsid, found: false}. Maximum 100 rsids per call."
    )]
    async fn lookup_rsids(
        &self,
        Parameters(args): Parameters<LookupRsidsParams>,
    ) -> Result<CallToolResult, McpError> {
        into_tool_result(
            vcf::lookup_rsids(
                self.registry.clone(),
                self.rsid_cache.clone(),
                LookupRsidsArgs {
                    sample: args.sample,
                    rsids: args.rsids,
                },
            )
            .await,
        )
    }
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct QueryRegionParams {
    /// Sample name as registered on the server.
    pub sample: String,
    /// Chromosome, with or without 'chr' prefix (e.g. "17" or "chr17").
    pub chrom: String,
    /// 1-based inclusive start position.
    pub start: u32,
    /// 1-based inclusive end position.
    pub end: u32,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct LookupRsidsParams {
    /// Sample name as registered on the server.
    pub sample: String,
    /// dbSNP rsids to look up (e.g. ["rs429358", "rs7412"]). Max 100 per call.
    pub rsids: Vec<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct QueryGeneParams {
    /// Sample name as registered on the server.
    pub sample: String,
    /// HGNC gene symbol, case-insensitive (e.g. "COL1A1").
    pub gene: String,
    /// Optional flank in base pairs added to each side of the gene's coords.
    pub flank_bp: Option<u32>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct CompareSamplesParams {
    /// Sample names to compare; minimum 2.
    pub samples: Vec<String>,
    /// Either a list of rsids or a chromosome region.
    pub query: CompareQueryParams,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum CompareQueryParams {
    Rsids { rsids: Vec<String> },
    Region { chrom: String, start: u32, end: u32 },
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct AddSampleParams {
    /// Absolute path to a bgzipped, tabix-indexed VCF file.
    pub path: String,
    /// Optional name; auto-derived from filename if omitted.
    pub name: Option<String>,
    /// Optional explicit genome build ("GRCh37" or "GRCh38"); auto-detected
    /// from the VCF header if omitted.
    pub build: Option<String>,
    /// Optional human-readable description.
    pub description: Option<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct AddSamplesFromFolderParams {
    /// Absolute path to a folder containing .vcf.gz files.
    pub folder: String,
    /// Walk subdirectories? Default false.
    pub recursive: Option<bool>,
    /// Maximum number of files to register in this call. Default 50, hard cap 200.
    pub max_files: Option<usize>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct RemoveSampleParams {
    /// Registered sample name to remove.
    pub name: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct ResetSamplesParams {
    /// Required acknowledgement. Must be `true`; the tool refuses otherwise.
    /// Drops every registered sample, all per-sample caches, and the
    /// persistent state file.
    pub confirm: bool,
}

/// Helper: a successful tool result wrapping any Serialize value as its
/// JSON text content. Internal serialization failures (our own bug) still
/// flow back as a JSON-RPC McpError::internal_error since they're truly
/// transport-level, not a domain outcome the LLM should reason about.
fn tool_ok<T: Serialize>(value: &T) -> Result<CallToolResult, McpError> {
    let payload =
        serde_json::to_string(value).map_err(|e| McpError::internal_error(e.to_string(), None))?;
    Ok(CallToolResult::success(vec![Content::text(payload)]))
}

/// Helper: collapse a domain Result into the right CallToolResult shape —
/// success → tool_ok; Err → tool-result with isError=true via
/// `domain_error_to_result`. Errors never go back through the JSON-RPC
/// error channel where clients tend to mask them as "failed to call tool".
fn into_tool_result<T: Serialize>(r: Result<T, Error>) -> Result<CallToolResult, McpError> {
    match r {
        Ok(v) => tool_ok(&v),
        Err(e) => Ok(domain_error_to_result(e)),
    }
}

/// Convert a domain Error into a tool result with `isError: true`. The
/// payload includes the kind, category, message, and an optional hint so
/// the LLM has everything it needs to either coach the user (`user_input`),
/// help them diagnose a data/env issue (`user_data`), or offer to file a
/// bug (`unexpected`).
///
/// This returns the error in the tool-result channel (not the JSON-RPC
/// error channel) because clients — including Claude Desktop — often
/// surface protocol errors as a generic "failed to call tool" while
/// faithfully passing tool-result content to the model.
fn domain_error_to_result(e: Error) -> CallToolResult {
    let category = e.category();
    let kind = e.kind();
    let message = e.to_string();
    let hint = e.hint();

    // Log every error that crosses the MCP boundary, tagged with kind +
    // category. Greppable via `error_category=unexpected` in
    // mcp-server-vcf-mcp.log. Unexpected errors warn so they stand out.
    match category {
        Category::Unexpected => tracing::warn!(
            error_kind = kind,
            error_category = category.as_str(),
            error_message = %message,
            "tool returned an unexpected error"
        ),
        _ => tracing::info!(
            error_kind = kind,
            error_category = category.as_str(),
            error_message = %message,
            "tool returned an error"
        ),
    }

    let mut payload = serde_json::json!({
        "kind": kind,
        "category": category.as_str(),
        "message": message,
    });
    if let Some(h) = hint {
        payload["hint"] = serde_json::Value::String(h.to_string());
    }
    CallToolResult::structured_error(payload)
}

#[tool_handler]
impl ServerHandler for VcfServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo::new(ServerCapabilities::builder().enable_tools().build())
            .with_server_info(Implementation::from_build_env())
            .with_protocol_version(ProtocolVersion::V_2024_11_05)
            .with_instructions(
                "VCF query MCP server. To get started, call `add_sample` with a full file path, or `add_samples_from_folder` with a directory containing .vcf.gz files. Genome build is auto-detected from the VCF header and the sample name is derived from the filename. Then query the registered samples via list_samples, query_region, lookup_rsids, query_gene, or compare_samples. Use remove_sample to unregister, or reset_samples (confirm=true required) to wipe the entire registry for a clean slate — useful as a workflow step.\n\nWhen a tool returns an error, the JSON-RPC error `data` field carries a triage hint: `category` is one of `user_input` (the user typed something the tool can't accept — help them fix their args), `user_data` (their VCF file or environment looks problematic — help them diagnose, e.g. re-download, check md5, recheck path), or `unexpected` (this looks like a vcf-mcp bug — consider offering the user to file an issue). `kind` is the stable error variant name (e.g. InvalidVcfFile, QueryTimeout).",
            )
    }
}
