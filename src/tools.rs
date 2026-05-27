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

use crate::config::{Config, Sample};
use crate::error::Error;
use crate::vcf::{
    self, CompareQuery, CompareSamplesArgs, LookupRsidsArgs, QueryGeneArgs, QueryRegionArgs,
    RsidCache,
};

#[derive(Clone)]
pub struct VcfServer {
    config: Arc<Config>,
    rsid_cache: Arc<RsidCache>,
    // Populated and read by the #[tool_router] / #[tool_handler] proc-macros.
    #[allow(dead_code)]
    tool_router: ToolRouter<VcfServer>,
}

#[derive(Serialize)]
struct SampleSummary<'a> {
    name: &'a str,
    build: &'a str,
    description: &'a str,
}

impl<'a> From<&'a Sample> for SampleSummary<'a> {
    fn from(s: &'a Sample) -> Self {
        SampleSummary {
            name: &s.name,
            build: &s.build,
            description: &s.description,
        }
    }
}

#[tool_router]
impl VcfServer {
    pub fn new(config: Arc<Config>) -> Self {
        Self {
            config,
            rsid_cache: Arc::new(RsidCache::default()),
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "List the VCF samples configured on this server.")]
    async fn list_samples(&self) -> Result<CallToolResult, McpError> {
        let summaries: Vec<SampleSummary<'_>> = self
            .config
            .samples
            .iter()
            .map(SampleSummary::from)
            .collect();
        let payload = serde_json::to_string(&summaries)
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        Ok(CallToolResult::success(vec![Content::text(payload)]))
    }

    #[tool(
        description = "Return variant calls overlapping a chromosome region. Coordinates are 1-based inclusive. Chromosome may be given with or without the 'chr' prefix; the server normalizes to the file's convention. Multi-allelic ALTs are comma-separated. The result is capped at 500 records; when more would match, `truncated` is true."
    )]
    async fn query_region(
        &self,
        Parameters(args): Parameters<QueryRegionParams>,
    ) -> Result<CallToolResult, McpError> {
        let result = vcf::query_region(
            self.config.clone(),
            QueryRegionArgs {
                sample: args.sample,
                chrom: args.chrom,
                start: args.start,
                end: args.end,
            },
        )
        .await
        .map_err(map_domain_error)?;
        let payload = serde_json::to_string(&result)
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        Ok(CallToolResult::success(vec![Content::text(payload)]))
    }

    #[tool(
        description = "Return variants in a gene's coordinates using the embedded Ensembl 115 (GRCh38) / 87 (GRCh37) gene table. Gene symbol is case-insensitive HGNC (e.g. COL1A1, brca1). Optional `flank_bp` extends the queried window on each side. If the gene is unknown, the error suggests close matches."
    )]
    async fn query_gene(
        &self,
        Parameters(args): Parameters<QueryGeneParams>,
    ) -> Result<CallToolResult, McpError> {
        let result = vcf::query_gene(
            self.config.clone(),
            QueryGeneArgs {
                sample: args.sample,
                gene: args.gene,
                flank_bp: args.flank_bp,
            },
        )
        .await
        .map_err(map_domain_error)?;
        let payload = serde_json::to_string(&result)
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        Ok(CallToolResult::success(vec![Content::text(payload)]))
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
        let result = vcf::compare_samples(
            self.config.clone(),
            self.rsid_cache.clone(),
            CompareSamplesArgs {
                samples: args.samples,
                query: inner_query,
            },
        )
        .await
        .map_err(map_domain_error)?;
        let payload = serde_json::to_string(&result)
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        Ok(CallToolResult::success(vec![Content::text(payload)]))
    }

    #[tool(
        description = "Look up variants by dbSNP rsid. First call against a sample triggers a one-time cache build (streams the full VCF, takes seconds-to-minutes depending on file size). Returns one entry per input rsid in input order; entries not found yield {rsid, found: false}. Maximum 100 rsids per call."
    )]
    async fn lookup_rsids(
        &self,
        Parameters(args): Parameters<LookupRsidsParams>,
    ) -> Result<CallToolResult, McpError> {
        let result = vcf::lookup_rsids(
            self.config.clone(),
            self.rsid_cache.clone(),
            LookupRsidsArgs {
                sample: args.sample,
                rsids: args.rsids,
            },
        )
        .await
        .map_err(map_domain_error)?;
        let payload = serde_json::to_string(&result)
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        Ok(CallToolResult::success(vec![Content::text(payload)]))
    }
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct QueryRegionParams {
    /// Sample name as configured on the server.
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
    /// Sample name as configured on the server.
    pub sample: String,
    /// dbSNP rsids to look up (e.g. ["rs429358", "rs7412"]). Max 100 per call.
    pub rsids: Vec<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct QueryGeneParams {
    /// Sample name as configured on the server.
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

fn map_domain_error(e: Error) -> McpError {
    match e {
        Error::SampleNotFound(_)
        | Error::InvalidChromosome { .. }
        | Error::InvalidRange { .. }
        | Error::RegionTooLarge { .. }
        | Error::EmptyRsidList
        | Error::TooManyRsids { .. }
        | Error::GeneNotFound(_)
        | Error::TooFewSamples { .. } => McpError::invalid_params(e.to_string(), None),
        _ => McpError::internal_error(e.to_string(), None),
    }
}

#[tool_handler]
impl ServerHandler for VcfServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo::new(ServerCapabilities::builder().enable_tools().build())
            .with_server_info(Implementation::from_build_env())
            .with_protocol_version(ProtocolVersion::V_2024_11_05)
            .with_instructions(
                "VCF query MCP server. Use list_samples to discover available samples.",
            )
    }
}
