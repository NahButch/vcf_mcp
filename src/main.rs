mod config;
mod error;
mod genes;
mod registry;
mod tools;
mod vcf;

use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Arc;

use clap::{Parser, Subcommand};
use rmcp::ServiceExt;
use rmcp::transport::stdio;
use tracing_subscriber::EnvFilter;

use crate::config::Config;
use crate::registry::{SampleRegistry, default_state_path};
use crate::tools::VcfServer;

// Composite version shown by `--version`. Includes Cargo's semver + the git
// commit count (build number) + short SHA so the running binary is uniquely
// identifiable from the command line alone.
const FULL_VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    "+build.",
    env!("VCF_MCP_BUILD"),
    ".",
    env!("VCF_MCP_COMMIT"),
);

#[derive(Parser)]
#[command(name = "vcf-mcp", version = FULL_VERSION, about = "MCP server for querying VCF files")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Start the MCP server, speaking JSON-RPC over stdio.
    Serve {
        /// Optional TOML config to bootstrap registered samples on startup.
        /// Without this, samples come from the state file (or empty if
        /// `--ephemeral`).
        #[arg(short, long)]
        config: Option<PathBuf>,

        /// Validate config / startup and exit without serving.
        #[arg(long)]
        check: bool,

        /// Disable disk-backed state persistence. Each server start is blank;
        /// the LLM has to add_sample on every fresh session.
        #[arg(long)]
        ephemeral: bool,

        /// Path to the state file. Defaults to the OS-standard data location.
        #[arg(long)]
        state_file: Option<PathBuf>,

        /// Restrict add_sample paths to directories under one of these roots.
        /// Repeatable. Default: no restriction (Claude reads anywhere the
        /// user account can).
        #[arg(long = "allowed-root")]
        allowed_roots: Vec<PathBuf>,
    },
}

fn main() -> ExitCode {
    // Tracing MUST go to stderr — stdout is reserved for the MCP protocol.
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_writer(std::io::stderr)
        .with_ansi(false)
        .init();

    let cli = Cli::parse();

    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            tracing::error!("{e}");
            let mut source = std::error::Error::source(&*e);
            while let Some(s) = source {
                tracing::error!("caused by: {s}");
                source = s.source();
            }
            ExitCode::from(1)
        }
    }
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    match cli.command {
        Command::Serve {
            config,
            check,
            ephemeral,
            state_file,
            allowed_roots,
        } => {
            // Resolve state file path: explicit > default. Disabled if --ephemeral.
            let state_path = if ephemeral {
                None
            } else {
                state_file.or_else(default_state_path)
            };
            if let Some(ref p) = state_path {
                tracing::info!(state = %p.display(), "using state file");
            } else {
                tracing::info!("ephemeral mode: no state file");
            }

            let registry = Arc::new(SampleRegistry::new(state_path));
            // Load from disk (no-op if missing). Errors here are fatal because
            // they indicate corrupt state.
            registry.load_from_disk()?;
            tracing::info!(samples = registry.len(), "state loaded");

            // Optional bootstrap from a TOML config (import existing samples
            // into the registry on first run).
            if let Some(toml_path) = &config {
                tracing::info!(config = %toml_path.display(), "importing samples from TOML");
                let cfg = Config::load(toml_path)?;
                cfg.validate()?;
                let added = registry.import(cfg.samples)?;
                tracing::info!(imported = added, "TOML bootstrap complete");
            }

            if !allowed_roots.is_empty() {
                tracing::info!(
                    roots = ?allowed_roots,
                    "restricting add_sample to listed roots"
                );
            }

            if check {
                tracing::info!(
                    samples = registry.len(),
                    "--check requested; exiting before server startup"
                );
                return Ok(());
            }

            let rt = tokio::runtime::Builder::new_multi_thread()
                .enable_all()
                .build()?;
            rt.block_on(serve_stdio(registry, allowed_roots))?;
            Ok(())
        }
    }
}

async fn serve_stdio(
    registry: Arc<SampleRegistry>,
    allowed_roots: Vec<PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::sync::Arc;
    let rsid_cache = Arc::new(vcf::RsidCache::default());

    // Kick off background warmup BEFORE serving so the cold-cache cost
    // doesn't land on the first tool call (which would otherwise risk
    // exceeding Claude Desktop's stdio request timeout → pipe teardown
    // → respawn → cold caches again). See warmup_samples_in_background
    // doc comment for the full story.
    vcf::warmup_samples_in_background(registry.clone(), rsid_cache.clone());

    let server = VcfServer::new(registry, rsid_cache, allowed_roots);
    tracing::info!("starting MCP server on stdio");
    let service = server.serve(stdio()).await?;
    service.waiting().await?;
    tracing::info!("MCP server stopped");
    Ok(())
}
