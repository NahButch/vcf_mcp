mod config;
mod error;
mod tools;
mod vcf;

use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Arc;

use clap::{Parser, Subcommand};
use rmcp::ServiceExt;
use rmcp::transport::stdio;
use tracing_subscriber::EnvFilter;

use crate::config::{Config, default_config_path};
use crate::tools::VcfServer;

#[derive(Parser)]
#[command(name = "vcf-mcp", version, about = "MCP server for querying VCF files")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Start the MCP server, speaking JSON-RPC over stdio.
    Serve {
        /// Path to the config file. Defaults to the OS-standard config location.
        #[arg(short, long)]
        config: Option<PathBuf>,

        /// Validate config and exit without starting the server.
        #[arg(long)]
        check: bool,
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
            // Print full error chain to stderr for the operator.
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
        Command::Serve { config, check } => {
            let path = match config {
                Some(p) => p,
                None => default_config_path().ok_or_else(|| {
                    "could not resolve a default config path on this platform; pass --config"
                        .to_string()
                })?,
            };
            tracing::info!(config = %path.display(), "loading config");
            let cfg = Config::load(&path)?;
            cfg.validate()?;
            tracing::info!(samples = cfg.samples.len(), "config validated");

            if check {
                tracing::info!("--check requested; exiting before server startup");
                return Ok(());
            }

            let rt = tokio::runtime::Builder::new_multi_thread()
                .enable_all()
                .build()?;
            rt.block_on(serve_stdio(Arc::new(cfg)))?;
            Ok(())
        }
    }
}

async fn serve_stdio(cfg: Arc<Config>) -> Result<(), Box<dyn std::error::Error>> {
    let server = VcfServer::new(cfg);
    tracing::info!("starting MCP server on stdio");
    let service = server.serve(stdio()).await?;
    service.waiting().await?;
    tracing::info!("MCP server stopped");
    Ok(())
}
