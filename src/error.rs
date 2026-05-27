use std::path::PathBuf;

use thiserror::Error;

// Several variants are wired in but not yet constructed; they cover the
// surface that the next pass (real VCF queries + remaining tools) will use.
#[allow(dead_code)]
#[derive(Debug, Error)]
pub enum Error {
    #[error("failed to read config from {path}: {source}")]
    ConfigRead {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse config at {path}: {source}")]
    ConfigParse {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },

    #[error(
        "invalid build {build:?} for sample {sample:?}: expected one of \"GRCh37\", \"GRCh38\""
    )]
    InvalidBuild { sample: String, build: String },

    #[error("VCF path does not exist for sample {sample:?}: {path}")]
    VcfMissing { sample: String, path: PathBuf },

    #[error("tabix index not found for sample {sample:?}: expected {expected} alongside {vcf}")]
    IndexMissing {
        sample: String,
        vcf: PathBuf,
        expected: PathBuf,
    },

    #[error("duplicate sample name {0:?} in config")]
    DuplicateSample(String),

    #[error("no samples configured")]
    NoSamples,

    #[error("sample {0:?} not found")]
    SampleNotFound(String),

    #[error("region length {length} exceeds maximum {max} bp; narrow the query")]
    RegionTooLarge { length: u64, max: u64 },

    #[error("invalid range: start ({start}) must be <= end ({end})")]
    InvalidRange { start: u32, end: u32 },

    #[error("chromosome {chrom:?} not found in sample {sample:?}. Available: {available:?}")]
    InvalidChromosome {
        sample: String,
        chrom: String,
        available: Vec<String>,
    },

    #[error("query exceeded {secs}s timeout")]
    QueryTimeout { secs: u64 },

    #[error("rsids list has {count} entries; max {max}")]
    TooManyRsids { count: usize, max: usize },

    #[error("rsids list is empty")]
    EmptyRsidList,

    #[error("{0}")]
    GeneNotFound(String),

    #[error("compare_samples needs at least 2 samples; got {count}")]
    TooFewSamples { count: usize },

    #[error("invalid path {path}: {reason}")]
    PathInvalid { path: PathBuf, reason: String },

    #[error(
        "path {path} is not under any allowed root. Allowed: {allowed:?}. Restart the server without --allowed-root to remove this restriction."
    )]
    PathNotAllowed {
        path: PathBuf,
        allowed: Vec<PathBuf>,
    },

    #[error("file at {path} is not a valid bgzipped VCF: {reason}")]
    InvalidVcfFile { path: PathBuf, reason: String },

    #[error(
        "could not detect genome build from {path}. Pass `build` explicitly as one of \"GRCh37\", \"GRCh38\"."
    )]
    BuildNotDetectable { path: PathBuf },

    #[error(
        "folder scan found {found} VCF files; current limit is {max}. Re-call with max_files={found} (server hard cap is {hard_cap})."
    )]
    TooManyFiles {
        found: usize,
        max: usize,
        hard_cap: usize,
    },

    #[error("state file error at {path}: {message}")]
    StateFile { path: PathBuf, message: String },

    #[error("failed to open VCF {path}: {source}")]
    VcfOpen {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("VCF read error in {path}: {message}")]
    VcfRead { path: PathBuf, message: String },
}

pub type Result<T> = std::result::Result<T, Error>;
