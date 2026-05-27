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

    #[error(
        "reset_samples is destructive and requires `confirm: true`; refusing to drop the registry without explicit consent"
    )]
    ResetNotConfirmed,
}

pub type Result<T> = std::result::Result<T, Error>;

/// Triage hint for the calling LLM. Determines whether Claude should help the
/// user fix their input, help them diagnose a data / environment problem, or
/// offer to file a bug report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    /// User-side input mistake — bad args, exceeded limits, unknown sample,
    /// unknown gene symbol. Don't offer bug reports for these.
    UserInput,
    /// User-side data or environment problem — corrupt VCF, missing file,
    /// unreadable index, unparseable TOML. Tell the user what to check;
    /// don't offer a bug report.
    UserData,
    /// Something internal vcf-mcp didn't expect — panic, timeout, unhandled
    /// edge case in validation, IO error in a path that should have worked.
    /// Worth offering to file a bug.
    Unexpected,
}

impl Category {
    pub fn as_str(&self) -> &'static str {
        match self {
            Category::UserInput => "user_input",
            Category::UserData => "user_data",
            Category::Unexpected => "unexpected",
        }
    }
}

impl Error {
    /// Tag this error for the LLM's triage heuristic. See the module-level
    /// docs / README for what each category means; new variants must opt in
    /// explicitly here so the categorization stays auditable.
    pub fn category(&self) -> Category {
        match self {
            // User typed something wrong / asked for too much / referenced
            // something that isn't registered. No bug to report.
            Error::SampleNotFound(_)
            | Error::InvalidChromosome { .. }
            | Error::InvalidRange { .. }
            | Error::RegionTooLarge { .. }
            | Error::EmptyRsidList
            | Error::TooManyRsids { .. }
            | Error::GeneNotFound(_)
            | Error::TooFewSamples { .. }
            | Error::TooManyFiles { .. }
            | Error::DuplicateSample(_)
            | Error::NoSamples
            | Error::InvalidBuild { .. }
            | Error::PathNotAllowed { .. }
            | Error::ResetNotConfirmed => Category::UserInput,

            // The user's filesystem / file content is the source of the
            // problem — corrupt download, missing file, malformed config.
            // Claude should help them diagnose.
            Error::ConfigRead { .. }
            | Error::ConfigParse { .. }
            | Error::VcfMissing { .. }
            | Error::IndexMissing { .. }
            | Error::PathInvalid { .. }
            | Error::InvalidVcfFile { .. }
            | Error::BuildNotDetectable { .. }
            | Error::VcfOpen { .. } => Category::UserData,

            // Smells like a bug in vcf-mcp itself, or an environment problem
            // we don't have a clean diagnostic for. Worth a bug report.
            Error::QueryTimeout { .. } | Error::VcfRead { .. } | Error::StateFile { .. } => {
                Category::Unexpected
            }
        }
    }

    /// Stable identifier for the error variant — useful for log filtering
    /// and (eventually) for stable bug-report titles.
    pub fn kind(&self) -> &'static str {
        match self {
            Error::ConfigRead { .. } => "ConfigRead",
            Error::ConfigParse { .. } => "ConfigParse",
            Error::InvalidBuild { .. } => "InvalidBuild",
            Error::VcfMissing { .. } => "VcfMissing",
            Error::IndexMissing { .. } => "IndexMissing",
            Error::DuplicateSample(_) => "DuplicateSample",
            Error::NoSamples => "NoSamples",
            Error::SampleNotFound(_) => "SampleNotFound",
            Error::RegionTooLarge { .. } => "RegionTooLarge",
            Error::InvalidRange { .. } => "InvalidRange",
            Error::InvalidChromosome { .. } => "InvalidChromosome",
            Error::QueryTimeout { .. } => "QueryTimeout",
            Error::TooManyRsids { .. } => "TooManyRsids",
            Error::EmptyRsidList => "EmptyRsidList",
            Error::GeneNotFound(_) => "GeneNotFound",
            Error::TooFewSamples { .. } => "TooFewSamples",
            Error::PathInvalid { .. } => "PathInvalid",
            Error::PathNotAllowed { .. } => "PathNotAllowed",
            Error::InvalidVcfFile { .. } => "InvalidVcfFile",
            Error::BuildNotDetectable { .. } => "BuildNotDetectable",
            Error::TooManyFiles { .. } => "TooManyFiles",
            Error::StateFile { .. } => "StateFile",
            Error::VcfOpen { .. } => "VcfOpen",
            Error::VcfRead { .. } => "VcfRead",
            Error::ResetNotConfirmed => "ResetNotConfirmed",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn input_errors_are_user_input() {
        assert_eq!(
            Error::SampleNotFound("x".into()).category(),
            Category::UserInput
        );
        assert_eq!(
            Error::InvalidRange { start: 5, end: 3 }.category(),
            Category::UserInput
        );
        assert_eq!(
            Error::RegionTooLarge {
                length: 20_000_000,
                max: 10_000_000
            }
            .category(),
            Category::UserInput
        );
        assert_eq!(Error::EmptyRsidList.category(), Category::UserInput);
        assert_eq!(
            Error::GeneNotFound("nope".into()).category(),
            Category::UserInput
        );
    }

    #[test]
    fn data_errors_are_user_data() {
        assert_eq!(
            Error::InvalidVcfFile {
                path: PathBuf::from("/x"),
                reason: "magic".into()
            }
            .category(),
            Category::UserData
        );
        assert_eq!(
            Error::BuildNotDetectable {
                path: PathBuf::from("/x")
            }
            .category(),
            Category::UserData
        );
        assert_eq!(
            Error::VcfMissing {
                sample: "x".into(),
                path: PathBuf::from("/missing")
            }
            .category(),
            Category::UserData
        );
    }

    #[test]
    fn internal_errors_are_unexpected() {
        assert_eq!(
            Error::QueryTimeout { secs: 30 }.category(),
            Category::Unexpected
        );
        assert_eq!(
            Error::StateFile {
                path: PathBuf::from("/x"),
                message: "y".into()
            }
            .category(),
            Category::Unexpected
        );
    }

    #[test]
    fn kind_matches_variant_name() {
        assert_eq!(Error::EmptyRsidList.kind(), "EmptyRsidList");
        assert_eq!(Error::QueryTimeout { secs: 30 }.kind(), "QueryTimeout");
    }
}
