use std::collections::HashSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sample {
    pub name: String,
    pub vcf_path: PathBuf,
    pub build: String,
    #[serde(default)]
    pub description: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Config {
    #[serde(default)]
    pub samples: Vec<Sample>,
}

impl Config {
    pub fn load(path: &Path) -> Result<Self> {
        let raw = std::fs::read_to_string(path).map_err(|source| Error::ConfigRead {
            path: path.to_path_buf(),
            source,
        })?;
        let mut cfg: Config = toml::from_str(&raw).map_err(|source| Error::ConfigParse {
            path: path.to_path_buf(),
            source,
        })?;
        // Resolve relative vcf_paths against the directory containing the config file.
        if let Some(parent) = path.parent() {
            for s in &mut cfg.samples {
                if s.vcf_path.is_relative() {
                    s.vcf_path = parent.join(&s.vcf_path);
                }
            }
        }
        Ok(cfg)
    }

    pub fn validate(&self) -> Result<()> {
        if self.samples.is_empty() {
            return Err(Error::NoSamples);
        }
        // Pass 1: name uniqueness and build validity (purely metadata).
        let mut seen: HashSet<&str> = HashSet::new();
        for s in &self.samples {
            if !seen.insert(&s.name) {
                return Err(Error::DuplicateSample(s.name.clone()));
            }
            match s.build.as_str() {
                "GRCh37" | "GRCh38" => {}
                other => {
                    return Err(Error::InvalidBuild {
                        sample: s.name.clone(),
                        build: other.to_string(),
                    });
                }
            }
        }
        // Pass 2: filesystem checks.
        for s in &self.samples {
            if !s.vcf_path.exists() {
                return Err(Error::VcfMissing {
                    sample: s.name.clone(),
                    path: s.vcf_path.clone(),
                });
            }
            let tbi = tabix_path_for(&s.vcf_path);
            if !tbi.exists() {
                return Err(Error::IndexMissing {
                    sample: s.name.clone(),
                    vcf: s.vcf_path.clone(),
                    expected: tbi,
                });
            }
        }
        Ok(())
    }
}

fn tabix_path_for(vcf: &Path) -> PathBuf {
    let mut tbi = vcf.as_os_str().to_owned();
    tbi.push(".tbi");
    PathBuf::from(tbi)
}

pub fn default_config_path() -> Option<PathBuf> {
    directories::ProjectDirs::from("", "", "vcf-mcp").map(|p| p.config_dir().join("config.toml"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_build() {
        let cfg = Config {
            samples: vec![Sample {
                name: "x".into(),
                vcf_path: PathBuf::from("/nonexistent.vcf.gz"),
                build: "hg19".into(),
                description: String::new(),
            }],
        };
        match cfg.validate() {
            Err(Error::InvalidBuild { sample, build }) => {
                assert_eq!(sample, "x");
                assert_eq!(build, "hg19");
            }
            other => panic!("expected InvalidBuild, got {other:?}"),
        }
    }

    #[test]
    fn rejects_empty_samples() {
        let cfg = Config { samples: vec![] };
        assert!(matches!(cfg.validate(), Err(Error::NoSamples)));
    }

    #[test]
    fn rejects_duplicate_names() {
        let cfg = Config {
            samples: vec![
                Sample {
                    name: "dup".into(),
                    vcf_path: PathBuf::from("/x.vcf.gz"),
                    build: "GRCh38".into(),
                    description: String::new(),
                },
                Sample {
                    name: "dup".into(),
                    vcf_path: PathBuf::from("/y.vcf.gz"),
                    build: "GRCh38".into(),
                    description: String::new(),
                },
            ],
        };
        assert!(matches!(cfg.validate(), Err(Error::DuplicateSample(_))));
    }

    #[test]
    fn tabix_path_appends_tbi() {
        assert_eq!(
            tabix_path_for(Path::new("/data/x.vcf.gz")),
            PathBuf::from("/data/x.vcf.gz.tbi")
        );
    }
}
