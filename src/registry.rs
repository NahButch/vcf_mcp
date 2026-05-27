//! Mutable in-memory sample registry with optional disk-backed persistence.
//!
//! Replaces the old static `Config.samples` list. The registry is the source
//! of truth at runtime; the state file is a write-through cache so that
//! `add_sample` / `remove_sample` survive server restarts.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, RwLock};

use noodles_tabix as tabix;
use serde::{Deserialize, Serialize};

use crate::config::Sample;
use crate::error::{Error, Result};

#[derive(Serialize, Deserialize, Default)]
struct StateFile {
    #[serde(default)]
    samples: Vec<Sample>,
}

pub struct SampleRegistry {
    inner: Mutex<Inner>,
    state_path: Option<PathBuf>,
    // Ephemeral per-sample tabix-index cache. Not persisted to state.json;
    // rebuilt on demand. Lives on the registry so queries don't have to
    // thread a separate cache through every call site.
    tabix_cache: RwLock<HashMap<String, Arc<tabix::Index>>>,
}

struct Inner {
    by_name: HashMap<String, Sample>,
    // canonical absolute path → name, for idempotent re-add by path.
    by_path: HashMap<PathBuf, String>,
}

impl SampleRegistry {
    pub fn new(state_path: Option<PathBuf>) -> Self {
        Self {
            inner: Mutex::new(Inner {
                by_name: HashMap::new(),
                by_path: HashMap::new(),
            }),
            state_path,
            tabix_cache: RwLock::new(HashMap::new()),
        }
    }

    pub fn get_tabix_index(&self, sample: &str) -> Option<Arc<tabix::Index>> {
        self.tabix_cache.read().ok()?.get(sample).cloned()
    }

    pub fn set_tabix_index(&self, sample: String, index: Arc<tabix::Index>) {
        if let Ok(mut w) = self.tabix_cache.write() {
            w.insert(sample, index);
        }
    }

    fn drop_tabix_index(&self, sample: &str) {
        if let Ok(mut w) = self.tabix_cache.write() {
            w.remove(sample);
        }
    }

    /// Load samples from the configured state file, if any. Missing file is OK
    /// (fresh install).
    pub fn load_from_disk(&self) -> Result<()> {
        let Some(path) = &self.state_path else {
            return Ok(());
        };
        if !path.exists() {
            return Ok(());
        }
        let raw = std::fs::read_to_string(path).map_err(|e| Error::StateFile {
            path: path.clone(),
            message: e.to_string(),
        })?;
        let state: StateFile = serde_json::from_str(&raw).map_err(|e| Error::StateFile {
            path: path.clone(),
            message: format!("parse: {e}"),
        })?;
        let mut inner = self.inner.lock().unwrap();
        inner.by_name.clear();
        inner.by_path.clear();
        for s in state.samples {
            inner.by_path.insert(s.vcf_path.clone(), s.name.clone());
            inner.by_name.insert(s.name.clone(), s);
        }
        Ok(())
    }

    /// Import samples (typically from a TOML config). Entries with a name
    /// collision get a random suffix appended. Same path under different
    /// names is allowed — TOML can intentionally register a file under
    /// multiple views.
    pub fn import(&self, samples: Vec<Sample>) -> Result<usize> {
        let mut added = 0;
        for mut s in samples {
            let canonical =
                std::fs::canonicalize(&s.vcf_path).unwrap_or_else(|_| s.vcf_path.clone());
            s.vcf_path = canonical.clone();

            let mut inner = self.inner.lock().unwrap();
            let final_name = resolve_collision(&inner, &s.name);
            s.name = final_name;
            // First-name-wins for the by_path lookup; we don't overwrite if
            // another name is already registered against this canonical path.
            inner
                .by_path
                .entry(canonical)
                .or_insert_with(|| s.name.clone());
            inner.by_name.insert(s.name.clone(), s);
            added += 1;
        }
        self.persist()?;
        Ok(added)
    }

    pub fn get(&self, name: &str) -> Option<Sample> {
        let inner = self.inner.lock().unwrap();
        inner.by_name.get(name).cloned()
    }

    pub fn get_by_path(&self, canonical_path: &Path) -> Option<Sample> {
        let inner = self.inner.lock().unwrap();
        let name = inner.by_path.get(canonical_path)?;
        inner.by_name.get(name).cloned()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.inner.lock().unwrap().by_name.contains_key(name)
    }

    pub fn list(&self) -> Vec<Sample> {
        let inner = self.inner.lock().unwrap();
        let mut v: Vec<Sample> = inner.by_name.values().cloned().collect();
        v.sort_by(|a, b| a.name.cmp(&b.name));
        v
    }

    pub fn len(&self) -> usize {
        self.inner.lock().unwrap().by_name.len()
    }

    /// Register a fully-validated sample. The caller is responsible for the
    /// file validation chain. Returns the registered sample (whose name may
    /// carry a random suffix if there was a collision). Path-idempotency is
    /// NOT enforced here — same path under different names is allowed. The
    /// `add_sample` tool path applies path-idempotency at a higher level
    /// only when no explicit name was provided.
    pub fn add_validated(&self, mut sample: Sample) -> Result<Sample> {
        let canonical = sample.vcf_path.clone();
        let mut inner = self.inner.lock().unwrap();
        let final_name = resolve_collision(&inner, &sample.name);
        sample.name = final_name;
        inner
            .by_path
            .entry(canonical)
            .or_insert_with(|| sample.name.clone());
        inner.by_name.insert(sample.name.clone(), sample.clone());
        drop(inner);
        self.persist()?;
        Ok(sample)
    }

    pub fn remove(&self, name: &str) -> Result<Option<Sample>> {
        let mut inner = self.inner.lock().unwrap();
        let Some(sample) = inner.by_name.remove(name) else {
            return Ok(None);
        };
        inner.by_path.remove(&sample.vcf_path);
        drop(inner);
        self.drop_tabix_index(name);
        self.persist()?;
        Ok(Some(sample))
    }

    fn persist(&self) -> Result<()> {
        let Some(path) = &self.state_path else {
            return Ok(());
        };
        let samples = self.list();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| Error::StateFile {
                path: path.clone(),
                message: format!("mkdir parent: {e}"),
            })?;
        }
        let payload = StateFile { samples };
        let json = serde_json::to_string_pretty(&payload).map_err(|e| Error::StateFile {
            path: path.clone(),
            message: format!("serialize: {e}"),
        })?;
        // Atomic write: tmp + rename.
        let tmp_path = path.with_extension("tmp");
        std::fs::write(&tmp_path, &json).map_err(|e| Error::StateFile {
            path: tmp_path.clone(),
            message: format!("write tmp: {e}"),
        })?;
        std::fs::rename(&tmp_path, path).map_err(|e| Error::StateFile {
            path: path.clone(),
            message: format!("rename: {e}"),
        })?;
        Ok(())
    }
}

fn resolve_collision(inner: &Inner, base: &str) -> String {
    if !inner.by_name.contains_key(base) {
        return base.to_string();
    }
    // Append _<4 random alnum> until we find a free slot. ~1.6M values; one
    // collision-resolution attempt should always be enough.
    for _ in 0..16 {
        let candidate = format!("{base}_{}", random_suffix(4));
        if !inner.by_name.contains_key(&candidate) {
            return candidate;
        }
    }
    // Vanishingly unlikely; fall back to time-based suffix.
    format!("{base}_{}", random_suffix(8))
}

fn random_suffix(n: usize) -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let mut x = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(0xC0FFEE)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let alphabet = b"abcdefghijklmnopqrstuvwxyz0123456789";
    let mut out = String::with_capacity(n);
    for _ in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        out.push(alphabet[((x >> 33) as usize) % alphabet.len()] as char);
    }
    out
}

/// Default location for the state file.
pub fn default_state_path() -> Option<PathBuf> {
    directories::ProjectDirs::from("", "", "vcf-mcp").map(|p| p.data_dir().join("state.json"))
}

/// Strip `.vcf.gz` (case-insensitive) and any preceding `.snp-indel`-style
/// sub-extensions from a filename, sanitize the result to a name suitable for
/// a sample id.
pub fn derive_name_from_path(path: &Path) -> String {
    let filename = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("sample");
    // Strip .vcf.gz case-insensitively.
    let stem = match filename
        .rfind('.')
        .filter(|&i| filename[i..].eq_ignore_ascii_case(".gz"))
        .and_then(|gz_idx| {
            let before_gz = &filename[..gz_idx];
            before_gz
                .rfind('.')
                .filter(|&j| before_gz[j..].eq_ignore_ascii_case(".vcf"))
                .map(|vcf_idx| &filename[..vcf_idx])
        }) {
        Some(s) => s,
        None => filename,
    };
    let mut out: String = stem
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if out.len() > 32 {
        out.truncate(32);
    }
    // Strip trailing underscores left over from punctuation collapsing.
    while out.ends_with('_') {
        out.pop();
    }
    if out.is_empty() {
        out = "sample".to_string();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derive_name_strips_extensions() {
        assert_eq!(derive_name_from_path(Path::new("foo.vcf.gz")), "foo");
        assert_eq!(derive_name_from_path(Path::new("foo.VCF.GZ")), "foo");
        assert_eq!(
            derive_name_from_path(Path::new("/x/y/patient_001.snp-indel.vcf.gz")),
            "patient_001_snp-indel"
        );
        assert_eq!(
            derive_name_from_path(Path::new("Jane Doe (WGS).vcf.gz")),
            "Jane_Doe__WGS"
        );
    }

    #[test]
    fn derive_name_truncates_to_32() {
        let long = "a".repeat(80) + ".vcf.gz";
        let n = derive_name_from_path(Path::new(&long));
        assert_eq!(n.len(), 32);
        assert!(n.chars().all(|c| c == 'a'));
    }

    #[test]
    fn collision_appends_suffix() {
        let reg = SampleRegistry::new(None);
        let s1 = Sample {
            name: "patient".into(),
            vcf_path: PathBuf::from("/tmp/a.vcf.gz"),
            build: "GRCh38".into(),
            description: String::new(),
        };
        let s2 = Sample {
            name: "patient".into(),
            vcf_path: PathBuf::from("/tmp/b.vcf.gz"),
            build: "GRCh38".into(),
            description: String::new(),
        };
        let r1 = reg.add_validated(s1).unwrap();
        let r2 = reg.add_validated(s2).unwrap();
        assert_eq!(r1.name, "patient");
        assert_ne!(r2.name, "patient");
        assert!(r2.name.starts_with("patient_"));
        assert_eq!(reg.len(), 2);
    }

    #[test]
    fn same_path_under_different_names_is_two_entries() {
        // Allowing the same canonical path under multiple names is
        // intentional — useful for comparing the same file under different
        // metadata views, or for TOML bootstrap clarity.
        let reg = SampleRegistry::new(None);
        let s = Sample {
            name: "patient".into(),
            vcf_path: PathBuf::from("/tmp/a.vcf.gz"),
            build: "GRCh38".into(),
            description: "first".into(),
        };
        reg.add_validated(s.clone()).unwrap();
        reg.add_validated(Sample {
            name: "patient_alias".into(),
            ..s
        })
        .unwrap();
        assert_eq!(reg.len(), 2);
    }

    #[test]
    fn get_by_path_returns_first_registered_name() {
        let reg = SampleRegistry::new(None);
        let p = PathBuf::from("/tmp/a.vcf.gz");
        reg.add_validated(Sample {
            name: "first".into(),
            vcf_path: p.clone(),
            build: "GRCh38".into(),
            description: String::new(),
        })
        .unwrap();
        reg.add_validated(Sample {
            name: "second".into(),
            vcf_path: p.clone(),
            build: "GRCh38".into(),
            description: String::new(),
        })
        .unwrap();
        // by_path keeps the first-registered name, so the LLM can reuse the
        // canonical name when re-adding by path without an explicit name.
        assert_eq!(reg.get_by_path(&p).unwrap().name, "first");
    }

    #[test]
    fn remove_works() {
        let reg = SampleRegistry::new(None);
        reg.add_validated(Sample {
            name: "x".into(),
            vcf_path: PathBuf::from("/tmp/x.vcf.gz"),
            build: "GRCh38".into(),
            description: String::new(),
        })
        .unwrap();
        assert!(reg.contains("x"));
        let removed = reg.remove("x").unwrap();
        assert!(removed.is_some());
        assert!(!reg.contains("x"));
        // Re-add should now succeed
        reg.add_validated(Sample {
            name: "x".into(),
            vcf_path: PathBuf::from("/tmp/x.vcf.gz"),
            build: "GRCh38".into(),
            description: String::new(),
        })
        .unwrap();
        assert!(reg.contains("x"));
    }
}
