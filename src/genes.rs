//! Embedded gene coordinate tables.
//!
//! The TSVs in `data/genes_grch3{7,8}.tsv` are preprocessed from Ensembl
//! release 115 (GRCh38) and the GRCh37 archive release 87 GTFs by the
//! `generate_gene_tables` test (see `tests/integration.rs`). Only
//! protein-coding genes are included, and gene symbols mapping to more
//! than one Ensembl gene_id are dropped to avoid ambiguity.
//!
//! Source URLs (verified 2026-05-26):
//!   - https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
//!   - https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

use std::collections::HashMap;
use std::sync::OnceLock;

#[allow(dead_code)] // for documentation / future use in tool descriptions
pub const ENSEMBL_RELEASE_GRCH38: u32 = 115;
#[allow(dead_code)]
pub const ENSEMBL_RELEASE_GRCH37: u32 = 87;

const TSV_GRCH38: &str = include_str!("../data/genes_grch38.tsv");
const TSV_GRCH37: &str = include_str!("../data/genes_grch37.tsv");

#[derive(Debug, Clone)]
pub struct GeneRegion {
    pub symbol: String,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    #[allow(dead_code)] // present in the TSV; not currently read by callers
    pub strand: char,
    pub ensembl_id: String,
}

pub struct GeneTable {
    by_symbol_upper: HashMap<String, GeneRegion>,
}

impl GeneTable {
    fn parse(tsv: &str) -> Self {
        let mut by_symbol_upper = HashMap::with_capacity(20_000);
        for (i, line) in tsv.lines().enumerate() {
            if i == 0 && line.starts_with("gene_symbol") {
                continue; // header
            }
            if line.is_empty() {
                continue;
            }
            let mut it = line.split('\t');
            let Some(symbol) = it.next() else { continue };
            let Some(chrom) = it.next() else { continue };
            let Some(start) = it.next().and_then(|s| s.parse::<u32>().ok()) else {
                continue;
            };
            let Some(end) = it.next().and_then(|s| s.parse::<u32>().ok()) else {
                continue;
            };
            let strand = it.next().and_then(|s| s.chars().next()).unwrap_or('.');
            let ensembl_id = it.next().unwrap_or("").to_string();
            by_symbol_upper.insert(
                symbol.to_ascii_uppercase(),
                GeneRegion {
                    symbol: symbol.to_string(),
                    chrom: chrom.to_string(),
                    start,
                    end,
                    strand,
                    ensembl_id,
                },
            );
        }
        Self { by_symbol_upper }
    }

    pub fn lookup(&self, symbol: &str) -> Option<&GeneRegion> {
        self.by_symbol_upper.get(&symbol.to_ascii_uppercase())
    }

    /// Up to `limit` symbols within Levenshtein distance 1 of `symbol`,
    /// excluding `symbol` itself. Returned in deterministic alphabetical order.
    pub fn suggestions(&self, symbol: &str, limit: usize) -> Vec<String> {
        let needle = symbol.to_ascii_uppercase();
        let mut hits: Vec<&str> = self
            .by_symbol_upper
            .iter()
            .filter(|(k, _)| **k != needle && levenshtein_le_1(k, &needle))
            .map(|(_, v)| v.symbol.as_str())
            .collect();
        hits.sort_unstable();
        hits.truncate(limit);
        hits.into_iter().map(|s| s.to_string()).collect()
    }

    #[allow(dead_code)] // used by tests; useful for diagnostics
    pub fn len(&self) -> usize {
        self.by_symbol_upper.len()
    }
}

static GRCH38: OnceLock<GeneTable> = OnceLock::new();
static GRCH37: OnceLock<GeneTable> = OnceLock::new();

pub fn table_for(build: &str) -> Option<&'static GeneTable> {
    match build {
        "GRCh38" => Some(GRCH38.get_or_init(|| GeneTable::parse(TSV_GRCH38))),
        "GRCh37" => Some(GRCH37.get_or_init(|| GeneTable::parse(TSV_GRCH37))),
        _ => None,
    }
}

/// True iff `a` and `b` are within Levenshtein distance 1 (insertion,
/// deletion, or substitution of a single character). Distance 0 (equal
/// strings) also returns true — callers filter that case explicitly.
fn levenshtein_le_1(a: &str, b: &str) -> bool {
    let av: Vec<char> = a.chars().collect();
    let bv: Vec<char> = b.chars().collect();
    let la = av.len();
    let lb = bv.len();
    if la.abs_diff(lb) > 1 {
        return false;
    }
    if la == lb {
        return av.iter().zip(bv.iter()).filter(|(x, y)| x != y).count() <= 1;
    }
    let (shorter, longer) = if la < lb { (&av, &bv) } else { (&bv, &av) };
    let mut i = 0usize;
    let mut j = 0usize;
    let mut skipped = false;
    while i < shorter.len() && j < longer.len() {
        if shorter[i] == longer[j] {
            i += 1;
            j += 1;
        } else if !skipped {
            j += 1;
            skipped = true;
        } else {
            return false;
        }
    }
    // After the loop: shorter is fully consumed (i == shorter.len()).
    // Either we already used our skip, or longer has exactly one trailing char
    // we'll implicitly skip. Both states satisfy distance == 1.
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grch38_table_loads_and_has_known_genes() {
        let t = table_for("GRCh38").expect("GRCh38 table");
        assert!(
            t.len() > 15_000,
            "expected ~19k protein-coding genes, got {}",
            t.len()
        );

        let col1a1 = t.lookup("COL1A1").expect("COL1A1");
        assert_eq!(col1a1.chrom, "17");
        // Stable coords for COL1A1 in Ensembl 115 GRCh38
        assert!(col1a1.start >= 50_180_000 && col1a1.end <= 50_205_000);

        let brca1 = t.lookup("BRCA1").expect("BRCA1");
        assert_eq!(brca1.chrom, "17");

        let apoe = t.lookup("APOE").expect("APOE");
        assert_eq!(apoe.chrom, "19");
    }

    #[test]
    fn grch37_table_loads() {
        let t = table_for("GRCh37").expect("GRCh37 table");
        assert!(t.len() > 15_000);
        let col1a1 = t.lookup("COL1A1").expect("COL1A1");
        assert_eq!(col1a1.chrom, "17");
        // GRCh37 coords are different from GRCh38
        assert!(col1a1.start >= 48_000_000 && col1a1.end <= 48_500_000);
    }

    #[test]
    fn lookup_is_case_insensitive() {
        let t = table_for("GRCh38").unwrap();
        assert!(t.lookup("col1a1").is_some());
        assert!(t.lookup("Col1A1").is_some());
        assert!(t.lookup("COL1A1").is_some());
    }

    #[test]
    fn unknown_build_returns_none() {
        assert!(table_for("hg19").is_none());
        assert!(table_for("").is_none());
    }

    #[test]
    fn distance_one_basics() {
        // distance 0
        assert!(levenshtein_le_1("BRCA1", "BRCA1"));
        // substitution
        assert!(levenshtein_le_1("BRCA1", "BRCA2"));
        // insertion
        assert!(levenshtein_le_1("BRCA", "BRCA1"));
        // deletion
        assert!(levenshtein_le_1("BRCA1", "BRCA"));
        // distance 2
        assert!(!levenshtein_le_1("BRCA1", "BRCA22"));
        assert!(!levenshtein_le_1("COL1A1", "COL3A1Z"));
        // empty cases
        assert!(levenshtein_le_1("", ""));
        assert!(levenshtein_le_1("", "A"));
        assert!(!levenshtein_le_1("", "AB"));
    }

    #[test]
    fn suggestions_for_typo() {
        let t = table_for("GRCh38").unwrap();
        // BRCA1 typo: BRCA1 → BRCA2 is distance 1
        let s = t.suggestions("BRCA2X", 5);
        assert!(s.contains(&"BRCA2".to_string()), "got {s:?}");
    }

    #[test]
    fn suggestions_exclude_exact_match() {
        let t = table_for("GRCh38").unwrap();
        let s = t.suggestions("BRCA1", 5);
        assert!(!s.contains(&"BRCA1".to_string()));
    }
}
