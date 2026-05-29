//! Embedded "Fun Genetics Starter Panel" dataset.
//!
//! A static, curated table of well-studied trait-associated markers (mostly
//! dbSNP rsIDs, a few gene regions) used to offer users a friendly starter
//! analysis of their registered samples. The server only supplies the data
//! points (marker + genotype + the curator's direction note); all scoring,
//! orientation resolution, interpretation, and rendering are left to the
//! assistant. This is curiosity tooling, not medical advice.
//!
//! Source: `data/fun_panel.tsv` (category, trait, gene, marker, effect_allele,
//! direction). rsIDs are extracted from the `marker` field at parse time; a row
//! with no rsID (e.g. opsins, GNPTAB/NAGPA) is a gene-region row the assistant
//! can follow up on with `query_gene`.

use std::sync::OnceLock;

const PANEL_TSV: &str = include_str!("../data/fun_panel.tsv");

#[derive(Debug, Clone)]
pub struct PanelEntry {
    pub category: String,
    pub trait_name: String,
    pub gene: String,
    pub marker: String,
    pub effect_allele: String,
    pub direction: String,
    /// rsIDs parsed out of `marker` (lowercased). Empty for gene-region rows.
    pub rsids: Vec<String>,
}

impl PanelEntry {
    pub fn is_gene_region(&self) -> bool {
        self.rsids.is_empty()
    }
}

fn parse(tsv: &str) -> Vec<PanelEntry> {
    let mut out = Vec::new();
    for (i, line) in tsv.lines().enumerate() {
        if i == 0 && line.starts_with("category\t") {
            continue; // header
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        let Some(category) = it.next() else { continue };
        let Some(trait_name) = it.next() else {
            continue;
        };
        let Some(gene) = it.next() else { continue };
        let Some(marker) = it.next() else { continue };
        let effect_allele = it.next().unwrap_or("");
        let direction = it.next().unwrap_or("");
        out.push(PanelEntry {
            category: category.to_string(),
            trait_name: trait_name.to_string(),
            gene: gene.to_string(),
            marker: marker.to_string(),
            effect_allele: effect_allele.to_string(),
            direction: direction.to_string(),
            rsids: extract_rsids(marker),
        });
    }
    out
}

/// Pull every `rs<digits>` token out of a marker string (case-insensitive).
/// A single row may reference several (e.g. the CYP2C19 *2/*3/*17 row).
fn extract_rsids(marker: &str) -> Vec<String> {
    let chars: Vec<char> = marker.chars().collect();
    let n = chars.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        let is_rs = (chars[i] == 'r' || chars[i] == 'R')
            && i + 2 < n
            && (chars[i + 1] == 's' || chars[i + 1] == 'S')
            && chars[i + 2].is_ascii_digit();
        if is_rs {
            let mut s = String::from("rs");
            let mut j = i + 2;
            while j < n && chars[j].is_ascii_digit() {
                s.push(chars[j]);
                j += 1;
            }
            out.push(s);
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

/// All panel entries, parsed once.
pub fn entries() -> &'static [PanelEntry] {
    static ENTRIES: OnceLock<Vec<PanelEntry>> = OnceLock::new();
    ENTRIES.get_or_init(|| parse(PANEL_TSV))
}

/// Unique rsIDs across the whole panel, in first-seen order. These are the IDs
/// to look up per sample (one batched pass).
pub fn unique_rsids() -> &'static [String] {
    static RSIDS: OnceLock<Vec<String>> = OnceLock::new();
    RSIDS.get_or_init(|| {
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for e in entries() {
            for r in &e.rsids {
                if seen.insert(r.clone()) {
                    out.push(r.clone());
                }
            }
        }
        out
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_loads_and_has_rows() {
        let e = entries();
        assert!(e.len() > 80, "expected the full panel, got {}", e.len());
    }

    #[test]
    fn rsids_are_extracted_and_deduped() {
        let rsids = unique_rsids();
        // Well under the 100-per-lookup cap, so one batched pass suffices.
        assert!(
            rsids.len() <= 100,
            "rsid count {} exceeds lookup cap",
            rsids.len()
        );
        assert!(rsids.contains(&"rs713598".to_string()));
        // rs1800497 and the CYP2C9 SNPs appear in multiple rows but must dedupe.
        let count = rsids.iter().filter(|r| *r == "rs1800497").count();
        assert_eq!(count, 1, "rs1800497 should be deduped");
    }

    #[test]
    fn multi_rsid_row_parses_all() {
        // The CYP2C19 row carries three rsIDs in one marker field.
        let cyp = entries()
            .iter()
            .find(|e| e.gene == "CYP2C19")
            .expect("CYP2C19 row");
        assert_eq!(cyp.rsids.len(), 3, "got {:?}", cyp.rsids);
        assert!(cyp.rsids.contains(&"rs4244285".to_string()));
    }

    #[test]
    fn gene_region_rows_have_no_rsid() {
        let opsin = entries()
            .iter()
            .find(|e| e.gene.contains("OPN1LW"))
            .expect("opsin row");
        assert!(opsin.is_gene_region());
    }
}
