use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;
use std::sync::{Arc, RwLock};
use std::time::{Duration, Instant};

use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
use noodles_vcf::{
    self as vcf,
    variant::record::{
        AlternateBases as _, Filters as _, Ids as _,
        samples::{
            Sample as _,
            keys::key,
            series::{Value, value::genotype::Phasing},
        },
    },
};
use serde::Serialize;

use crate::config::{Config, Sample};
use crate::error::{Error, Result};

const MAX_REGION_BP: u64 = 10_000_000;
const MAX_RECORDS: usize = 500;
const QUERY_TIMEOUT_SECS: u64 = 30;
const MAX_RSIDS_PER_LOOKUP: usize = 100;

#[derive(Debug, Serialize)]
pub struct VariantRecord {
    pub chrom: String,
    pub pos: u32,
    #[serde(rename = "ref")]
    pub ref_bases: String,
    pub alt: String,
    pub rsid: Option<String>,
    pub genotype: Option<String>,
    pub genotype_alleles: Option<String>,
    pub depth: Option<i32>,
    pub gq: Option<i32>,
    pub filter: String,
}

#[derive(Debug, Serialize)]
pub struct QueryRegionResponse {
    pub sample: String,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub count: usize,
    pub truncated: bool,
    pub variants: Vec<VariantRecord>,
}

#[derive(Debug, Clone)]
pub struct QueryRegionArgs {
    pub sample: String,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

pub async fn query_region(cfg: Arc<Config>, args: QueryRegionArgs) -> Result<QueryRegionResponse> {
    if args.start > args.end {
        return Err(Error::InvalidRange {
            start: args.start,
            end: args.end,
        });
    }
    let length = u64::from(args.end) - u64::from(args.start) + 1;
    if length > MAX_REGION_BP {
        return Err(Error::RegionTooLarge {
            length,
            max: MAX_REGION_BP,
        });
    }

    let sample = cfg
        .samples
        .iter()
        .find(|s| s.name == args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?
        .clone();
    let vcf_path = sample.vcf_path.clone();

    let timeout = Duration::from_secs(QUERY_TIMEOUT_SECS);
    let work = tokio::task::spawn_blocking(move || query_region_blocking(&sample, args));

    match tokio::time::timeout(timeout, work).await {
        Err(_) => Err(Error::QueryTimeout {
            secs: QUERY_TIMEOUT_SECS,
        }),
        Ok(Err(join_err)) => Err(Error::VcfRead {
            path: vcf_path,
            message: format!("blocking task panicked: {join_err}"),
        }),
        Ok(Ok(r)) => r,
    }
}

fn query_region_blocking(sample: &Sample, args: QueryRegionArgs) -> Result<QueryRegionResponse> {
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(&sample.vcf_path)
        .map_err(|source| Error::VcfOpen {
            path: sample.vcf_path.clone(),
            source,
        })?;
    let header = reader.read_header().map_err(|e| Error::VcfRead {
        path: sample.vcf_path.clone(),
        message: e.to_string(),
    })?;

    let available_contigs: Vec<String> = header.contigs().keys().cloned().collect();
    let chrom = normalize_chrom(&args.chrom, &available_contigs).ok_or_else(|| {
        Error::InvalidChromosome {
            sample: args.sample.clone(),
            chrom: args.chrom.clone(),
            available: available_contigs.clone(),
        }
    })?;

    let start = Position::try_from(args.start as usize).map_err(|_| Error::InvalidRange {
        start: args.start,
        end: args.end,
    })?;
    let end = Position::try_from(args.end as usize).map_err(|_| Error::InvalidRange {
        start: args.start,
        end: args.end,
    })?;
    let region = Region::new(chrom.as_ref(), start..=end);

    let query = reader.query(&header, &region).map_err(|e| Error::VcfRead {
        path: sample.vcf_path.clone(),
        message: e.to_string(),
    })?;

    let mut variants: Vec<VariantRecord> = Vec::new();
    let mut truncated = false;

    for r in query.records() {
        if variants.len() >= MAX_RECORDS {
            truncated = true;
            break;
        }
        let rec = r.map_err(|e| Error::VcfRead {
            path: sample.vcf_path.clone(),
            message: e.to_string(),
        })?;
        let v = extract_variant(&header, &rec).map_err(|e| Error::VcfRead {
            path: sample.vcf_path.clone(),
            message: e,
        })?;
        variants.push(v);
    }

    Ok(QueryRegionResponse {
        sample: args.sample,
        chrom: chrom.into_owned(),
        start: args.start,
        end: args.end,
        count: variants.len(),
        truncated,
        variants,
    })
}

fn normalize_chrom<'a>(input: &'a str, available: &[String]) -> Option<Cow<'a, str>> {
    if available.iter().any(|c| c == input) {
        return Some(Cow::Borrowed(input));
    }
    if let Some(stripped) = input.strip_prefix("chr") {
        if available.iter().any(|c| c == stripped) {
            return Some(Cow::Owned(stripped.to_string()));
        }
    } else {
        let prefixed = format!("chr{input}");
        if available.iter().any(|c| c == &prefixed) {
            return Some(Cow::Owned(prefixed));
        }
    }
    None
}

fn extract_variant(
    header: &vcf::Header,
    rec: &vcf::Record,
) -> std::result::Result<VariantRecord, String> {
    let chrom = rec.reference_sequence_name().to_string();
    let pos = rec
        .variant_start()
        .transpose()
        .map_err(|e| e.to_string())?
        .map(usize::from)
        .ok_or_else(|| "record missing variant_start".to_string())? as u32;
    let ref_bases = rec.reference_bases().to_string();

    // gVCF reference-call records have no ALT (file shows "."), which noodles
    // surfaces as an empty iterator — emit "." rather than an ambiguous "".
    let alt = {
        let mut alts: Vec<String> = Vec::new();
        for r in rec.alternate_bases().iter() {
            alts.push(r.map_err(|e| e.to_string())?.to_string());
        }
        if alts.is_empty() {
            ".".to_string()
        } else {
            alts.join(",")
        }
    };

    let rsid = match rec.ids().iter().next() {
        None | Some(".") => None,
        Some(s) => Some(s.to_string()),
    };

    let filter = {
        let mut parts: Vec<String> = Vec::new();
        for r in rec.filters().iter(header) {
            parts.push(r.map_err(|e| e.to_string())?.to_string());
        }
        if parts.is_empty() {
            ".".to_string()
        } else {
            parts.join(";")
        }
    };

    let mut genotype: Option<String> = None;
    let mut depth: Option<i32> = None;
    let mut gq: Option<i32> = None;
    if let Some(s0) = rec.samples().iter().next() {
        // GT is parsed by noodles into a structured Genotype value, not a raw
        // string — rebuild the canonical "0/1" / "0|1" / "./." text from the
        // (allele, phasing) pairs.
        if let Some(Ok(Some(Value::Genotype(g)))) = s0.get(header, key::GENOTYPE) {
            genotype = Some(render_genotype(g.as_ref()).map_err(|e| e.to_string())?);
        }
        if let Some(Ok(Some(Value::Integer(v)))) = s0.get(header, key::READ_DEPTH) {
            depth = Some(v);
        }
        if let Some(Ok(Some(Value::Integer(v)))) = s0.get(header, key::CONDITIONAL_GENOTYPE_QUALITY)
        {
            gq = Some(v);
        }
    }

    let genotype_alleles = genotype
        .as_deref()
        .map(|gt| compute_genotype_alleles(gt, &ref_bases, &alt));

    Ok(VariantRecord {
        chrom,
        pos,
        ref_bases,
        alt,
        rsid,
        genotype,
        genotype_alleles,
        depth,
        gq,
        filter,
    })
}

fn render_genotype(
    g: &dyn noodles_vcf::variant::record::samples::series::value::genotype::Genotype,
) -> std::io::Result<String> {
    let mut out = String::new();
    for (i, result) in g.iter().enumerate() {
        let (allele, phasing) = result?;
        if i > 0 {
            out.push(match phasing {
                Phasing::Phased => '|',
                Phasing::Unphased => '/',
            });
        }
        match allele {
            Some(idx) => out.push_str(&idx.to_string()),
            None => out.push('.'),
        }
    }
    Ok(out)
}

// Translate a VCF GT string ("0/1", "1|2", "./.") into a concatenated allele
// string ("CT", "AG", "./."). Phasing separators ('/', '|') are stripped.
// Missing alleles ('.') or any parse failure cause us to return the raw GT
// unchanged — which preserves "./." and any unexpected forms losslessly.
fn compute_genotype_alleles(gt: &str, ref_bases: &str, alt_csv: &str) -> String {
    let alts: Vec<&str> = if alt_csv.is_empty() {
        Vec::new()
    } else {
        alt_csv.split(',').collect()
    };
    let alleles: Vec<&str> = std::iter::once(ref_bases)
        .chain(alts.iter().copied())
        .collect();

    let parts: Vec<&str> = gt.split(['/', '|']).collect();
    if parts.iter().any(|p| *p == "." || p.is_empty()) {
        return gt.to_string();
    }
    let mut out = String::new();
    for p in parts {
        let Ok(idx) = p.parse::<usize>() else {
            return gt.to_string();
        };
        let Some(a) = alleles.get(idx) else {
            return gt.to_string();
        };
        out.push_str(a);
    }
    out
}

// ---------- query_gene ----------

#[derive(Debug, Clone)]
pub struct QueryGeneArgs {
    pub sample: String,
    pub gene: String,
    pub flank_bp: Option<u32>,
}

#[derive(Debug, Serialize)]
pub struct QueryGeneResponse {
    pub gene: String,
    pub ensembl_id: String,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub build: String,
    pub flank_bp: u32,
    pub count: usize,
    pub truncated: bool,
    pub variants: Vec<VariantRecord>,
}

pub async fn query_gene(cfg: Arc<Config>, args: QueryGeneArgs) -> Result<QueryGeneResponse> {
    let sample = cfg
        .samples
        .iter()
        .find(|s| s.name == args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?
        .clone();

    let table = crate::genes::table_for(&sample.build).ok_or_else(|| Error::InvalidBuild {
        sample: sample.name.clone(),
        build: sample.build.clone(),
    })?;

    let gene = match table.lookup(&args.gene) {
        Some(g) => g.clone(),
        None => {
            let suggestions = table.suggestions(&args.gene, 3);
            let msg = if suggestions.is_empty() {
                format!("gene {:?} not found in {}", args.gene, sample.build)
            } else {
                format!(
                    "gene {:?} not found in {}. Did you mean: {}?",
                    args.gene,
                    sample.build,
                    suggestions.join(", ")
                )
            };
            return Err(Error::GeneNotFound(msg));
        }
    };

    let flank = args.flank_bp.unwrap_or(0);
    let start = gene.start.saturating_sub(flank).max(1);
    let end = gene.end.saturating_add(flank);

    let region_response = query_region(
        cfg.clone(),
        QueryRegionArgs {
            sample: args.sample.clone(),
            chrom: gene.chrom.clone(),
            start,
            end,
        },
    )
    .await?;

    Ok(QueryGeneResponse {
        gene: gene.symbol,
        ensembl_id: gene.ensembl_id,
        chrom: region_response.chrom,
        start: gene.start,
        end: gene.end,
        build: sample.build,
        flank_bp: flank,
        count: region_response.count,
        truncated: region_response.truncated,
        variants: region_response.variants,
    })
}

// ---------- lookup_rsids ----------

#[derive(Debug, Clone)]
pub struct RsidLocation {
    pub chrom: String,
    pub pos: u32,
}

pub type RsidIndex = HashMap<String, RsidLocation>;

#[derive(Default)]
pub struct RsidCache {
    inner: RwLock<HashMap<String, Arc<RsidIndex>>>,
}

impl RsidCache {
    pub fn get(&self, sample: &str) -> Option<Arc<RsidIndex>> {
        self.inner.read().ok()?.get(sample).cloned()
    }
    pub fn set(&self, sample: String, idx: Arc<RsidIndex>) {
        if let Ok(mut w) = self.inner.write() {
            w.insert(sample, idx);
        }
    }
    #[allow(dead_code)] // used by tests; useful for diagnostics
    pub fn has(&self, sample: &str) -> bool {
        self.inner
            .read()
            .map(|r| r.contains_key(sample))
            .unwrap_or(false)
    }
}

#[derive(Debug, Clone)]
pub struct LookupRsidsArgs {
    pub sample: String,
    pub rsids: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct LookupRsidEntry {
    pub rsid: String,
    pub found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub variant: Option<VariantRecord>,
}

#[derive(Debug, Serialize)]
pub struct LookupRsidsResponse {
    pub sample: String,
    pub count: usize,
    pub found_count: usize,
    pub results: Vec<LookupRsidEntry>,
}

pub async fn lookup_rsids(
    cfg: Arc<Config>,
    cache: Arc<RsidCache>,
    args: LookupRsidsArgs,
) -> Result<LookupRsidsResponse> {
    if args.rsids.is_empty() {
        return Err(Error::EmptyRsidList);
    }
    if args.rsids.len() > MAX_RSIDS_PER_LOOKUP {
        return Err(Error::TooManyRsids {
            count: args.rsids.len(),
            max: MAX_RSIDS_PER_LOOKUP,
        });
    }

    let sample = cfg
        .samples
        .iter()
        .find(|s| s.name == args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?
        .clone();

    let idx = ensure_rsid_index(&cache, &sample).await?;

    let vcf_path = sample.vcf_path.clone();
    let rsids = args.rsids.clone();
    let work =
        tokio::task::spawn_blocking(move || lookup_rsids_blocking(&sample.vcf_path, &idx, &rsids));
    let results = match tokio::time::timeout(Duration::from_secs(QUERY_TIMEOUT_SECS), work).await {
        Err(_) => {
            return Err(Error::QueryTimeout {
                secs: QUERY_TIMEOUT_SECS,
            });
        }
        Ok(Err(e)) => {
            return Err(Error::VcfRead {
                path: vcf_path,
                message: format!("blocking task panicked: {e}"),
            });
        }
        Ok(Ok(r)) => r?,
    };

    let found_count = results.iter().filter(|e| e.found).count();
    Ok(LookupRsidsResponse {
        sample: args.sample,
        count: results.len(),
        found_count,
        results,
    })
}

async fn ensure_rsid_index(cache: &RsidCache, sample: &Sample) -> Result<Arc<RsidIndex>> {
    if let Some(idx) = cache.get(&sample.name) {
        return Ok(idx);
    }
    let path = sample.vcf_path.clone();
    let started = Instant::now();
    let result = tokio::task::spawn_blocking(move || build_rsid_index(&path))
        .await
        .map_err(|e| Error::VcfRead {
            path: sample.vcf_path.clone(),
            message: format!("cache build panicked: {e}"),
        })?
        .map_err(|e| Error::VcfRead {
            path: sample.vcf_path.clone(),
            message: format!("cache build failed: {e}"),
        })?;
    let entries = result.len();
    let elapsed = started.elapsed();
    tracing::info!(
        sample = %sample.name,
        entries = entries,
        elapsed_ms = elapsed.as_millis() as u64,
        "built rsid cache"
    );
    let arc = Arc::new(result);
    cache.set(sample.name.clone(), arc.clone());
    Ok(arc)
}

// VCFs have no rsid index, so the only way to map rs→position is to stream the
// whole file once. Bgzf reads + line-level parsing of CHROM/POS/ID is fast: a
// 12M-record WGS file completes in a few seconds. We deliberately skip the
// full noodles record parser here — only the first three columns matter.
fn build_rsid_index(path: &Path) -> std::io::Result<RsidIndex> {
    let mut reader = bgzf::io::Reader::new(File::open(path)?);
    let mut idx = HashMap::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        if reader.read_line(&mut buf)? == 0 {
            break;
        }
        if buf.starts_with('#') {
            continue;
        }
        let mut it = buf.splitn(4, '\t');
        let chrom = match it.next() {
            Some(s) => s,
            None => continue,
        };
        let pos_str = match it.next() {
            Some(s) => s,
            None => continue,
        };
        let id_field = match it.next() {
            Some(s) => s,
            None => continue,
        };
        if id_field == "." || id_field.is_empty() {
            continue;
        }
        let pos: u32 = match pos_str.trim().parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        // The ID column may contain multiple ;-separated IDs.
        for id in id_field.split(';') {
            if id.starts_with("rs") {
                idx.insert(
                    id.to_string(),
                    RsidLocation {
                        chrom: chrom.to_string(),
                        pos,
                    },
                );
            }
        }
    }
    Ok(idx)
}

fn lookup_rsids_blocking(
    vcf_path: &Path,
    idx: &RsidIndex,
    rsids: &[String],
) -> Result<Vec<LookupRsidEntry>> {
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(vcf_path)
        .map_err(|source| Error::VcfOpen {
            path: vcf_path.to_path_buf(),
            source,
        })?;
    let header = reader.read_header().map_err(|e| Error::VcfRead {
        path: vcf_path.to_path_buf(),
        message: e.to_string(),
    })?;

    let mut out = Vec::with_capacity(rsids.len());
    for rsid in rsids {
        let Some(loc) = idx.get(rsid) else {
            out.push(LookupRsidEntry {
                rsid: rsid.clone(),
                found: false,
                variant: None,
            });
            continue;
        };

        let pos = Position::try_from(loc.pos as usize).map_err(|_| Error::VcfRead {
            path: vcf_path.to_path_buf(),
            message: format!("invalid cached position for {rsid}"),
        })?;
        let region = Region::new(loc.chrom.as_str(), pos..=pos);
        let query = reader.query(&header, &region).map_err(|e| Error::VcfRead {
            path: vcf_path.to_path_buf(),
            message: e.to_string(),
        })?;
        let mut matched: Option<VariantRecord> = None;
        for r in query.records() {
            let rec = r.map_err(|e| Error::VcfRead {
                path: vcf_path.to_path_buf(),
                message: e.to_string(),
            })?;
            let v = extract_variant(&header, &rec).map_err(|e| Error::VcfRead {
                path: vcf_path.to_path_buf(),
                message: e,
            })?;
            if v.rsid.as_deref() == Some(rsid.as_str()) {
                matched = Some(v);
                break;
            }
        }
        out.push(LookupRsidEntry {
            rsid: rsid.clone(),
            found: matched.is_some(),
            variant: matched,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_chrom_matches_chr_prefix() {
        let avail = vec!["chr17".to_string(), "chr1".to_string()];
        assert_eq!(normalize_chrom("17", &avail).as_deref(), Some("chr17"));
        assert_eq!(normalize_chrom("chr17", &avail).as_deref(), Some("chr17"));
        assert_eq!(normalize_chrom("chr99", &avail), None);
    }

    #[test]
    fn normalize_chrom_matches_numeric() {
        let avail = vec!["17".to_string(), "1".to_string()];
        assert_eq!(normalize_chrom("chr17", &avail).as_deref(), Some("17"));
        assert_eq!(normalize_chrom("17", &avail).as_deref(), Some("17"));
    }

    #[test]
    fn alleles_diploid_het() {
        assert_eq!(compute_genotype_alleles("0/1", "C", "T"), "CT");
        assert_eq!(compute_genotype_alleles("1|0", "C", "T"), "TC");
        assert_eq!(compute_genotype_alleles("1/1", "C", "T"), "TT");
    }

    #[test]
    fn alleles_multi_allelic() {
        assert_eq!(compute_genotype_alleles("1/2", "C", "T,G"), "TG");
    }

    #[test]
    fn alleles_missing_returns_raw() {
        assert_eq!(compute_genotype_alleles("./.", "C", "T"), "./.");
        assert_eq!(compute_genotype_alleles("0/.", "C", "T"), "0/.");
    }

    #[test]
    fn rsid_index_builds_against_slice() {
        // GIAB benchmark slice has all IDs = ".", so the resulting index is
        // expected to be empty. This is a smoke test that build_rsid_index
        // streams the entire file without erroring out.
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/na12878_chr17_slice.vcf.gz");
        let idx = build_rsid_index(&path).expect("build cache");
        assert_eq!(idx.len(), 0, "slice has no rsids in ID column");
    }

    #[test]
    fn rsid_cache_get_set() {
        let cache = RsidCache::default();
        assert!(!cache.has("x"));
        let mut idx: RsidIndex = HashMap::new();
        idx.insert(
            "rs1".to_string(),
            RsidLocation {
                chrom: "chr1".to_string(),
                pos: 100,
            },
        );
        cache.set("x".to_string(), Arc::new(idx));
        assert!(cache.has("x"));
        let got = cache.get("x").expect("present");
        assert_eq!(got.get("rs1").unwrap().pos, 100);
    }
}
