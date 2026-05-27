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

use crate::config::Sample;
use crate::error::{Error, Result};
use crate::registry::{SampleRegistry, derive_name_from_path};

const MAX_REGION_BP: u64 = 10_000_000;
const MAX_RECORDS: usize = 500;
const QUERY_TIMEOUT_SECS: u64 = 30;
const MAX_RSIDS_PER_LOOKUP: usize = 100;

#[derive(Debug, Clone, Serialize)]
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

pub async fn query_region(
    registry: Arc<SampleRegistry>,
    args: QueryRegionArgs,
) -> Result<QueryRegionResponse> {
    if args.start < 1 {
        return Err(Error::InvalidStart { start: args.start });
    }
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

    let sample = registry
        .get(&args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?;
    let vcf_path = sample.vcf_path.clone();

    let timeout = Duration::from_secs(QUERY_TIMEOUT_SECS);
    let registry_for_blocking = registry.clone();
    let work = tokio::task::spawn_blocking(move || {
        query_region_blocking(&sample, &registry_for_blocking, args)
    });

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

fn query_region_blocking(
    sample: &Sample,
    registry: &SampleRegistry,
    args: QueryRegionArgs,
) -> Result<QueryRegionResponse> {
    let started = Instant::now();
    let mut reader = open_indexed_reader_with_cache(registry, sample)?;
    let header = reader.read_header().map_err(|e| Error::VcfRead {
        path: sample.vcf_path.clone(),
        message: e.to_string(),
    })?;
    tracing::info!(
        perf = true,
        phase = "query_region:setup",
        sample = %sample.name,
        elapsed_ms = started.elapsed().as_millis() as u64,
        "indexed reader open + header read"
    );

    let available_contigs: Vec<String> = header.contigs().keys().cloned().collect();
    let chrom = normalize_chrom(&args.chrom, &available_contigs).ok_or_else(|| {
        // Truncate the available list — GRCh38+alt+decoy builds have ~2500
        // contigs and dumping them all bloats the error to hundreds of KB.
        // The standard chromosomes come first in a well-formed header, so the
        // first 25 are the useful ones.
        const SHOW: usize = 25;
        let total = available_contigs.len();
        Error::InvalidChromosome {
            sample: args.sample.clone(),
            chrom: args.chrom.clone(),
            available: available_contigs.iter().take(SHOW).cloned().collect(),
            shown: total.min(SHOW),
            total,
        }
    })?;

    // Bounds check against the contig's declared length (Issue 3): a start
    // past the end of the chromosome otherwise surfaces as a noodles
    // "invalid start bound" read error miscategorized as `unexpected`.
    if let Some(contig) = header.contigs().get(chrom.as_ref()) {
        if let Some(len) = contig.length() {
            if args.start as usize > len {
                return Err(Error::PositionOutOfBounds {
                    chrom: chrom.into_owned(),
                    pos: args.start,
                    length: len as u32,
                });
            }
        }
    }

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
    tracing::info!(
        perf = true,
        phase = "query_region:query_issued",
        sample = %sample.name,
        chrom = %chrom,
        start = args.start,
        end = args.end,
        elapsed_ms = started.elapsed().as_millis() as u64,
        "tabix query issued; iterating records"
    );

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

    tracing::info!(
        perf = true,
        phase = "query_region:done",
        sample = %sample.name,
        records = variants.len(),
        truncated = truncated,
        elapsed_ms = started.elapsed().as_millis() as u64,
        "query_region complete"
    );

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

pub async fn query_gene(
    registry: Arc<SampleRegistry>,
    args: QueryGeneArgs,
) -> Result<QueryGeneResponse> {
    let sample = registry
        .get(&args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?;

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
        registry.clone(),
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

// ---------- compare_samples ----------

#[derive(Debug, Clone)]
pub enum CompareQuery {
    Rsids(Vec<String>),
    Region { chrom: String, start: u32, end: u32 },
}

#[derive(Debug, Clone)]
pub struct CompareSamplesArgs {
    pub samples: Vec<String>,
    pub query: CompareQuery,
}

#[derive(Debug, Serialize)]
pub struct CompareSampleResult {
    pub sample: String,
    pub found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype_alleles: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gq: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct CompareVariantEntry {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rsid: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chrom: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pos: Option<u32>,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    pub ref_bases: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alt: Option<String>,
    pub results: Vec<CompareSampleResult>,
}

#[derive(Debug, Serialize)]
pub struct CompareSamplesResponse {
    pub samples: Vec<String>,
    pub query_type: &'static str,
    pub count: usize,
    pub truncated: bool,
    pub results: Vec<CompareVariantEntry>,
}

pub async fn compare_samples(
    registry: Arc<SampleRegistry>,
    cache: Arc<RsidCache>,
    args: CompareSamplesArgs,
) -> Result<CompareSamplesResponse> {
    if args.samples.len() < 2 {
        return Err(Error::TooFewSamples {
            count: args.samples.len(),
        });
    }
    // Up-front sample existence check so we fail before doing any IO.
    for name in &args.samples {
        if !registry.contains(name) {
            return Err(Error::SampleNotFound(name.clone()));
        }
    }
    match args.query.clone() {
        CompareQuery::Region { chrom, start, end } => {
            compare_region(registry, &args.samples, &chrom, start, end).await
        }
        CompareQuery::Rsids(rsids) => {
            compare_by_rsids(registry, cache, &args.samples, &rsids).await
        }
    }
}

async fn compare_region(
    registry: Arc<SampleRegistry>,
    samples: &[String],
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<CompareSamplesResponse> {
    use std::collections::HashMap;

    let mut per_sample: HashMap<String, Vec<VariantRecord>> = HashMap::new();
    let mut any_truncated = false;
    for name in samples {
        let resp = query_region(
            registry.clone(),
            QueryRegionArgs {
                sample: name.clone(),
                chrom: chrom.to_string(),
                start,
                end,
            },
        )
        .await?;
        if resp.truncated {
            any_truncated = true;
        }
        per_sample.insert(name.clone(), resp.variants);
    }

    // Key by (chrom, pos, ref, alt). Iterate samples in input order so the
    // first-time-seen variants get inserted deterministically (only matters
    // for which rsid wins when two samples disagree on rsid annotation —
    // first sample wins).
    type Key = (String, u32, String, String);
    let mut by_key: HashMap<Key, CompareVariantEntry> = HashMap::new();
    let mut key_order: Vec<Key> = Vec::new();
    for name in samples {
        for v in per_sample.get(name).into_iter().flatten() {
            let key = (v.chrom.clone(), v.pos, v.ref_bases.clone(), v.alt.clone());
            if let std::collections::hash_map::Entry::Vacant(slot) = by_key.entry(key.clone()) {
                key_order.push(key);
                slot.insert(CompareVariantEntry {
                    rsid: v.rsid.clone(),
                    chrom: Some(v.chrom.clone()),
                    pos: Some(v.pos),
                    ref_bases: Some(v.ref_bases.clone()),
                    alt: Some(v.alt.clone()),
                    results: Vec::with_capacity(samples.len()),
                });
            }
        }
    }

    // For each entry, build per-sample results in input order.
    for entry in by_key.values_mut() {
        for name in samples {
            let matching = per_sample.get(name).and_then(|vs| {
                vs.iter().find(|v| {
                    Some(&v.chrom) == entry.chrom.as_ref()
                        && Some(v.pos) == entry.pos
                        && Some(&v.ref_bases) == entry.ref_bases.as_ref()
                        && Some(&v.alt) == entry.alt.as_ref()
                })
            });
            entry
                .results
                .push(make_sample_result(name, matching.cloned()));
        }
    }

    // Emit in (chrom, pos) sort order for stable output.
    let mut results: Vec<CompareVariantEntry> = key_order
        .into_iter()
        .map(|k| by_key.remove(&k).unwrap())
        .collect();
    results.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.pos.cmp(&b.pos))
            .then(a.ref_bases.cmp(&b.ref_bases))
            .then(a.alt.cmp(&b.alt))
    });

    Ok(CompareSamplesResponse {
        samples: samples.to_vec(),
        query_type: "region",
        count: results.len(),
        truncated: any_truncated,
        results,
    })
}

async fn compare_by_rsids(
    registry: Arc<SampleRegistry>,
    cache: Arc<RsidCache>,
    samples: &[String],
    rsids: &[String],
) -> Result<CompareSamplesResponse> {
    use std::collections::HashMap;

    let mut per_sample: HashMap<String, Vec<LookupRsidEntry>> = HashMap::new();
    for name in samples {
        let sub_started = Instant::now();
        tracing::info!(
            perf = true,
            phase = "compare:subcall",
            event = "start",
            sample = %name,
            rsid_count = rsids.len(),
            "compare_samples per-sample subcall begin"
        );
        let resp = lookup_rsids(
            registry.clone(),
            cache.clone(),
            LookupRsidsArgs {
                sample: name.clone(),
                rsids: rsids.to_vec(),
            },
        )
        .await?;
        tracing::info!(
            perf = true,
            phase = "compare:subcall",
            event = "end",
            sample = %name,
            rsid_count = rsids.len(),
            found_count = resp.found_count,
            elapsed_ms = sub_started.elapsed().as_millis() as u64,
            "compare_samples per-sample subcall complete"
        );
        per_sample.insert(name.clone(), resp.results);
    }

    let mut results: Vec<CompareVariantEntry> = Vec::with_capacity(rsids.len());
    for rsid in rsids {
        // Use the first sample with a hit (in input sample order) to populate
        // chrom/pos/ref/alt metadata.
        let mut meta: Option<&VariantRecord> = None;
        for name in samples {
            let entries = per_sample.get(name).unwrap();
            if let Some(e) = entries
                .iter()
                .find(|e| e.rsid == *rsid && e.found && e.variant.is_some())
            {
                meta = e.variant.as_ref();
                break;
            }
        }

        let mut sample_results = Vec::with_capacity(samples.len());
        for name in samples {
            let entries = per_sample.get(name).unwrap();
            let variant = entries
                .iter()
                .find(|e| e.rsid == *rsid && e.found)
                .and_then(|e| e.variant.clone());
            sample_results.push(make_sample_result(name, variant));
        }

        results.push(CompareVariantEntry {
            rsid: Some(rsid.clone()),
            chrom: meta.map(|v| v.chrom.clone()),
            pos: meta.map(|v| v.pos),
            ref_bases: meta.map(|v| v.ref_bases.clone()),
            alt: meta.map(|v| v.alt.clone()),
            results: sample_results,
        });
    }

    Ok(CompareSamplesResponse {
        samples: samples.to_vec(),
        query_type: "rsids",
        count: results.len(),
        truncated: false,
        results,
    })
}

fn make_sample_result(name: &str, variant: Option<VariantRecord>) -> CompareSampleResult {
    match variant {
        Some(v) => CompareSampleResult {
            sample: name.to_string(),
            found: true,
            genotype: v.genotype,
            genotype_alleles: v.genotype_alleles,
            depth: v.depth,
            gq: v.gq,
            filter: Some(v.filter),
        },
        None => CompareSampleResult {
            sample: name.to_string(),
            found: false,
            genotype: None,
            genotype_alleles: None,
            depth: None,
            gq: None,
            filter: None,
        },
    }
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

    pub fn clear(&self) {
        if let Ok(mut w) = self.inner.write() {
            w.clear();
        }
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
    registry: Arc<SampleRegistry>,
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

    let sample = registry
        .get(&args.sample)
        .ok_or_else(|| Error::SampleNotFound(args.sample.clone()))?;

    let idx = ensure_rsid_index(&cache, &sample).await?;

    let vcf_path = sample.vcf_path.clone();
    let rsids = args.rsids.clone();
    let registry_for_blocking = registry.clone();
    let work = tokio::task::spawn_blocking(move || {
        lookup_rsids_blocking(&sample, &registry_for_blocking, &idx, &rsids)
    });
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
        tracing::info!(
            perf = true,
            phase = "rsid_cache",
            event = "hit",
            sample = %sample.name,
            entries = idx.len(),
            "rsid cache hit"
        );
        return Ok(idx);
    }
    let path = sample.vcf_path.clone();
    tracing::info!(
        perf = true,
        phase = "rsid_cache",
        event = "build_start",
        sample = %sample.name,
        vcf = %path.display(),
        "rsid cache miss; streaming VCF to build index"
    );
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
    let elapsed_secs = elapsed.as_secs_f64().max(0.001);
    tracing::info!(
        perf = true,
        phase = "rsid_cache",
        event = "build_done",
        sample = %sample.name,
        entries = entries,
        elapsed_ms = elapsed.as_millis() as u64,
        rate_per_sec = (entries as f64 / elapsed_secs) as u64,
        "rsid cache built"
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
    let started = Instant::now();
    let mut records: u64 = 0;
    let mut last_checkpoint = started;
    loop {
        buf.clear();
        if reader.read_line(&mut buf)? == 0 {
            break;
        }
        if buf.starts_with('#') {
            continue;
        }
        records += 1;
        // Heartbeat every ~2 seconds — if you see one of these in the log
        // followed by silence, the next records that were being read are
        // where the hang is.
        if records % 500_000 == 0 && last_checkpoint.elapsed().as_secs() >= 2 {
            let secs = started.elapsed().as_secs_f64().max(0.001);
            tracing::info!(
                perf = true,
                phase = "rsid_cache",
                event = "build_progress",
                records = records,
                rsids_so_far = idx.len(),
                rate_per_sec = (records as f64 / secs) as u64,
                "still building rsid cache"
            );
            last_checkpoint = Instant::now();
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
    sample: &Sample,
    registry: &SampleRegistry,
    idx: &RsidIndex,
    rsids: &[String],
) -> Result<Vec<LookupRsidEntry>> {
    let started = Instant::now();
    let vcf_path = sample.vcf_path.as_path();
    let mut reader = open_indexed_reader_with_cache(registry, sample)?;
    let header = reader.read_header().map_err(|e| Error::VcfRead {
        path: vcf_path.to_path_buf(),
        message: e.to_string(),
    })?;
    let setup_ms = started.elapsed().as_millis() as u64;
    tracing::info!(
        perf = true,
        phase = "lookup_rsids:setup",
        sample = %sample.name,
        rsid_count = rsids.len(),
        elapsed_ms = setup_ms,
        "indexed reader open + header read"
    );

    let queries_started = Instant::now();
    let mut hits = 0usize;
    let mut misses = 0usize;
    let mut out = Vec::with_capacity(rsids.len());
    for rsid in rsids {
        let Some(loc) = idx.get(rsid) else {
            misses += 1;
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
        if matched.is_some() {
            hits += 1;
        } else {
            misses += 1;
        }
        out.push(LookupRsidEntry {
            rsid: rsid.clone(),
            found: matched.is_some(),
            variant: matched,
        });
    }
    let queries_ms = queries_started.elapsed().as_millis() as u64;
    tracing::info!(
        perf = true,
        phase = "lookup_rsids:queries",
        sample = %sample.name,
        rsid_count = rsids.len(),
        hits = hits,
        misses = misses,
        elapsed_ms = queries_ms,
        avg_ms_per_rsid = (queries_ms as f64 / rsids.len().max(1) as f64) as u64,
        "per-rsid tabix queries complete"
    );
    Ok(out)
}

// ---------- add_sample / add_samples_from_folder / remove_sample ----------

const BGZF_MAGIC: [u8; 4] = [0x1F, 0x8B, 0x08, 0x04];
const FOLDER_DEFAULT_MAX: usize = 50;
const FOLDER_HARD_CAP: usize = 200;

#[derive(Debug, Clone)]
pub struct AddSampleArgs {
    pub path: String,
    pub name: Option<String>,
    pub build: Option<String>,
    pub description: Option<String>,
}

#[derive(Debug, Clone)]
pub struct AddSamplesFromFolderArgs {
    pub folder: String,
    pub recursive: Option<bool>,
    pub max_files: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct SkippedFile {
    pub path: String,
    pub reason: String,
}

#[derive(Debug, Serialize)]
pub struct FolderScanResponse {
    pub folder: String,
    pub scanned: usize,
    pub registered: Vec<Sample>,
    pub skipped: Vec<SkippedFile>,
}

pub async fn add_sample(
    registry: Arc<SampleRegistry>,
    allowed_roots: Arc<Vec<std::path::PathBuf>>,
    args: AddSampleArgs,
) -> Result<Sample> {
    tokio::task::spawn_blocking(move || add_sample_blocking(&registry, &allowed_roots, args))
        .await
        .map_err(|e| Error::PathInvalid {
            path: std::path::PathBuf::new(),
            reason: format!("join: {e}"),
        })?
}

pub async fn remove_sample(registry: Arc<SampleRegistry>, name: String) -> Result<Option<Sample>> {
    tokio::task::spawn_blocking(move || registry.remove(&name))
        .await
        .map_err(|e| Error::PathInvalid {
            path: std::path::PathBuf::new(),
            reason: format!("join: {e}"),
        })?
}

#[derive(Debug, Serialize)]
pub struct ResetSamplesResponse {
    pub removed_count: usize,
    pub removed: Vec<Sample>,
}

/// Kick off background warmup for every currently-registered sample.
///
/// Returns immediately. Spawns one task per sample which:
/// 1. Acquires the tabix index (loads `.tbi` from disk if present, builds in
///    memory otherwise) and caches it on the registry.
/// 2. Builds the rsid → position cache and stores it on the rsid_cache.
///
/// Why this matters: on a cold start (server respawn after Claude Desktop's
/// stdio timeout, or first run after `add_sample`), the first tool call that
/// needs a cache blocks for ~4 s per sample while the cache builds. For 2-3
/// samples that's >10 s, which can exceed the MCP client's request timeout
/// → client kills the pipe → respawn → cold caches again → death spiral.
///
/// By doing the work in the background as soon as the server starts, the
/// first user-driven tool call (typically list_samples / server_info / a
/// region query) returns instantly; by the time the user issues a big
/// `compare_samples` or `lookup_rsids`, the caches are usually warm.
///
/// Warmup errors are logged but don't propagate — startup must always
/// succeed even if e.g. a previously-registered VCF was deleted off disk.
pub fn warmup_samples_in_background(registry: Arc<SampleRegistry>, rsid_cache: Arc<RsidCache>) {
    tokio::spawn(async move {
        let samples = registry.list();
        let total = samples.len();
        if total == 0 {
            tracing::info!(
                perf = true,
                phase = "warmup",
                event = "skip",
                "no registered samples to warm up"
            );
            return;
        }
        tracing::info!(
            perf = true,
            phase = "warmup",
            event = "all_start",
            sample_count = total,
            "background sample warmup begin"
        );
        let started = Instant::now();
        for sample in samples {
            warmup_one_sample(&registry, &rsid_cache, &sample).await;
        }
        tracing::info!(
            perf = true,
            phase = "warmup",
            event = "all_end",
            sample_count = total,
            elapsed_ms = started.elapsed().as_millis() as u64,
            "background sample warmup complete"
        );
    });
}

async fn warmup_one_sample(registry: &SampleRegistry, rsid_cache: &RsidCache, sample: &Sample) {
    let started = Instant::now();
    tracing::info!(
        perf = true,
        phase = "warmup",
        event = "sample_start",
        sample = %sample.name,
        "warming caches for sample"
    );

    // Tabix index: skip if already cached. Errors warn-log but don't abort
    // — a missing or unreadable VCF here just means the user gets the real
    // error the next time they query this sample.
    if registry.get_tabix_index(&sample.name).is_none() {
        let path = sample.vcf_path.clone();
        let name = sample.name.clone();
        let result = tokio::task::spawn_blocking(move || acquire_tabix_index(&path)).await;
        match result {
            Ok(Ok(idx)) => registry.set_tabix_index(name, Arc::new(idx)),
            Ok(Err(e)) => tracing::warn!(
                perf = true,
                phase = "warmup",
                event = "tabix_failed",
                sample = %sample.name,
                error = %e,
                "tabix index warmup failed (continuing — query path will retry)"
            ),
            Err(e) => tracing::warn!(
                perf = true,
                phase = "warmup",
                event = "tabix_panicked",
                sample = %sample.name,
                error = %e,
                "tabix index warmup task panicked"
            ),
        }
    }

    // RSID cache: skip if already populated. ensure_rsid_index emits its own
    // perf events and logs errors; we just want the side-effect of warming.
    if rsid_cache.get(&sample.name).is_none() {
        let _ = ensure_rsid_index(rsid_cache, sample).await;
    }

    tracing::info!(
        perf = true,
        phase = "warmup",
        event = "sample_end",
        sample = %sample.name,
        elapsed_ms = started.elapsed().as_millis() as u64,
        "sample warmup complete"
    );
}

#[derive(Debug, Serialize)]
pub struct ServerLimits {
    pub max_region_bp: u64,
    pub max_records_per_query: usize,
    pub max_rsids_per_lookup: usize,
    pub folder_default_max_files: usize,
    pub folder_hard_cap_files: usize,
    pub query_timeout_secs: u64,
}

#[derive(Debug, Serialize)]
pub struct ServerInfoResponse {
    pub name: &'static str,
    /// Semver from Cargo.toml.
    pub version: &'static str,
    /// Monotonic build counter — git commit count when built from a git
    /// checkout, otherwise 0.
    pub build: &'static str,
    /// Short git SHA of the commit this binary was built from, or "unknown".
    pub commit: &'static str,
    /// Human-friendly composite, e.g. "vcf-mcp 0.1.0+build.42.83f78b1".
    pub version_string: String,
    pub ensembl_release_grch38: u32,
    pub ensembl_release_grch37: u32,
    pub registered_samples: usize,
    pub limits: ServerLimits,
}

/// Snapshot of the running server's identity and state, exposed as a tool so
/// the LLM can answer "what version are you?" / "which gene table?" / "what
/// are the limits?" without guessing.
pub fn server_info(registry: &SampleRegistry) -> ServerInfoResponse {
    let name = env!("CARGO_PKG_NAME");
    let version = env!("CARGO_PKG_VERSION");
    let build = env!("VCF_MCP_BUILD");
    let commit = env!("VCF_MCP_COMMIT");
    let version_string = format!("{name} {version}+build.{build}.{commit}");
    ServerInfoResponse {
        name,
        version,
        build,
        commit,
        version_string,
        ensembl_release_grch38: crate::genes::ENSEMBL_RELEASE_GRCH38,
        ensembl_release_grch37: crate::genes::ENSEMBL_RELEASE_GRCH37,
        registered_samples: registry.len(),
        limits: ServerLimits {
            max_region_bp: MAX_REGION_BP,
            max_records_per_query: MAX_RECORDS,
            max_rsids_per_lookup: MAX_RSIDS_PER_LOOKUP,
            folder_default_max_files: FOLDER_DEFAULT_MAX,
            folder_hard_cap_files: FOLDER_HARD_CAP,
            query_timeout_secs: QUERY_TIMEOUT_SECS,
        },
    }
}

/// Wipe the registry: drop all samples, clear caches, and rewrite the state
/// file as empty. Useful for workflow / cowork steps ("clean slate before
/// re-registering"). Destructive — caller must pass `confirm=true` upstream;
/// this function trusts that the tool layer enforced it.
pub async fn reset_samples(
    registry: Arc<SampleRegistry>,
    rsid_cache: Arc<RsidCache>,
) -> Result<ResetSamplesResponse> {
    let removed = tokio::task::spawn_blocking(move || {
        let r = registry.clear()?;
        rsid_cache.clear();
        Ok::<Vec<Sample>, Error>(r)
    })
    .await
    .map_err(|e| Error::StateFile {
        path: std::path::PathBuf::new(),
        message: format!("reset join: {e}"),
    })??;
    tracing::info!(
        removed = removed.len(),
        "registry reset (all samples and caches cleared)"
    );
    Ok(ResetSamplesResponse {
        removed_count: removed.len(),
        removed,
    })
}

pub async fn add_samples_from_folder(
    registry: Arc<SampleRegistry>,
    allowed_roots: Arc<Vec<std::path::PathBuf>>,
    args: AddSamplesFromFolderArgs,
) -> Result<FolderScanResponse> {
    tokio::task::spawn_blocking(move || folder_scan_blocking(&registry, &allowed_roots, args))
        .await
        .map_err(|e| Error::PathInvalid {
            path: std::path::PathBuf::new(),
            reason: format!("join: {e}"),
        })?
}

fn add_sample_blocking(
    registry: &SampleRegistry,
    allowed_roots: &[std::path::PathBuf],
    args: AddSampleArgs,
) -> Result<Sample> {
    let raw = std::path::PathBuf::from(&args.path);
    let canonical = std::fs::canonicalize(&raw).map_err(|e| Error::PathInvalid {
        path: raw.clone(),
        reason: e.to_string(),
    })?;

    // Allowlist check (optional).
    if !allowed_roots.is_empty() {
        let canonical_roots: Vec<_> = allowed_roots
            .iter()
            .filter_map(|r| std::fs::canonicalize(r).ok())
            .collect();
        if !canonical_roots
            .iter()
            .any(|root| canonical.starts_with(root))
        {
            return Err(Error::PathNotAllowed {
                path: canonical,
                allowed: canonical_roots,
            });
        }
    }

    validate_vcf_path(&canonical)?;

    // Acquire the tabix index (load .tbi if present, otherwise build in
    // memory — never write). The Arc-wrapped index gets cached on the
    // registry after we successfully complete validation.
    let index = acquire_tabix_index(&canonical)?;

    let file = std::fs::File::open(&canonical).map_err(|source| Error::VcfOpen {
        path: canonical.clone(),
        source,
    })?;
    // Clone the freshly-built index for the reader; we still hold the
    // original to cache against the eventual sample name once we know it.
    let mut reader = vcf::io::IndexedReader::new(file, index.clone());
    let header = reader.read_header().map_err(|e| Error::InvalidVcfFile {
        path: canonical.clone(),
        reason: format!("header parse: {e}"),
    })?;

    if header.contigs().is_empty() {
        return Err(Error::InvalidVcfFile {
            path: canonical,
            reason: "no ##contig lines in header".to_string(),
        });
    }
    if header.sample_names().is_empty() {
        return Err(Error::InvalidVcfFile {
            path: canonical,
            reason: "no sample columns in header (sites-only VCF not supported)".to_string(),
        });
    }

    // Functional probe: pick a contig the tabix INDEX actually knows about
    // (not just one from the VCF ##contig list — a sliced fixture often has
    // a full set of header contigs but data for only a subset). Query a wide
    // range; zero records is fine, what matters is that the operation
    // completes, proving index ↔ data consistency.
    let probe_contig = reader
        .index()
        .header()
        .and_then(|h| h.reference_sequence_names().iter().next().cloned())
        .ok_or_else(|| Error::InvalidVcfFile {
            path: canonical.clone(),
            reason: "tabix index has no indexed reference sequences".to_string(),
        })?;
    // 250M comfortably covers the longest human chromosome (chr1 ≈ 249 Mb).
    let probe_start = Position::try_from(1usize).unwrap();
    let probe_end = Position::try_from(250_000_000usize).unwrap();
    let region = Region::new(probe_contig.clone(), probe_start..=probe_end);
    let probe = reader
        .query(&header, &region)
        .map_err(|e| Error::InvalidVcfFile {
            path: canonical.clone(),
            reason: format!("functional probe on {probe_contig}: {e}"),
        })?;
    drop(probe);

    // Detect or validate the build.
    let build = match args.build {
        Some(b) if b == "GRCh37" || b == "GRCh38" => b,
        Some(b) => {
            return Err(Error::InvalidBuild {
                sample: args
                    .name
                    .clone()
                    .unwrap_or_else(|| derive_name_from_path(&canonical)),
                build: b,
            });
        }
        None => detect_build(&header).ok_or_else(|| Error::BuildNotDetectable {
            path: canonical.clone(),
        })?,
    };

    // Path-idempotency only when the caller didn't ask for a specific name.
    // If the user explicitly says "register this file under name X", honor
    // it even if the same path is already registered under name Y.
    let provided_name = args.name.filter(|n| !n.is_empty());
    if provided_name.is_none() {
        if let Some(existing) = registry.get_by_path(&canonical) {
            return Ok(existing);
        }
    }
    let name = provided_name.unwrap_or_else(|| derive_name_from_path(&canonical));

    let sample = Sample {
        name,
        vcf_path: canonical,
        build,
        description: args.description.unwrap_or_default(),
    };

    let registered = registry.add_validated(sample)?;
    // Cache the index we just loaded/built so the first query against this
    // sample doesn't re-acquire it.
    registry.set_tabix_index(registered.name.clone(), Arc::new(index));
    Ok(registered)
}

fn folder_scan_blocking(
    registry: &SampleRegistry,
    allowed_roots: &[std::path::PathBuf],
    args: AddSamplesFromFolderArgs,
) -> Result<FolderScanResponse> {
    let recursive = args.recursive.unwrap_or(false);
    let max_files = args
        .max_files
        .unwrap_or(FOLDER_DEFAULT_MAX)
        .min(FOLDER_HARD_CAP);

    let folder = std::path::PathBuf::from(&args.folder);
    let folder_canonical = std::fs::canonicalize(&folder).map_err(|e| Error::PathInvalid {
        path: folder.clone(),
        reason: e.to_string(),
    })?;
    let meta = std::fs::metadata(&folder_canonical).map_err(|e| Error::PathInvalid {
        path: folder_canonical.clone(),
        reason: e.to_string(),
    })?;
    if !meta.is_dir() {
        return Err(Error::PathInvalid {
            path: folder_canonical,
            reason: "not a directory".to_string(),
        });
    }

    let mut candidates: Vec<std::path::PathBuf> = Vec::new();
    collect_vcf_paths(&folder_canonical, recursive, &mut candidates).map_err(|e| {
        Error::PathInvalid {
            path: folder_canonical.clone(),
            reason: format!("walk: {e}"),
        }
    })?;
    candidates.sort();

    if candidates.len() > max_files {
        return Err(Error::TooManyFiles {
            found: candidates.len(),
            max: max_files,
            hard_cap: FOLDER_HARD_CAP,
        });
    }

    let mut registered = Vec::new();
    let mut skipped = Vec::new();
    for path in candidates {
        let attempt = add_sample_blocking(
            registry,
            allowed_roots,
            AddSampleArgs {
                path: path.to_string_lossy().into_owned(),
                name: None,
                build: None,
                description: None,
            },
        );
        match attempt {
            Ok(s) => registered.push(s),
            Err(e) => skipped.push(SkippedFile {
                path: path.to_string_lossy().into_owned(),
                reason: e.to_string(),
            }),
        }
    }

    Ok(FolderScanResponse {
        folder: folder_canonical.to_string_lossy().into_owned(),
        scanned: registered.len() + skipped.len(),
        registered,
        skipped,
    })
}

fn collect_vcf_paths(
    dir: &std::path::Path,
    recursive: bool,
    out: &mut Vec<std::path::PathBuf>,
) -> std::io::Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let ft = entry.file_type()?;
        if ft.is_dir() && recursive {
            collect_vcf_paths(&entry.path(), recursive, out)?;
        } else if ft.is_file() {
            let n = entry.file_name();
            let ns = n.to_string_lossy().to_ascii_lowercase();
            if ns.ends_with(".vcf.gz") {
                out.push(entry.path());
            }
        }
    }
    Ok(())
}

fn validate_vcf_path(canonical: &std::path::Path) -> Result<()> {
    // Existence + regular file.
    let meta = std::fs::metadata(canonical).map_err(|e| Error::PathInvalid {
        path: canonical.to_path_buf(),
        reason: e.to_string(),
    })?;
    if !meta.is_file() {
        return Err(Error::PathInvalid {
            path: canonical.to_path_buf(),
            reason: "not a regular file".to_string(),
        });
    }

    // Filename extension.
    let fname = canonical.file_name().and_then(|s| s.to_str()).unwrap_or("");
    if !fname.to_ascii_lowercase().ends_with(".vcf.gz") {
        return Err(Error::InvalidVcfFile {
            path: canonical.to_path_buf(),
            reason: format!("filename {fname:?} does not end in .vcf.gz"),
        });
    }

    // BGZF magic bytes.
    use std::io::Read as _;
    let mut f = std::fs::File::open(canonical).map_err(|e| Error::PathInvalid {
        path: canonical.to_path_buf(),
        reason: e.to_string(),
    })?;
    let mut magic = [0u8; 4];
    let n = f.read(&mut magic).map_err(|e| Error::InvalidVcfFile {
        path: canonical.to_path_buf(),
        reason: format!("read magic bytes: {e}"),
    })?;
    if n < 4 || magic != BGZF_MAGIC {
        return Err(Error::InvalidVcfFile {
            path: canonical.to_path_buf(),
            reason: format!(
                "BGZF magic mismatch (got {:02X} {:02X} {:02X} {:02X}; expected 1F 8B 08 04)",
                magic[0], magic[1], magic[2], magic[3]
            ),
        });
    }

    // Note: .tbi presence is no longer required here. The acquire_tabix_index
    // helper in add_sample_blocking either loads an existing .tbi or builds
    // the index entirely in memory — no write to the source folder.
    Ok(())
}

/// Acquire a tabix index for the given VCF without touching the filesystem.
/// If `<vcf>.tbi` exists, deserialize it (cheap). If not, stream the bgzipped
/// VCF and build the index in memory. The returned index is owned by the
/// caller and typically gets cached on the `SampleRegistry`.
fn acquire_tabix_index(canonical: &std::path::Path) -> Result<noodles_tabix::Index> {
    use noodles_tabix as tabix;

    let mut tbi = canonical.as_os_str().to_owned();
    tbi.push(".tbi");
    let tbi_path = std::path::PathBuf::from(tbi);

    // The tabix index always lives in memory. If a `.tbi` happens to exist on
    // disk (e.g. produced by a sequencing pipeline or our index_vcf example),
    // we load it as a fast shortcut; otherwise we build the same index in
    // memory by streaming the file. Both are normal — neither is an error or
    // a fallback worth flagging.
    if tbi_path.exists() {
        let started = Instant::now();
        let result = tabix::fs::read(&tbi_path).map_err(|e| Error::InvalidVcfFile {
            path: canonical.to_path_buf(),
            reason: format!("loading existing .tbi: {e}"),
        });
        tracing::info!(
            perf = true,
            phase = "tabix_index",
            event = "ready",
            source = "disk",
            vcf = %canonical.display(),
            elapsed_ms = started.elapsed().as_millis() as u64,
            "tabix index ready (loaded a .tbi shortcut from disk)"
        );
        return result;
    }

    let started = Instant::now();
    tracing::info!(
        perf = true,
        phase = "tabix_index",
        event = "build_start",
        source = "memory",
        vcf = %canonical.display(),
        "building in-memory tabix index"
    );
    let index = build_tabix_index_in_memory(canonical).map_err(|e| Error::InvalidVcfFile {
        path: canonical.to_path_buf(),
        reason: format!("tabix index build failed: {e}"),
    })?;
    tracing::info!(
        perf = true,
        phase = "tabix_index",
        event = "ready",
        source = "memory",
        vcf = %canonical.display(),
        elapsed_ms = started.elapsed().as_millis() as u64,
        "tabix index ready (built in memory)"
    );
    Ok(index)
}

/// Stream a bgzipped VCF and produce a tabix index in memory. Same logic as
/// `examples/index_vcf.rs` but returns the index instead of writing it.
fn build_tabix_index_in_memory(
    vcf_path: &std::path::Path,
) -> std::io::Result<noodles_tabix::Index> {
    use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
    use noodles_tabix as tabix;

    let file = std::fs::File::open(vcf_path)?;
    let mut reader = bgzf::io::Reader::new(file);

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut line_start = reader.virtual_position();
    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        let line_end = reader.virtual_position();
        if buf.starts_with('#') {
            line_start = line_end;
            continue;
        }
        let mut it = buf.split('\t');
        let chrom = it
            .next()
            .ok_or_else(|| std::io::Error::other("missing CHROM"))?;
        let pos_str = it
            .next()
            .ok_or_else(|| std::io::Error::other("missing POS"))?;
        let _id = it.next();
        let ref_bases = it
            .next()
            .ok_or_else(|| std::io::Error::other("missing REF"))?;

        let pos: usize = pos_str
            .trim()
            .parse()
            .map_err(|e| std::io::Error::other(format!("bad POS {pos_str:?}: {e}")))?;
        let start = Position::try_from(pos)
            .map_err(|e| std::io::Error::other(format!("start out of range: {e}")))?;
        let end_pos = pos + ref_bases.len().saturating_sub(1);
        let end = Position::try_from(end_pos.max(pos))
            .map_err(|e| std::io::Error::other(format!("end out of range: {e}")))?;
        indexer
            .add_record(chrom, start, end, Chunk::new(line_start, line_end))
            .map_err(|e| std::io::Error::other(format!("indexer.add_record: {e}")))?;
        line_start = line_end;
    }
    Ok(indexer.build())
}

/// Open an `IndexedReader` for the given sample using a cached or
/// freshly-built tabix index — never reads `<vcf>.tbi` from disk twice and
/// never writes one. The returned reader's lifetime is bound to the caller.
fn open_indexed_reader_with_cache(
    registry: &SampleRegistry,
    sample: &Sample,
) -> Result<vcf::io::IndexedReader<bgzf::io::Reader<std::fs::File>>> {
    let index = match registry.get_tabix_index(&sample.name) {
        Some(idx) => (*idx).clone(),
        None => {
            let fresh = acquire_tabix_index(&sample.vcf_path)?;
            registry.set_tabix_index(sample.name.clone(), Arc::new(fresh.clone()));
            fresh
        }
    };
    let file = std::fs::File::open(&sample.vcf_path).map_err(|source| Error::VcfOpen {
        path: sample.vcf_path.clone(),
        source,
    })?;
    Ok(vcf::io::IndexedReader::new(file, index))
}

/// Detect genome build from a VCF header. Tries, in order:
/// 1. `##reference=...` substring match for known build names.
/// 2. `##contig=<...assembly=...>` field on any contig.
/// 3. chr1 length heuristic (GRCh38: 248,956,422; GRCh37: 249,250,621).
fn detect_build(header: &vcf::Header) -> Option<String> {
    use noodles_vcf::header::record::value::Collection;

    fn classify(s: &str) -> Option<&'static str> {
        let l = s.to_ascii_lowercase();
        if l.contains("grch38") || l.contains("hg38") {
            return Some("GRCh38");
        }
        if l.contains("grch37") || l.contains("hg19") {
            return Some("GRCh37");
        }
        None
    }

    // Method 1: ##reference line (or any other unstructured key carrying a build hint).
    for (key, records) in header.other_records() {
        if !key.as_ref().eq_ignore_ascii_case("reference") {
            continue;
        }
        if let Collection::Unstructured(values) = records {
            for s in values {
                if let Some(b) = classify(s) {
                    return Some(b.to_string());
                }
            }
        }
    }

    // Method 2: contig with assembly field (non-standard tag in other_fields).
    for (_name, contig) in header.contigs() {
        for (k, v) in contig.other_fields() {
            if k.as_ref().eq_ignore_ascii_case("assembly") {
                if let Some(b) = classify(v) {
                    return Some(b.to_string());
                }
            }
        }
    }

    // Method 3: chr1 / 1 length heuristic.
    for (name, contig) in header.contigs() {
        let n: &str = name;
        if n == "chr1" || n == "1" {
            if let Some(len) = contig.length() {
                let len = len as i64;
                if (len - 248_956_422).abs() < 1000 {
                    return Some("GRCh38".to_string());
                }
                if (len - 249_250_621).abs() < 1000 {
                    return Some("GRCh37".to_string());
                }
            }
        }
    }

    None
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
