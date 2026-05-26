use std::borrow::Cow;
use std::sync::Arc;
use std::time::Duration;

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

    let alt = {
        let mut alts: Vec<String> = Vec::new();
        for r in rec.alternate_bases().iter() {
            alts.push(r.map_err(|e| e.to_string())?.to_string());
        }
        alts.join(",")
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
}
