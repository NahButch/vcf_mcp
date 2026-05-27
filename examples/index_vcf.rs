//! Build a tabix index (.tbi) for an existing bgzipped VCF.
//!
//! Usage:
//!     cargo run --release --example index_vcf -- <path.vcf.gz>
//!
//! Writes `<path.vcf.gz>.tbi` next to the input. Useful when the VCF you
//! were handed by a sequencing vendor / pipeline doesn't ship with a
//! tabix index and you don't have `tabix` / `bcftools` on your PATH.
//!
//! Assumes the input is BGZF-compressed (the standard for tabix-indexable
//! VCFs) and sorted by chromosome, then position. A plain-gzip VCF will
//! fail to read; recompress with bgzip first.

use std::fs::File;
use std::io::BufRead;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() != 2 {
        bail!("usage: index_vcf <path.vcf.gz>");
    }
    let src = PathBuf::from(&argv[1]);
    let dst = PathBuf::from(format!("{}.tbi", argv[1]));

    let file = File::open(&src).with_context(|| format!("opening {}", src.display()))?;
    let mut reader = bgzf::io::Reader::new(file);

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut line_start = reader.virtual_position();
    let mut buf = String::new();
    let mut records: u64 = 0;
    let started = Instant::now();
    let mut last_report = started;

    loop {
        buf.clear();
        let n = reader
            .read_line(&mut buf)
            .context("reading line from bgzf")?;
        if n == 0 {
            break;
        }
        let line_end = reader.virtual_position();

        // Headers don't get indexed; skip.
        if buf.starts_with('#') {
            line_start = line_end;
            continue;
        }

        // Minimal field parse: CHROM \t POS \t ID \t REF \t ...
        let mut it = buf.split('\t');
        let chrom = it.next().context("missing CHROM column")?;
        let pos_str = it.next().context("missing POS column")?;
        let _id = it.next();
        let ref_bases = it.next().context("missing REF column")?;

        let pos: usize = pos_str
            .trim()
            .parse()
            .with_context(|| format!("bad POS value {pos_str:?} on record {}", records + 1))?;
        let start = Position::try_from(pos).context("POS out of valid range")?;
        // End is pos + len(REF) - 1 (1-based inclusive).
        let end_pos = pos + ref_bases.len().saturating_sub(1);
        let end = Position::try_from(end_pos.max(pos)).context("end position out of range")?;

        indexer
            .add_record(chrom, start, end, Chunk::new(line_start, line_end))
            .with_context(|| format!("indexer.add_record at record {}", records + 1))?;

        line_start = line_end;
        records += 1;

        if last_report.elapsed().as_secs() >= 2 {
            let rate = records as f64 / started.elapsed().as_secs_f64();
            eprintln!("  {} records indexed ({:.0}/s)", fmt_count(records), rate);
            last_report = Instant::now();
        }
    }

    let elapsed = started.elapsed();
    eprintln!(
        "Indexed {} records in {:.1}s",
        fmt_count(records),
        elapsed.as_secs_f64()
    );

    let index = indexer.build();
    let f = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
    let mut writer = tabix::io::Writer::new(f);
    writer.write_index(&index).context("writing tabix index")?;
    eprintln!("Wrote {}", dst.display());
    Ok(())
}

fn fmt_count(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}
