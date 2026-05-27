//! Rewrite a bgzipped VCF's ID column with deterministic synthetic rsids,
//! emit a new bgzipped VCF + tabix index.
//!
//! Used to build a lookup_rsids test fixture from our existing NA12878 slice,
//! which doesn't carry dbSNP IDs of its own. The Nth data record gets ID
//! `rs<1_000_000+N>` — so a 279-record input produces rs1000001..rs1000279.
//!
//! Usage:
//!     cargo run --release --example annotate_rsids -- <input.vcf.gz> <output.vcf.gz>

use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() != 3 {
        bail!("usage: annotate_rsids <input.vcf.gz> <output.vcf.gz>");
    }
    let src = PathBuf::from(&argv[1]);
    let dst = PathBuf::from(&argv[2]);
    let idx_dst = PathBuf::from(format!("{}.tbi", &argv[2]));

    let in_file = File::open(&src).with_context(|| format!("opening {}", src.display()))?;
    let mut reader = bgzf::io::Reader::new(in_file);

    let out_file = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
    let mut writer = bgzf::io::Writer::new(out_file);

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut buf = String::new();
    let mut record_idx: u32 = 0;
    let mut line_start = writer.virtual_position();

    loop {
        buf.clear();
        let n = reader.read_line(&mut buf).context("reading line")?;
        if n == 0 {
            break;
        }

        if buf.starts_with('#') {
            writer
                .write_all(buf.as_bytes())
                .context("writing header line")?;
            line_start = writer.virtual_position();
            continue;
        }

        // Split on tabs only for the first 4 columns; rest stays as a single tail.
        let trimmed = buf.trim_end_matches('\n');
        let mut parts = trimmed.splitn(5, '\t');
        let chrom = parts.next().context("missing CHROM")?;
        let pos_str = parts.next().context("missing POS")?;
        let _orig_id = parts.next().context("missing ID")?;
        let ref_bases = parts.next().context("missing REF")?;
        let rest = parts.next().unwrap_or(""); // ALT \t QUAL \t FILTER \t INFO ...

        record_idx += 1;
        let new_id = format!("rs{}", 1_000_000 + record_idx);

        let line_out = format!("{chrom}\t{pos_str}\t{new_id}\t{ref_bases}\t{rest}\n");
        writer
            .write_all(line_out.as_bytes())
            .context("writing data line")?;

        let line_end = writer.virtual_position();

        let pos: usize = pos_str
            .trim()
            .parse()
            .with_context(|| format!("bad POS {pos_str:?}"))?;
        let start = Position::try_from(pos).context("POS out of range")?;
        let end_pos = pos + ref_bases.len().saturating_sub(1);
        let end = Position::try_from(end_pos.max(pos)).context("end out of range")?;

        indexer
            .add_record(chrom, start, end, Chunk::new(line_start, line_end))
            .with_context(|| format!("indexer.add_record at record {record_idx}"))?;

        line_start = line_end;
    }

    writer.finish().context("finishing bgzf stream")?;

    let index = indexer.build();
    let idx_file =
        File::create(&idx_dst).with_context(|| format!("creating {}", idx_dst.display()))?;
    let mut idx_writer = tabix::io::Writer::new(idx_file);
    idx_writer.write_index(&index).context("writing index")?;

    eprintln!(
        "Wrote {} records (rs1000001..rs{}) to {}",
        record_idx,
        1_000_000 + record_idx,
        dst.display()
    );
    eprintln!("Index: {}", idx_dst.display());
    Ok(())
}
