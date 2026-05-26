//! Slice a region out of a bgzipped, tabix-indexed VCF and write a new
//! bgzipped slice with its own tabix index.
//!
//! Usage:
//!     cargo run --release --example slice_vcf -- <input.vcf.gz> <output.vcf.gz> <region>
//!
//! Region format follows tabix conventions, e.g. `17:50100000-50300000`.

use std::fs::File;

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use noodles_vcf::{self as vcf, variant::Record as _, variant::io::Write as _};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        bail!("usage: slice_vcf <input.vcf.gz> <output.vcf.gz> <region>");
    }
    let src = &args[1];
    let dst = &args[2];
    let region_str = &args[3];

    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(src)
        .with_context(|| format!("opening indexed reader for {src}"))?;
    let header = reader.read_header().context("reading header")?;
    let region: Region = region_str
        .parse()
        .with_context(|| format!("parsing region {region_str:?}"))?;

    let out_file = File::create(dst).with_context(|| format!("creating {dst}"))?;
    let mut writer = vcf::io::Writer::new(bgzf::io::Writer::new(out_file));
    writer.write_header(&header).context("writing header")?;

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut start_position = writer.get_ref().virtual_position();
    let mut count = 0usize;

    let query = reader
        .query(&header, &region)
        .context("starting region query")?;
    for result in query.records() {
        let record = result.context("reading record")?;
        writer
            .write_variant_record(&header, &record)
            .context("writing record")?;
        let end_position = writer.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let chrom = record.reference_sequence_name();
        let start = record
            .variant_start()
            .context("record missing variant_start")??;
        let end = record
            .variant_end(&header)
            .context("computing variant_end")?;
        indexer
            .add_record(chrom, start, end, chunk)
            .context("indexer add_record")?;

        start_position = end_position;
        count += 1;
    }

    let bgzf_writer = writer.into_inner();
    bgzf_writer.finish().context("finishing bgzf stream")?;

    let index = indexer.build();
    let idx_path = format!("{dst}.tbi");
    let idx_file = File::create(&idx_path).with_context(|| format!("creating {idx_path}"))?;
    let mut idx_writer = tabix::io::Writer::new(idx_file);
    idx_writer
        .write_index(&index)
        .context("writing tabix index")?;

    eprintln!("Wrote {count} records from {region_str} to {dst}");
    eprintln!("Index: {idx_path}");
    Ok(())
}
