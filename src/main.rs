// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

extern crate clap;
extern crate csv;
extern crate data_encoding;
extern crate failure;
extern crate rayon;
extern crate ring;
extern crate rust_htslib;
extern crate simplelog;
extern crate tempfile;
extern crate terminal_size;
#[macro_use]
extern crate log;
extern crate faccess;
extern crate human_panic;

use clap::{App, Arg};
use faccess::{AccessMode, PathExt};
use failure::Error;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use simplelog::*;
use std::cmp;
use std::collections::HashSet;
use std::fs;
use std::io::prelude::*;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;
use terminal_size::{terminal_size, Width};

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("subset-bam")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("1.1.0")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com>, Wyatt McDonnell <wyatt.mcdonnell@10xgenomics.com>")
        .about("Subsetting 10x Genomics BAM files")
        .arg(Arg::with_name("bam")
             .short("b")
             .long("bam")
             .value_name("FILE")
             .help("Cellranger BAM file.")
             .required(true))
        .arg(Arg::with_name("cell_barcodes")
             .short("c")
             .long("cell-barcodes")
             .value_name("FILE")
             .help("File with cell barcodes to be extracted.")
             .required(true))
        .arg(Arg::with_name("out_bam")
             .short("o")
             .long("out-bam")
             .value_name("OUTPUT_FILE")
             .help("Output BAM.")
             .required(true))
        .arg(Arg::with_name("log_level")
             .long("log-level")
             .possible_values(&["info", "debug", "error"])
             .default_value("error")
             .help("Logging level."))
        .arg(Arg::with_name("cores")
             .long("cores")
             .default_value("1")
             .value_name("INTEGER")
             .help("Number of cores to use. If larger than 1, will write BAM subsets to temporary files before merging."))
        .arg(Arg::with_name("bam_tag")
             .long("bam-tag")
             .default_value("CB")
             .help("Change from default value (CB) to subset alignments based on alternative tags."));
    args
}

pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

pub struct Metrics {
    pub total_reads: usize,
    pub barcoded_reads: usize,
    pub kept_reads: usize,
}

pub struct ChunkArgs<'a> {
    cell_barcodes: &'a HashSet<Vec<u8>>,
    i: usize,
    bam_file: &'a str,
    tmp_dir: &'a Path,
    bam_tag: String,
    virtual_start: Option<i64>,
    virtual_stop: Option<i64>,
}

pub struct ChunkOuts {
    metrics: Metrics,
    out_bam_file: PathBuf,
}

fn main() {
    //setup_panic!();  // pretty panics for users
    let mut cli_args = Vec::new();
    for arg in std::env::args_os() {
        cli_args.push(arg.into_string().unwrap());
    }
    _main(cli_args);
}

fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args
        .value_of("cell_barcodes")
        .expect("You must provide a cell barcodes file");
    let out_bam_file = args
        .value_of("out_bam")
        .expect("You must provide a path to write the new BAM file");
    let ll = args.value_of("log_level").unwrap();
    let cores = args
        .value_of("cores")
        .unwrap_or_default()
        .parse::<u64>()
        .expect("Failed to convert cores to integer");
    let bam_tag = args.value_of("bam_tag").unwrap_or_default().to_string();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => {
            println!("Log level not valid");
            process::exit(1);
        }
    };
    let _ = SimpleLogger::init(ll, Config::default());

    check_inputs_exist(bam_file, cell_barcodes, out_bam_file);
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();
    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(&bam_file, &cores).unwrap();

    let mut chunks = Vec::new();
    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        let c = ChunkArgs {
            cell_barcodes: &cell_barcodes,
            i: i,
            bam_file: &bam_file,
            tmp_dir: tmp_dir.path(),
            bam_tag: bam_tag.clone(),
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
        };
        chunks.push(c);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();
    let results: Vec<_> = pool.install(|| {
        chunks
            .par_iter()
            .map(|chunk| slice_bam_chunk(chunk))
            .collect()
    });

    // combine metrics
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };

    fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
        metrics.total_reads += m.total_reads;
        metrics.barcoded_reads += m.barcoded_reads;
        metrics.kept_reads += m.kept_reads;
    }

    let mut tmp_bams = Vec::new();
    for c in results.iter() {
        add_metrics(&mut metrics, &c.metrics);
        tmp_bams.push(&c.out_bam_file);
    }

    if metrics.kept_reads == 0 {
        error!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
        process::exit(2);
    }

    // just copy the temp file over
    if cores == 1 {
        fs::copy(tmp_bams[0], out_bam_file).unwrap();
    } else {
        info!("Merging {} BAM chunks into final output", cores);
        merge_bams(tmp_bams, Path::new(out_bam_file));
    }

    info!("Done!");
    info!(
        "Visited {} alignments, found {} with barcodes and kept {}",
        metrics.total_reads, metrics.barcoded_reads, metrics.kept_reads
    );
}

pub fn check_inputs_exist(bam_file: &str, cell_barcodes: &str, out_bam_path: &str) {
    for path in [bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            error!("File {} does not exist", path);
            process::exit(1);
        }
    }
    let path = Path::new(out_bam_path);
    if path.exists() {
        error!("Output path already exists");
        process::exit(1);
    }
    if path.is_dir() {
        error!("Output path is a directory");
        process::exit(1);
    }
    let _parent_dir = path.parent();
    if _parent_dir.is_none() {
        error!("Unable to parse directory from {}", out_bam_path);
        process::exit(1);
    }
    let parent_dir = _parent_dir.unwrap();
    if (parent_dir.to_str().unwrap().len() > 0) & !parent_dir.exists() {
        error!("Output directory {:?} does not exist", parent_dir);
        process::exit(1);
    }

    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = bam_file.to_owned() + ".bai";
            if !Path::new(&bai).exists() {
                error!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = bam_file.to_owned() + ".crai";
            if !Path::new(&crai).exists() {
                error!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        &_ => {
            error!("BAM file did not end in .bam or .cram. Unable to validate");
            process::exit(1);
        }
    }
}

pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashSet<Vec<u8>>, Error> {
    let r = fs::File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashSet::new();

    for l in reader.lines() {
        let seq = l?.into_bytes();
        bc_set.insert(seq);
    }
    let num_bcs = bc_set.len();
    if num_bcs == 0 {
        error!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    debug!("Loaded {} barcodes", num_bcs);
    Ok(bc_set)
}

pub fn get_cell_barcode(rec: &Record, bam_tag: &str) -> Option<Vec<u8>> {
    //println!("{:?}", rec.aux(bam_tag.as_bytes()));
    match rec.aux(bam_tag.as_bytes()) {
        Some(Aux::String(hp)) => {
            let cb = hp.to_vec();
            Some(cb)
        }
        _ => None,
    }
}

pub fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr)?;
    Ok(out_handle)
}

pub fn bgzf_noffsets(
    bam_path: &str,
    num_chunks: &u64,
) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);

        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    // if we only have one thread, this is easy
    if *num_chunks == 1 as u64 {
        let final_offsets = vec![(None, None)];
        return Ok(final_offsets);
    }

    let bam_bytes = fs::metadata(bam_path)?.len();
    let mut initial_offsets = Vec::new();
    let step_size = bam_bytes / num_chunks;
    for n in 1..*num_chunks {
        initial_offsets.push((step_size * n) as u64);
    }

    let num_bytes = if initial_offsets.len() > 1 {
        let diff = vec_diff(&initial_offsets);
        let m = diff.iter().max().unwrap();
        cmp::min(1 << 16, *m)
    } else {
        1 << 16
    };

    // linear search to the right of each possible offset until
    // a valid virtual offset is found
    let mut adjusted_offsets = Vec::new();
    let mut fp = fs::File::open(bam_path)?;
    for offset in initial_offsets {
        fp.seek(io::SeekFrom::Start(offset))?;
        let mut buffer = [0; 2 << 16];
        fp.read(&mut buffer)?;
        for i in 0..num_bytes {
            if is_valid_bgzf_block(&buffer[i as usize..]) {
                adjusted_offsets.push(offset + i);
                break;
            }
        }
    }
    // bit-shift and produce start/stop intervals
    let mut final_offsets = Vec::new();

    // handle special case where we only found one offset
    if adjusted_offsets.len() == 1 {
        final_offsets.push((None, None));
        return Ok(final_offsets);
    }

    final_offsets.push((None, Some(((adjusted_offsets[1]) as i64) << 16)));
    for n in 2..num_chunks - 1 {
        let n = n as usize;
        final_offsets.push((
            Some((adjusted_offsets[n - 1] as i64) << 16),
            Some((adjusted_offsets[n] as i64) << 16),
        ));
    }
    final_offsets.push((
        Some(((adjusted_offsets[adjusted_offsets.len() - 1]) as i64) << 16),
        None,
    ));
    Ok(final_offsets)
}

pub fn is_valid_bgzf_block(block: &[u8]) -> bool {
    // look for the bgzip magic characters \x1f\x8b\x08\x04
    // TODO: is this sufficient?
    if block.len() < 18 {
        return false;
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false;
    }
    true
}

pub fn slice_bam_chunk(args: &ChunkArgs) -> ChunkOuts {
    let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
    let out_bam_file = args.tmp_dir.join(format!("{}.bam", args.i));
    let mut out_bam = load_writer(&bam, &out_bam_file).unwrap();
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };
    for r in bam.iter_chunk(args.virtual_start, args.virtual_stop) {
        let rec = r.unwrap();
        metrics.total_reads += 1;
        let barcode = get_cell_barcode(&rec, &args.bam_tag);
        if barcode.is_some() {
            metrics.barcoded_reads += 1;
            let barcode = barcode.unwrap();
            if args.cell_barcodes.contains(&barcode) {
                metrics.kept_reads += 1;
                out_bam.write(&rec).unwrap();
            }
        }
    }
    let r = ChunkOuts {
        metrics: metrics,
        out_bam_file: out_bam_file,
    };
    info!("Chunk {} is done", args.i);
    r
}

pub fn merge_bams(tmp_bams: Vec<&PathBuf>, out_bam_file: &Path) {
    use rust_htslib::bam::Read; // collides with fs::Read
    let bam = bam::Reader::from_path(tmp_bams[0]).unwrap();
    let mut out_bam = load_writer(&bam, out_bam_file).unwrap();
    for b in tmp_bams.iter() {
        let mut rdr = bam::Reader::from_path(b).unwrap();
        for _rec in rdr.records() {
            let rec = _rec.unwrap();
            out_bam.write(&rec).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use data_encoding::HEXUPPER;
    use ring::digest::{Context, Digest, SHA256};
    use tempfile::tempdir;

    /// Compute digest value for given `Reader` and print it
    /// This is taken from the Rust cookbook
    /// https://rust-lang-nursery.github.io/rust-cookbook/cryptography/hashing.html
    fn sha256_digest<R: Read>(mut reader: R) -> Result<Digest, Error> {
        let mut context = Context::new(&SHA256);
        let mut buffer = [0; 1024];

        loop {
            let count = reader.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            context.update(&buffer[..count]);
        }

        Ok(context.finish())
    }

    #[test]
    fn test_bam_single_core() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.bam");
        let out_file = out_file.to_str().unwrap();
        for l in &[
            "subset-bam",
            "-b",
            "test/test.bam",
            "-c",
            "test/barcodes.csv",
            "-o",
            out_file,
            "--cores",
            "1",
        ] {
            cmds.push(l.to_string());
        }
        _main(cmds);
        let fh = fs::File::open(&out_file).unwrap();
        let d = sha256_digest(fh).unwrap();
        let d = HEXUPPER.encode(d.as_ref());
        assert_eq!(
            d,
            "65061704E9C15BFC8FECF07D1DE527AF666E7623525262334C3FDC62F366A69E"
        );
    }
}
