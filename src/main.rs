// Copyright (c) 2020 10X Genomics, Inc.
// Modified 2025: add split-per-barcode mode and --dry-run; update to clap v4 and rust-htslib 0.46 APIs.

use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Format;
use rust_htslib::bam::Record;
use simplelog::*;
use log::{info, error, debug};
use std::cmp;
use std::collections::HashSet;
use std::fs;
use std::io::prelude::*;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;
use terminal_size::{terminal_size, Width};
use failure::Error;

fn get_args() -> Command {
    Command::new("subset-bam")
        .term_width(terminal_size().map(|(Width(w), _)| w as usize).unwrap_or(120))
        .version("1.2.0")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com>, Wyatt McDonnell <wyatt.mcdonnell@10xgenomics.com>, Modified by Ben Johnson")
        .about("Subsetting 10x Genomics BAM files")
        .arg(
            Arg::new("bam")
                .short('b')
                .long("bam")
                .value_name("FILE")
                .help("Cellranger BAM/CRAM file.")
                .required(true),
        )
        .arg(
            Arg::new("cell_barcodes")
                .short('c')
                .long("cell-barcodes")
                .value_name("FILE")
                .help("File with cell barcodes to be extracted (one per line).")
                .required(true),
        )
        .arg(
            Arg::new("out_bam")
                .short('o')
                .long("out-bam")
                .value_name("OUTPUT_FILE")
                .help("Output BAM (single-output mode)."),
        )
        .arg(
            Arg::new("split_per_barcode")
                .long("split-per-barcode")
                .help("If set, write one BAM per barcode by looping over the BAM once per barcode.")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("out_dir")
                .long("out-dir")
                .value_name("DIR")
                .help("Directory to receive <barcode>.bam files (used only with --split-per-barcode)."),
        )
        .arg(
            Arg::new("dry_run")
                .long("dry-run")
                .help("If set with --split-per-barcode, only print which BAMs would be written (no processing).")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("log_level")
                .long("log-level")
                .value_parser(["info", "debug", "error"])
                .default_value("error")
                .help("Logging level."),
        )
        .arg(
            Arg::new("cores")
                .long("cores")
                .default_value("1")
                .value_name("INTEGER")
                .help("Number of cores to use. If larger than 1, will write BAM subsets to temporary files before merging."),
        )
        .arg(
            Arg::new("bam_tag")
                .long("bam-tag")
                .default_value("CB")
                .help("Change from default value (CB) to subset alignments based on alternative tags."),
        )
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
    let cli_args: Vec<String> = std::env::args().collect();
    _main(cli_args);
}

fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let bam_file = args.get_one::<String>("bam").expect("You must provide a BAM file");
    let cell_barcodes_path = args
        .get_one::<String>("cell_barcodes")
        .expect("You must provide a cell barcodes file");

    let ll = args.get_one::<String>("log_level").unwrap().as_str();
    let cores = args
        .get_one::<String>("cores")
        .unwrap()
        .parse::<u64>()

        .expect("Failed to convert cores to integer");
    let bam_tag = args.get_one::<String>("bam_tag").unwrap().to_string();

    let split = *args.get_one::<bool>("split_per_barcode").unwrap_or(&false);
    let dry_run = *args.get_one::<bool>("dry_run").unwrap_or(&false);

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

    if split {
        // New minimal-change mode: loop over barcodes and emit one BAM per barcode.
        let out_dir = args
            .get_one::<String>("out_dir")
            .map(|s| s.as_str())
            .expect("--out-dir is required with --split-per-barcode");
        ensure_bam_and_index(bam_file);
        ensure_dir_exists(out_dir);

        let barcodes_vec = load_barcodes_vec(&cell_barcodes_path).unwrap();
        info!("Preparing to split into {} per-barcode BAMs", barcodes_vec.len());

        if dry_run {
            println!("--dry-run enabled; the following BAMs would be created in '{}':", out_dir);
            for bc in &barcodes_vec {
                println!("  {}", Path::new(out_dir).join(format!("{}.bam", bc)).display());
            }
            println!("Dry run complete; no BAMs were written.");
            return;
        }

        for bc in barcodes_vec {
            let mut set = HashSet::new();
            set.insert(bc.as_bytes().to_vec());
            let out_path = Path::new(out_dir).join(format!("{}.bam", bc));
            if out_path.exists() {
                eprintln!("Output already exists: {:?}", out_path);
                process::exit(1);
            }
            let m = subset_for_barcode_set(bam_file, &set, &out_path, cores, &bam_tag);
            info!("Wrote {:?} (kept {} reads)", out_path, m.kept_reads);
        }
        info!("Done.");
        return;
    }

    // Original single-output behavior
    let out_bam_file = args
        .get_one::<String>("out_bam")
        .map(|s| s.as_str())
        .expect("You must provide a path to write the new BAM file");

    check_inputs_exist(bam_file, cell_barcodes_path, out_bam_file);
    let cell_barcodes = load_barcodes(&cell_barcodes_path).unwrap();
    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(&bam_file, &cores).unwrap();

    let mut chunks = Vec::new();
    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        let c = ChunkArgs {
            cell_barcodes: &cell_barcodes,
            i,
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
        eprintln!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
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

// NEW: ensure BAM exists and has index (.bai/.crai)
fn ensure_bam_and_index(bam_file: &str) {
    if !Path::new(bam_file).exists() {
        eprintln!("BAM file {} does not exist", bam_file);
        process::exit(1);
    }
    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = format!("{}.bai", bam_file);
            if !Path::new(&bai).exists() {
                eprintln!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = format!("{}.crai", bam_file);
            if !Path::new(&crai).exists() {
                eprintln!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        _ => {
            eprintln!("BAM file did not end in .bam or .cram. Unable to validate");
            process::exit(1);
        }
    }
}

fn ensure_dir_exists(dir: &str) {
    let p = Path::new(dir);
    if !p.exists() {
        fs::create_dir_all(p).unwrap();
    }
    if !p.is_dir() {
        eprintln!("--out-dir must be a directory");
        process::exit(1);
    }
}

// NEW: minimal wrapper to reuse existing chunking per-barcode-set
fn subset_for_barcode_set(
    bam_file: &str,
    barcodes: &HashSet<Vec<u8>>,
    out_bam_file: &Path,
    cores: u64,
    bam_tag: &str,
) -> Metrics {
    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(&bam_file, &cores).unwrap();

    let mut chunks = Vec::new();
    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        chunks.push(ChunkArgs {
            cell_barcodes: barcodes,
            i,
            bam_file: &bam_file,
            tmp_dir: tmp_dir.path(),
            bam_tag: bam_tag.to_string(),
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
        });
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();
    let results: Vec<_> = pool.install(|| {
        chunks.par_iter().map(|chunk| slice_bam_chunk(chunk)).collect()
    });

    // combine metrics and temp outputs
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };
    let mut tmp_bams = Vec::new();
    for c in results.iter() {
        metrics.total_reads += c.metrics.total_reads;
        metrics.barcoded_reads += c.metrics.barcoded_reads;
        metrics.kept_reads += c.metrics.kept_reads;
        tmp_bams.push(&c.out_bam_file);
    }

    if cores == 1 {
        fs::copy(tmp_bams[0], out_bam_file).unwrap();
    } else {
        merge_bams(tmp_bams, out_bam_file);
    }

    metrics
}

pub fn check_inputs_exist(bam_file: &str, cell_barcodes: &str, out_bam_path: &str) {
    for path in [bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            eprintln!("File {} does not exist", path);
            process::exit(1);
        }
    }
    let path = Path::new(out_bam_path);
    if path.exists() {
        eprintln!("Output path already exists");
        process::exit(1);
    }
    if path.is_dir() {
        eprintln!("Output path is a directory");
        process::exit(1);
    }
    let _parent_dir = path.parent();
    if _parent_dir.is_none() {
        eprintln!("Unable to parse directory from {}", out_bam_path);
        process::exit(1);
    }
    let parent_dir = _parent_dir.unwrap();
    if (parent_dir.to_str().unwrap().len() > 0) & !parent_dir.exists() {
        eprintln!("Output directory {:?} does not exist", parent_dir);
        process::exit(1);
    }

    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = format!("{}.bai", bam_file);
            if !Path::new(&bai).exists() {
                eprintln!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = format!("{}.crai", bam_file);
            if !Path::new(&crai).exists() {
                eprintln!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        _ => {
            eprintln!("BAM file did not end in .bam or .cram. Unable to validate");
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
        eprintln!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    Ok(bc_set)
}

// NEW: load barcodes preserving order as Vec<String>
fn load_barcodes_vec(filename: impl AsRef<Path>) -> Result<Vec<String>, Error> {
    let r = fs::File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut v = Vec::new();
    for l in reader.lines() {
        let s = l?;
        let s = s.trim().to_string();
        if !s.is_empty() {
            v.push(s);
        }
    }
    if v.is_empty() {
        eprintln!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    Ok(v)
}

pub fn get_cell_barcode(rec: &Record, bam_tag: &str) -> Option<Vec<u8>> {
    match rec.aux(bam_tag.as_bytes()) {
        Ok(Aux::String(s)) => {
            Some(s.as_bytes().to_vec())
        }
        Ok(_) => None,
        Err(_) => None,
    }
}

pub fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr, Format::Bam)?;
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
        if let Some(barcode) = barcode {
            metrics.barcoded_reads += 1;
            if args.cell_barcodes.contains(&barcode) {
                metrics.kept_reads += 1;
                out_bam.write(&rec).unwrap();
            }
        }
    }
    let r = ChunkOuts {
        metrics,
        out_bam_file,
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
    use std::io::Read;
    use tempfile::tempdir;

    /// Compute digest value for given `Reader` and print it
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
