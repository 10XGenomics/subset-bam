extern crate rust_htslib;
extern crate clap;
extern crate csv;
extern crate rayon;
extern crate terminal_size;
extern crate tempfile;
extern crate simplelog;
#[macro_use]
extern crate log;
#[macro_use]
extern crate failure;
#[macro_use]
extern crate human_panic;

use std::io;
use std::process;
use std::cmp;
use std::fs::{self, File};
use std::path::Path;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::collections::HashSet;
use clap::{Arg, App};
use terminal_size::{Width, terminal_size};
use simplelog::*;
use failure::Error;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use rust_htslib::bam;
use rayon::prelude::*;

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("subset-bam")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("0")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com>")
        .about("Subsetting 10x Genomics BAM files")
        .arg(Arg::with_name("bam")
             .short("b")
             .long("bam")
             .value_name("FILE")
             .help("Cellranger BAM file")
             .required(true))
        .arg(Arg::with_name("cell_barcodes")
             .short("c")
             .long("cell-barcodes")
             .value_name("FILE")
             .help("File with cell barcodes to be extracted")
             .required(true))
        .arg(Arg::with_name("out_bam")
             .short("o")
             .long("out-bam")
             .value_name("OUTPUT_FILE")
             .help("Output BAM")
             .required(true))
        .arg(Arg::with_name("log_level")
             .long("log-level")
             .possible_values(&["info", "debug", "error"])
             .default_value("error")
             .help("Logging level"))
        .arg(Arg::with_name("threads")
             .long("threads")
             .short("t")
             .value_name("INTEGER")
             .default_value("1")
             .help("Number of parallel threads to use"))
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

struct BamHolder {
    bam: bam::Reader
}

impl BamHolder {
    fn new(bam_path: &str) -> BamHolder {
        let bam = bam::Reader::from_path(bam_path).unwrap();
        BamHolder {bam: bam}
    }
}

struct ReadFilter<'a> {
    cell_barcodes: &'a HashSet<Vec<u8>>,
    bam_iter: rust_htslib::bam::ChunkIterator<'a, rust_htslib::bam::Reader>,
    bam_tag: &'a String,
}

impl<'a> ReadFilter<'a> {
    fn new(offset: (Option<i64>, Option<i64>), 
           cell_barcodes: &'a HashSet<Vec<u8>>,
           bam_holder: &'a mut BamHolder,
           bam_tag: &'a String) -> ReadFilter<'a> {
        let bam_iter = bam_holder.bam.iter_chunk(offset.0, offset.1);
        ReadFilter { cell_barcodes: cell_barcodes,
                     bam_iter: bam_iter,
                     bam_tag: bam_tag}
    }
}

impl<'a> Iterator for ReadFilter<'a> {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let read = self.bam_iter.next();
            // catch if iteration is done
            let read = match read {
                None => return None,
                Some(_) => read.unwrap()
            };
            let read = read.unwrap();
            let barcode = get_cell_barcode(&read, self.bam_tag);
            if barcode.is_some() {
                let barcode = barcode.unwrap();
                if self.cell_barcodes.contains(&barcode) {
                    return Some(read)
                }
            }
        }
    }
}


fn main() {
    setup_panic!();  // pretty panics for users
    let mut cli_args = Vec::new();
    for arg in std::env::args_os() {
        cli_args.push(arg.into_string().unwrap());
    }
    _main(cli_args);
}

fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args.value_of("cell_barcodes").expect("You must provide a cell barcodes file");
    let out_bam_file = args.value_of("out_bam").expect("You must provide a path to write the new BAM file");
    let threads = args.value_of("threads").unwrap_or_default()
                                          .parse::<usize>()
                                          .expect("Failed to convert threads to integer");
    let ll = args.value_of("log_level").unwrap();
    let bam_tag = args.value_of("bam_tag").unwrap_or_default().to_string();
    check_inputs_exist(bam_file, cell_barcodes, out_bam_file);
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => { println!("Log level not valid"); process::exit(1); }
    };

    let _ = SimpleLogger::init(ll, Config::default());

    let virtual_offsets = bgzf_noffsets(&bam_file, &(threads as u64)).unwrap();
    let mut out_fh = load_writer(&bam_file, &out_bam_file).unwrap();

    for vo in virtual_offsets.iter() {
        let mut bh = BamHolder::new(bam_file);
        let i = ReadFilter::new(*vo, &cell_barcodes, &mut bh, &bam_tag);
        for rec in i {
            out_fh.write(&rec).unwrap();
        }
    }

    //let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    //debug!("Initialized a thread pool with {} threads", threads);

    //let results: Vec<_> = pool.install(|| virtual_offsets.iter().par_bridge()
    //                                .map_with(out_fh, |out_fh, offsets| { 
    //                                    extract_reads(offsets, out_fh, &cell_barcodes).unwrap() 
    //                                    } )
    //                                .collect());
}

pub fn check_inputs_exist(bam_file: &str, cell_barcodes: &str, out_bam_path: &str) {
    for path in [bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            error!("File {} does not exist", path);
            process::exit(1);
        }
    }

    let _dir = Path::new(out_bam_path).parent();
    if _dir.is_none() {
        error!("Unable to parse directory from {}", out_bam_path);
        process::exit(1);
    }
    let dir = _dir.unwrap();
    if (dir.to_str().unwrap().len() > 0) & !dir.exists() {
        error!("Output directory {:?} does not exist", dir);
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

pub fn get_cell_barcode(rec: &Record, bam_tag: &String) -> Option<Vec<u8>> {
    match rec.aux(bam_tag.as_bytes()) {
        Some(Aux::String(hp)) => {
            let cb = hp.to_vec();
            Some(cb)
        },
        _ => None,
    }
}

pub fn bgzf_noffsets(bam_path: &str, num_chunks: &u64) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);

        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    // if we only have one thread, this is easy
    if *num_chunks == 1 as u64 {
        let final_offsets = vec![(None, None)];
        return Ok(final_offsets)
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
    let mut fp = File::open(bam_path)?;
    for offset in initial_offsets {
        fp.seek(io::SeekFrom::Start(offset))?;
        let mut buffer = [0; 2 << 16];
        fp.read(&mut buffer)?;
        for i in 0..num_bytes {
            if is_valid_bgzf_block(&buffer[i as usize..]) {
                adjusted_offsets.push(offset + i);
                break
            }
        }
    }
    // bit-shift and produce start/stop intervals
    let mut final_offsets = Vec::new();
    final_offsets.push((None, 
                        Some(((adjusted_offsets[1]) as i64) << 16)));
    for n in 2..num_chunks-1 {
        let n = n as usize;
        final_offsets.push((Some((adjusted_offsets[n-1] as i64) << 16), 
                            Some((adjusted_offsets[n] as i64) << 16)));
    }
    final_offsets.push((Some(((adjusted_offsets[adjusted_offsets.len() - 1]) as i64) << 16),
                        None));
    Ok(final_offsets)
}


pub fn is_valid_bgzf_block(block: &[u8]) -> bool {
    // look for the bgzip magic characters \x1f\x8b\x08\x04
    // TODO: is this sufficient?
    if block.len() < 18 {
        return false
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false
    }
    true
}

pub fn load_writer(bam_path: &str, out_bam_path: &str) -> Result<bam::Writer, Error> {
    let bam = bam::IndexedReader::from_path(bam_path)?;
    use rust_htslib::bam::Read;
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr)?;
    Ok(out_handle)
}
