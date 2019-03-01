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

    let bam = bam::IndexedReader::from_path(bam_file).unwrap();
    let mut out_bam = load_writer(&bam, &out_bam_file).unwrap();
    out_bam.set_threads(2).unwrap();
    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    use rust_htslib::bam::Read; // collides with fs::Read
    bam.set_threads(2).unwrap();
    for r in bam.iter_chunk(None, None) {
        let rec = r.unwrap();
        let barcode = get_cell_barcode(&rec, &bam_tag);
        if barcode.is_some() {
            let barcode = barcode.unwrap();
            if cell_barcodes.contains(&barcode) {
                out_bam.write(&rec).unwrap();
            }
        }
    }
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

pub fn get_cell_barcode(rec: &Record, bam_tag: &str) -> Option<Vec<u8>> {
    //println!("{:?}", rec.aux(bam_tag.as_bytes()));
    match rec.aux(bam_tag.as_bytes()) {
        Some(Aux::String(hp)) => {
            let cb = hp.to_vec();
            Some(cb)
        },
        _ => None,
    }
}

pub fn load_writer(bam: &bam::IndexedReader, out_bam_path: &str) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr)?;
    Ok(out_handle)
}
