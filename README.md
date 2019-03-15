# subset-bam

`subset-bam` is a tool to subset a 10x Genomics BAM file based on a tag, most commonly the cell barcode tag. 

**This tool is currently under development and is not yet ready for general use.**

## Overview of how it works
`subset-bam` is a simple tool implemented in Rust that takes a 10x Genomics BAM file, a CSV file defining the subset of cells you want to isolate, and produces a new BAM file with only alignments associated with those cells.

In the subsetting process, temporary BAM files will be written to your temporary file (`$TMPDIR`) location before a final concatenation step. Please make sure this location has enough space to support this operation.

## Support
This tool is not officially supported. If you have any comments, please submit a GitHub issue.

## Installation

`subset-bam` has automatically generated downloadable binaries for generic linux and Mac OSX under the [releases page](https://github.com/10XGenomics/vartrix/releases). The linux binaries are expected to work on [our supported Operating Systems](https://support.10xgenomics.com/os-support). 

## Compiling from source
`subset-bam` is standard Rust executable project, that works with stable Rust >=1.13. Install Rust through the standard channels, then type `cargo build --release`. The executable will appear at `target/release/subset-bam`. As usual it's important to use a release build to get good performance.

## Usage

`--bam (-b)`: Input 10x Genomics BAM. This BAM must have the `CB` tag to define the barcodes of cell barcodes (or the tag defined by `--bam-tag`). Must also have an index (.bai) file. REQUIRED.

`--cell-barcodes (-c)`: A cell barcodes file as produced by Cell Ranger that defines which barcodes were called as cells. One barcode per line. In Cell Ranger runs, this can be found in the sub-folder `outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv` where `${refGenome}` is the name of the reference genome used in your Cell Ranger run. This file can be used as column labels for the output matrix. REQUIRED.

`--out-bam (-o)`: A path to write the subsetted BAM file to. REQUIRED.

`--cores`: Number of parallel cores to use. DEFAULT: 1.

`--log-level`: One of `info`, `error` or `debug`. Increasing levels of logging. DEFAULT: error.

`--bam-tag`: Change this to use an alternative tag to the default `CB` tag. This can be useful for subsetting BAMs from LongRanger.

## License
`subset-bam` is licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.
