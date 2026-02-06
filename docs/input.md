# Spl-IsoQuant input

To run Spl-IsoQuant on single-cell or spatial, you should provide:

* Reference genome in FASTA format (can be gzipped);
* Reference gene annotation in gffutils database or GTF/GFF format (can be gzipped);
* Long RNA reads (PacBio or Oxford Nanopore) in one of the following formats:
  * FASTA/FASTQ (can be gzipped);
  * Unmapped BAM files (typical for PacBio CCS reads);
  * Sorted and indexed BAM;
* Barcode information in either of the following formats:
  * A barcode whitelist file in TSV format (can be gzipped);
  * A TSV file with mapping between reads and barcodes (can be gzipped);
  * An MDF file with custom molecule description;
* Optionally, you can provide a TSV file with mapping between barcodes and cell types / spatial spots.

Spl-IsoQuant is also capable of using short Illumina reads to correct long-read alignments.

The ways of providing input files are described below.

## Specifying input reads via command line

The main options are `--fastq`, `--unmapped_bam` and `--bam` (see description below). 
These options accept one or multiple files separated by space.
All provided files are treated as a single experiment, which means a single combined GTF will
be generated. If multiple files are provided, Spl-IsoQuant will compute tables with each column
corresponding to an individual file (per-sample counts).
To set a specific label for each sample use the `--label` option. Number of labels must be equal to the number of files.
To a set a prefix for the output files use the `--prefix` option.

This pipeline is typical for the cases when a user is
interested in comparing expression between different replicas/conditions within the same experiment.

### Barcodes and single-cell/spatial data

To enable single-cell or spatial analysis, you should:
* Specify mode via `--mode` option;
* Provide barcode information using one of the following options:
  * A barcode whitelist file via `--barcode_whitelist`;
  * A TSV able with mapping between reads and barcodes via `--barcoded_reads`;
  * An MDF file with custom molecule description via `--molecule`;
* Additionally, you can provide a table with mapping between barcodes and cell types / spatial spots via `--barcode2spot`.

See details [here](single_cell.md) and [here](cmd.md#single-cell-and-spatial-transcriptomics-options).

### Short reads for alignment correction

A BAM file with Illumina reads can be provided via `--illumina_bam`. It cannot be the only input, but may only be used with either `--bam` or `--fastq`.
The option accepts one or multiple bam files separated by space. All files will be combined and used to correct offsets between introns in long and short reads as well as skipped exons.


## Specifying input reads via yaml file

To provide all input files in a single description file, you can use a [YAML](https://www.redhat.com/en/topics/automation/what-is-yaml) file via `--yaml` (see description below).
You can provide multiple experiments in a single YAML file with each experiment containing an arbitrary number of smaples/replicas.
A distinct output folder with individual GTFs and abundance tables will be generated for each experiment.
In this option, BAM files with short reads for correction can be provided for each experiment.

The YAML file contains a list of experiments (e.g. in square brackets).
The first entry in the list should be the type of files the experiments contain, written as `data format: `
followed by the type in quotation marks. The type can be either `fastq`, `bam` or `unmapped_bam`.

Each experiment is represented as set of parameters (e.g. in curly brackets).
Each experiment must have a name and a list of long-read files in the specified format.
Additionally, it may contain one or multiple BAM files with short reads.
The name is provided as `name: ` followed by the experiment name in quotation marks.
Both short and long read files are provided as a list of file paths in quotation marks,
following `long read files: ` and `illumina bam: ` respectively.
Labels for the files can also be set with `labels: `.
The number of labels needs to be the same as the number of files with long reads.
All paths should be either absolute or relative to the YAML file.

For example:

```
[
  data format: "fastq",
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/FILE1.fastq",
      "/PATH/TO/FILE2.fastq"
    ],
    labels: [
      "Sample1",
      "Sample2"
    ],
    illumina bam: ["PATH/TO/ILLUMINA1.bam"]
  },
  {
    name: "Experiment2",
    long read files: [
      "/PATH/TO/FILE3.fastq"
    ],
    illumina bam: ["PATH/TO/ILLUMINA2.bam"]
  }
]

```


Output sub-folders will be named `Experiment1` and `Experiment2`.
Both sub-folders will contain predicted transcript models and abundance tables.
Abundance table for `Experiment2` with have columns "Sample1" and "Sample2".

Note, that  `--bam`, `--unmapped_bam`, `--fastq` and `--label` options are not compatible with `--yaml`.
See more in [examples](examples.md).

