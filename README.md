[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)


# Spl-IsoQuant manual

Spl-IsoQuant is a forked version of [IsoQuant](https://github.com/ablab/IsoQuant) developed for spatial long-read data analysis.


**Quick start:**  

*   Spl-IsoQuant can be downloaded from [https://github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant)


*   You will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.


*   Verify your installation by running:

        isoquant.py --test

*   To run IsoQuant on raw FASTQ/FASTA files use the following command

        isoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/data.fastq.gz \
        --barcode_whitelist /PATH/TO/barcodes.tsv \
        --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER


For example, using the toy data provided within this repository,

        ./isoquant.py --reference tests/toy_data/MAPT.Mouse.reference.fasta \
        --genedb tests/toy_data/MAPT.Mouse.genedb.gtf \
        --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq \
        --data_type nanopore -o toy_data_out


<a name="sec1"></a>
# About Spl-IsoQuant

Spl-IsoQuant is a tool for the genome-based analysis of spatial long RNA reads obtained with Spl-ISO-Seq protocol.
It supports both PacBio and Oxford Nanopore sequencing data. 
Spl-IsoQuant performs barcode detection, read-to-isoform assignment, quantification and PCR deduplication. 

Spl-IsoQuant is released under GPLv2 on June 12th, 2024 and can be downloaded from [https://github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant).

<a name="sec1.1"></a>
## Supported data types

Spl-IsoQuant supports long-read RNA sequencing data (PacBio CCS, ONT cDNA/dRNA) obtained with Spl-Iso-Seq protocol.
It also supports bulk long-read data.

Reads must be provided in FASTQ or FASTA format (can be gzipped). If you have already aligned your reads to the reference genome,
simply provide sorted and indexed BAM files.


<a name="sec1.2"></a>
## Supported reference data

Reference genome should be provided in multi-FASTA format (can be gzipped).
Reference genome is mandatory even when BAM files are provided.

Reference gene annotation is not mandatory, but is likely to increase precision and recall.
It can be provided in GFF/GTF format (can be gzipped).
In this case it will be converted to [gffutils](https://pythonhosted.org/gffutils/installation.html) database. Information on converted databases will be stored in your `~/.config/IsoQuant/db_config.json` to increase speed of future runs. You can also provide gffutils database manually. Make sure that chromosome/scaffold names are identical in FASTA file and gene annotation.
Note, that gffutils databases may not work correctly on NFS shares. It is possible to set a designated folder for 
the database with `--genedb_output` (different from the output directory).

Pre-constructed aligner index can also be provided to increase mapping time.

<a name="sec2"></a>
# Installation
IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.8 and higher) to be pre-installed on it.
You will also need
* [gffutils](https://pythonhosted.org/gffutils/installation.html)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [biopython](https://biopython.org/)
* [pybedtools](https://daler.github.io/pybedtools/)
* [pyfaidx](https://pypi.org/project/pyfaidx/)
* [pandas](https://pandas.pydata.org/)
* [pyyaml](https://pypi.org/project/PyYAML/)
* [editdistance](https://pypi.org/project/editdistance/)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/download/)
 

<a name="sec2.1"></a>
## Installation and requirements
To obtain Spl-IsoQuant you can download repository and install requirements.  
Clone Spl-IsoQuant repository and switch to the latest release:
```bash
git clone https://github.com/algbio/spl-IsoQuant
cd spl-IsoQuant
```
Install requirements:
```bash
pip install -r requirements.txt
```
You also need [samtools](http://www.htslib.org/download/) and [minimap2](https://github.com/lh3/minimap2) to be in the `$PATH` variable.

Alternatively, you can install all dependencies via conda by installing original IsoQuant:
```
conda create -c conda-forge -c bioconda -n isoquant python=3.8 isoquant
```

Typical installation should take no more than a few minutes.

<a name="sec2.3"></a>
## Verifying your installation
To verify IsoQuant installation type
```bash
splisoquant.py --test
```
to run on toy dataset.  
If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```

Running on the test data should not require more than a few minutes.


<a name="sec3"></a>
# Running IsoQuant
<a name="sec3.1"></a>
## IsoQuant input
To run IsoQuant, you should provide:
* Long RNA reads (PacBio or Oxford Nanopore) in one of the following formats:
  * FASTA/FASTQ (can be gzipped);
  * Sorted and indexed BAM;
* Reference sequence in FASTA format (can be gzipped);
* Reference gene annotation in gffutils database or GTF/GFF format (can be gzipped).



### Specifying input data via command line

Two main options are `--fastq` and `--bam` (see description below). Both options accept one or multiple files separated by space.
All provided files are treated as a single experiment, which means a single combined GTF will
be generated. If multiple files are provided, IsoQuant will compute tables with each column
corresponding to an individual file (per-sample counts).
To set a specific label for each sample use the `--label` option. Number of labels must be equal to the number of files.
To a set a prefix for the output files use the `--prefix` option.

This pipeline is typical for the cases when a user is
interested in comparing expression between different replicas/conditions within the same experiment.



<a name="sec3.2"></a>
## Spl-IsoQuant command line options


### Basic options
`--output` (or `-o`)
    Output folder, will be created automatically.

Note: if your output folder is located on a shared disk, use `--genedb_output` for storing
reference annotation database.

`--help` (or `-h`)
    Prints help message.

`--full_help`
    Prints all available options (including hidden ones).

`--test`
    Runs IsoQuant on the toy data set.   


### Input options
`--data_type` or `-d`
    Type of data to process, supported values are:  `pacbio_ccs` (same as `pacbio`), `nanopore` (same as `ont`)
and  `assembly` (same as `transcripts`). This option affects the algorithm parameters.

`--reference` or `-r`
    Reference genome in FASTA format (can be gzipped), required even when BAM files are provided.

`--index`
    Reference genome index for the specified aligner (`minimap2` by default),
can be provided only when raw reads are used as an input (constructed automatically if not set).

`--genedb` or `-g`
    Gene database in gffutils database format or GTF/GFF format (can be gzipped).
If you use official gene annotations we recommend to set `--complete_genedb` option.

`--complete_genedb`
    Set this flag if gene annotation contains transcript and gene meta-features.
Use this flag when providing official annotations, e.g. GENCODE.
This option will set `disable_infer_transcripts` and `disable_infer_genes` gffutils options,
which dramatically speeds up gene database conversion (see more [here](https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html)).

#### Providing input reads via command line option:

`--fastq`
    Input FASTQ/FASTA file(s), can be gzipped;  a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](#sec3.1).


`--bam`
    Sorted and indexed BAM file(s); a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](#sec3.1).


#### Other input options:
`--stranded`
    Reads strandness type, supported values are: `forward`, `reverse`, `none`.

`--fl_data`
    Input sequences represent full-length transcripts; both ends of the sequence are considered to be reliable.

`--prefix` or `-p`
    Prefix for all output files and sub-folder name. `OUT` if not set.

`--labels` or `-l`
    Sets space-separated sample names. Make sure that the number of labels is equal to the number of files.
Input file names are used as labels if not set.



### Pipeline options

`--resume`
    Resume a previously unfinished run. Output folder with previous run must be specified.
    Allowed options are `--threads` and `--debug`, other options cannot be changed.
    IsoQuant will run from the beginning if the output folder does not contain the previous run.

`--force`
    force to overwrite the folder with previous run.

`--threads` or `-t`
    Number of threads to use, 16 by default.

`--clean_start`
    Do not use previously generated gene database, genome indices or BAM files, run pipeline from the very beginning (will take more time).





### Hidden options
<a name="hidden"></a>
Options below are shown only with `--full_help` option.
We recommend _not_ to modify these options unless you are clearly aware of their effect.

`--no_gzip`
    Do not compress large output files.

`--no_gtf_check`
    Do not perform input GTF checks.

`--no_secondary`
    Ignore secondary alignments.

`--aligner`
    Force to use this alignment method, can be `starlong` or `minimap2`; `minimap2` is currently used as default. Make sure the specified aligner is in the `$PATH` variable.

`--no_junc_bed`
    Do not use gene annotation for read mapping.

`--junc_bed_file`
    Annotation in BED12 format produced by `paftools.js gff2bed` (can be found in `minimap2`), will be created automatically if not given.

`--delta`
    Delta for inexact splice junction comparison, chosen automatically based on data type (e.g. 4bp for PacBio, 6pb for ONT).

`--genedb_output`
    If your output folder is located on a shared storage (e.g. NFS share), use this option to set another path
    for storing the annotation database, because SQLite database cannot be created on a shared disks.
    The folder will be created automatically.

`--high_memory`
    Cache read alignments instead for making several passes over a BAM file, noticeably increases RAM usage, 
but may improve running time when disk I/O is relatively slow.


### Examples
<a name="examples"></a>

* Mapped PacBio CCS reads in BAM format; pre-converted gene annotation:

```bash
splisoquant.py -d pacbio_ccs --bam mapped_reads.bam \
 --genedb annotation.db --output output_dir
```

* Nanopore dRNA stranded reads; official annotation in GTF format, use custon prefix for output:
```bash
splisoquant.py -d nanopore --stranded forward --fastq ONT.raw.fastq.gz \
 --reference reference.fasta --genedb annotation.gtf --complete_genedb \
 --output output_dir --prefix My_ONT
```

* Nanopore cDNA reads; no reference annotation:
```bash
splisoquant.py -d nanopore --fastq ONT.cDNA.raw.fastq.gz \
 --reference reference.fasta --output output_dir --prefix My_ONT_cDNA
```

* PacBio FL reads; custom annotation in GTF format, which contains only exon features:
```bash
splisoquant.py -d pacbio_ccs --fl_data --fastq CCS.fastq \
 --reference reference.fasta --genedb genes.gtf --output output_dir
```

* Nanopore cDNA reads, multiple samples/replicas within a single experiment; official annotation in GTF format:
```bash
splisoquant.py -d nanopore --bam ONT.cDNA_1.bam ONT.cDNA_2.bam ONT.cDNA_3.bam \
 --reference reference.fasta --genedb annotation.gtf --complete_genedb --output output_dir
 --predix ONT_3samples --labels A1 A2 A3
```

* ONT cDNA reads; 2 experiments with 3 replicates; official annotation in GTF format:
```bash
splisoquant.py -d nanopore --yaml dataset.yaml  \
 --complete_genedb --genedb genes.gtf \
 --reference reference.fasta --output output_dir
```

dataset.yaml file :

```
[
  data format: "fastq",
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE1/file1.fastq",
      "/PATH/TO/SAMPLE1/file2.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate2",
      "Replicate3"
    ]
  },
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE2/file1.fastq",
      "/PATH/TO/SAMPLE2/file2.fastq",
      "/PATH/TO/SAMPLE2/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate2",
      "Replicate3"
    ]
  }
]

```


IsoQuant will produce 2 sets of resulting files (including annotations and expression tables), one for each experiment.
Output sub-folder will be named `Experiment1` and `Experiment2`.
Expression tables will have columns "Replicate1", "Replicate2" and "Replicate3".


* ONT cDNA reads; 1 experiment with 2 replicates, each replicate has 2 files; official annotation in GTF format:
```bash
splisoquant.py -d nanopore --yaml dataset.yaml  \
  --complete_genedb --genedb genes.gtf \
 --reference reference.fasta --prefix MY_SAMPLE \
 --output output_dir  
```

dataset.yaml file :


```
[
  data format: "fastq",
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE1/file1.fastq",
      "/PATH/TO/SAMPLE1/file2.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate1",
      "Replicate2",
      "Replicate2"
    ]
  }
]

```


IsoQuant will produce one output sub-folder `Experiment1`.
Expression tables will have columns "Replicate1" and "Replicate2".
Files having identical labels will be treated as a single replica (and thus the counts will be combined).


<a name="sec3.3"></a>
## IsoQuant output

### Output files

IsoQuant output files will be stored in `<output_dir>`, which is set by the user.
If the output directory was not specified the files are stored in `isoquant_output`.

IsoQuant consists of two stages, which generate its own output:
1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isofrom assignment,
splice site correction and abundance quantification for reference genes/transcripts.
2. Transcript discovery. Reconstructs transcript models and performs abundance quantification for discovered isoforms.

#### Reference-based analysis output

_Will be produced only if a reference gene annotation is provided._

* `SAMPLE_ID.read_assignments.tsv.gz` - TSV file with read to isoform assignments (gzipped by default);
* `SAMPLE_ID.corrected_reads.bed.gz` - BED file with corrected read alignments (gzipped by default);
* `SAMPLE_ID.transcript_tpm.tsv` - TSV file with reference transcript expression in TPM;
* `SAMPLE_ID.transcript_counts.tsv` - TSV file with raw read counts for reference transcript;
* `SAMPLE_ID.gene_tpm.tsv` - TSV file with reference gene expression in TPM;
* `SAMPLE_ID.gene_counts.tsv` - TSV file with raw read counts for reference genes;

If `--sqanti_output` is set, IsoQuant will produce output in [SQANTI](https://github.com/ConesaLab/SQANTI3)-like format:
* `SAMPLE_ID.novel_vs_known.SQANTI-like.tsv` - discovered novel transcripts vs reference transcripts (similar, but not identical to SQANTI `classification.txt`);

If `--count_exons` is set, exon and intron counts will be produced:
* `SAMPLE_ID.exon_counts.tsv` - reference exon inclusion/exclusion read counts;
* `SAMPLE_ID.intron_counts.tsv` - reference intron inclusion/exclusion read counts;

If `--read_group` is set, the per-group expression values for reference features will be also computed:
* `SAMPLE_ID.gene_grouped_tpm.tsv`
* `SAMPLE_ID.transcript_grouped_tpm.tsv`
* `SAMPLE_ID.gene_grouped_counts.tsv`
* `SAMPLE_ID.transcript_grouped_counts.tsv`
* `SAMPLE_ID.exon_grouped_counts.tsv`
* `SAMPLE_ID.intron_grouped_counts.tsv`

#### Transcript discovery output

_Will not be produced if `--no_model_construction` is set._

File names typically contain `transcript_model` in their name.

* `SAMPLE_ID.transcript_models.gtf` - GTF file with discovered expressed transcript (both known and novel transcripts);
* `SAMPLE_ID.transcript_model_reads.tsv.gz` - TSV file indicating which reads contributed to transcript models (gzipped by default);
* `SAMPLE_ID.transcript_model_tpm.tsv` - expression of discovered transcripts models in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.transcript_model_counts.tsv` - raw read counts for discovered transcript models (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.extended_annotation.gtf` - GTF file with the entire reference annotation plus all discovered novel transcripts;


If `--read_group` is set, the per-group counts for discovered transcripts will be also computed:
* `SAMPLE_ID.transcript_model_grouped_counts.tsv`
* `SAMPLE_ID.transcript_model_grouped_tpm.tsv`


If multiple experiments are provided, aggregated expression matrices will be placed in `<output_dir>`:
* `combined_gene_counts.tsv`
* `combined_gene_tpm.tsv`
* `combined_transcript_counts.tsv`
* `combined_transcript_tpm.tsv`

Additionally, a log file will be saved to the directory.  
* <output_dir>/isoquant.log   

If raw reads were provided, BAM file(s) will be stored in `<output_dir>/<SAMPLE_ID>/aux/`.  
In case `--keep_tmp` option was specified this directory will also contain temporary files.

### Output file formats

Although most output files include headers that describe the data, a brief explanation of the output files is provided below.

#### Read to isoform assignment

Tab-separated values, the columns are:

* `read_id` - read id;
* `chr` - chromosome id;
* `strand` - strand of the assigned isoform (not to be confused with read mapping strand);
* `isoform_id` - isoform id to which the read was assigned;
* `gene_id` - gene id to which the read was assigned;
* `assignment_type` - assignment type, can be:
    - `unique` - reads was unambiguously assigned to a single known isoform;
    - `unique_minor_difference` - read was assigned uniquely but has alignment artifacts;
    - `inconsistent` - read was matched with inconsistencies, closest match(es) are reported;
    - `ambiguous` - read was assigned to multiple isoforms equally well;
    - `noninfomative` - reads is intronic/intergenic.
* `assignment_events` - list of detected inconsistencies; for each assigned isoform a list of detected inconsistencies relative to the respective isoform is stored; values in each list are separated by `+` symbol, lists are separated by comma, the number of lists equals to the number of assigned isoforms; possible events are (see graphical representation below):
    - consistent events:
        - `none` / `.` / `undefined` - no special event detected;
        - `mono_exon_match` mono-exonic read matched to mono-exonic transcript;
        - `fsm` - full splice match;
        - `ism_5/3` - incomplete splice match, truncated on 5'/3' side;
        - `ism_internal` - incomplete splice match, truncated on both sides;
        - `mono_exonic` - mono-exonic read matching spliced isoform;
        - `tss_match` / `tss_match_precise` - 5' read is located less than 50 / `delta` bases from the TSS of the assigned isoform
        - `tes_match` / `tes_match_precise` - 3' read is located less than 50 / `delta` bases from the TES of the assigned isoform (can be reported without detecting polyA sites)
    - alignment artifacts:
        - `intron_shift` - intron that seems to be shifted due to misalignment (typical for Nanopores);
        - `exon_misalignment` - short exon that seems to be missed due to misalignment  (typical for Nanopores);
        - `fake_terminal_exon_5/3` - short terminal exon at 5'/3' end that looks like an alignment artifact (typical for Nanopores);  
        - `terminal_exon_misalignment_5/3` - missed reference short terminal exon;
        - `exon_elongation_5/3` - minor exon extension at 5'/3' end (not exceeding 30bp);
        - `fake_micro_intron_retention` - short annotated introns are often missed by the aligners and thus are not considered as intron retention;
    - intron retentions:
        - `intron_retention` - intron retention;
        - `unspliced_intron_retention`  - intron retention by mono-exonic read;
        - `incomplete_intron_retention_5/3` - terminal exon at 5'/3' end partially covers adjacent intron;
    - significant inconsistencies (each type end with `_known` if _all_ resulting read introns are annotated and `_novel` otherwise):
        - `major_exon_elongation_5/3` - significant exon extension at 5'/3' end (exceeding 30bp);
        - `extra_intron_5/3` - additional intron on the 5'/3' end of the isoform;
        - `extra_intron` - read contains additional intron in the middle of exon;
        - `alt_donor_site` - read contains alternative donor site;
        - `alt_acceptor_site` - read contains alternative annotated acceptor site;
        - `intron_migration` - read contains alternative annotated intron of approximately the same length as in the isoform;
        - `intron_alternation` - read contains alternative intron, which doesn't fall intro any of the categories above;
        - `mutually_exclusive_exons` - read contains different exon(s) of the same total length comparing to the isoform;
        - `exon_skipping` - read skips exon(s) comparing to the isoform;
        - `exon_merge` - read skips exon(s) comparing to the isoform, but a sequence of a similar length is attached to a neighboring exon;
        - `exon_gain` - read contains additional exon(s) comparing to the isoform;
        - `exon_detach` - read contains additional exon(s) comparing to the isoform, but a neighboring exon looses a sequnce of a similar length;
        - `terminal_exon_shift` - read has alternative terminal exon;   
        - `alternative_structure` - reads has different intron chain that does not fall into any of categories above;
    - alternative transcription start / end (reported when poly-A tails are present):
        - `alternative_polya_site` - read has alternative polyadenylation site;
        - `internal_polya_site` - poly-A tail detected but seems to be originated from A-rich intronic region;
        - `correct_polya_site` - poly-A site matches reference transcript end;
        - `aligned_polya_tail` - poly-A tail aligns to the reference;  
        - `alternative_tss` - alternative transcription start site.
* `exons` - list of coordinates for normalized read exons (1-based, indels and polyA exons are excluded);
* `additional` - field for supplementary information, which may include:
    - `PolyA` - True if poly-A tail is detected;
    - `Canonical` - True if all read introns are canonical, Unspliced is used for mono-exon reads; (use `--check_canonical`)

Note, that a single read may occur more than once if assigned ambiguously.

#### Expression table format

Tab-separated values, the columns are:

* `feature_id` - genomic feature ID;
* `TPM` or `count` - expression value (float).

For grouped counts, each column contains expression values of a respective group.
In the number of groups exceeds 10, file will contain 3 columns:

* `feature_id` - genomic feature ID;
* `group_id` - name of the assigned group;
* `TPM` or `count` - expression value (float).



<a name="sec4"></a>
## Citation
The manuscript is currently under preparation.


<a name="sec5"></a>
## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve IsoQuant. If you have any troubles running IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/algbio/spl-IsoQuant/).


