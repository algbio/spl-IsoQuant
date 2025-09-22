[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)


# Spl-IsoQuant manual

Spl-IsoQuant is a forked version of [IsoQuant](https://github.com/ablab/IsoQuant) developed for upstream analysis of 
single-cell and spatial long-read data analysis.

Current version supports the following protocols:
- 10x 3' v3;
- 10x Visium 3' and Visium HD;
- Spl-Iso-Seq (Curio);
- StereoSeq.


**Quick start:**  

*   Spl-IsoQuant can be downloaded from [https://github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant)


*   You will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.


*   Verify your installation by running:

        splisoquant.py --test

*   To run Spl-IsoQuant on raw FASTQ/FASTA files use the following command

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/data.fastq.gz \
        --barcode_whitelist /PATH/TO/barcodes.tsv \
        --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER


For example, using the toy data provided within this repository

        splisoquant.py --threads 2 --mode stereoseq \
        --fastq', tests/stereo/S1.4K.subsample.fq.gz \
        --barcode_whitelist tests/stereo/barcodes.tsv \
        --reference tests/stereo/GRCm39.chrX.fa.gz \
        --genedb tests/stereo/gencode.chrX.ENSMUSG00000031153.gtf \
        --clean_start --data_type nanopore --complete_genedb \
        -p TEST_DATA --output splisoquant2_test 


<a name="sec1"></a>
# About Spl-IsoQuant

Spl-IsoQuant is a tool for the genome-based analysis of spatial long RNA reads obtained with Spl-ISO-Seq protocol.
It supports both PacBio and Oxford Nanopore sequencing data. 
Spl-IsoQuant performs barcode detection, read-to-isoform assignment, quantification and PCR deduplication. 

Similarly to IsoQuant, it can also discover novel isoforms. 
However,  read-to-isoform assignment, quantification and PCR deduplication are perfromed using reference annotaton. 
Currently, to incorporate novel gene and isoforms into your anaysis, you will need to run Spl-IsoQuant twice.

Spl-IsoQuant is released under GPLv2 on June 12th, 2024 and can be downloaded from [https://github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant).

<a name="sec1.1"></a>
## Supported data types

Spl-IsoQuant supports long-read RNA sequencing data (PacBio CCS, ONT cDNA/dRNA) obtained with one of the following protocols:
- 10x 3' v3;
- 10x Visium 3' and Visium HD;
- Spl-Iso-Seq (Curio);
- StereoSeq.

It also supports bulk long-read data. When bulk data is provided, the behaviour will identical to the original [IsoQuant](https://github.com/ablab/IsoQuant).

Reads must be provided in FASTQ or FASTA format (can be gzipped). If you have already aligned your reads to the reference genome,
simply provide sorted and indexed BAM files.

Spl-IsoQuant can also take aligned Illumina reads to correct long-read spliced alignments. However, short reads are _not_
used to discover transcript models or compute abundances.

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
Spl-IsoQuant requires a 64-bit Linux system or MacOS and Python (3.8 and higher) to be pre-installed on it.
There are no specific hardware requirements.

You will also need
* [gffutils](https://pythonhosted.org/gffutils/installation.html)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [biopython](https://biopython.org/)
* [pyfaidx](https://pypi.org/project/pyfaidx/)
* [pandas](https://pandas.pydata.org/)
* [pyyaml](https://pypi.org/project/PyYAML/)
* [editdistance](https://pypi.org/project/editdistance/)
* [ssw-py](https://pypi.org/project/ssw-py/)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/download/)
 

<a name="sec2.1"></a>
## Installation and requirements
To obtain Spl-IsoQuant you can download the repository and install requirements.  
Clone Spl-IsoQuant repository and switch to the latest release:
```bash
git clone https://github.com/algbio/spl-IsoQuant
cd spl-IsoQuant
```
Install requirements:
```bash
pip install -r requirements.txt
```
You also need [samtools 1.14+](http://www.htslib.org/download/) and [minimap2 2.18+](https://github.com/lh3/minimap2) to be in the `$PATH` variable.

Alternatively, you can install most of the dependencies via conda by installing original IsoQuant:

`conda create -c conda-forge -c bioconda -n isoquant python=3.8 isoquant`

Spl-IsoQuant was tested with version given in the `requirements.txt`.

Spl-IsoQuant will be available on bioconda soon.


<a name="sec2.3"></a>
## Verifying your installation
To verify Spl-IsoQuant installation type
```bash
splisoquant.py --test
```
to run on toy dataset.  
If the installation is successful, you will find the following information at the end of the log:
```bash
=== Spl-IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```

Running on the test data requires no more than a few minutes.


<a name="sec3"></a>
# Running Spl-IsoQuant
<a name="sec3.1"></a>
## Spl-IsoQuant input
To run Spl-IsoQuant, you should provide:
* Long RNA reads (PacBio or Oxford Nanopore) in one of the following formats:
  * FASTA/FASTQ (can be gzipped);
  * Sorted and indexed BAM;
* Reference sequence in FASTA format (can be gzipped);
* Reference gene annotation in gffutils database or GTF/GFF format (can be gzipped);
* Whitelist barcodes or a TSV with read ids and their respective barcodes and UMIs.



### Specifying input data via command line

Two main options are `--fastq` and `--bam` (see description below). Both options accept one or multiple files separated by space.
To a set a prefix for the output files use the `--prefix` option.


<a name="sec3.2"></a>
## Spl-IsoQuant command line options


### Basic options
`--output` (or `-o`)
    Output folder, will be created automatically.

Note: if your output folder is located on a shared disk, use `--genedb_output` for storing
reference annotation database.

`--help` (or `-h`)
    Prints help message.

`--test`
    Runs Spl-IsoQuant on the toy data set.   


### Input options
`--data_type` or `-d`
    Type of data to process, supported values are:  `pacbio_ccs` (same as `pacbio`), `nanopore` (same as `ont`)
and  `assembly` (same as `transcripts`). This option affects the algorithm parameters.

`--mode`
    Spl-IsoQuant mode to be used, can be one of: `bulk`, `tenX_v3`, `curio`, `stereoseq_nosplit`, `stereoseq`, `visium_hd` , `visium_5prime`.

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

#### Providing input data:

`--fastq`
    Input FASTQ/FASTA file(s), can be gzipped;  a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](#sec3.1).


`--bam`
    Sorted and indexed BAM file(s); a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](#sec3.1).


`--barcode_whitelist`
    A file with a list of possible barcodes (one per line).

`--barcoded_reads`
    A TSV file where each line contains 5 columns:
    read id, barcode, UMI sequence, barcode score and a flag indicating whether UMI is reliable (True/False).
    If you use a custom barcode calling software, barcode score can be set to the length of the barcode
    and the last column can be set to `True`.


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
    Spl-IsoQuant will run from the beginning if the output folder does not contain the previous run.

`--force`
    force to overwrite the folder with previous run.

`--threads` or `-t`
    Number of threads to use, 16 by default.

`--clean_start`
    Do not use previously generated gene database, genome indices or BAM files, run pipeline from the very beginning (will take more time).



### Hidden options
<a name="hidden"></a>
We recommend _not_ to modify these options unless you are clearly aware of their effect.

`--no_gzip`
    Do not compress large output files.

`--no_gtf_check`
    Do not perform input GTF checks.

`--no_secondary`
    Ignore secondary alignments.

`--genedb_output`
    If your output folder is located on a shared storage (e.g. NFS share), use this option to set another path
    for storing the annotation database, because SQLite database cannot be created on a shared disks.
    The folder will be created automatically.

`--high_memory`
    Cache read alignments instead for making several passes over a BAM file, noticeably increases RAM usage, 
but may improve running time when disk I/O is relatively slow.


### Example

To run Spl-IsoQuant on a simple Spl-Iso-Seq dataset provided within the package, run

        splisoquant.py --reference tests/splisoseq/GRCh38.chrX.fa.gz \
        --genedb tests/splisoseq/ref.gtf --complete_genedb \
        --fastq tests/splisoseq/ONT.fasta.gz --barcode_whitelist tests/splisoseq/barcodes.tsv.gz \
        -d ont -p TEST -o splisoseq_test

Typically, this test takes no more than a few minutes.

As a result, in the `splisoseq_test/TEST` folder you will get the following files:

Reference-based analysis:
* `TEST_DATA.split_reads_0.fa` - extracted cDNAs;
* `TEST_DATA.barcoded_reads_0.tsv` - detected barcodes;
* `TEST_DATA.UMI_filtered.ED4.allinfo` - PCR deduplicated reads with barcodes and gene assignments;
* `TEST_DATA.UMI_filtered.ED4.stats.tsv` - brief stats for PCR deduplication and read assignments;
* `TEST_DATA.read_assignments.tsv.gz` - read-to-isoform assignments (after PCR deduplication); 
* `TEST_DATA.corrected_reads.bed.gz` - read alignments after splice site correction;
* `TEST_DATA.gene_counts.tsv, TEST_DATA.gene_tpm.tsv` - bulk gene counts/TPMs after PCR deduplication;
* `TEST_DATA.transcript_counts.tsv, TEST_DATA.transcript_tpm.tsv` - bulk gene counts/TPMs after PCR deduplication;

Discovered novel transcripts and genes (reads after PCR deduplication are used):
* `TEST_DATA.extended_annotation.gtf` - extended gene annotation:original input annotation plus novel genes and transcripts;
* `TEST_DATA.transcript_models.gtf` - expressed known and novel transcripts;* 
* `TEST_DATA.discovered_gene_counts.tsv`, `TEST_DATA.discovered_gene_tpm.tsv` - bulk gene counts/TPMs including known and novel genes;
* `TEST_DATA.discovered_transcript_counts.tsv`, `TEST_DATA.discovered_transcript_tpm.tsv` - bulk gene counts/TPMs including known and novel transcripts;
* `TEST_DATA.transcript_model_reads.tsv.gz` - reads assigned to novel and known transcripts.


<a name="sec3.3"></a>
## Spl-IsoQuant output

### Output files

Spl-IsoQuant output files will be stored in `<output_dir>`, which is set by the user.
If the output directory was not specified the files are stored in `isoquant_output`.

Output files are:

* `SAMPLE_ID.split_reads_#.fa` - FASTA file with cDNAs extracted from raw reads, one file per every input file (only when `--mode stereo` is set);.
* `SAMPLE_ID.barcoded_reads_#.tsv` - TSV files containing read ids, their respective barcodes and UMIs, one file per every input file.
* `SAMPLE_ID.UMI_filtered.ED#.allinfo` - a TSV files reads kept after PCR deduplication.
   It contains all necessary information for downstream data analysis, such assigned gene, transcript, barcode, UMI etc.
* `SAMPLE_ID.UMI_filtered.ED#.stats.tsv` - brief statistics for the PCD deduplication and the resulting TSV files.
* `SAMPLE_ID.read_assignments.tsv.gz` - TSV file with read to isoform assignments (gzipped by default);
* `SAMPLE_ID.corrected_reads.bed.gz` - BED file with corrected read alignments (gzipped by default);
* `SAMPLE_ID.transcript_tpm.tsv` - TSV file with reference transcript expression in TPM;
* `SAMPLE_ID.transcript_counts.tsv` - TSV file with raw read counts for reference transcript;
* `SAMPLE_ID.gene_tpm.tsv` - TSV file with reference gene expression in TPM;
* `SAMPLE_ID.gene_counts.tsv` - TSV file with raw read counts for reference genes;

Discovered novel transcripts and genes (reads after PCR deduplication are used):
* `SAMPLE_ID.extended_annotation.gtf` - extended gene annotation:original input annotation plus novel genes and transcripts;
* `SAMPLE_ID.transcript_models.gtf` - expressed known and novel transcripts;* 
* `SAMPLE_ID.discovered_gene_counts.tsv`, `SAMPLE_ID.discovered_gene_tpm.tsv` - bulk gene counts/TPMs including known and novel genes;
* `SAMPLE_ID.discovered_transcript_counts.tsv`, `SAMPLE_ID.discovered_transcript_tpm.tsv` - bulk gene counts/TPMs including known and novel transcripts;
* `SAMPLE_ID.transcript_model_reads.tsv.gz` - reads assigned to novel and known transcripts.


Additionally, a log file will be saved to the directory.  
* <output_dir>/isoquant.log   

If raw reads were provided, BAM file(s) will be stored in `<output_dir>/<SAMPLE_ID>/aux/`.  

### Output file formats

Although most output files include headers that describe the data, a brief explanation of the output files is provided below.

#### Detected barcodes
* `read_id` - read id;
* `barcode` - detected barcode sequence, `*` if not detected;
* `UMI` - detected UMI sequence, `*` if not detected;
* `BC_score` - barcode alignemnt score, `-1` if not detected;
* `valid_UMI` - indicates whether UMI has the length similar to expected (True/False);
* `strand` - read strand (+/-/.)

Other columns contain information on positions of various primers and linkers, depending on the protocol being used.

#### UMI-filtered reads

* `read_id` - read id;
* `gene_id` - assigned gene id;
* `cell_type` - assigned cell type (None if not assigned);
* `barcode` - detected barcode sequence;
* `UMI` - detected UMI sequence;
* `exons` - exon coordinates;
* `TSS` - transcript start position;
* `TES` - transcript end position;
* `intons` - intron coordinates;
* `transcript_type` - transcript type (known/novel);
* `exon_num` - number of exons;
* `transcript_id` - assigned transcript id;
* `transcript_type` - transcript type derived from the gene annotation;

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


<a name="sec4"></a>
## Citation

If you use Spl-IsoQuant in your research, please, cite the following papers.

- Bulk analysis: [Prjibelski at al.](https://doi.org/10.1038/s41587-022-01565-y);
- Spl-Iso-Seq (Curio): [Foord, Prjibelski, Hu et al.](https://www.nature.com/articles/s41467-025-63301-9);
- StereoSeq, 10x Visium or Visium HD: [Michielsen, Prjibelski, Foord et al.](https://www.biorxiv.org/content/10.1101/2025.06.25.661563v1.abstract).

If you use it for 10x single-cell data, feel free to cite any of the last two papers or both of them.

<a name="sec5"></a>
## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve Spl-IsoQuant. 
If you have any troubles running Spl-IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/algbio/spl-IsoQuant/).


