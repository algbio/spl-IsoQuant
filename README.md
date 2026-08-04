[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/algbio/spl-IsoQuant)](https://github.com/algbio/spl-IsoQuant/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/algbio/spl-IsoQuant/total.svg?style=social&logo=github&label=Download)](https://github.com/algbio/spl-IsoQuant/releases)
[![UnitTests](https://github.com/algbio/spl-IsoQuant/actions/workflows/Unit_tests.yml/badge.svg)](https://github.com/algbio/spl-IsoQuant/actions/workflows/Unit_tests.yml)
[![User manual](https://github.com/algbio/spl-IsoQuant/actions/workflows/docs.yml/badge.svg)](https://algbio.github.io/spl-IsoQuant/)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/isoquant.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/isoquant)
[![PyPI Downloads](https://img.shields.io/pypi/v/isoquant)](https://pypi.org/project/isoquant/)


<img src="docs/splisoquant_logo.png" width="400" alt="Spl-IsoQuant">

### This is a fork of the original [IsoQuant](https://github.com/ablab/IsoQuant) project. The same group maintains both repositories equally.

[Full Spl-IsoQuant documentation can be found here](https://algbio.github.io/spl-IsoQuant/).
Information in this README is given only for convenience and is not a full user manual.

Current version: see `VERSION` file.

* [Citation information](#citation)
* [Feedback and bug reports](#feedback-and-bug-reports)
* [Quick start examples](#quick-start)


## About Spl-IsoQuant

Spl-IsoQuant is a tool for the genome-based analysis of long RNA reads, such as PacBio or
Oxford Nanopores. Spl-IsoQuant can perform the following analysis:
- **Quantification**:
  - Gene and transcript quantification;
  - Exon, splice junction, and intron retention quantification;
  - Detecting PolyA site and TSS;
  - Per-condition quantification (counts grouping);
- **Isoform and gene discovery**:
  - Annotation-based isoform discovery;
  - _De novo_ transcript and gene discovery;
  - Fusion gene discovery;
- **Single-cell and spatial transcriptomic analysis**:
  - Barcode calling for 10x single-cell, 10x Visium and Visium HD, Curio, Stereo-seq;
  - Universal barcode calling for user-specified protocols;
  - PCR deduplication;
  - Per-barcode and per-cell-type quantification;

All algorithms implemented within Spl-IsoQuant have been rigorously tested on real and simulated data and exhibit strong performance.

The latest Spl-IsoQuant version can be downloaded from [github.com/algbio/spl-IsoQuant/releases/latest](https://github.com/algbio/spl-IsoQuant/releases/latest).

Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

## Supported sequencing data

Spl-IsoQuant supports all kinds of long RNA data:
* PacBio CCS;
* ONT dRNA / ONT cDNA;
* Assembled / corrected transcript sequences.

Spl-IsoQuant can process bulk, single-cell and spatial long-read data. 
It supports 10x single-cell v2/v3/v4, 10x Visium/VisiumHD, Curio, Stereo-seq and arbitrary user-defined protocols. 
See [the user manual](https://algbio.github.io/spl-IsoQuant/single_cell.html) for details.

Reads must be provided in FASTQ/FASTA format (can be gzipped) or unmapped BAM format. 
If you have already aligned your reads to the reference genome, simply provide sorted and indexed BAM files.
Spl-IsoQuant expects reads to contain polyA tails. For more reliable transcript model construction do not trim polyA tails.

Spl-IsoQuant can also take aligned Illumina reads to correct long-read spliced alignments. However, short reads are _not_
used to discover transcript models or compute abundances.


## Supported reference data

Reference genome is mandatory and should be provided in multi-FASTA format (can be gzipped).

Reference gene annotation is not mandatory but is likely to increase precision and recall.
It can be provided in GFF/GTF format (can be gzipped).

Pre-constructed `minimap2` index can also be provided to reduce mapping time.


## Citation

If you use Spl-IsoQuant's transcript discovery and quantification in your research, please cite 
[Prjibelski, Mikheenko et al., 2023](https://doi.org/10.1038/s41587-022-01565-y).
Simulated data used in this paper is available here [zenodo.org/record/7611877](https://zenodo.org/record/7611877).

If you use Spl-IsoQuant for analysing Curio spatial data, please cite
[Foord, Prjibelski, Hu et al., 2025](https://doi.org/10.1038/s41467-025-63301-9).

If you use Spl-IsoQuant for processing Stereo-seq, 10x singlce-cell or 10x Visium/VisiumHD data, please cite 
[Michielsen, Prjibelski, Foord et al., 2026](https://www.biorxiv.org/content/10.1101/2025.06.25.661563v3).


## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve Spl-IsoQuant. If you have any troubles running Spl-IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/algbio/spl-IsoQuant/issues) or send them via email: isoquant.rna@gmail.com.



## Quick start

*   Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

*   Spl-IsoQuant is installed from GitHub:

        git clone https://github.com/algbio/spl-IsoQuant.git 
        cd spl-IsoQuant
        git checkout latest
        pip install -e .

Installation typically takes no more than a few minutes.

*   Spl-IsoQuant is **not** distributed via PyPI or bioconda and will not be published there.
    The original [IsoQuant](https://github.com/ablab/IsoQuant) provides the same functionality
    and is available from both:

        pip install isoquant

        conda create -c conda-forge -c bioconda -n isoquant python=3.12 isoquant

    These install the `isoquant` executable, not `splisoquant`. If you install this way,
    use `isoquant` in place of `splisoquant` in all commands below.

* If running simply from [the source archive](https://github.com/algbio/spl-IsoQuant/releases/) you will need:
  - Python3 (3.8 or higher)
  - [minimap2](https://github.com/lh3/minimap2)
  - [samtools](http://www.htslib.org/download/)
  - Python libraries from `requirements.txt` to be installed: 
  
        pip install -r requirements.txt

*   Verify your installation by running (typically takes less than 1 minute):

        splisoquant --test

*   To run Spl-IsoQuant on raw FASTQ/FASTA files, use the following command

        splisoquant --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        splisoquant --fastq isoquant_tests/simple_data/chr9.4M.ont.sim.fq.gz \
        --reference isoquant_tests/simple_data/chr9.4M.fa.gz \
        --genedb isoquant_tests/simple_data/chr9.4M.gtf.gz \
        --data_type nanopore --complete_genedb -p TEST_DATA --output splisoquant_test 


* To run Spl-IsoQuant on aligned reads (make sure your BAM is sorted and indexed), use the following command:

        splisoquant --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --bam /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

* If using official annotations containing `gene` and `transcript` features use `--complete_genedb` to save time.

* Using reference annotation is optional since version 3.0, you may preform de novo transcript discovery without providing `--genedb` option:

        splisoquant --reference /PATH/TO/reference_genome.fasta \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

* If multiple files are provided, Spl-IsoQuant will create a single output annotation and a single set of gene/transcript expression tables.

* To perform different types of analysis, check out `--analysis` option.

* For single-cell and spatial transcriptomic data, check out `--mode` option.

