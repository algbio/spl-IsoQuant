[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/algbio/spl-IsoQuant)](https://github.com/algbio/spl-IsoQuant/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/algbio/spl-IsoQuant/total.svg?style=social&logo=github&label=Download)](https://github.com/algbio/spl-IsoQuant/releases)


<img src="docs/splisoquant_logo.png" width="400" alt="IsoQuant">

This is a fork of the original [IsoQuant](https://github.com/ablab/IsoQuant) project.
The same group maintains both repositories equally.

[Full Spl-IsoQuant documentation can be found here](https://algbio.github.io/spl-IsoQuant/).
Information in this README is given only for convenience and is not a full user manual.

Current version: see `VERSION` file.

* [Citation information](#citation)
* [Feedback and bug reports](#feedback-and-bug-reports)
* [Quick start examples](#quick-start)


## About Spl-IsoQuant

Spl-IsoQuant is a tool for **single-cell and spatial long-read transcriptomics** analysis.
It performs genome-based analysis of long RNA reads from platforms such as PacBio or
Oxford Nanopore, with specialized support for single-cell and spatial protocols.
Spl-IsoQuant is capable of perfroming barcode and UMI detection for various sequencing protocols, 
UMI deduplication and barcode-aware quantification of reads, where reads are grouped
(e.g. according to cell types or spatial location), counts are reported according to the provided grouping.
We recommend providing smaller barcode whitelists (e.g. obtained from short-read sequencing) to achieve higher accuracy. 

Similarly to IsoQuant, it can also perform novel transcript discovery. 
However, in single-cell/spatial mode, Spl-IsoQuant will only discover 
novel isoforms for known genes, as reads that are not assigned to any 
known gene are discarded during the PCR deduplication step. 
To achieve full transcript discovery, run Spl-IsoQuant in bulk mode.

The latest Spl-IsoQuant version can be downloaded from [github.com/algbio/spl-IsoQuant/releases/latest](https://github.com/algbio/spl-IsoQuant/releases/latest).

Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

## Supported sequencing data

Spl-IsoQuant supports all kinds of long RNA data:
* PacBio CCS
* ONT dRNA / ONT cDNA
* Assembled / corrected transcript sequences

Reads must be provided in FASTQ/FASTA format (can be gzipped) or unmapped BAM format.
If you have already aligned your reads to the reference genome, simply provide sorted and indexed BAM files.

Spl-IsoQuant supports the following protocols:

* 10x 3' v3 single-cell;
* 10x 3' Visium spatial data;
* 10x Visium HD;
* Curio Biosciences spatial data;
* Stereo-seq spatial data;
* Any other single-cell or spatial protocol with barcode and UMI sequences (see more about [custom molecule description](https://algbio.github.io/spl-IsoQuant/single_cell.html#molecule-description-format-mdf)).


## Supported reference data

Reference genome is mandatory and should be provided in multi-FASTA format (can be gzipped).

Reference gene annotation is also mandatory for single-cell / spatial analysis.
It can be provided in GFF/GTF format (can be gzipped).

Pre-constructed `minimap2` index can also be provided to reduce mapping time.

## Citation

- [Spatial isoform sequencing at sub-micrometer single-cell resolution reveals novel patterns of spatial isoform variability in brain cell types
L Michielsen, AD Prjibelski, C Foord et al., bioRxiv, 2025.06. 25.661563. ](https://www.biorxiv.org/content/10.1101/2025.06.25.661563.abstract)

- [A spatial long-read approach at near-single-cell resolution reveals developmental regulation of splicing and polyadenylation sites in distinct cortical layers and cell types
C Foord, AD Prjibelski, W Hu et al., Nature Communications 16 (1), 8093](https://www.nature.com/articles/s41467-025-63301-9)

## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve Spl-IsoQuant. If you have any troubles running Spl-IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/algbio/spl-IsoQuant/issues).



## Quick start

*   Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

*   Spl-IsoQuant can be downloaded from [github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant):

        git clone https://github.com/algbio/spl-IsoQuant.git
        cd spl-IsoQuant
        pip install -r requirements.txt

*   If installing manually, you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

*   Verify your installation by running:

        splisoquant.py --test

*  To run Spl-IsoQuant on 10x single-cell data use the following command:

```
splisoquant.py --reference /PATH/TO/reference_genome.fasta \
--genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
--fastq /PATH/TO/10x.fastq.gz --barcode_whitelist /PATH/TO/barcodes.tsv \
--barcode2spot /PATH/TO/barcodes_to_celltype.tsv
--mode tenX_v3 --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER
```

*  Or provide your own table with barcoded reads:

```
splisoquant.py --reference /PATH/TO/reference_genome.fasta \
--genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
--fastq /PATH/TO/10x.fastq.gz --barcoded_reads /PATH/TO/barcoded_reads.tsv \
--barcode2spot /PATH/TO/barcodes_to_celltype.tsv
--mode tenX_v3 --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER
```

*  For example, using the toy Stereo-seq data provided within this repository:

```
./splisoquant.py --data_type nanopore --mode stereoseq_nosplit  \
--fastq /home/andreyp/ablab/spl-IsoQuant/tests/stereo/S1.4K.subsample.fq.gz \
--barcode_whitelist /home/andreyp/ablab/spl-IsoQuant/tests/stereo/barcodes.tsv \
--reference /home/andreyp/ablab/spl-IsoQuant/tests/stereo/GRCm39.chrX.7.fa.gz \
--genedb /home/andreyp/ablab/spl-IsoQuant/tests/stereo/gencode.chrX.ENSMUSG00000031153.gtf \
--complete_genedb --output splisoquant_test  -p TEST_DATA
```

*  You can also define your own molecule structure using the [molecule description format (MDF)](https://algbio.githu.io/spl-IsoQuant/single_cell.html#molecule-description-format-mdf) 
and provided to Spl-IsoQuant via `--molecule` option:

```
splisoquant.py --reference /PATH/TO/reference_genome.fasta \
--genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
--fastq /PATH/TO/10x.fastq.gz --molecule /PATH/TO/my_protocol.mdf \
--mode custom_sc --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER
```