[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/algbio/spl-IsoQuant)](https://github.com/algbio/spl-IsoQuant/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/algbio/spl-IsoQuant/total.svg?style=social&logo=github&label=Download)](https://github.com/algbio/spl-IsoQuant/releases)


# Spl-IsoQuant

This is a fork of the original [IsoQuant](https://github.com/ablab/IsoQuant) project.
The same group maintains both repositories equally well.

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
Spl-IsoQuant allows to reconstruct and quantify transcript models with
high precision and decent recall. If the reference annotation is given, Spl-IsoQuant also
assigns reads to the annotated isoforms based on their intron and exon structure.
Spl-IsoQuant further performs annotated gene, isoform, exon and intron quantification.
If reads are grouped (e.g. according to cell type or spatial location), counts are reported according to the provided grouping.

Latest Spl-IsoQuant version can be downloaded from [github.com/algbio/spl-IsoQuant/releases/latest](https://github.com/algbio/spl-IsoQuant/releases/latest).

Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

## Supported sequencing data

Spl-IsoQuant supports all kinds of long RNA data:
* PacBio CCS
* ONT dRNA / ONT cDNA
* Assembled / corrected transcript sequences

Reads must be provided in FASTQ/FASTA format (can be gzipped) or unmapped BAM format.
If you have already aligned your reads to the reference genome, simply provide sorted and indexed BAM files.
Spl-IsoQuant expects reads to contain polyA tails. For more reliable transcript model construction do not trim polyA tails.

Spl-IsoQuant can also take aligned Illumina reads to correct long-read spliced alignments. However, short reads are _not_
used to discover transcript models or compute abundances.


## Supported reference data

Reference genome is mandatory and should be provided in multi-FASTA format (can be gzipped).

Reference gene annotation is not mandatory, but is likely to increase precision and recall.
It can be provided in GFF/GTF format (can be gzipped).

Pre-constructed `minimap2` index can also be provided to reduce mapping time.


## Citation
The paper describing IsoQuant algorithms and benchmarking is available at [10.1038/s41587-022-01565-y](https://doi.org/10.1038/s41587-022-01565-y).

To try Spl-IsoQuant you can use the data that was used in the publication [zenodo.org/record/7611877](https://zenodo.org/record/7611877).


## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve Spl-IsoQuant. If you have any troubles running Spl-IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/algbio/spl-IsoQuant/issues) or send them via email: isoquant.rna@gmail.com.



## Quick start

*   Full Spl-IsoQuant documentation is available at [algbio.github.io/spl-IsoQuant](https://algbio.github.io/spl-IsoQuant/).

*   Spl-IsoQuant can be downloaded from [github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant):

        git clone https://github.com/algbio/spl-IsoQuant.git
        cd spl-IsoQuant
        pip install -r requirements.txt

*   If installing manually, you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

*   Verify your installation by running:

        splisoquant.py --test

*   To run Spl-IsoQuant on raw FASTQ/FASTA files use the following command

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./splisoquant.py --reference tests/toy_data/MAPT.Mouse.reference.fasta \
        --genedb tests/toy_data/MAPT.Mouse.genedb.gtf \
        --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq \
        --data_type nanopore -o toy_data_out


* To run Spl-IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --bam /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./splisoquant.py --reference tests/toy_data/MAPT.Mouse.reference.fasta \
        --genedb tests/toy_data/MAPT.Mouse.genedb.gtf \
        --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq \
        --data_type nanopore -o toy_data_out

* If using official annotations containing `gene` and `transcript` features use `--complete_genedb` to save time.

* Using reference annotation is optional since version 3.0, you may perform de novo transcript discovery without providing `--genedb` option:

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio|nanopore) -o OUTPUT_FOLDER

* If multiple files are provided, Spl-IsoQuant will create a single output annotation and a single set of gene/transcript expression tables.

