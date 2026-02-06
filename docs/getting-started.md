# Quick start

*   Spl-IsoQuant can be downloaded from [https://github.com/algbio/spl-IsoQuant](https://github.com/algbio/spl-IsoQuant):

        git clone https://github.com/algbio/spl-IsoQuant.git
        cd spl-IsoQuant
        pip install -r requirements.txt

*   If installing manually, you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pyfaidx](https://pypi.org/project/pyfaidx/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

*   Verify your installation by running:

        splisoquant.py --test

*   To run Spl-IsoQuant on 10x single-cell data use the following command:

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
        --fastq /PATH/TO/10x.fastq.gz --barcode_whitelist /PATH/TO/barcodes.tsv \
        --barcode2spot /PATH/TO/barcodes_to_celltype.tsv
        --mode tenX_v3 --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER

*    Or provide your own table with barcoded reads:


        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
        --fastq /PATH/TO/10x.fastq.gz --barcoded_reads /PATH/TO/barcoded_reads.tsv \
        --barcode2spot /PATH/TO/barcodes_to_celltype.tsv
        --mode tenX_v3 --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER

*    For example, using the toy Stereo-seq data provided within this repository:


        ./splisoquant.py --data_type nanopore --mode stereoseq_nosplit  \
        --fastq /home/andreyp/ablab/spl-IsoQuant/tests/stereo/S1.4K.subsample.fq.gz \
        --barcode_whitelist /home/andreyp/ablab/spl-IsoQuant/tests/stereo/barcodes.tsv \
        --reference /home/andreyp/ablab/spl-IsoQuant/tests/stereo/GRCm39.chrX.7.fa.gz \
        --genedb /home/andreyp/ablab/spl-IsoQuant/tests/stereo/gencode.chrX.ENSMUSG00000031153.gtf \
        --complete_genedb --output splisoquant_test  -p TEST_DATA


* You can also define your own molecule structure using the [molecule description format (MDF)](https://algbio.githu.io/spl-IsoQuant/single_cell.html#molecule-description-format-mdf) 
and provided to Spl-IsoQuant via `--molecule` option:

        splisoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf --complete_genedb \
        --fastq /PATH/TO/10x.fastq.gz --molecule /PATH/TO/my_protocol.mdf \
        --mode custom_sc --data_type (pacbio_ccs|nanopore) -o OUTPUT_FOLDER