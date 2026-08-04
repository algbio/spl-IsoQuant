# Quick start

* Spl-IsoQuant is installed from GitHub:

        git clone https://github.com/algbio/spl-IsoQuant.git 
        cd spl-IsoQuant
        git checkout latest
        pip install -e .

Installation typically takes no more than a few minutes.

* Spl-IsoQuant is **not** distributed via PyPI or bioconda and will not be published there.
  The original [IsoQuant](https://github.com/ablab/IsoQuant) provides the same functionality
  and is available from both:

        pip install isoquant

        conda create -c conda-forge -c bioconda -n isoquant python=3.12 isoquant

  These install the `isoquant` executable, not `splisoquant`. If you install this way,
  use `isoquant` in place of `splisoquant` in all commands below.

* If running simply from [the source archive](https://github.com/algbio/spl-IsoQuant/releases/) you will need:
  - Python3 (3.8 or higher)
  - [gffutils](https://pythonhosted.org/gffutils/installation.html)
  - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  - [biopython](https://biopython.org/)
  - [pyfaidx](https://pypi.org/project/pyfaidx/)
  - [ssw-py](https://pypi.org/project/ssw-py/)
  - [editdistance](https://pypi.org/project/editdistance/)
  - [numba](https://numba.pydata.org/)
  - [mappy](https://pypi.org/project/mappy/)
  - [intervaltree](https://pypi.org/project/intervaltree/)
  - [xgboost](https://xgboost.readthedocs.io/)
  - [scikit-learn](https://scikit-learn.org/)

and some other common Python libraries to be installed. See `requirements.txt` for details. 
  
You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

* All required Python libraries can be installed via: 

      pip install -r requirements.txt

* Verify your installation by running (typically takes less than 1 minute):

        splisoquant --test

* To run Spl-IsoQuant on raw FASTQ/FASTA files, use the following command

        splisoquant --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        splisoquant --fastq isoquant_tests/simple_data/chr9.4M.ont.sim.fq.gz \
        --reference isoquant_tests/simple_data/chr9.4M.fa.gz \
        --genedb isoquant_tests/simple_data/chr9.4M.gtf.gz \
        --data_type nanopore --complete_genedb -p TEST_DATA --output splisoquant_test 


* To run Spl-IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

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