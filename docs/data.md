# Supported data types

## Sequencing data

Spl-IsoQuant support all kinds of long RNA data:

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

Reference genome should be provided in multi-FASTA format (can be gzipped).
The reference genome is mandatory even when BAM files are provided.

Reference gene annotation is also mandatory for single-cell / spatial analysis
It should be provided in GFF/GTF format (can be gzipped).
It will be converted to [gffutils](https://pythonhosted.org/gffutils/installation.html) database. Information on converted databases will be stored in your `~/.config/Spl-IsoQuant/db_config.json` to increase speed of future runs. You can also provide gffutils database manually. 
Make sure that chromosome/scaffold names are identical in FASTA file and gene annotation.
Note that gffutils databases may not work correctly on NFS shares. It is possible to set a designated folder for 
the database with `--genedb_output` (different from the output directory).

Pre-constructed aligner index can also be provided to reduce mapping time.
