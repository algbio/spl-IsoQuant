

<img src="splisoquant_logo.png" width="400" alt="Spl-IsoQuant">

### This is a fork of the original [IsoQuant](https://github.com/ablab/IsoQuant) project. The same group maintains both repositories equally.

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

The latest Spl-IsoQuant version can be downloaded from [https://github.com/algbio/spl-IsoQuant/releases/latest](https://github.com/algbio/spl-IsoQuant/releases/latest).


### Spl-IsoQuant pipeline for single-cell and spatial transcriptomic data
![Pipeline](isoquant_pipeline_sc.png)

### Spl-IsoQuant pipeline for bulk data
![Pipeline](isoquant_pipeline.png)



