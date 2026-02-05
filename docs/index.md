# About Spl-IsoQuant

Spl-IsoQuant is a tool for **single-cell and spatial long-read transcriptomics** analysis.
It performs genome-based analysis of long RNA reads from platforms such as PacBio or
Oxford Nanopore, with specialized support for single-cell and spatial protocols.
Spl-IsoQuant allows to reconstruct and quantify transcript models with
high precision and decent recall. If the reference annotation is given, Spl-IsoQuant also
assigns reads to the annotated isoforms based on their intron and exon structure.
Spl-IsoQuant further performs annotated gene, isoform, exon and intron quantification.
If reads are grouped (e.g. according to cell type or spatial location), counts are reported according to the provided grouping.

Spl-IsoQuant consists of two stages, which generate its own output:

1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isoform assignment,
splice site correction and abundance quantification for reference genes/transcripts.
2. Transcript discovery. Reconstructs transcript models and performs abundance quantification for discovered isoforms.

Latest Spl-IsoQuant version can be downloaded from [https://github.com/algbio/spl-IsoQuant/releases/latest](https://github.com/algbio/spl-IsoQuant/releases/latest).

### Spl-IsoQuant pipeline
![Pipeline](isoquant_pipeline.png)
