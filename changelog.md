# Spl-IsoQuant changelog

## Spl-IsoQuant 4.1.0, 28 October 2026

- New options for cell barcode calling:
  - `--n_cells` to select cell barcodes from a large whitelist (e.g. a stock 10x whitelist), which significantly improves barcode calling recall;
  - `--barcode_whitelist auto` to detect cell barcodes without a whitelist;
  - `--split_molecules` to control splitting of concatenated cDNAs, replaces `tenX_v3_split`, `tenX_v2_split` and `stereoseq_nosplit` modes;
- New `--polya_trimmed` values to use external poly-A tail information: `tag:<TAG>` (e.g. Dorado `pt` tag), `list:<FILE>` (list of read IDs) and `flnc:<FILE>` (`isoseq refine` report);
- Spliced/unspliced counts for RNA velocity analysis in single-cell and spatial modes, including `.loom` output;
- HTML report with run summary and QC statistics for every experiment;
- New exon counts with independent exon quantification ([#426](https://github.com/ablab/IsoQuant/discussions/426));
- Fixed intron retention counts, exclusion counts now correspond to reads with the intron spliced ([#426](https://github.com/ablab/IsoQuant/discussions/426));
- Fixed a minor bug in UMI deduplication;
- New BAM outputs for single-cell and spatial modes, set via `--large_output`:
  - `tagged_bam` with barcode and UMI tags;
  - `deduplicated_bam` containing only UMI-deduplicated alignments;
- New `pacbio_all` preset for `--model_construction_strategy` for very sensitive transcript discovery (e.g. FLNCs);
- Fixed structural categories in the SQANTI-like output ([#265](https://github.com/ablab/IsoQuant/issues/265));
- Fixed polyA/TSS prediction in rare cases of overlapping genes;
- Fixed crash with `--check_canonical` ([#425](https://github.com/ablab/IsoQuant/issues/425), [#427](https://github.com/ablab/IsoQuant/issues/427));
- Reworked `--resume` option allowing resuming crashed runs (partial runs produced by previous versions cannot be resumed).

## Spl-IsoQuant 4.0.0, 5 August 2026

Since Spl-IsoQuant 4.0.0 versions duplicate IsoQuant version. 

- Fusion gene detection (developed by @solatar);
- PolyA site/TSS discovery (developed by @arina-c);
- Significantly improved accuracy of PolyA and TSS positions for discovered transcripts;
- Improved splice-site accuracy for discovered transcripts (developed by @heidi-holappa);
- Reworked read assignment and quantification for discovered transcripts to improve consistency with reference-based analysis;
- New exon counting strategies (two new output files);
- Added intron retention counts;
- Added `fsm_only` transcript/gene counting strategy ([#338](https://github.com/ablab/IsoQuant/issues/338));
- New universal read assignment format for all pipelines, possible to convert to old formats; 
- New key options:
  - `--mode` to switch between bulk, single-cell and spatial protocols;
  - `--analysis` to control the output produced;
- Revised and updated the documentation;
- Multiple minor bug-fixes.

## Spl-IsoQuant 2.4.0, 26 April 2026

- Improved recall for 10x single-cell barcode calling;

- Added cDNA splitting mode for 10x barcode calling;

- Improved exon and splice junction counts output (TSV matrix and mtx format);

- Added `no_auto` read grouping mode to avoid automatic grouping.

## Spl-IsoQuant 2.3.0, 13 March 2026

- TSV output files do not have a `#` in headers for easier downstream processing (thanks to @yinshiyi).

- Fixed Python 3.14 compatibility.

- Added `--barcoded_bam` flag for reading barcodes and UMIs directly from a input BAM.

- Added `--barcode2barcode` option for Visium HD data processing.

- Dramatic speed-up for 10x and universal barcode calling on large whitelists.

- Per-barcode counts grouping is now performed automatically if the barcodes are set.

## Spl-IsoQuant 2.2.1, 11 February 2026

- Fixed Python 3.8 compatibility.

- Fixed MTX conversion.

- Fixed`--test` message.

## Spl-IsoQuant 2.2.0, 5 February 2026

Support for 10x single-cell v3 and 10x Visium 3' data.

Universal barcode detection via user-defined molecule structure.

`--read_group` now supports multiple read grouping strategies. 
You can now simultaneously group counts by samples, BAM tags, barcode and barcode attributesm, 
read attributes provided in separate TSV files or within read ids themselves.

New `--large_output` option to control which large output files are generated.

New `--barcode2spot` option that allows to group counts by barcode attributes, such as cell type or spatial location.

Major performance improvements, especially related to Stereo-seq barocode detection.

## Spl-IsoQuant 2.1.0, 22 September 2025

Support for 10x Visium HD and Visium 3' data.

Significant performance and usability improvements.

## Spl-IsoQuant 2.0.0, 19 June 2025

Spl-IsoQuant 2 for Stereo-seq long read data.

## Spl-IsoQuant 1.0.0, 12 December 2024

Improve poly(A) processing and treatment of ambiguously assigned reads.

## Spl-IsoQuant 0.1.0, 12 June 2024

A initial release for Spl-IsoQuant, based on IsoQuant 3.4.1
