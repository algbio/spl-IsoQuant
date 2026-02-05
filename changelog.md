# Spl-IsoQuant changelog

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