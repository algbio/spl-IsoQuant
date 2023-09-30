#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import logging
import editdistance

from umi_filtering import UMIFilter, filter_bam


logger = logging.getLogger('IsoQuant')


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", type=str, help="read - barcode - UMI table", required=True)
    parser.add_argument("--read_assignments", "-r", type=str, help="IsoQuant read assignments", required=True)
    parser.add_argument("--bam", type=str, help="original BAM file, provide only if you need a BAM file"
                                                " with UMI-filtered alignments")

    parser.add_argument("--min_distance", type=int, help="minimal edit distance for UMIs to be considered distinct;"
                                                         "length difference is added to this values by default;"
                                                         "0 for equal UMIs, -1 for keeping only a single gene-barcode "
                                                         "read (default: 3)", default=3)
    parser.add_argument("--untrusted_umis", action="store_true", help="allow untrusted UMIs to be used", default=False)
    parser.add_argument("--only_spliced", action="store_true", help="keep only spliced reads", default=False)
    parser.add_argument("--only_unique", action="store_true", help="keep only non-ambiguous reads", default=False)
    parser.add_argument("--disregard_length_diff", action="store_true", help="do not account for length difference "
                                                                             "when comparing UMIs", default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    umi_filter = UMIFilter(args)
    umi_filter.process(args.read_assignments, args.output)
    if args.bam:
        filter_bam(args.bam, args.output  + ".UMI_filtered.reads.bam", umi_filter.selected_reads)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
