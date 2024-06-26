#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc

from common import *


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", default="gtf_stats")
    parser.add_argument("--genedb", "-d", type=str, help="prefix to reduced gene db")
    parser.add_argument("--gtf", "-g", type=str, help="gtf to assess")
    parser.add_argument("--tool", type=str, choices=SEPARATE_FUNCTORS.keys(),
                        help="tool used for generating GTF")

    args = parser.parse_args()
    if not args.genedb or not args.gtf or not args.tool:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    out_full_path = os.path.join(args.output, args.tool + ".full.gtf")
    out_known_path = os.path.join(args.output, args.tool + ".known.gtf")
    out_novel_path= os.path.join(args.output, args.tool + ".novel.gtf")
    print("Seprating known and novel transcripts")
    separator = SEPARATE_FUNCTORS[args.tool](args.gtf)
    split_gtf(args.gtf, separator, out_full_path, out_known_path, out_novel_path)
    print("Running gffcompare for entire GTF")
    expressed_gtf = args.genedb + ".expressed.gtf"
    run_gff_compare(expressed_gtf, out_full_path, os.path.join(args.output, args.tool + ".full.stats"))
    print("Running gffcompare for known transcripts")
    expressed_gtf = args.genedb + ".expressed_kept.gtf"
    run_gff_compare(expressed_gtf, out_known_path, os.path.join(args.output, args.tool + ".known.stats"))
    print("Running gffcompare for novel transcripts")
    expressed_gtf = args.genedb + ".excluded.gtf"
    run_gff_compare(expressed_gtf, out_novel_path, os.path.join(args.output, args.tool + ".novel.stats"))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

