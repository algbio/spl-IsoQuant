#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import sys
import os
from traceback import print_exc
import logging
import argparse
import subprocess
import pysam
import random
import numpy as np
from collections import defaultdict

logger = logging.getLogger('LRGASP')


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output file with counts")
    parser.add_argument("--fastq", "-f", help="long reads in FASTQ format (PacBio/ONT)", type=str)
    parser.add_argument("--threads", "-t", help="number of CPU threads for minimap [16]", default=16, type=int)
    parser.add_argument("--reference_transcripts", "-r", help="reference transcripts in FASTA format", type=str)
    parser.add_argument("--mandatory", "-m", help="file with a list of mandatory transcripts to be included, "
                                                  "counts are assigned randomly, "
                                                  "can be a TSV with transcript ids in the first column", type=str)
    parser.add_argument("--seed", "-s", help="randomizer seed [11]", default=11, type=int)

    args = parser.parse_args(args, namespace)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


def check_params(args):
    return args.fastq is not None and args.reference_transcripts is not None


def set_logger(args, logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info(" === LRGASP quantification pipeline started === ")
    random.seed(args.seed)
    logger.info("Mapping reads with minimap2...")
    samfile_name = args.output + ".sam"
    result = subprocess.run(["minimap2", args.reference_transcripts, args.fastq, "-x", "map-ont",
                             "-t", str(args.threads), "-a", "--secondary=no", "-o", samfile_name])
    if result.returncode != 0:
        logger.error("minimap2 failed with code %d, make sure it is in you $PATH variable" % result.returncode)
        exit(-1)

    logger.info("Quantifying transcripts...")
    transcript_counts = defaultdict(int)
    with pysam.AlignmentFile(samfile_name, "r") as samfile_in:
        for alignment in samfile_in:
            transcript_id = alignment.reference_name
            if alignment.reference_id == -1 or alignment.is_supplementary or alignment.is_secondary:
                continue
            transcript_counts[transcript_id] += 1
    os.remove(samfile_name)

    # reading mandatory transcripts
    mandatory_transcripts = []
    if args.mandatory is not None:
        for l in open(args.mandatory):
            l = l.strip()
            if not l or l.startswith("#"):
                continue
            mandatory_transcripts.append(l.split()[0])

    # adding "novel" transcript that must be in the simulated data
    for transcript_id in mandatory_transcripts:
        # random counts with the possibility to have a few very low-expressed and high-expressed isoforms
        transcript_counts[transcript_id] = 1 + int(np.random.gamma(1.4, 50))

    count_sum = 0.0
    for count in transcript_counts.values():
        # we don't devide by gene length here as we assume each long read is a single transcripts (unlike in short reads)
        count_sum += count
    scaling_factor = count_sum / 1000000.0

    outf = open(args.output, "w")
    outf.write("#transcript_id\tcounts\ttpm\n")
    for transcript_id in transcript_counts.keys():
        outf.write("%s\t%d\t%0.6f\n" % (transcript_id, transcript_counts[transcript_id],
                                        transcript_counts[transcript_id] / scaling_factor))
    outf.close()
    logger.info("Done. Output counts are stored in " + args.output)
    logger.info(" === LRGASP quantification pipeline finished === ")


def main(args):
    args = parse_args(args)
    set_logger(args, logger)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
