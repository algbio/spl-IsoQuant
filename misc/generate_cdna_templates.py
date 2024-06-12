############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import math
import os
import sys
import argparse
from Bio import SeqIO
import random
from common import *
from traceback import print_exc


POLYA_LEN = 30
LINKER = "TCTTCAGCGTTCCCGAGA"
UPS_PRIMER = "AAGCAGTGGTATCAACGCAGAGT"
PCR_PRIMER = "CTACACGACGCTCTTCCGATCT"
UPS_PRIMER_REV = reverese_complement(UPS_PRIMER)
LEFT_BC_LENGTH = 8
RIGHT_BC_LENGTH = 6
BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
UMI_LENGTH = 9

TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
R1 = "CTACACGACGCTCTTCCGATCT" # R1
BARCODE_LEN_10X = 16
UMI_LEN_10X = 12

NUCLS = ['A', 'C', 'G', 'T']
VNUCLS =  ['A', 'C', 'G']


def get_random_seq(length, non_t_tail=0):
    seq = ""
    for i in range(length-non_t_tail):
        index = random.randint(0, len(NUCLS) - 1)
        seq += NUCLS[index]
    for i in range(non_t_tail):
        index = random.randint(0, len(VNUCLS) - 1)
        seq += VNUCLS[index]
    return seq


def create_template_spatial(sequence, barcode, umi):
    assert len(barcode) == BC_LENGTH
    barcoded_part = PCR_PRIMER + barcode[:LEFT_BC_LENGTH] + LINKER + barcode[LEFT_BC_LENGTH:] + umi
    template_seq = barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + UPS_PRIMER_REV
    return reverese_complement(template_seq)


#AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT—--DNA------CCCATGTACTCTGCGTTGATACCACTGCTT

# p5="AATGATACGGCGACCACCGAGATCTACAC"
# i5="NNNNNNNNNN"
# truseqread1="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# 10xbc="NNNNNNNNNNNNNNNN"
# umi="NNNNNNNNNNNN"
# polyt="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
# tso="CCCATGTACTCTGCGTTGATACCACTGCTT"
def create_template_10x(sequence, barcode, umi):
    assert len(barcode) == BARCODE_LEN_10X
    barcoded_part = R1 + barcode + umi
    template_seq =  barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + TSO
    return reverese_complement(template_seq)


def load_counts(inf, total_count):
    count_dict = {}
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split("\t")
        count_dict[v[0]] = math.ceil(float(v[2]) * total_count / 1000000)
    return count_dict


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--transcriptome", "-t", help="reference transcriptome", required=True, type=str)
    parser.add_argument("--counts", "-c", help="expression table", type=str, required=True)
    parser.add_argument("--template_count", help="total number of templates to generate", type=int, default=100000)
    parser.add_argument("--barcodes", "-b", help="barcode list (random if not set)", type=str)
    parser.add_argument("--umis", "-u", help="UMI list (random if not set)", type=str)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--mode", help="[spatial | 10x]", type=str, default="spatial")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    isoforms = SeqIO.to_dict(SeqIO.parse(args.transcriptome, "fasta"))
    print("Loaded %d reference sequences from %s" % (len(isoforms), args.transcriptome))
    barcodes = load_barcodes(args.barcodes) if args.barcodes else []
    print("Using %d reference barcodes from %s" % (len(barcodes), args.barcodes))
    umis = load_barcodes(args.umis) if args.umis else []
    print("Using %d reference UMIs from %s" % (len(umis), args.umis))
    count_dict = load_counts(args.counts, args.template_count)
    print("Loaded %d counts from %s" % (len(count_dict), args.counts))

    bc_len = BC_LENGTH if args.mode == "spatial" else BARCODE_LEN_10X
    umi_len = UMI_LENGTH if args.mode == "spatial" else UMI_LEN_10X
    template_func = create_template_spatial if args.mode == "spatial" else create_template_10x

    scale_factor = 1000000.0 / sum(count_dict.values())
    counts_file = args.output + ".counts.tsv"
    with open(counts_file, "w") as expr_table:
        expr_table.write("#transcript_id\tcounts\ttpm\n")
        for t_id in count_dict:
            expr_table.write("%s\t%d\t%.4f\n" % (t_id, count_dict[t_id], scale_factor * count_dict[t_id]))
    print("Saved templates counts to %s" % counts_file)

    template_fname = args.output + ".templates.fasta"
    output_fasta = open(template_fname, "w")
    output_counts = open(args.output + ".template_counts.tsv", "w")
    for t_id in count_dict:
        transcript_seq = isoforms[t_id]
        template_list = []
        for i in range(count_dict[t_id]):
            bc = barcodes[random.randint(0, len(barcodes) - 1)] if barcodes else get_random_seq(bc_len)
            umi = barcodes[random.randint(0, len(umis) - 1)] if umis else get_random_seq(umi_len, 2)
            template = template_func(transcript_seq, bc, umi)
            seq_id = "READ_%d_%s_%s_%s" % (i, t_id, bc, umi)
            output_counts.write("%s\t%d\t%.4f\n" % (seq_id, 1, scale_factor))
            template_list.append(SeqIO.SeqRecord(seq=Seq.Seq(template), id=seq_id,  description="", name=""))
        SeqIO.write(template_list, output_fasta, "fasta")
    output_fasta.close()
    output_counts.close()
    print("Saved templates counts to %s" % template_fname)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
