############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import logging
from collections import defaultdict
import pickle

import pandas as pd

from .common import proper_plural_form


logger = logging.getLogger('IsoQuant')


class EnumStats:
    def __init__(self, inf=None):
        self.stats_dict = defaultdict(int)
        if inf is not None:
            self.load(inf)

    # element must be Enum
    def add(self, element, count=1):
        self.stats_dict[element] += count

    def print_start(self, header_string=""):
        if header_string:
            logger.info(header_string)
        for e in sorted(self.stats_dict.keys(), key=lambda x: x.name):
            logger.info("%s: %d" % (e.name, self.stats_dict[e]))

    def dump(self, out_file):
        pickler = pickle.Pickler(open(out_file, "wb"), -1)
        pickler.dump(self.stats_dict)

    def load(self, in_file):
        unpickler = pickle.Unpickler(open(in_file, "rb"), fix_imports=False)
        self.stats_dict = unpickler.load()

    def merge(self, other):
        for e in other.stats_dict.keys():
            self.stats_dict[e] += other.stats_dict[e]


def transform_counts(path_to_csv, label, column_name='count', full=False):
    df = pd.read_csv(path_to_csv, sep='\t')
    df_features = df.copy() if full else df[:-3].copy()
    df_features.rename(columns={column_name: label}, inplace=True)
    return df_features


def combine_table(input_data, output, get_file_name, output_file_name, column_name='count', full=False):
    sample_0 = input_data.samples[0]
    combined_table = transform_counts(get_file_name(sample_0), sample_0.prefix, column_name, full)
    for sample in input_data.samples[1:]:
        combined_table = pd.merge(combined_table,
                                  transform_counts(get_file_name(sample), sample.prefix, column_name, full),
                                  on='#feature_id', how='outer')

    combined_table.to_csv(os.path.join(output, output_file_name), sep='\t', index=False)


def combine_counts(input_data, output):
    logger.info("Aggregating counts for " + proper_plural_form("sample", len(input_data.samples)))
    combine_table(input_data, output, lambda x: x.out_gene_counts_tsv + "_counts.tsv", "combined_gene_counts.tsv")
    combine_table(input_data, output, lambda x: x.out_gene_counts_tsv + "_tpm.tsv", "combined_gene_tpm.tsv",
                  column_name='TPM', full=True)
    combine_table(input_data, output, lambda x: x.out_transcript_counts_tsv + "_counts.tsv",
                  "combined_transcript_counts.tsv")
    combine_table(input_data, output, lambda x: x.out_transcript_counts_tsv + "_tpm.tsv",
                  "combined_transcript_tpm.tsv", column_name='TPM', full=True)
    logger.info("Aggregation finished")
