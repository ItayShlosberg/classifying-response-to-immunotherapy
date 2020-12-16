"""
This script is useful when there is a new sample (given in FASTQ format) that underwent CellRanger_count pipeline
and now we've got the output of CellRanger (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz).
After running it we'll have droplet_RNAseq object (python object).
"""

from os.path import join
import csv
import gzip
import os
from scipy.io import mmread
import numpy as np
import pickle
from utilities.droplet_dataset import *

# features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz should be in that path.
SAMPLE_PATH = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\milestones\milestone 3 - 2.12.20\M145_M146\M146_filtered_feature_bc_matrix'
OUTPUT_PATH = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples'
SAMPLE_NAME = 'M146'


def convert_sample():
    # Extract matrix
    mat = np.array(mmread(join(SAMPLE_PATH, "matrix.mtx.gz")).todense()).astype(np.uint16).T

    # Extract genes names and features
    features_path = join(SAMPLE_PATH, "features.tsv.gz")
    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]

    # Not in ues
    feature_types = [row[2] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]

    # Extract cells' barcodes
    barcodes_path = join(SAMPLE_PATH, "barcodes.tsv.gz")
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, 'rt'), delimiter="\t")]

    # Create sample object
    rna_sample = RNAseq_Sample(mat, gene_names, barcodes, feature_ids)
    return rna_sample


if __name__ == '__main__':
    rna_sample = convert_sample()

    # Save sample
    output_folder = join(OUTPUT_PATH, SAMPLE_NAME)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    pickle.dump((rna_sample), open(os.path.join(output_folder, "RNA_sample.pkl"), "wb"))