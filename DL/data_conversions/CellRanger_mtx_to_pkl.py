"""
This script is useful when there is a new sample (given in FASTQ format) that underwent CellRanger_count pipeline
and now we've got the output of CellRanger (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz).
After running it we'll have droplet_RNAseq object (python object).


dir hierarchy input required:
SAMPLES >> M1 >> (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz)
           M2 >> (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz)
           M3 >> (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz)
           M4 >> (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz)

"""
import sys
# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts'
import sys
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
from os.path import join
import csv
import gzip
import os
from scipy.io import mmread
import numpy as np
import pickle
from utilities.droplet_dataset import *

# features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz should be in that path.
SAMPLES_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\FASTAQ_OUTPUTS'
OUTPUT_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\PYTHON_OBJECTS'
# SAMPLE_NAME = 'M146'


def convert_sample(sample_name):
    sample_path = join(SAMPLES_PATH, sample_name)
    # Extract matrix
    mat = np.array(mmread(join(sample_path, "matrix.mtx.gz")).todense()).astype(np.uint16).T

    # Extract genes names and features
    features_path = join(sample_path, "features.tsv.gz")
    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]

    # Not in ues
    feature_types = [row[2] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]

    # Extract cells' barcodes
    barcodes_path = join(sample_path, "barcodes.tsv.gz")
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, 'rt'), delimiter="\t")]

    # Create sample object
    rna_sample = RNAseq_Sample(mat, gene_names, barcodes, feature_ids)
    return rna_sample


if __name__ == '__main__':
    if not os.path.isdir(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    samples = sorted([subfolder for subfolder in os.listdir(SAMPLES_PATH)])
    for i, sample_name in enumerate(samples):
        print(f'Working on {sample_name} ({i+1}/{len(samples)})')
        rna_sample = convert_sample(sample_name)
        pickle.dump((rna_sample), open(os.path.join(OUTPUT_PATH, f"{sample_name}.pkl"), "wb"))