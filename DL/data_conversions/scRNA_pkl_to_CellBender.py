"""
Use that script to create input to CellBender.
You should define the input path of the samples in PKL format.
"""

import sklearn
from utilities.droplet_dataset import *
from utilities import *
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib.pyplot as plt
import pickle
import random
from scipy.stats import pearsonr
from matplotlib.pyplot import figure
import pandas as pd
import os.path as path
from DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from scipy.io import mmwrite

OUTPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\CellBender\dummy'
INPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\10.12.20'




def create_gene_tsv(rna_sample, sample_id):
    path = join(OUTPUT_DIR, sample_id, 'genes.tsv')
    with open(path, 'wb') as writer:
        for idx in range(rna_sample.number_of_genes):
            line = rna_sample.features[idx]+'\t'
            line += rna_sample.gene_names[idx] + '\t'
            line += 'Gene Expression\n'
            writer.write(str.encode(line))


def create_barcode_tsv(rna_sample, sample_id):
    path = join(OUTPUT_DIR, sample_id, 'barcodes.tsv')
    with open(path, 'wb') as writer:
        for idx in range(rna_sample.number_of_cells):
            line = rna_sample.barcodes[idx] + '\n'
            writer.write(str.encode(line))


def create_mtx(rna_sample, sample_id):
    path = join(OUTPUT_DIR, sample_id, 'matrix.mtx')
    with open(path, 'wb') as writer:
        writer.write(str.encode('%%MatrixMarket matrix coordinate integer general\n'))
        writer.write(str.encode('%metadata_json: {"format_version": 2, "software_version": "3.1.0"}\n'))
        writer.write(str.encode(f'{rna_sample.number_of_genes} {rna_sample.number_of_cells} {len(np.where(rna_sample.counts!=0)[0])}\n'))
        for cell_idx in range(rna_sample.number_of_cells):
            for gene_idx in range(rna_sample.number_of_genes - 1, -1, -1):
                if rna_sample.counts[cell_idx, gene_idx]:
                    line = f'{gene_idx+1} {cell_idx+1} {rna_sample.counts[cell_idx, gene_idx]}\n'
                    writer.write(str.encode(line))


def convert_all_samples():
    samples = [subfolder for subfolder in os.listdir(INPUT_DIR)]

    create_folder(OUTPUT_DIR)
    for sample_id in samples:
        print(sample_id)
        create_folder(join(OUTPUT_DIR, sample_id))
        # Extracts one of the samples from PC
        sample_path = join(INPUT_DIR, sample_id, f'{sample_id}.pkl')
        rna_sample = extract_droplet_data_from_pickle(sample_path)

        # create matrix file.
        create_mtx(rna_sample, sample_id)

        # create gene.tsv file.
        create_gene_tsv(rna_sample, sample_id)

        # create barcodes.tsv file.
        create_barcode_tsv(rna_sample, sample_id)


if __name__ == '__main__':
    convert_all_samples()

    # path = r'D:\Technion studies\Keren Laboratory\Data\CellBender_data\barcodes.tsv'
    # with open(path, 'r') as f:
    #     lines = f.readlines()
    # path = join(OUTPUT_DIR, r'barcodes.tsv')
    # with open(path, 'r') as f:
    #     lines2 = f.readlines()
    #
    # _breakpoint = 0



