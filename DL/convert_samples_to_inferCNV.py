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


OUTPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\converted_data'
INPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\5.12.20'
GENOME_PATH = r'D:\Technion studies\Keren Laboratory\Data\inferCNV_data\gencode_v19_gene_pos.txt'


def convert_matrix(rna_sample, sample_id):

    gene_duplications = {g:False for g in rna_sample.gene_names}

    writer_path = join(OUTPUT_DIR, sample_id)
    create_folder(writer_path)
    with open(join(writer_path, 'matrix.matrix'), "wb") as writer:
        header = str.encode('\t'.join(rna_sample.barcodes)+'\n')
        writer.write(header)
        for idx, gene_name in enumerate(rna_sample.gene_names):
            if gene_duplications[gene_name]:
                continue
            gene_duplications[gene_name] = True
            print(f'{idx}/{len(rna_sample.gene_names)}')
            values = '\t'.join(rna_sample.counts[:, idx].astype(str).tolist())
            line = gene_name+'\t'+values + '\n'
            writer.write(str.encode(line))


    _breakpoint = 0


def create_annotations(rna_sample, sample_id):
    immune_indexes = [idx for idx, ci in enumerate(rna_sample.cells_information) if ci.is_immune]
    cancer_indexes = [idx for idx, ci in enumerate(rna_sample.cells_information) if not ci.is_immune]

    immune_barcodes = [rna_sample.barcodes[i] for i in immune_indexes]
    cancer_barcodes = [rna_sample.barcodes[i] for i in cancer_indexes]

    writer_path = join(OUTPUT_DIR, sample_id)
    create_folder(writer_path)
    with open(join(writer_path, 'annotation.txt'), "wb") as writer:
        for barcode in immune_barcodes:
            st = str.encode(barcode + '\t' + 'immune' + '\n')
            writer.write(st)

        for barcode in cancer_barcodes:
            st = str.encode(barcode + '\t' + 'cancer' + '\n')
            writer.write(st)


def convert_all_samples():
    samples = [subfolder for subfolder in os.listdir(INPUT_DIR)]

    # Extract ImmuneCellsMarkersUpdated Excel file from PC and load it into DataFrame.

    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    for sample_id in samples:
        print(sample_id)
        # Extracts one of the samples from PC
        sample_path = join(INPUT_DIR, sample_id, f'{sample_id}.pkl')
        rna_sample = extract_droplet_data_from_pickle(sample_path)

        convert_matrix(rna_sample, sample_id)
        create_annotations(rna_sample, sample_id)


if __name__ == '__main__':

    convert_all_samples()
    # with open(GENOME_PATH, 'rb') as f:
    #     genome = f.readlines()
    #     genome = [str(gene).split('\\t')[0][2:] for gene in genome]
    _breakpoint = 0
