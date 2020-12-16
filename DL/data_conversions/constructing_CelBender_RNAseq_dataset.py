import scanpy as sc
from  os.path import join
import sys
import os
import os
from os.path import join
import sklearn
from utilities.droplet_dataset import *
from utilities import *
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import random
from scipy.stats import pearsonr
from matplotlib.pyplot import figure
from termcolor import colored
import tables
from utilities.general_helpers import *

CELLBENDER_OUTPUT_PATH = r'D:\PycharmProjects\CellBender\server_output\parameterized_with_repairs'
RAW_SAMPLES_OUTPUTS = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples_10.12.20'
OUTPUT_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\CellBender\PARAMETERIZED_RUN_RNASEQ'



def copy_rna_object(rna_sample):
    """
    Creating a copy of RNAseq object
    :param rna_sample: Original sample
    :return: a copy.
    """
    copy_rna_sample = RNAseq_Sample(rna_sample.counts, rna_sample.gene_names, rna_sample.barcodes, rna_sample.features)
    # for c_inf1, c_inf2 in zip(copy_rna_sample.cells_information.cells_information_list,
    #                           rna_sample.cells_information.cells_information_list):
    #     for k, v in c_inf2.__dict__.items():
    #         if k in c_inf1.__dict__.keys():
    #             c_inf1.__setattr__(k, v)
    return copy_rna_sample


def constructing_scRNAseq_Dataset_keep_empty():
    """
    Here we create a droplet RNAseq dataset of samples after CellBender takes place as a preprocess.
    We use CellBender output as a get-go. Clean droplet'll be used while droplets that were classified as empty
    will be used too with the original values.
    PKLs will be saved in OUTPUT_PATH + empty_kept Dir.
    """
    samples = [subfolder for subfolder in os.listdir(RAW_SAMPLES_OUTPUTS) if not 'xlsx' in subfolder]
    KEEP_EMPTY_OUTPUT_PATH = join(OUTPUT_PATH, 'empty_kept')
    create_folder(KEEP_EMPTY_OUTPUT_PATH)
    for sample_id in samples:
        print(colored(sample_id, 'blue'))
        create_folder(join(KEEP_EMPTY_OUTPUT_PATH, sample_id))
        cb_path = join(CELLBENDER_OUTPUT_PATH, sample_id, 'out_filtered.h5')
        sample_path = join(RAW_SAMPLES_OUTPUTS, sample_id, f'{sample_id}.pkl')
        rna_sample = pickle.load(open(sample_path, 'rb'))

        # If CellBender process failed on this sample there is no cleaning matrix
        # we'll stay with old and non-clean background rna_seq sample.
        if not os.path.isfile(cb_path):
            print(colored(f'CellBender output for {sample_id} doesn\'t exist', 'red'))
            CB_rna_sample = copy_rna_object(rna_sample)
            pickle.dump((CB_rna_sample), open(os.path.join(KEEP_EMPTY_OUTPUT_PATH, sample_id, f"{sample_id}.pkl"), "wb"))
            continue

        # Retrieving samples and CellBender output.
        number_of_droplets = rna_sample.number_of_cells
        filterd_cb = sc.read_10x_h5(cb_path, genome='background_removed')

        # Creating rna copy with background removal changes
        cb_counts = filterd_cb.to_df().values
        cb_barcodes = list(filterd_cb.to_df().index)
        rna_sample.counts[[rna_sample.barcodes.index(b) for b in cb_barcodes], :] = cb_counts
        CB_rna_sample = copy_rna_object(rna_sample)

        # Marking empty droplets.
        empty_droplets_idxs = [idx for idx, _b in enumerate(rna_sample.barcodes) if not _b in cb_barcodes]
        for empty_droplet_idx in empty_droplets_idxs:
            CB_rna_sample.cells_information[empty_droplet_idx].is_CelBender_empty = True

        # Saving CellBender sample in RNAseq pkl object
        pickle.dump((CB_rna_sample ), open(os.path.join(KEEP_EMPTY_OUTPUT_PATH, sample_id, f"{sample_id}.pkl"), "wb"))



def constructing_scRNAseq_Dataset_drop_empty():
    """
    Here we create a droplet RNAseq dataset of samples after CellBender takes place as a preprocess.
    We use CellBender output as a get-go. Clean droplets'll be the used while droplets that were classified as empty
    will be removed.
    PKLs will be saved in OUTPUT_PATH + empty_kept Dir.
    """
    samples = [subfolder for subfolder in os.listdir(RAW_SAMPLES_OUTPUTS) if not 'xlsx' in subfolder]
    REMOVE_EMPTY_OUTPUT_PATH = join(OUTPUT_PATH, 'empty_removed')
    create_folder(REMOVE_EMPTY_OUTPUT_PATH)
    for sample_id in samples:
        print(colored(sample_id, 'blue'))
        create_folder(join(REMOVE_EMPTY_OUTPUT_PATH, sample_id))
        cb_path = join(CELLBENDER_OUTPUT_PATH, sample_id, 'out_filtered.h5')
        sample_path = join(RAW_SAMPLES_OUTPUTS, sample_id, f'{sample_id}.pkl')
        rna_sample = pickle.load(open(sample_path, 'rb'))

        # If CellBender process failed on this sample there is no cleaning matrix
        # we'll stay with old and non-clean background rna_seq sample.
        if not os.path.isfile(cb_path):
            print(colored(f'CellBender output for {sample_id} doesn\'t exist', 'red'))
            CB_rna_sample = copy_rna_object(rna_sample)

            pickle.dump((CB_rna_sample), open(os.path.join(REMOVE_EMPTY_OUTPUT_PATH, sample_id, f"{sample_id}.pkl"), "wb"))
            continue

        # Retrieving samples and CellBender output.
        number_of_droplets = rna_sample.number_of_cells
        filterd_cb = sc.read_10x_h5(cb_path, genome='background_removed')

        # Creating rna copy with background removal changes
        cb_counts = filterd_cb.to_df().values
        cb_barcodes = list(filterd_cb.to_df().index)
        # non_empty_droplets_idxs = [rna_sample.barcodes.index(b) for b in cb_barcodes]
        CB_rna_sample = RNAseq_Sample(cb_counts, rna_sample.gene_names, cb_barcodes, rna_sample.features)

        # Saving CellBender sample in RNAseq pkl object
        pickle.dump((CB_rna_sample), open(os.path.join(REMOVE_EMPTY_OUTPUT_PATH, sample_id, f"{sample_id}.pkl"), "wb"))


if __name__ == '__main__':
    constructing_scRNAseq_Dataset_keep_empty()
    constructing_scRNAseq_Dataset_drop_empty()