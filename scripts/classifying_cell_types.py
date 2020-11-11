"""
For droplet_seq csRNAseq 2020 data. Classifying cell_types by markers table in ImmuneCellsMarkersUpdated.xlsx.

There are Positive markers so that the cells that are suspected to be classified into them should have value greater
than 1 in order to satisfy the conditions and there are negative markers so that the cells that are suspected to be
classified into them should have value 1.
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

PROJECT_PATH = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy'
MARKERS_PATH = path.join(PROJECT_PATH, r'Data\ImmuneCellsMarkersUpdated.xlsx')


def dict_append(dic, plc, val):
    dic[plc] = dic[plc] + [val]


def find_indexes_of_markers_in_sample(_sample_genes, markers):
    return [idx1 for idx1, g1 in enumerate(_sample_genes) for idx2, g2 in enumerate(markers) if g1 == g2]


def get_all_possible_combinations_of_markers(cell_type_name, celltype_markers):
    def recursive_comnibation_building(leftover, all_combinations=[[]]):
        def meiosis_of_combination(and_markers, or_markers):
            new_combinations = []
            for or_marker in or_markers:
                for and_marker in and_markers:
                    new_combinations.append(and_marker + [or_marker])
            return new_combinations

        if len(leftover) == 0:
            return all_combinations
        all_currrent_markers = leftover[0]
        del leftover[0]
        all_combination = meiosis_of_combination(all_combinations, all_currrent_markers)
        return recursive_comnibation_building(leftover, all_combination)

    initial_marker_structure = [s.split(';') for s in celltype_markers]
    combinations = recursive_comnibation_building(initial_marker_structure)
    return combinations


def builds_cell_type_markers_table(df):
    """
    Using Excel file, can  be used both with positive and negative markers df.
    :return: Dictionary containing for each cell-type (keys) the corresponding markers (values).
    """

    # Reorganize the table, the column names actually accouns the first row of the DataFrame
    columns = [v if not 'Unnamed' in v else np.NAN for v in list(df.columns)]
    df.loc[28] = columns
    df = df.reindex([28] + list(range(0, 28))).reset_index(drop=True)

    # convert df to dictionary markers table
    rm_nan = lambda lst: list([vv for vv in lst if not vv is np.NAN])
    df_rows = df.values.tolist()
    markers_table = {v[0]: rm_nan(v[1:]) for v in df_rows}

    return markers_table


def maps_positive_cell_types_to_count_table(positive_markers_table, rna_sample):
    """
    Builds mapping-table of sample's count-table. TODO: check counts shape to see is cells are the rows or the cols.
    :param positive_markers_table: positive markers dictionary.
    :param rna_sample: in which you are interested to classify its counts table.
    :return: positive_cell_types_mapping_table, dictionary mapping for each cell its possible cell-types.
    """

    positive_cell_types_mapping_table = {v: [] for v in range(rna_sample.number_of_cells)}
    for curr_cell_type, curr_markers in positive_markers_table.items():
        print(curr_cell_type)
        print(curr_markers)

        possible_combinations_of_markers = get_all_possible_combinations_of_markers(curr_cell_type, curr_markers)
        for curr_combination_markers in possible_combinations_of_markers:
            curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names, curr_combination_markers)
            counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
            cells_satisfying_markers = np.all(counts_only_interesting_genes > 0, axis=1)
            cells_satisfying_markers_indexs = np.where(cells_satisfying_markers)
            print(cells_satisfying_markers_indexs[0].shape)

        [dict_append(positive_cell_types_mapping_table, cell, curr_cell_type) for cell in
         cells_satisfying_markers_indexs[0]]
    return positive_cell_types_mapping_table


def maps_negative_cell_types_to_count_table(negative_markers_table, rna_sample):
    """
    Builds mapping-table of sample's count-table. TODO: check counts shape to see is cells are the rows or the cols.
    :param positive_markers_table: positive markers dictionary.
    :param rna_sample: in which you are interested to classify its counts table.
    :return: positive_cell_types_mapping_table, dictionary mapping for each cell its possible cell-types.
    """

    negative_cell_types_mapping_table = {v: [] for v in range(rna_sample.number_of_cells)}
    for curr_cell_type, curr_markers in negative_markers_table.items():
        print(curr_cell_type)
        print(curr_markers)

        curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names, curr_markers)
        counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
        cells_satisfying_markers = np.any(counts_only_interesting_genes > 0, axis=1)
        cells_satisfying_markers_indexs = np.where(cells_satisfying_markers)
        print(cells_satisfying_markers_indexs[0].shape)

        [dict_append(negative_cell_types_mapping_table, cell, curr_cell_type) for cell in
         cells_satisfying_markers_indexs[0]]
    return negative_cell_types_mapping_table


def cross_check_mapping(rna_sample, positive_cell_types_mapping_table, negative_cell_types_mapping_table):
    n_overlapping = 0
    number_of_classified_cells = 0

    for i in range(rna_sample.number_of_cells):
        pos_celltypes = positive_cell_types_mapping_table[i]
        neg_celltypes = negative_cell_types_mapping_table[i]
        overlapping = intersection_of_lists(pos_celltypes, neg_celltypes)
        if len(pos_celltypes) > 0:
            number_of_classified_cells += 1
        if len(overlapping) > 0:
            n_overlapping += 1
    print(f'Number of cell with overlap between pos and neg markers: {n_overlapping}')
    print(f'number of classified cells: {number_of_classified_cells}')
    print(f'number of cells (in sample): {rna_sample.number_of_cells}')
    print(f'portion of classified cells in sample: {number_of_classified_cells / rna_sample.number_of_cells}')


def classifying_cell_type():
    # Step 1: Extracts one of the samples from PC
    sample_id = 'M123'
    sample_path = fr'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples\{sample_id}\RNA_sample.pkl'
    rna_sample = extract_droplet_data_from_pickle(sample_path)
    rna_sample.counts = rna_sample.counts.T # TODO: fix that in a way that the RNA sample would be in the right shape in advance (cells are rows, genes are cols).
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')

    # Step 2: Builds positive/negative cell type marker table.
    # Extract ImmuneCellsMarkersUpdated Excel file from PC and load it into DataFrame.
    xls = pd.ExcelFile(MARKERS_PATH)
    positive_markers_df = pd.read_excel(xls, 'and_or')
    negative_markers_df = pd.read_excel(xls, 'none')
    positive_markers_table = builds_cell_type_markers_table(positive_markers_df)
    negative_markers_table = builds_cell_type_markers_table(negative_markers_df)

    # Step 3: maps cell-type markers to sample's count-table cells.
    positive_cell_types_mapping_table = maps_positive_cell_types_to_count_table(positive_markers_table, rna_sample)
    negative_cell_types_mapping_table = maps_negative_cell_types_to_count_table(negative_markers_table, rna_sample)

    # Step 4: Cross-checks determinations between neg-markers and pos-markers and finds conflicts.
    cross_check_mapping(rna_sample, positive_cell_types_mapping_table, negative_cell_types_mapping_table)


if __name__ == '__main__':
    classifying_cell_type()

