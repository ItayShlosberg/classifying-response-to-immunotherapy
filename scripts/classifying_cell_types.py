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
from os.path import join

PROJECT_PATH = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy'
MARKERS_PATH = path.join(PROJECT_PATH, r'Data\ImmuneCellsMarkersUpdated.xlsx')
CANCER_MARKERS = {'tumor':  ['MLANA', 'PMEL', 'TYR', 'MITF', 'AXL']}
LYMPHOID = ['T cells', 'CD4 helper T cells', 'CD8 Cytotoxic T cells', 'Regulatory T cells', 'Regulatory CD4 T cells', 'Regulatory CD8 T cells', 'Regulatory CD4_CD8 T cells', 'NKT cells', 'NK cells', 'B cells', 'Activated T cells', 'Senescence T cells', 'Terminal effector', 'Exhausted T cells', 'Stem_like T cells', 'Memory T cells', 'Memory CD4 T cells', 'Memory CD8 T cells', 'Memory CD4_CD8 T cells']
MYELOID = [['Macrophage_immature', 'Macrophage_mature', 'Monocyte_immature', 'Monocyte_mature', 'cDCs_dendritic_cells', 'pDCs', 'myeloid cells_general_immature', 'myeloid cells_general_mature', 'Neutrophils', 'Granolocytes',]]
SAMPLES_PATH = fr'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples'
PKL_NAME = r'RNA_sample.pkl'
OUT_FOLDER = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\output files\runs'

def get_length_of_mapping_table(mapping_table):
    accumulated_length = 0
    for key, val in mapping_table.items():
        if len(val)>0:
            accumulated_length += 1
    return accumulated_length


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
        # print(curr_cell_type)
        # print(curr_markers)
        possible_combinations_of_markers = get_all_possible_combinations_of_markers(curr_cell_type, curr_markers)
        for curr_combination_markers in possible_combinations_of_markers:
            curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names, curr_combination_markers)
            counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
            if counts_only_interesting_genes.sum() > 0:
                cells_satisfying_markers = np.all(counts_only_interesting_genes > 0, axis=1)
            else:
                cells_satisfying_markers = np.zeros(rna_sample.counts.shape[0], dtype=bool)
            cells_satisfying_markers_indexs = np.where(cells_satisfying_markers)
            # print(cells_satisfying_markers_indexs[0].shape)

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
        # print(curr_cell_type)
        # print(curr_markers)

        curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names, curr_markers)
        counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
        cells_satisfying_markers = np.any(counts_only_interesting_genes > 0, axis=1)
        cells_satisfying_markers_indexs = np.where(cells_satisfying_markers)
        # print(cells_satisfying_markers_indexs[0].shape)

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


def finds_pos_neg_conflicts(positive_cell_types_mapping_table, negative_cell_types_mapping_table):
    conflict_indicators = np.zeros((len(positive_cell_types_mapping_table)),
                                   dtype=bool)  # True if there is a conflict in place i, False otherwise
    for idx, pair in enumerate(
            zip(positive_cell_types_mapping_table.values(), negative_cell_types_mapping_table.values())):
        pos_classification = pair[0]
        neg_classification = pair[1]
        if len([cls for cls in pos_classification if cls in neg_classification]) > 0:
            conflict_indicators[idx] = True
    return conflict_indicators


def finds_cancer_conflicts(positive_cell_types_mapping_table, cancer_mapping_table):
    """
    Should be run after cleaning negative markers classification.
    """
    conflict_indicators = np.zeros((len(positive_cell_types_mapping_table)),
                                   dtype=bool)  # True if there is a conflict in place i, False otherwise
    for idx, pair in enumerate(zip(positive_cell_types_mapping_table.values(), cancer_mapping_table.values())):
        pos_classification = pair[0]
        cancer_classification = pair[1]
        if len(pos_classification) > 0 and len(cancer_classification) > 0:
            conflict_indicators[idx] = True
    return conflict_indicators


def cleans_by_indicators(mapping_table, conflict_indicators):
    conflict_indexes = np.where(conflict_indicators)
    for idx in conflict_indexes[0]:
        mapping_table[idx] = []


def classifying_cell_type(sample_id):
    # Step 1: Extracts one of the samples from PC
    sample_path = join(SAMPLES_PATH, sample_id, PKL_NAME)
    rna_sample = extract_droplet_data_from_pickle(sample_path)
    rna_sample.counts = rna_sample.counts.T # TODO: fix that in a way that the RNA sample would be in the right shape in advance (cells are rows, genes are cols).
    print(f'sample id {sample_id}')
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')

    # Step 2: Builds positive/negative cell type marker table.
    xls = pd.ExcelFile(MARKERS_PATH) # Extract ImmuneCellsMarkersUpdated Excel file from PC and load it into DataFrame.
    positive_markers_df = pd.read_excel(xls, 'and_or')
    negative_markers_df = pd.read_excel(xls, 'none')
    positive_markers_table = builds_cell_type_markers_table(positive_markers_df)
    negative_markers_table = builds_cell_type_markers_table(negative_markers_df)

    # Step 3: maps cell-type markers to sample's count-table cells.
    positive_cell_types_mapping_table = maps_positive_cell_types_to_count_table(positive_markers_table, rna_sample)
    negative_cell_types_mapping_table = maps_negative_cell_types_to_count_table(negative_markers_table, rna_sample)
    cancer_cell_mapping_table = maps_negative_cell_types_to_count_table(CANCER_MARKERS, rna_sample)

    # Step 4: Cross-checks determinations and finds conflicts.
    pos_neg_conflicts = finds_pos_neg_conflicts(positive_cell_types_mapping_table, negative_cell_types_mapping_table)
    cleans_by_indicators(positive_cell_types_mapping_table, pos_neg_conflicts)
    # TODO: add myeloid<>lymphoid conflict here
    cancers_conflicts = finds_cancer_conflicts(positive_cell_types_mapping_table, cancer_cell_mapping_table)
    cleans_by_indicators(positive_cell_types_mapping_table, cancers_conflicts)
    cleans_by_indicators(cancer_cell_mapping_table, cancers_conflicts)

    return {'positive_cell_types_mapping_table': positive_cell_types_mapping_table,
            'cancer_cell_mapping_table': cancer_cell_mapping_table,
            'pos_neg_conflicts': pos_neg_conflicts,
            'cancers_conflicts': cancers_conflicts,
            'number_of_cells': rna_sample.number_of_cells}


def save_classification_summary_sample(sample,
                                       cell_mapping_table,
                                       cancer_cell_mapping_table,
                                       pos_neg_conflicts,
                                       cancers_conflicts):

    folder = join(OUT_FOLDER, sample)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    pickle.dump(cell_mapping_table, open(join(folder, 'cell_mapping_table.pkl'), "wb"))
    pickle.dump(cancer_cell_mapping_table, open(join(folder, 'cancer_cell_mapping_table.pkl'), "wb"))
    pickle.dump(pos_neg_conflicts, open(join(folder, 'pos_neg_conflicts.pkl'), "wb"))
    pickle.dump(cancers_conflicts, open(join(folder, 'cancers_conflicts.pkl'), "wb"))


def summary():
    summary_df = pd.DataFrame(columns=['sample name',
                                       'number of cells',
                                       'total classified cells',
                                       'number of cells classified immune',
                                       'number of cells classified cancer',
                                       'number of pos-neg markers conflicts',
                                       'number of cancer-immune conflicts'])
    samples = [subfolder for subfolder in os.listdir(SAMPLES_PATH)]
    for sample in samples:
        # Classify sample's cells
        result = classifying_cell_type(sample)
        cells_mapping_table = result['positive_cell_types_mapping_table']
        cancer_cell_mapping_table = result['cancer_cell_mapping_table']
        pos_neg_conflicts = result['pos_neg_conflicts']
        cancers_conflicts = result['cancers_conflicts']
        number_of_cells = result['number_of_cells']

        # Extract summary of classification.
        number_of_cells_classified_immune = get_length_of_mapping_table(cells_mapping_table)
        number_of_cells_classified_cancer = get_length_of_mapping_table(cancer_cell_mapping_table)
        number_of_posneg_markers_conflicts = np.sum(pos_neg_conflicts)
        number_of_cancer_immune_conflicts = np.sum(cancers_conflicts)

        # Append to DF
        summary_df = summary_df.append(pd.DataFrame([[sample,
                                                      number_of_cells,
                                                      (number_of_cells_classified_immune+number_of_cells_classified_cancer)/number_of_cells,
                                                      number_of_cells_classified_immune,
                                                      number_of_cells_classified_cancer,
                                                      number_of_posneg_markers_conflicts,
                                                      number_of_cancer_immune_conflicts]],
                                       columns=summary_df.columns))

        # Save result
        save_classification_summary_sample(sample,
                                           cells_mapping_table,
                                           cancer_cell_mapping_table,
                                           pos_neg_conflicts,
                                           cancers_conflicts)
    summary_df.to_csv(join(OUT_FOLDER, 'summary.csv'))
    return summary_df


if __name__ == '__main__':
    summary_df = summary()
