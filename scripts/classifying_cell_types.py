"""
For droplet_seq csRNAseq 2020 data.
Classifying cell_types by markers table in ImmuneCellsMarkersUpdated.xlsx, and by cancer markers.

There are Positive markers so that the cells suspected to be classified into them should have value greater
than zero in order to satisfy the conditions and there are negative markers so that the cells that are suspected to be
classified into them should have value zero.

Outputs:
After running this script each sample the following files will be saved:
cell_mapping_table.pkl - dictionary mapping between cells and their associated cell-types after removing cell-type with conflicts.
cancer_cell_mapping_table.pkl - dictionary plays a rule as list of all cells classified as cancer.
pos_neg_conflicts.pkl - indicators array for all cells, which classified cells had a conflict with negative markers.
cancers_conflicts.pkl - indicators array for all cells, which cells had a conflict with cancer.
cells_with_neg_pos_conflict.csv - table mapping between cells with a negative markers conflict and the this list of negative markers.
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
MARKERS_PATH = path.join(PROJECT_PATH, r'Data\ImmuneCellsMarkersUpdated_12.11.20.xlsx')
CANCER_MARKERS = {'tumor':  ['MLANA', 'PMEL', 'TYR', 'MITF', 'AXL']}
LYMPHOID = ['T cells', 'CD4 helper T cells', 'CD8 Cytotoxic T cells', 'Regulatory T cells', 'Regulatory CD4 T cells', 'Regulatory CD8 T cells', 'Regulatory CD4_CD8 T cells', 'NKT cells', 'NK cells', 'B cells', 'Activated T cells', 'Senescence T cells', 'Terminal effector', 'Exhausted T cells', 'Stem_like T cells', 'Memory T cells', 'Memory CD4 T cells', 'Memory CD8 T cells', 'Memory CD4_CD8 T cells']
MYELOID = [['Macrophage_immature', 'Macrophage_mature', 'Monocyte_immature', 'Monocyte_mature', 'cDCs_dendritic_cells', 'pDCs', 'myeloid cells_general_immature', 'myeloid cells_general_mature', 'Neutrophils', 'Granolocytes',]]
SAMPLES_PATH = fr'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples'
PKL_NAME = r'RNA_sample.pkl'
OUT_FOLDER = r'D:\Technion\output files\runs\1.12.20'
MHC2_GENES = ['HLA-DM', 'HLA-DMA', 'HLA-DMB', 'HLA-DO',
             'HLA-DOA', 'HLA-DOB', 'HLA-DP', 'HLA-DPA1',
             'HLA-DPB1', 'HLA-DQ', 'HLA-DQA1', 'HLA-DQA2',
             'HLA-DQB1', 'HLA-DQB2', 'HLA-DR', 'HLA-DRA',
             'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']
# The genes that appear in the samples (optimization):
REDUCED_MHC2_GENES = ['HLA-DMA',  'HLA-DMB',  'HLA-DOA',  'HLA-DOB',
                      'HLA-DPA1',  'HLA-DPB1',  'HLA-DQA1',  'HLA-DQB1',
                      'HLA-DQB2',  'HLA-DRA',  'HLA-DRB1',  'HLA-DRB5']


def get_length_of_mapping_table(mapping_table):
    """
    Mapping table is a dictionary, where each key represents a index in sample count-table. That is, a cell.
    The function returns the number of classified cells, those whose value in the dictionary isn't empty.
    :param mapping_table:
    :return:
    """
    accumulated_length = 0
    for key, val in mapping_table.items():
        if len(val)>0:
            accumulated_length += 1
    return accumulated_length


def convert_MHC2_markers_tostr(markers):
    """
    In case where the marker "MHCII' is in the list of markers we'll replace it with the list of genes that
    make up MHCII. This list will be added as one string.
    :param markers: list of markers
    :return: The new markers list after converting it.
    """
    if 'MHCII' in markers:
        markers = [m for m in markers if m != 'MHCII'] + [';'.join(REDUCED_MHC2_GENES)]
    return markers


def convert_MHC2_markers_tolist(markers):
    """
    In case where the marker "MHCII' is in the list of markers we'll replace it with the list of genes that
    make up MHCII. This list will be added as list of string.
    :param markers: list of markers
    :return: The new markers list after converting it.
    """
    if 'MHCII' in markers:
        markers = [m for m in markers if m != 'MHCII'] + REDUCED_MHC2_GENES
    return markers


def dict_append(dic, plc, val):
    """
    Dictionary of lists. Appends value to one of the lists in given place in the dictionary (by key).
    :param dic: Dictionary of lists.
    :param plc: place in the dictionary (key).
    :param val: the values you want to aadd to the list.
    """
    dic[plc] = dic[plc] + [val]


def find_indexes_of_markers_in_sample(_sample_genes, markers):
    """
    find the indexes of markers in gene list.
    :param _sample_genes: list of gene names ordered according to the count-table order.
    :param markers: whom you want to know their indexes in the gene list.
    :return: list of indexes corresponding to given markers.
    """
    return [idx for g1 in markers for idx, g2 in enumerate(_sample_genes) if g1 == g2]


def get_all_possible_combinations_of_markers(celltype_markers):
    """
    While there are markers of a cell-type that have to be shown in order to classify a cell, the are other markers
    separated with ';', so that only having one of them can be enough in order to be sure  that cell corresponding to
    this cell-type. This function gives you all possible combinations of lists of markers, so that you can check each
    of those list separately and if even in one of them you find a match you can define classify this cell.
    :param celltype_markers: list of this cell-type markers. part of them can be a number of markers separated by ';'.
    :return: All combinations.
    """
    def recursive_combination_building(leftover, all_combinations=[[]]):
        """
        building recursively the combination-list.
        :param leftover: The remaining part of the list that you should add to the final output.
        :param all_combinations: the output that are being built recursively.
        :return: all_combinations when we reach the stop-condition.
        """
        def meiosis_of_combination(and_markers, or_markers):
            """
            Divides and_markers is actually the all_combination list being built. or_markers is the current marker
            (of list of markers separated by ';') and we want to add it for all combination lists.
            :return: new_combinations - updated and_markers list.
            """
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
        return recursive_combination_building(leftover, all_combination)

    initial_marker_structure = [s.split(';') for s in convert_MHC2_markers_tostr(celltype_markers)]
    combinations = recursive_combination_building(initial_marker_structure)
    return combinations


def builds_cell_type_markers_table(df):
    """
    Extract cell-types markers from Excel file. the function can be used both with positive and negative markers df.
    :return: Dictionary containing for each cell-type (keys) the corresponding markers (values).
    """

    # Reorganize the table, the column names actually accouns the first row of the DataFrame
    columns = [v if not 'Unnamed' in v else np.NAN for v in list(df.columns)]
    df.loc[len(df)] = columns
    df = df.reindex([len(df)-1] + list(range(0, len(df)-1))).reset_index(drop=True)

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
        possible_combinations_of_markers = get_all_possible_combinations_of_markers(curr_markers)
        for curr_combination_markers in possible_combinations_of_markers:
            curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names, curr_combination_markers)
            counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
            if counts_only_interesting_genes.sum() > 0: # In case of there is not appearances of those markers in sample's genes
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

        curr_genes_indexes = find_indexes_of_markers_in_sample(rna_sample.gene_names,
                                                               convert_MHC2_markers_tolist(curr_markers))
        counts_only_interesting_genes = np.take(rna_sample.counts, curr_genes_indexes, axis=1)
        cells_satisfying_markers = np.any(counts_only_interesting_genes > 0, axis=1)
        cells_satisfying_markers_indexs = np.where(cells_satisfying_markers)
        # print(cells_satisfying_markers_indexs[0].shape)

        [dict_append(negative_cell_types_mapping_table, cell, curr_cell_type) for cell in
         cells_satisfying_markers_indexs[0]]
    return negative_cell_types_mapping_table


def finds_pos_neg_conflicts(positive_cell_types_mapping_table, negative_cell_types_mapping_table):
    """
    Check if one of cells classified to one of the cell-type but also negatively classified to that cell-type.
    :param positive_cell_types_mapping_table: Positive mapping table - by markers from 'and_or' tab.
    :param negative_cell_types_mapping_table: Negative Mapping table - by markers from 'none' tab.
    :return: conflict_indicators - numpy array each place is an indicator says if the corresponding cell
    in the count table has a pos-neg conflict.
    conflict_list - in addition to the indicator list we want the names of the problematic cell-type with their markers
    are responsible to the conflict.
    """
    conflict_indicators = np.zeros((len(positive_cell_types_mapping_table)),
                                   dtype=bool)  # True if there is a conflict in place i, False otherwise
    conflict_list = {}
    for idx, pair in enumerate(
            zip(positive_cell_types_mapping_table.values(), negative_cell_types_mapping_table.values())):
        pos_classification = pair[0]
        neg_classification = pair[1]
        overlapping_list = list(set([cls for cls in pos_classification if cls in neg_classification]))
        if len(overlapping_list) > 0:
            conflict_indicators[idx] = True
            conflict_list[idx] = overlapping_list
    return conflict_indicators, conflict_list


def finds_cancer_conflicts(positive_cell_types_mapping_table, cancer_mapping_table):
    """
    Check if one of cells classified to one of the cell-type but also classified as cancer-cell.
    Note: Should be run after cleaning negative markers classification.
    :param positive_cell_types_mapping_table: Positive mapping table - by markers from 'and_or' tab. After removing
    cells having a conflict with negative markers.
    :param cancer_mapping_table:
    :return: conflict_indicators.
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
    """
    Remove completely the classification from cells having a conflict.
    Note: Should be very cautious using this function. not always you want to remove all the classification having
    only a conflict corresponding to one of the classifications.
    :param mapping_table: between cells and their cell-type classifications (can be tumor).
    :param conflict_indicators: boolean numpy array indicates if it's a conflict in each place.
    """
    conflict_indexes = np.where(conflict_indicators)
    for idx in conflict_indexes[0]:
        mapping_table[idx] = []


def cleans_by_overlapping_cell_type_names(mapping_table, conflict_list):
    """
    Remove the classification from cells having a conflict. Only the cell-types that are responsible to the conflict.
    :param mapping_table: between cells and their cell-type classifications (can be tumor).
    :param conflict_list: dict mapping between cells and their list of cell-types with conflict.
    :return: cell_types_removed - a dictionary that contains for each cell which had a cell-type that removed
    its list of cell-types were removed.
    """
    cell_types_removed = {}
    for key, values in conflict_list.items():
        cell_types_removed[key] = [cell_type for cell_type in mapping_table[key] if cell_type in values]
        mapping_table[key] = list(set([cell_type for cell_type in mapping_table[key] if not cell_type in values]))
    return cell_types_removed


def extract_sample(sample_id):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    sample_path = join(SAMPLES_PATH, sample_id, PKL_NAME)
    rna_sample = extract_droplet_data_from_pickle(sample_path)
    print(f'sample id {sample_id}')
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')
    return rna_sample


def classifying_cell_type(rna_sample, positive_markers_df, negative_markers_df):
    """
    Classifies all cells in one sample.
    function steps:
    1. Builds positive/negative cell type marker table.
    2. maps cell-type markers to sample's count-table cells.
    3. Cross-checks positive and negative markers tables and finds conflicts.
    4. Cleans conflicts of positive markers with negative markers.
    5. Cross-checks positive and negative markers tables and finds conflicts.
    6. Cleans cancer conflicts from both of the tables.
    :param sample_id: The is of the sample you want to classify its cells.
    :return: positive cell-types mapping table - removing classification if cell classified also as cancer, and removed
                                                 cell-types classification with a negative conflict.
            cancer cell mapping table - boolean indicator numpy array.
            pos_neg_conflicts - boolean indicator numpy array.
            pos_neg_conflict_list - dictionary mapping cell index and list of cell-types with conflict.
            cancers_conflicts - boolean indicator numpy array.
            number of cells in sample - int.
    """
    # Step 1: Builds positive/negative cell type marker table.
    positive_markers_table = builds_cell_type_markers_table(positive_markers_df)
    negative_markers_table = builds_cell_type_markers_table(negative_markers_df)

    # Step 2: maps cell-type markers to sample's count-table cells.
    positive_cell_types_mapping_table = maps_positive_cell_types_to_count_table(positive_markers_table, rna_sample)
    negative_cell_types_mapping_table = maps_negative_cell_types_to_count_table(negative_markers_table, rna_sample)
    cancer_cell_mapping_table = maps_negative_cell_types_to_count_table(CANCER_MARKERS, rna_sample)

    # Step 3: Cross-checks positive and negative markers tables and finds conflicts.
    pos_neg_conflict_indicators, pos_neg_conflict_list = finds_pos_neg_conflicts(positive_cell_types_mapping_table, negative_cell_types_mapping_table)

    # Step 4: Cleans conflicts of positive markers with negative markers.
    cell_types_removed = cleans_by_overlapping_cell_type_names(positive_cell_types_mapping_table, pos_neg_conflict_list)
    # TODO: add myeloid<>lymphoid conflict here
    # Step 5: Cross-checks positive and negative markers tables and finds conflicts.
    cancers_conflicts = finds_cancer_conflicts(positive_cell_types_mapping_table, cancer_cell_mapping_table)
    # Step 6: Cleans cancer conflicts from both of the tables.
    cleans_by_indicators(positive_cell_types_mapping_table, cancers_conflicts)
    cleans_by_indicators(cancer_cell_mapping_table, cancers_conflicts)

    return {'positive_cell_types_mapping_table': positive_cell_types_mapping_table,
            'cancer_cell_mapping_table': cancer_cell_mapping_table,
            'pos_neg_conflicts': pos_neg_conflict_indicators,
            'pos_neg_conflict_list': pos_neg_conflict_list,
            'cancers_conflicts': cancers_conflicts,
            'number_of_cells': rna_sample.number_of_cells,
            'cell_types_removed': cell_types_removed}


def save_classification_summary_sample(sample_id,
                                       number_of_cells,
                                       cell_mapping_table,
                                       cancer_cell_mapping_table,
                                       pos_neg_conflict_indicators,
                                       cancers_conflict_indicators,
                                       pos_neg_conflict_list,
                                       rna_sample):
    """
    Makes folder for current sample in OUT_FOLDER, saving PKLs of cell_mapping_table, cancer_cell_mapping_table,
    pos_neg_conflicts, cancers_conflicts. And CSV file for pos_neg_conflict_list.
    :param: All outputs from classifying_cell_type function.
    """
    folder = join(OUT_FOLDER, sample_id)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    pickle.dump(cell_mapping_table, open(join(folder, 'cell_mapping_table.pkl'), "wb"))
    pickle.dump(cancer_cell_mapping_table, open(join(folder, 'cancer_cell_mapping_table.pkl'), "wb"))
    pickle.dump(pos_neg_conflict_indicators, open(join(folder, 'pos_neg_conflicts.pkl'), "wb"))
    pickle.dump(cancers_conflict_indicators, open(join(folder, 'cancers_conflicts.pkl'), "wb"))

    conflict_df = pd.DataFrame([[x] + [';'.join(y)] for x, y in list(pos_neg_conflict_list.items())],
                               columns=['cell idx', 'problematic cell-types'])

    # mapping_table_df = pd.DataFrame([[key, ';'.join(val)] for key, val in cell_mapping_table.items()],
    #                                 columns=['cell index', 'values'])
    conflict_df.to_csv(join(folder, 'cells_with_neg_pos_conflict.csv'))

    # save also the updated rna object that now contains all the classification information.
    pickle.dump(rna_sample, open(join(folder, f'{sample_id}.pkl'), 'wb'))


def update_rna_sample_cells_findings(rna_sample,
                                     cells_mapping_table,
                                     cell_types_removed,
                                     cancer_cell_mapping_table,
                                     cancers_conflicts):
    cells_information = [Cell_information() for _ in range(rna_sample.number_of_cells)]

    # update cells' cell-type list
    for k, v in cells_mapping_table.items():
        cells_information[k].cell_type_list = v
        if len(v):
            cells_information[k].is_classified = True

    # update cell-types related conflict, cell with conflict and without classification
    # means it could have been classified if that conflict hadn't happened.
    for k, v in cell_types_removed.items():
        cells_information[k].conflict_related_cell_types = v
        if not cells_information[k].is_classified:
            cells_information[k].could_have_been_classified = True

    for k, v in cancer_cell_mapping_table.items():
        if 'tumor' in v:
            cells_information[k].is_cancer = True

    for idx in np.where(cancers_conflicts)[0]:
        cells_information[idx].cancer_immune_conflict = True

    rna_sample.cells_information = cells_information
    return cells_information


def summary_over_all_samples():
    """
    Loop over all sample and classify all cells. Saves results and summary.
    :return: Summary.
    """
    summary_df = pd.DataFrame(columns=['sample name',
                                       'number of cells',
                                       'total classified cells',
                                       'number of cells classified as immune',
                                       'number of cells classified as cancer',
                                       'number of pos-neg markers conflicts',
                                       'number of cancer-immune conflicts',
                                       'number of cancer-immune conflicts without tag'])
    samples = [subfolder for subfolder in os.listdir(SAMPLES_PATH)]

    # Extract ImmuneCellsMarkersUpdated Excel file from PC and load it into DataFrame.
    xls = pd.ExcelFile(MARKERS_PATH)
    positive_markers_df = pd.read_excel(xls, 'and_or')
    negative_markers_df = pd.read_excel(xls, 'none')
    all_samples_cells_information = {}
    if not os.path.isdir(OUT_FOLDER):
        os.mkdir(OUT_FOLDER)
    for sample in samples:
        # Extracts one of the samples from PC
        rna_sample = extract_sample(sample)

        # Classify sample's cells
        result = classifying_cell_type(rna_sample, positive_markers_df, negative_markers_df)
        cells_mapping_table = result['positive_cell_types_mapping_table']
        cancer_cell_mapping_table = result['cancer_cell_mapping_table']
        pos_neg_conflicts = result['pos_neg_conflicts']
        pos_neg_conflict_list = result['pos_neg_conflict_list']
        cancers_conflicts = result['cancers_conflicts']
        number_of_cells = result['number_of_cells']
        cell_types_removed = result['cell_types_removed']



        # Update RNA sample object: cells meta-data for each cells its findings
        all_samples_cells_information[sample] = update_rna_sample_cells_findings(rna_sample,
                                                                                 cells_mapping_table,
                                                                                 cell_types_removed,
                                                                                 cancer_cell_mapping_table,
                                                                                 cancers_conflicts)

        # Extract summary of classification.
        number_of_cells_classified_immune = get_length_of_mapping_table(cells_mapping_table)
        number_of_cells_classified_cancer = get_length_of_mapping_table(cancer_cell_mapping_table)
        number_of_posneg_markers_conflicts = len(cell_types_removed)
        number_of_cancer_immune_conflicts = np.sum(cancers_conflicts)
        # Could have been classified if there was not a conflict.
        number_of_cells_could_have_been_classified = sum(
            [c.could_have_been_classified for c in all_samples_cells_information[sample]])


        # Append to DF
        summary_df = summary_df.append(pd.DataFrame([[sample,
                                                      number_of_cells,
                                                      (number_of_cells_classified_immune+number_of_cells_classified_cancer)/number_of_cells,
                                                      number_of_cells_classified_immune,
                                                      number_of_cells_classified_cancer,
                                                      number_of_posneg_markers_conflicts,
                                                      number_of_cancer_immune_conflicts,
                                                      number_of_cells_could_have_been_classified]],
                                       columns=summary_df.columns))

        # Save sample result
        save_classification_summary_sample(sample,
                                           number_of_cells,
                                           cells_mapping_table,
                                           cancer_cell_mapping_table,
                                           pos_neg_conflicts,
                                           cancers_conflicts,
                                           pos_neg_conflict_list,
                                           rna_sample)

    summary_conflict_related_cell_types(all_samples_cells_information).to_csv(join(OUT_FOLDER, 'conflict_related_cell_types_df.csv'))
    summary_df.to_csv(join(OUT_FOLDER, 'summary.csv'))
    return summary_df


def summary_conflict_related_cell_types(all_samples_cells_information):
    """
    Build a table of all conflict related cell types of cells that don't have any class after conflict removed.
    :param all_samples_cells_information: list of Cell_information python object.
    :return: conflict_related_cell_types_df
    """
    conflict_related_cell_types_df = pd.DataFrame(columns=['sample name',
                                       'cell index',
                                       'cell-types removed'])
    for sample_id, cells_information in all_samples_cells_information.items():
        for cell_index, cell_inf in enumerate(cells_information):
            if cell_inf.could_have_been_classified:
                conflict_related_cell_types_df = conflict_related_cell_types_df.append(pd.DataFrame([[sample_id,
                                                                                                      cell_index,
                                                                                                      ';'.join(cell_inf.conflict_related_cell_types)]],
                                                                                                    columns=[
                                                                                                        'sample name',
                                                                                                        'cell index',
                                                                                                        'cell-types removed']))
    return conflict_related_cell_types_df


if __name__ == '__main__':
     summary_df = summary_over_all_samples()