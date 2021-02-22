"""
After drawing conclusion with inferCNV and  classifying tumor cells and immune cells using the CNV maps
run this script to update the data and save a version of PKL of the data with the conclusion.


---------------------------------------------
For tumor and not immune cells (the lower-side of the cnv map):
dying cells will remain dying cells.
tumor cells will remain tumor.
cells that weren't tumor but found to be tumor now, will be marked as tumor.
cells with immune and cancer conflict will be removed and the cells that grouped together with them also.

Method (only for tumor):
# Reading xlsx table, for every sample we have a division of the cells into fragments. Each fragment (a group) of
cells got a classification of what we should do with the cells of that fragment (type table found in the excel also).

---------------------------------------------
For immune cells (the up-side of the cnv map):
cells that were found to have some signatures suspected to be related to tumor clusters will be removed.

Method (only for immune):
# We should have a folder that contains PKL files only for samples that should undergo removal of cells.
each PKL file contains all the barcode that **should be** removed! if a barcode isn't shown up it means that the cell is
an immune cell without problematic issues.

---------------------------------------------
########## Tumor cluster types ##########
0	stromal cells - conflicts will be classified as stromal, cancer cells remain cancer
1	stromal cells - conflicts will be removed, cancer cells remain cancer
2	tumor - conflicts will be classified as tumor
3	tumor - conflicts will be classified as stromal	otherwise become tumor
4	tumor - conflicts will be removed, otherwise become tumor
5	remove all - suspected to be doublets (in apoptosis case may be in use)


########## Immune cluster types ##########
0	Mannualy - some cluster will be remove	separated cluster - the pickle will contain cell for removing
1	Gaussian - some cluster will be removed	separated cluster - the pickle will contain cell for removing
2	No cluster will be removed	Valid - no pickle
"""

import pickle
import pandas as pd
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from utilities.general_helpers import *
from termcolor import colored


# path for samples which will be used to update. Important: taking the last-updated scrublet output. (after all QC processes).
ROW_SAMPLES_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\21.2.21'
INFERCNV_SAMPLES_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\executions\all_data_31.12.20'
# path of folder where all samples having cell needed be removed have PKL file containing all barcodes of the cell needed be removed.
IMMUNE_CELLS_REMOVAL_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\analysis_conclusions\immune_clustering'
# In that path all pkl of the updated properties will be saved.
OUTPUT_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'
# Tumor table contains row for each sample splioting the tumor cells (Not immune cells) into cluster ##sorted by InferCNV output##
TUMOR_TABLE_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\analysis_conclusions\tumor_classifying_clusters.xlsx'
# Immune table contains row for each sample that have processed, if there are cells needed to be removed it's indicated in cluster-type.
IMMUNE_TABLE_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\analysis_conclusions\immune_classifying_clusters.xlsx'


def extract_sample(sample_id):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    data_path = join(ROW_SAMPLES_PATH, sample_id, f'{sample_id}.pkl')
    rna_sample = extract_droplet_data_from_pickle(data_path)
    print(colored(f'sample id {sample_id}', 'blue'))
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')
    return rna_sample


def extract_tumor_table():
    tumor_df = pd.read_excel(TUMOR_TABLE_PATH)
    tumor_df = tumor_df[list(tumor_df.columns)[:15]]
    return tumor_df


def extract_cnv_tumor_barcodes(sample_id):
    tumor_file = r'infercnv.observations.txt'  # tumor cells

    with open(join(INFERCNV_SAMPLES_PATH, sample_id, tumor_file), 'r') as f:
        line = f.readline()
        barcodes_length = len(line[:-1].split(' '))
        print(f'number of barcodes {barcodes_length}')
        barcodes = [ii.replace('\"', '') for ii in line.split(' ')]
        barcodes[-1] = barcodes[-1][:-1]
    return barcodes


def extract_immune_table():
    immune_df = pd.read_excel(IMMUNE_TABLE_PATH)
    immune_df = immune_df[list(immune_df.columns)[:2]]
    return immune_df


def perform_tumor_clustering_procedure_of_one_sample():
    pass


def perform_immune_clustering_procedure_of_one_sample():
    pass


def rearrange_tumor_clusters(tumor_clusters, num_of_cells):
    rearranged_clusters = [(convert_cluster_str_to_int(tumor_clusters[i*2]), tumor_clusters[i*2+1]) for i in range(int(len(tumor_clusters)/2))
            if tumor_clusters[i*2]!='-']

    # if there is no clusters it means that all the cells are stromal
    # so we'll build one cluster of all indexes and mark it '1' (the code of stromal cells).
    if not len(rearranged_clusters):
        rearranged_clusters.append(((0, num_of_cells), 1))

    return rearranged_clusters


def convert_cluster_str_to_int(tumor_clusters):
    int_tumor_clusters = tuple([int(jj) for jj in ''.join(ii for ii in tumor_clusters if ii not in ['(', ')', ' ']).split(',')])
    return int_tumor_clusters


def update_tumor_cells(sample , tumor_clusters, tumor_barcodes_order):
    """

    Tumor cluster types:
    0	stromal cells - conflicts will be classified as stromal, cancer cells remain cancer
    1	stromal cells - conflicts will be removed, cancer cells remain cancer
    2	tumor - conflicts will be classified as tumor
    3	tumor - conflicts will be classified as stromal	otherwise become tumor
    4	tumor - conflicts will be removed, otherwise become tumor
    5	remove all - suspected to be doublets (in apoptosis case may be in use)
    :param sample:
    :param tumor_clusters:
    :param tumor_barcodes_order:
    :return:
    """
    sample.cells_information.setattr('is_stromal', None, False)
    subsets = [(sample.get_subset_by_barcodes(tumor_barcodes_order[cluster[0][0]: cluster[0][1]]), cluster[1])
               for cluster in tumor_clusters]

    for cluster, cluster_type in subsets:
        bool_is_cancer_indices = cluster.cells_information.getattr('is_cancer')
        bool_is_conflict_indices = cluster.cells_information.getattr('cancer_immune_conflict')
        bool_is_apoptosis = cluster.cells_information.getattr('is_apoptosis')

        if cluster_type == 1:
            # stromal cells are those which are not cancer and don't have conflicts.
            bool_stromal_indices = [not bool_is_cancer_indices[ii] and not bool_is_conflict_indices[ii] for ii
                                    in range(len(bool_is_cancer_indices))]
            stromal_cells_RNAseq = cluster[bool_stromal_indices]
            stromal_cells_RNAseq.cells_information.setattr('is_stromal', None, True)

        elif cluster_type == 2:
            bool_not_apoptosis_indices = [not val for val in bool_is_apoptosis]
            moving_to_be_cancer_cells = cluster[bool_not_apoptosis_indices]
            moving_to_be_cancer_cells.cells_information.setattr('is_cancer', None, True)

        elif cluster_type == 4:
            moving_to_be_cancer_cells_indices = [not bool_is_apoptosis[ii] and not bool_is_conflict_indices[ii] for ii
                                    in range(len(bool_is_apoptosis))]
            moving_to_be_cancer_cells = cluster[moving_to_be_cancer_cells_indices]
            moving_to_be_cancer_cells.cells_information.setattr('is_cancer', None, True)

        else:
            print(colored(f"There is use in cluster type: {cluster_type}, end there is no reference"))
            raise Exception(f"There is use in cluster type: {cluster_type}, end there is no reference")

        # mark cells for removal:
        bool_is_cancer_indices = cluster.cells_information.getattr('is_cancer') # Now we might have more cells marked tumor
        bool_is_stromal_indices = cluster.cells_information.getattr('is_stromal')

        # if there is no tumor or stromal flag we will remove that cell (it's dying or garbage)
        cells_for_removal_indices = [not bool_is_cancer_indices[ii] and not bool_is_stromal_indices[ii]
                                     for ii in range(len(bool_is_cancer_indices))]
        # print(1)
        if sum(cells_for_removal_indices):
            cells_for_removal = cluster[cells_for_removal_indices]
            cells_for_removal.cells_information.setattr('should_be_removed', None, True)


def update_immune_cells(sample_id, rna_sample, immune_df):
    """ ##### IMMUNE REMOVAL TYPE #####
    0	Mannualy - some cluster will be remove	seperated cluster - the pickle will contain cell for removing
    1	Gaussian - some cluster will be removed	seperated cluster - the pickle will contain cell for removing
    2	No cluster will be removed	Valid - no pickle

    :param sample_id:
    :param rna_sample:
    :param immune_df:
    :return:
    """
    cluster_type = list(immune_df[immune_df['sample'] == sample_id].iloc[0])[1]
    if cluster_type == 2:
        return
    immune_removal_barcodes = pickle.load(open(join(IMMUNE_CELLS_REMOVAL_PATH, f'{sample_id}.pkl'), 'rb'))
    rna_sample.get_subset_by_barcodes(immune_removal_barcodes).cells_information.setattr('should_be_removed', None,
                                                                                         True)


def go_over_all_samples(tumor_df, immune_df):
    samples = [subfolder for subfolder in os.listdir(ROW_SAMPLES_PATH)]
    create_folder(OUTPUT_PATH)
    for sample_id in [s for s in samples if (not 'csv' in s and not 'xlsx' in s)]:
        # Extracts one of the samples from PC
        rna_sample = extract_sample(sample_id)
        rna_sample.cells_information.setattr('should_be_removed', None, False)

        # Update tumor cells
        tumor_clusters = list(tumor_df[tumor_df['sample'] == sample_id].iloc[0])[1:]
        tumor_clusters = rearrange_tumor_clusters(tumor_clusters, rna_sample.number_of_cells)
        tumor_barcodes = extract_cnv_tumor_barcodes(sample_id)
        # perform_tumor_clustering_procedure_of_one_sample()
        update_tumor_cells(rna_sample, tumor_clusters, tumor_barcodes)


        # Update immune cells
        if sample_id in list(immune_df['sample']):
            update_immune_cells(sample_id, rna_sample, immune_df)

        # mark all apoptosis and doublets (scrublet) cells for removal
        rna_sample[rna_sample.cells_information.getattr('is_apoptosis')].cells_information.setattr('should_be_removed', None, True)
        doublets_indices = rna_sample.cells_information.getattr('is_doublet')
        if (sum(doublets_indices)/rna_sample.number_of_cells) < 0.1:
            rna_sample[doublets_indices].cells_information.setattr('should_be_removed', None, True)

        # Save an updated version of current sample_id with inferCNV changes.
        pickle.dump((rna_sample), open(join(OUTPUT_PATH, f'{sample_id}.pkl'), 'wb'))


if __name__ == '__main__':
    # if not os.path.isdir(OUTPUT_PATH):
    #     os.mkdir(OUTPUT_PATH)
    tumor_df = extract_tumor_table()
    immune_df = extract_immune_table()

    go_over_all_samples(tumor_df, immune_df)
    print("hello")
    _breakpoint = 0