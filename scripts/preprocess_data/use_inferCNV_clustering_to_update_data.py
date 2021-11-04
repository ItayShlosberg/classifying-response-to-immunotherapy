"""
After drawing conclusion with inferCNV and classifying tumor cells and immune cells using the CNV maps
run this script to update the data and save a version of PKL of the data with the conclusions.



This is the last step in QC and as such, we update here other properties:
- We mark for removal: Dying cells (found to be in apoptosis).
- We mark for removal: doublets (classified by Scrublet).
Update 6.5.21:
- We mark for removal: immune cells having myeloid & lymphoid conflicts.
- We mark for removal: we decided to remove all tumor cells with conflicts even those which have CNVs patterns.
- specify in cells_information whether a cell was removed in CellBender.
- specify contaminated cells in immune compartment (CSV)- moved to stromal compartment
Update 24.5.21:
- specify contaminated cells in immune compartment (CSV, kmeans - cluster #10, in k=10) - moved to epithelial compartment
update 4.11.21:
-  we usually remove samples with a low number of genes (<500), but only if they do not have neutrophil markers.
   Neutrophil are usually small and if we remove all cells with <500 genes we will lose them.


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
from utilities.droplet_dataset import loading_sample
from utilities.general_helpers import intersection_of_lists

# In that path all pkl of the updated properties will be saved.
# OUTPUT_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'
OUTPUT_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\inferCNV\update_runs\4.11.21'

# path for samples which will be used to update. Important: taking the last-updated scrublet output. (after all other QC processes).
ROW_SAMPLES_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\ROW_DATA'
SAMPLES_INFORMATION_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\scrublet\6.10.21'

INFERCNV_SAMPLES_PATH = r'D:\inferCNV\executions\6.10.21'
# path of folder where all samples having cell needed be removed have PKL file containing all barcodes of the cell needed be removed.
IMMUNE_CELLS_REMOVAL_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\inferCNV\analysis_conclusions\immune_clustering'

# Tumor table contains row for each sample splioting the tumor cells (Not immune cells) into cluster ##sorted by InferCNV output##
TUMOR_TABLE_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\inferCNV\analysis_conclusions\tumor_classifying_clusters.xlsx'
# Immune table contains row for each sample that have processed, if there are cells needed to be removed it's indicated in cluster-type.
IMMUNE_TABLE_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\inferCNV\analysis_conclusions\immune_classifying_clusters.xlsx'
# CellBender csv, barcodes of empty cells. should be marked as cellbender empty
EMPTY_BARCODES_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren Laboratory\python_playground\outputs\CellBender\empty_droplets_barcodes_v2.csv'

# Potential contaminated cells, we need to remove these cells from the immune compartment and move them to the stroma
# one, as they could definitely be related to this phenomenon where fibroblast “ate” the neutrophils as this is a common
# feature of removing dying ce  lls.
CONTAMINATED_FIBROBLAST_CELLS_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren Laboratory\Data\tables\kmeans_conclusion_tables\stroma_contaminated_cells_kmeans_k11_cluster7_expressing_neut.csv'

# Potential contaminated cells, we need to remove these cells from the immune compartment and move them to
# the epithelial one, as there were epithelial markers in cluster #10 in k=10. (kmeans 10.5.21)
CONTAMINATED_EPITHELIAL_CELLS_KMEANS10_5_21_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren Laboratory\Data\tables\kmeans_conclusion_tables\epithelial_contaminated_cells_kmeans10.5.21_k10_cluster10.csv'

# Potential contaminated cells, we need to remove these cells from the immune compartment and move them to
# the epithelial one, as there were epithelial markers in cluster #12 in k=15. (kmeans 24.5.21)
CONTAMINATED_EPITHELIAL_CELLS_KMEANS24_5_21_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren Laboratory\Data\tables\kmeans_conclusion_tables\epithelial_contaminated_cells_kmeans24.5.21_k15_cluster12.csv'

# update 10.6.21 - tumor cells can express some immune markers, we need to look at what markers they are allowed to express
TUMOR_IMMUNE_DOUBLETS_MARKER_POLICY_PATH =  r'C:\Users\KerenYlab\Desktop\Technion studies\Keren Laboratory\Data\tables\tumor_immune_doublets_marker_policy.xlsx'

NON_EMPTY_CELL_MIN_GENES = 500

def extract_sample(sample_id):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    data_path = join(SAMPLES_INFORMATION_PATH, sample_id, f'{sample_id}.pkl')
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


def remove_empty_cells(rna_sample):
    """
    update 4.11.21:
    we usually remove samples with a low number of genes (<500), but only if they do not have neutrophil markers.
    Neutrophil are usually small and if we remove all cells with <500 genes we will lose them.
    :param rna_sample:
    :return:
    """
    # First, find empty cells
    sample_num_exp_genes_in_cell = np.sum(rna_sample.counts > 0, axis=1)
    sample_empty_cells_bool = sample_num_exp_genes_in_cell < NON_EMPTY_CELL_MIN_GENES
    sample_empy_cells = rna_sample[sample_empty_cells_bool]

    # Find cells that are not neutrophils
    not_netr_indices = [not 'Neutrophils' in gg for gg in rna_sample.cells_information.getattr('cell_type_list')]
    sample_not_netr = rna_sample[not_netr_indices]

    # Take barcodes (identifiers) of overlap between  sample_empy_cells and sample_not_netr:
    overlap_barcodes = intersection_of_lists(sample_empy_cells.barcodes, sample_not_netr.barcodes)
    comment = f'Empty cell, less than 500 genes are expressed'
    cells_for_removal = rna_sample.get_subset_by_barcodes(overlap_barcodes)
    cells_for_removal.cells_information.setattr('should_be_removed', None, True)
    cells_for_removal.cells_information.setattr('comment', None, comment)


def rearrange_tumor_clusters(tumor_clusters, num_of_cells):
    rearranged_clusters = [(convert_cluster_str_to_int(tumor_clusters[i*2]), tumor_clusters[i*2+1]) for i in range(int(len(tumor_clusters)/2))
            if tumor_clusters[i*2]!='-']

    # if there is no clusters it means that all the cells are stromal
    # so we'll build one cluster of all indexes and mark it '1' (the code of stromal cells).
    if not len(rearranged_clusters):
        rearranged_clusters.append(((0, num_of_cells), 1))

    return rearranged_clusters


def treat_immune_cells_contaminated_with_fibroblast(rna_sample, sample_id):
    # Potential contaminated cells, we need to remove these cells from the immune compartment and move them to
    # the stroma one, as they could definitely be related to this phenomenon where fibroblast “ate” the neutrophils
    # as this is a common feature of removing dying cells.
    contaminated_fib_cells_df = pd.read_csv(CONTAMINATED_FIBROBLAST_CELLS_PATH)
    contaminated_fib_cells_list = list(
        contaminated_fib_cells_df[contaminated_fib_cells_df['sample'] == sample_id]['barcode'])
    rna_sample.get_subset_by_barcodes(contaminated_fib_cells_list).cells_information.setattr('is_immune',
                                                                                             None, False)
    rna_sample.get_subset_by_barcodes(contaminated_fib_cells_list).cells_information.setattr('is_lymphoid',
                                                                                             None, False)
    rna_sample.get_subset_by_barcodes(contaminated_fib_cells_list).cells_information.setattr('is_myeloid',
                                                                                             None, False)
    rna_sample.get_subset_by_barcodes(contaminated_fib_cells_list).cells_information.setattr('is_stromal',
                                                                                             None, True)
    contaminated_comment = f'Found as a potential contaminated cell in kmeans k=11, cluster 7 expressing neut markers'
    rna_sample.get_subset_by_barcodes(contaminated_fib_cells_list).cells_information.setattr('comment',
                                                                                             None, contaminated_comment)


def treat_immune_cells_contaminated_with_epithelial(rna_sample, sample_id, cells_path):
    # update 25.5.21 - Potential contaminated cells, we need to remove these cells from the immune compartment
    # and move them to the epithelial one, as there were epithelial markers in cluster #10 in k=10. (kmeans 10.5.21)
    contaminated_epith_cells_df = pd.read_csv(cells_path)
    contaminated_epith_cells_list = list(
        contaminated_epith_cells_df[contaminated_epith_cells_df['sample'] == sample_id]['barcode'])
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('is_immune',
                                                                                               None, False)
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('is_lymphoid',
                                                                                               None, False)
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('is_myeloid',
                                                                                               None, False)
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('is_epithelial',
                                                                                               None, True)
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('is_stromal',
                                                                                               None, True)
    contaminated_comment = f'10.5.21kmeans - Found as a potential contaminated cells in kmeans k=10, #cluster 10 epith markers'
    rna_sample.get_subset_by_barcodes(contaminated_epith_cells_list).cells_information.setattr('comment',
                                                                                               None,
                                                                                               contaminated_comment)


def treat_tumor_cells_with_immune_conflict(rna_sample):
    """
    update 10.6.21 - tumor cells can express some immune markers, we need to look at what markers they are allowed to
    express, we have this information in the table Moshe gave us.
    if a cell doesn't express any immune marker which it should not express we will remove the mark of 'conflict' and
    classify the cell back as tumor cell and use it from mow on in the analysis.
    :param rna_sample:
    :return:
    """
    # read table built after Moshe note was given. this is he list of all genes that should not be on tumor cells
    tumor_immune_doublets_marker_policy_df = pd.read_excel(TUMOR_IMMUNE_DOUBLETS_MARKER_POLICY_PATH)
    not_allowed_immune_markers = tumor_immune_doublets_marker_policy_df.iloc[:, 1].dropna().tolist()
    # allowed_immune_markers = tumor_immune_doublets_marker_policy_df.iloc[:, 0].dropna().tolist()

    # take cells marked as immune&tumor conflicts
    conflict_rna_sample = rna_sample.filter_cells_by_property('cancer_immune_conflict', True)
    not_allowed_gene_indices = [conflict_rna_sample.gene_names.index(marker) for marker in not_allowed_immune_markers if
                                marker in conflict_rna_sample.gene_names]

    # take indices of all cells that dont express even one of these markers.
    cells_of_tumor_compartment_indices = np.sum(conflict_rna_sample.counts[:, not_allowed_gene_indices] > 0, axis=1) == 0

    cells_moving = conflict_rna_sample[cells_of_tumor_compartment_indices]
    cells_moving.cells_information.setattr('is_cancer', None, True)
    cells_moving.cells_information.setattr('cancer_immune_conflict', None, False)
    contaminated_comment = f'10.6.21 tumor cells with immune markers that may be on tumor cells'
    cells_moving.cells_information.setattr('comment', None, contaminated_comment)


def convert_cluster_str_to_int(tumor_clusters):
    int_tumor_clusters = tuple([int(jj) for jj in ''.join(ii for ii in tumor_clusters if ii not in ['(', ')', ' ']).split(',')])
    return int_tumor_clusters


def update_tumor_cells(sample , tumor_clusters, tumor_barcodes_order):
    """

    Tumor cluster types:
    0	stromal cells - conflicts will be classified as stromal, cancer cells remain cancer
    1	stromal cells - conflicts will be removed, cancer cells remain cancer
    2	tumor - conflicts will be classified as tumor (originally, now (6.5.21 we always remove cancer&immune conflict cells)
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
        # 1	stromal cells - conflicts will be removed, cancer cells remain cancer
        if cluster_type == 1:
            # stromal cells are those which are not cancer and don't have conflicts.
            bool_stromal_indices = [not bool_is_cancer_indices[ii] and not bool_is_conflict_indices[ii] for ii
                                    in range(len(bool_is_cancer_indices))]
            stromal_cells_RNAseq = cluster[bool_stromal_indices]
            stromal_cells_RNAseq.cells_information.setattr('is_stromal', None, True)
        # 2 and 4	tumor - conflicts will be removed, otherwise become tumor
        elif cluster_type == 4 or cluster_type == 2:
            moving_to_be_cancer_cells_indices = [not bool_is_apoptosis[ii] and not bool_is_conflict_indices[ii] for ii
                                    in range(len(bool_is_apoptosis))]
            moving_to_be_cancer_cells = cluster[moving_to_be_cancer_cells_indices]
            moving_to_be_cancer_cells.cells_information.setattr('is_cancer', None, True)
        # Previous version (start-5.5.21) where we split 2 and 4 to 2 different cases and kept some
        # cancer&immune-conflict cells that carry CNVs patterns.
        # # 2	tumor - conflicts will be classified as tumor
        # elif cluster_type == 2:
        #     bool_not_apoptosis_indices = [not val for val in bool_is_apoptosis]
        #     moving_to_be_cancer_cells = cluster[bool_not_apoptosis_indices]
        #     moving_to_be_cancer_cells.cells_information.setattr('is_cancer', None, True)
        # # 4	tumor - conflicts will be removed, otherwise become tumor
        # elif cluster_type == 4:
        #     moving_to_be_cancer_cells_indices = [not bool_is_apoptosis[ii] and not bool_is_conflict_indices[ii] for ii
        #                             in range(len(bool_is_apoptosis))]
        #     moving_to_be_cancer_cells = cluster[moving_to_be_cancer_cells_indices]
        #     moving_to_be_cancer_cells.cells_information.setattr('is_cancer', None, True)
        elif cluster_type == 5:
            cluster.cells_information.setattr('should_be_removed', None, True)
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
    samples = [ss.replace(".pkl", "") for ss in os.listdir(SAMPLES_INFORMATION_PATH) if (not 'csv' in ss and not 'xlsx' in ss)]
    create_folder(OUTPUT_PATH)
    for iter_idx, sample_id in enumerate(samples):
        print(f'{sample_id};\t{iter_idx+1}/{len(samples)}')
        # Extracts a single sample from PC
        rna_sample = loading_sample(row_data_path=join(ROW_SAMPLES_PATH, f'{sample_id}.pkl'),
                                    cells_information_path=join(SAMPLES_INFORMATION_PATH, f'{sample_id}.pkl'))
        rna_sample.cells_information.setattr('should_be_removed', None, False)
        rna_sample.cells_information.setattr('is_epithelial', None, False)
        rna_sample.cells_information.setattr('comment', None, None)


        # update 10.6.21 - tumor cells can express some immune markers, we need to look at what markers they are
        # allowed to express
        treat_tumor_cells_with_immune_conflict(rna_sample)

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

        # (6.5.21) removes also all cells with lymphoid # myeloid markers conflicts.
        mye = np.array(rna_sample.cells_information.getattr('is_myeloid'))
        lym = np.array(rna_sample.cells_information.getattr('is_lymphoid'))
        cells_for_removal = rna_sample[lym & mye]
        cells_for_removal.cells_information.setattr('should_be_removed', None, True)
        comment = f'myeloid & lymphoid conflict'
        cells_for_removal.cells_information.setattr('comment', None, comment)

        # mark empty cells (by CellBender definition)
        empty_barcodes_df = pd.read_csv(EMPTY_BARCODES_PATH)
        empty_cells_barcodes_list = list(empty_barcodes_df[empty_barcodes_df['sample'] == sample_id]['barcode'])
        rna_sample.get_subset_by_barcodes(empty_cells_barcodes_list).cells_information.setattr('is_CelBender_empty', None, True)

        # Potential contaminated cells, we need to remove these cells from the immune compartment and move them to
        # the stroma one, as they could definitely be related to this phenomenon where fibroblast “ate” the neutrophils
        # as this is a common feature of removing dying cells.
        treat_immune_cells_contaminated_with_fibroblast(rna_sample, sample_id)

        # update 25.5.21 - Potential contaminated cells, we need to remove these cells from the immune compartment
        # and move them to the epithelial one, as there were epithelial markers in cluster #10 in k=10. (kmeans 10.5.21)
        treat_immune_cells_contaminated_with_epithelial(rna_sample, sample_id,
                                                        CONTAMINATED_EPITHELIAL_CELLS_KMEANS10_5_21_PATH)

        # update 26.5.21 - Potential contaminated cells, we need to remove these cells from the immune compartment
        # and move them to the epithelial one, as there were epithelial markers in cluster #12 in k=15. (kmeans 24.5.21)
        treat_immune_cells_contaminated_with_epithelial(rna_sample, sample_id,
                                                        CONTAMINATED_EPITHELIAL_CELLS_KMEANS24_5_21_PATH)

        # update 4.11.21:
        # we usually remove samples with a low number of genes (<500), but only if they do not have neutrophil markers.
        remove_empty_cells(rna_sample)

        # Save an updated version of current sample_id with inferCNV changes.
        # pickle.dump((rna_sample), open(join(OUTPUT_PATH, f'{sample_id}.pkl'), 'wb'))
        rna_sample.save_cells_information(join(OUTPUT_PATH, f'{sample_id}.pkl'))


if __name__ == '__main__':
    # if not os.path.isdir(OUTPUT_PATH):
    #     os.mkdir(OUTPUT_PATH)
    tumor_df = extract_tumor_table()
    immune_df = extract_immune_table()

    go_over_all_samples(tumor_df, immune_df)
    _breakpoint = 0