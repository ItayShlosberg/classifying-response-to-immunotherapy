"""
After kmeans was done, you can run this script to find markers.
after finding markers run list_markers_after_clustering_server.py on the output of this script to get more convenient
format for the output.

the script get as an input the index of k (the kmeans k) and save the markers gene names the pvals and logFC values.
"""

from os.path import join
import sys
# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts'
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
from utilities.droplet_dataset import *
from utilities.clustering_tools import find_satisfying_list_of_markers_in_clusters
from utilities.general_helpers import create_folder
import os

OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/cluster_analysis/subcohort_1.1.22_run_1.1.22' # cytotoxic_t_cells/myeloid

KMEANS_ROW_CLUSTERS_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/subcohort_1.1.22_run_1.1.22/row_kmeans'
KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'

# FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/cohort_normalized_26.6.21.pkl'
FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/subcohort/normalized/1.1.22/subcohort_immune_cells_normalized_1.1.22_protein_coding_genes.pkl'

SUBSET = 'MYELOIDS'    # None - all cells, MYELOIDS/CYTOTOXIC_T_CELLS

# use this path to temporarily store the subset of cohort if you can't run this script as a job on server.
# TEMPORARY_PATH = r'/storage/md_keren/shitay/garbage/myel_cohort.pkl' # cytox_cohort/myel_cohort


# def tempor_wrapper_loader():
#     if TEMPORARY_PATH and os.path.isfile(TEMPORARY_PATH):
#         filtered_cells = pickle.load(open(TEMPORARY_PATH, 'rb'))
#     else:
#         filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
#         filtered_cells = filtered_cells.filter_cells_by_property('is_immune', True)
#         if SUBSET:
#             filtered_cells = get_requested_subset(filtered_cells, SUBSET)
#         if TEMPORARY_PATH:
#             pickle.dump((filtered_cells), open(TEMPORARY_PATH, 'wb'), protocol=4)
#     return filtered_cells


def data_loader():
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
    filtered_cells = filtered_cells.filter_cells_by_property('is_immune', True)
    if SUBSET:
        filtered_cells = get_requested_subset(filtered_cells, SUBSET)
    return filtered_cells


if __name__ == '__main__':

    K = int(sys.argv[1])
    print(f'Running Analyze_kmeans_clusters with K={K}')
    print(f'Loading cohort file:\n{FILTERED_CELLS_PATH}')

    filtered_cells = data_loader()
    # filtered_cells = tempor_wrapper_loader()

    print(f'Number of cells: {filtered_cells.number_of_cells}')

    kmeans_file_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME + f'_k_{K}.pkl')
    print(f'Loading kmeans file:\n{kmeans_file_path}')
    clusters = pickle.load(open(kmeans_file_path, 'rb'))['clusters']

    markers = find_satisfying_list_of_markers_in_clusters(filtered_cells, clusters)
    create_folder(OUTPUT_PATH)
    pickle.dump((markers), open(os.path.join(OUTPUT_PATH, f'cluster_analysis_k_{K}.pkl'), 'wb'))