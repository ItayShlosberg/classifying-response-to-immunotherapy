from os.path import join
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
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
from scipy.stats import pearsonr


# COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/cohort.pkl'
# INPUT_PATH = r'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
INPUT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/immune_cells_26.6.21_4k_genes.pkl'
OUTPUT_PATH = r'/storage/md_keren/shitay/garbage/Debug/clustering'
FILE_NAME = r'kmeans_immune_cells_4k_genes'
SUBSET = 'CYTOTOXIC_T_CELLS' # None - all cells, MYELOIDS/CYTOTOXIC_T_CELLS


def pearson_distance_metric(point1, point2):
    return 1 - pearsonr(point1, point2)[0]


def get_requested_subset(cohort):
    if SUBSET is None:
        print(f'No subset request was given')
        return cohort
    elif SUBSET == 'MYELOIDS':
        print(f'Will be performed on myeloids')
        myeloid_indices = cohort.cells_information.getattr('is_myeloid')
        return cohort[myeloid_indices]
    elif SUBSET == 'CYTOTOXIC_T_CELLS':
        print(f'Will be performed on cytotoxic T cells')
        cytotoxic_T_cells_indices = ['CD8 Cytotoxic T cells' in types for types in
                                     cohort.cells_information.getattr('cell_type_list')]
        return cohort[cytotoxic_T_cells_indices]
    print(f'ERROR! No valid SUBSET value was passed!')
    raise Exception(f'ERROR! the SUBSET value is not valid!')


if __name__ == '__main__':

    file_num = int(sys.argv[1])
    K = 9
    create_folder(OUTPUT_PATH)

    # Loads cohort
    print(f'Loading cohort from:\n{INPUT_PATH}')
    cohort = pickle.load(open(INPUT_PATH, 'rb'))
    if SUBSET:
        cohort = get_requested_subset(cohort)

    sample = cohort.counts

    print(f'Running kmeans on {INPUT_PATH}, with K = {K}')
    # create K-Means algorithm with specific distance metric
    metric = distance_metric(type_metric.USER_DEFINED, func=pearson_distance_metric)
    initial_centers = kmeans_plusplus_initializer(sample, K).initialize()
    kmeans_instance = kmeans(sample, initial_centers, metric=metric)
    # run cluster analysis and obtain results
    kmeans_instance.process()

    # Saves results
    INSTANCES_DIRECTORY = r'kmeans_instances'
    ROW_CLUSTERS_DIRECTORY = r'row_kmeans'
    BARCODE_MAPPING_DIRECTORY = r'barcode_mapping'
    create_folder(join(OUTPUT_PATH, INSTANCES_DIRECTORY))
    create_folder(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY))
    create_folder(join(OUTPUT_PATH, BARCODE_MAPPING_DIRECTORY))

    pickle.dump((kmeans_instance), open(join(OUTPUT_PATH, INSTANCES_DIRECTORY, FILE_NAME+f'_k_{file_num}.pkl'), 'wb'))
    pickle.dump({'clusters': kmeans_instance.get_clusters(), 'centers': kmeans_instance.get_centers()}, open(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY, FILE_NAME+f'_k_{K}.pkl'), 'wb'))
    print(f'OUTPUT was saved in {OUTPUT_PATH}')

    # get DF of barcode division of clusters
    clusters_cells = [cohort[cluster_indices] for cluster_indices in kmeans_instance.get_clusters()]
    df = pd.DataFrame(columns=['Sample', 'Barcode', 'Cluster'])
    cluster_idx = 1
    for idx, cls in enumerate(clusters_cells):
        df = df.append(
            pd.DataFrame(transpose_list([cls.samples, cls.barcodes, [cluster_idx] * cls.number_of_cells]),
                         columns=df.columns))
        cluster_idx += 1
    df = df.reset_index()
    df.to_csv(join(OUTPUT_PATH, BARCODE_MAPPING_DIRECTORY, FILE_NAME+f'_k_{file_num}_clusters_mapping.csv'), index=False)
    print()