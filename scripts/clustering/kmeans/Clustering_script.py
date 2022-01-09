"""
Here we run kmeans a few times (max - 10 times) for requested K and take the result with the biggest number
of non_empty clusters.
"""

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
from utilities.droplet_dataset import build_cohort, get_requested_subset
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
from scipy.stats import pearsonr


# INPUT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/immune_cells_26.6.21_4k_genes.pkl'
# OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/cohort_26.6.21_run_11.8.21'
INPUT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/subcohort/normalized/1.1.22/subcohort_immune_cells_normalized_1.1.22_4k_genes.pkl'
OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/subcohort_1.1.22_run_1.1.22'
FILE_NAME = r'kmeans_immune_cells_4k_genes'
SUBSET = 'MYELOIDS'    # None - all cells, MYELOIDS/CYTOTOXIC_T_CELLS
MAX_ITERATION = 10  # Number of iteration to get non_empty clusters.
NON_EMPTY_CLUSTER_DEF = 20  # defines how many cells is a non_empty cluster


def pearson_distance_metric(point1, point2):
    return 1 - pearsonr(point1, point2)[0]


def save_result(kmeans_instance):
    create_folder(OUTPUT_PATH)
    INSTANCES_DIRECTORY = r'kmeans_instances'
    ROW_CLUSTERS_DIRECTORY = r'row_kmeans'
    BARCODE_MAPPING_DIRECTORY = r'barcode_mapping'
    create_folder(join(OUTPUT_PATH, INSTANCES_DIRECTORY))
    create_folder(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY))
    create_folder(join(OUTPUT_PATH, BARCODE_MAPPING_DIRECTORY))

    pickle.dump((kmeans_instance), open(join(OUTPUT_PATH, INSTANCES_DIRECTORY, FILE_NAME+f'_k_{K}.pkl'), 'wb'))
    pickle.dump({'clusters': kmeans_instance.get_clusters(), 'centers': kmeans_instance.get_centers()}, open(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY, FILE_NAME+f'_k_{K}.pkl'), 'wb'))
    print(f'OUTPUT was saved in {OUTPUT_PATH}')


def run_kmeans(K, gene_expression_matrix):
    # create K-Means algorithm with specific distance metric
    metric = distance_metric(type_metric.USER_DEFINED, func=pearson_distance_metric)
    initial_centers = kmeans_plusplus_initializer(gene_expression_matrix, K).initialize()
    kmeans_instance = kmeans(gene_expression_matrix, initial_centers, metric=metric)
    # run cluster analysis and obtain results
    kmeans_instance.process()
    return kmeans_instance


def save_mapping_df(kmeans_instance):
    BARCODE_MAPPING_DIRECTORY = r'barcode_mapping'
    clusters_cells = [cohort[cluster_indices] for cluster_indices in kmeans_instance.get_clusters()]
    df = pd.DataFrame(columns=['Sample', 'Barcode', 'Cluster'])
    cluster_idx = 1
    for idx, cls in enumerate(clusters_cells):
        df = df.append(
            pd.DataFrame(transpose_list([cls.samples, cls.barcodes, [cluster_idx] * cls.number_of_cells]),
                         columns=df.columns))
        cluster_idx += 1
    df = df.reset_index()
    df.to_csv(join(OUTPUT_PATH, BARCODE_MAPPING_DIRECTORY, FILE_NAME + f'_k_{K}_clusters_mapping.csv'), index=False)
    print()


if __name__ == '__main__':

    K = int(sys.argv[1])

    # Loads cohort
    print(f'Loading cohort from:\n{INPUT_PATH}')
    cohort = pickle.load(open(INPUT_PATH, 'rb'))
    if SUBSET:
        cohort = get_requested_subset(cohort, SUBSET)

    c_iteration = 1
    print(f'Running kmeans on {INPUT_PATH}, with K = {K}')
    kmeans_instance = run_kmeans(K, cohort.counts)

    # Supervise that there is non_empty clusters, runs until you get just non_empty clusters
    curr_most_clusters = sum([NON_EMPTY_CLUSTER_DEF <= len(cluster_indices) for cluster_indices in kmeans_instance.get_clusters()])
    instance_of_most_clusters = kmeans_instance
    while curr_most_clusters < K and c_iteration < 10:
        c_iteration += 1
        print(f'Cluster sizes:')
        print([len(cluster_indices) for cluster_indices in kmeans_instance.get_clusters()], end='\n\n')
        print(f'The number of non_empty clusters is {curr_most_clusters}. Running again kmeans. Attempt num: {c_iteration}')

        kmeans_instance = run_kmeans(K, cohort.counts)
        n_non_empty_clusters = sum([NON_EMPTY_CLUSTER_DEF <= len(cluster_indices) for cluster_indices in kmeans_instance.get_clusters()])

        if n_non_empty_clusters > curr_most_clusters:
            curr_most_clusters = n_non_empty_clusters
            instance_of_most_clusters = kmeans_instance

    # Saves results
    save_result(instance_of_most_clusters)

    # get DF of barcode division of clusters
    save_mapping_df(instance_of_most_clusters)

