"""
robustness analysis clustering

To determine the robustness of each clustering solution, we need to perform 100 iterations in which we randomly remove 10% of the cells,
and re-run the k-means algorithm and check the stability of the clustering solution. We then need to quantify the agreement of
a given solution with the original one as the number of pairs of cells that were either clustered together, or not clustered together,
in both solutions, divided by the total number pairs shared between the runs.

Here we have output of 100 runs of kmeans for each K.
and now we want to quantify the agreement of a given solution with the original one as the number of pairs of cells
that were either clustered together, or not clustered together, in both solutions,
divided by the total number pairs shared between the runs.

so we only calculate for each of the runs the score
"""

from os.path import join
# ------- SERVER EXTENSIONS ---------
lib = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
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
from utilities.general_helpers import create_folder
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
from scipy.stats import pearsonr
import random

# COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/cohort.pkl'

BASE_KMEANS_PATH = r''
ROBUSTNESS_RUNS_PATH = r'/storage/md_keren/shitay/outputs/clustering/robustness_analysis_clustering/runs'
# FILTERED_CELLS_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
NUMBER_OF_IMMUNE_CELLS = r'/storage/md_keren/shitay/outputs/clustering/number_of_immune_cells.txt'
ROW_KMEANS_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/row_kmeans'
RUNS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/robustness_analysis_clustering/runs'
OUTPUT = r'/storage/md_keren/shitay/outputs/clustering/robustness_analysis_clustering/run_scores'

Ks = list(range(2, 16))
runs_values = list(range(100))


if __name__ == '__main__':

    # job id, between 0-1499
    job_process = int(sys.argv[1]) + 1000

    k = Ks[int(job_process / 100)]
    run = runs_values[job_process % 100]

    print(f'Calculating score of k={k}, run={run}')

    original_clusters = pickle.load(open(join(ROW_KMEANS_PATH, fr'kmeans_immune_cells_var0.315_k_{k}.pkl'), 'rb'))
    current_it_clusters_dic = pickle.load(open(join(RUNS_PATH, f'{k}/{run}.pkl'), 'rb'))


    # transfer indices to the respective ones (because there was a filtering of 10% of the cells)
    # filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
    # n_cells = filtered_cells.number_of_cells
    # del filtered_cells
    with open(NUMBER_OF_IMMUNE_CELLS, 'r') as f:
        n_cells = int(f.read())


    mapping_orig_to_filtering = np.arange(n_cells)[current_it_clusters_dic['sample_indices']]
    mapping_orig_to_filtering_dic = {}
    for i in range(mapping_orig_to_filtering.shape[0]):
        mapping_orig_to_filtering_dic[mapping_orig_to_filtering[i]] = i

    # mapping orig indices to respective indices values after filtering
    base_clusters = [[mapping_orig_to_filtering_dic[jj] for jj in ii if jj in mapping_orig_to_filtering_dic.keys()] for
                     ii in original_clusters['clusters']]

    # retrieve iteration clusters.
    iteration_clusters = current_it_clusters_dic['clusters']
    filterd_indices = flatten_list(iteration_clusters)
    n_cells_filtered = len(filterd_indices)

    # build a map of pairs in base clustering
    pair_map_in_base_clusters = np.zeros((n_cells_filtered, n_cells_filtered), dtype=np.uint8)
    for cluster in base_clusters:
        for ii in cluster:
            for jj in cluster:
                if ii != jj:
                    pair_map_in_base_clusters[ii, jj] = 1

    # build a map of pairs in current iteration clustering
    pair_in_curr_iteration_clusters = np.zeros((n_cells_filtered, n_cells_filtered), dtype=np.uint8)
    for cluster in iteration_clusters:
        for ii in cluster:
            for jj in cluster:
                if ii != jj:
                    pair_in_curr_iteration_clusters[ii, jj] = 1

    # using the maps to find agreement of pairs. divided by 2 because each pair is counted twice.
    number_of_pairs_agreement = np.sum((pair_map_in_base_clusters + pair_in_curr_iteration_clusters) == 2) / 2

    # calculate the number of pairs in the original clusters:
    number_of_pairs_in_origin_clusters = np.sum(pair_map_in_base_clusters == 1)/2
    number_of_pairs_in_iter_clusters = np.sum(pair_in_curr_iteration_clusters == 1) / 2

    # the final score
    # result = number_of_pairs_agreement / number_of_pairs_in_origin_clusters

    # saving it
    with open(join(OUTPUT, f'K_{k}_run_{run}.txt'), "w") as text_file:
        text_file.write('number_of_pairs_agreement ' + str(number_of_pairs_agreement)+'\n')
        text_file.write('number_of_pairs_in_origin_clusters '+str(number_of_pairs_in_origin_clusters)+'\n')
        text_file.write('number_of_pairs_in_iter_clusters ' + str(number_of_pairs_in_iter_clusters) + '\n')