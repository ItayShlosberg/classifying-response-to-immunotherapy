"""
robustness analysis clustering

To determine the robustness of each clustering solution, we need to perform 100 iterations in which we randomly remove 10% of the cells,
and re-run the k-means algorithm and check the stability of the clustering solution. We then need to quantify the agreement of
a given solution with the original one as the number of pairs of cells that were either clustered together, or not clustered together,
in both solutions, divided by the total number pairs shared between the runs.

Here we just run kmeans 100 time for each cluster. running this script once will run kmeans for k=2-15 with a random sample of 90%
of the immune cells, and will save the indices of the sample and the clusters indices.
Important note: the cluster indices are after the filtering of the sample indices. should be taken into account when continue
to robustness analysis calculation.

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
INPUT_PATH = r'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
OUTPUT_DIR = r'/storage/md_keren/shitay/outputs/clustering/robustness_analysis_clustering/runs'


def pearson_distance(point1, point2):
    """
    Pearson distance.
    :param point1:
    :param point2:
    :return:
    """
    return 1 - pearsonr(point1, point2)[0]


if __name__ == '__main__':

    iteration_id = int(sys.argv[1])
    percent_of_cells = 0.9
    immune_cells_var0315 = pickle.load(open(INPUT_PATH, 'rb'))

    number_of_cells = immune_cells_var0315.number_of_cells
    sample_indices = random.sample(list(range(number_of_cells)), k=round(number_of_cells * percent_of_cells))

    sample = immune_cells_var0315.counts[sample_indices]

    create_folder(OUTPUT_DIR)
    for K in range(2, 16):
        create_folder(join(OUTPUT_DIR, str(K)))

        print(f'Running kmeans on {INPUT_PATH}, with K = {K}')

        metric = distance_metric(type_metric.USER_DEFINED, func=pearson_distance)

        # create K-Means algorithm with specific distance metric
        initial_centers = kmeans_plusplus_initializer(sample, K).initialize()

        kmeans_instance = kmeans(sample, initial_centers, metric=metric)

        # run cluster analysis and obtain results
        kmeans_instance.process()


        output_path = join(OUTPUT_DIR, str(K), f'{iteration_id}.pkl')
        pickle.dump(({'sample_indices': sample_indices, 'clusters': kmeans_instance.get_clusters()}), open(output_path, 'wb'))
        print(f'OUTPUT was saved in {OUTPUT_DIR}')

        del kmeans_instance, initial_centers