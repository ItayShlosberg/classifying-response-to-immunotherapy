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
INPUT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/5.21/immune_cells_4k_genes.pkl'
OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/24.5.21'
FILE_NAME = r'kmeans_immune_cells_4k_genes'


def pearson_distance_metric(point1, point2):
    return 1 - pearsonr(point1, point2)[0]


if __name__ == '__main__':

    K = int(sys.argv[1])
    create_folder(OUTPUT_PATH)

    # Loads cohort
    print(f'Loading cohort from:\n{INPUT_PATH}')
    cohort = pickle.load(open(INPUT_PATH, 'rb'))
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
    create_folder(join(OUTPUT_PATH, INSTANCES_DIRECTORY))
    create_folder(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY))

    pickle.dump((kmeans_instance), open(join(OUTPUT_PATH, INSTANCES_DIRECTORY, FILE_NAME+f'_k_{K}.pkl'), 'wb'))
    pickle.dump({'clusters': kmeans_instance.get_clusters(), 'centers': kmeans_instance.get_centers()}, open(join(OUTPUT_PATH, ROW_CLUSTERS_DIRECTORY, FILE_NAME+f'_k_{K}.pkl'), 'wb'))
    print(f'OUTPUT was saved in {OUTPUT_PATH}')
