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
INPUT_PATH = r'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/kmeans/'


def user_function(point1, point2):
    return 1 - pearsonr(point1, point2)[0]


if __name__ == '__main__':

    K = int(sys.argv[1])
    print(f'Running kmeans on {INPUT_PATH}, with K = {K}')

    metric = distance_metric(type_metric.USER_DEFINED, func=user_function)

    immune_cells_var0315 = pickle.load(open(INPUT_PATH, 'rb'))
    sample = immune_cells_var0315.counts


    # create K-Means algorithm with specific distance metric
    initial_centers = kmeans_plusplus_initializer(sample, K).initialize()

    kmeans_instance = kmeans(sample, initial_centers, metric=metric)

    # run cluster analysis and obtain results
    kmeans_instance.process()

    pickle.dump((kmeans_instance), open(join(OUTPUT_PATH, fr'kmeans_immune_cells_var0.315_k_{K}.pkl'), 'wb'))
    print(f'OUTPUT was saved in {OUTPUT_PATH}')
