"""
After kmeans was done, you can run this script to find markers.
after finding markers run list_markers_after_clustering_server.py on the output of this script to get more convenient
format for the output.

the script get as an input the index of k (the kmeans k) and save the markers gene names the pvals and logFC values.
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
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
from scipy.stats import pearsonr
from utilities.ML_environment import find_markers_in_clusters
from utilities.general_helpers import create_folder


OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/cluster_analysis/cluster_analysis_16.4.21'
FILTERED_CELLS_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
KMEANS_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/row_kmeans/'

if __name__ == '__main__':

    K = int(sys.argv[1])
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))

    kmeans_file_path = os.path.join(KMEANS_PATH, f'kmeans_immune_cells_var0.315_k_{K}.pkl')
    print(f'Loading kmeans file:\n{kmeans_file_path}')
    clusters = pickle.load(open(kmeans_file_path, 'rb'))['clusters']

    markers = find_markers_in_clusters(filtered_cells, clusters)
    create_folder(OUTPUT_PATH)
    pickle.dump((markers), open(os.path.join(OUTPUT_PATH, f'cluster_analysis_k_{K}.pkl'), 'wb'))