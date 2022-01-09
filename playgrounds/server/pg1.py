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
from scripts.clustering.kmeans.Clustering_script import *

import pickle
import pyclustering



for k in range(2, 16):
    path = fr'/storage/md_keren/shitay/outputs/kmeans/kmeans_immune_cells_var0.315_k_{k}.pkl'
    out_path = fr'/storage/md_keren/shitay/outputs/kmeans/row_kmeans/kmeans_immune_cells_var0.315_k_{k}.pkl'


    kmeans_instance = pickle.load(open(path, 'rb'))
    output = {'clusters': kmeans_instance.get_clusters(), 'centers': kmeans_instance.get_centers()}
    pickle.dump((output), open(out_path, 'wb'))