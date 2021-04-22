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


ROW_DATA_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/ROW_DATA/'
CELLS_INFORMATION = r'/storage/md_keren/shitay/Data/inferCNV_data/update_runs/4.3.21/'
OUTPUT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/5.4.21.pkl'


if __name__ == '__main__':
    gene_list = build_cohort_gene_list(samples_information_path = CELLS_INFORMATION)
    cohort = build_cohort(samples_row_data_path = ROW_DATA_PATH,
                          samples_cells_information_path = CELLS_INFORMATION,
                          gene_id_list = gene_list)
    pickle.dump((cohort), open(OUTPUT_PATH, 'wb'))