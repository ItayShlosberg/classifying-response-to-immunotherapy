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
# ------- SERVER EXTENSIONS ---------

import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
from scripts.clustering.Clustering_script import user_function


# SAMPLES = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'
#
# OUTPUT = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\Cohort\cohort_all_samples_3.2.21.pkl'
# # ss = r'C:\Users\itay\Desktop\New folder'
#
# samples = [subfolder for subfolder in os.listdir(SAMPLES) if not 'csv' in subfolder]
# rna_sample = extract_droplet_data_from_pickle(join(SAMPLES, samples[0]))
#
# _breakpoint = 0


import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.manifold import TSNE
import pickle
# from Bio.Cluster import kcluster
import os
import numpy as np
import yaml
import os
import pandas
from collections import Counter
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import sys
import pyclustering
from shutil import copyfile
import matplotlib as plt



for k in range(2, 16):
    path = r'/storage/md_keren/shitay/outputs/kmeans/kmeans_immune_cells_var0.315_k_2.pkl'
    kmeans_instance = pickle.load(open(path, 'rb'))
print(kmeans_instance.get_clusters())
