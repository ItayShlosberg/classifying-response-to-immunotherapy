"""
Running TSNE with different variance.
Conclusion after all attemp was to use variance 0.315 (~4K genes)
"""

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.manifold import TSNE
import pickle
from Bio.Cluster import kcluster
import os
import numpy as np
import yaml
import os
import pandas
from collections import Counter
import matplotlib.pyplot as plt
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import sys
from shutil import copyfile
import matplotlib as plt
import sys
from os.path import join
# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts'
import sys
import umap
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
# ------- SERVER EXTENSIONS ---------



OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/UMAP/UMAP_cohort_immune_o.315var.pkl'
COHORT_PATH = r'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'



if __name__ == '__main__':

    print("Running UMAP")

    print(f"Input path {COHORT_PATH}")
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    print(f"Input shape {cohort.counts.shape}")

    cells_embedded = umap.UMAP().fit_transform(cohort.counts)
    print(f"UMAP output size {cells_embedded.shape}")

    pickle.dump((cells_embedded), open(OUTPUT_PATH, 'wb'))





