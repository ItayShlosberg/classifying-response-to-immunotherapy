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
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
# ------- SERVER EXTENSIONS ---------




OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/TSNE/all_cells_variance_0.315.pkl'
COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/cohort.pkl'



if __name__ == '__main__':

    print("Running TSNE")
    variance = 0.315
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    cohort.filter_genes_by_variance(variance, in_place=True)
    print(f"Counts shape {cohort.counts.shape}")


    cells_embedded = TSNE(n_components=2, random_state=21).fit_transform(cohort.counts)
    print(f"TSNE output size {cells_embedded.shape}")
    pickle.dump((cells_embedded), open(OUTPUT_PATH, 'wb'))


import pyclustering.utils.metric