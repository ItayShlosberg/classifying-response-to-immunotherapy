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




# VARIANCE = [6, 4, 1, 0.75, 0.5, 0.325, 0.315, 0.25]
# perplexity
PERPLEXITY = [10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0]


OUTPUT_DIR = r'/storage/md_keren/shitay/outputs/TSNE/perplexity_euc/'
# COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/cohort.pkl'
COHORT_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'

if __name__ == '__main__':

    print("Running TSNE")
    print(f'ARG {sys.argv[1]}')
    # # print(type(int(sys.argv[1])))
    # variance = VARIANCE[int(sys.argv[1])]
    perplexity = PERPLEXITY[int(sys.argv[1])]

    # variance = 0.315
    # for variance in VARIANCE:
    print(f"Using perplexity {perplexity}")
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    # cohort = cohort.filter_cells_by_property('is_immune', True)
    # cohort.filter_genes_by_variance(variance, in_place=True)
    print(f"Counts shape {cohort.counts.shape}")
    # pickle.dump((cohort), open(join(OUTPUT_DIR, f'immune_cells_var{variance}.pkl'), 'wb'))

    # print(f'pickle has been saved in {join(OUTPUT_DIR, f"immune_cells_var{variance}.pkl")}')
    cells_embedded = TSNE(n_components=2, perplexity=perplexity, random_state=21).fit_transform(cohort.counts)
    print(f"TSNE output size {cells_embedded.shape}")
    pickle.dump((cells_embedded), open(join(OUTPUT_DIR, f'immune_TSNE_embedded_var0.315_perplexity{perplexity}.pkl'), 'wb'))


import pyclustering.utils.metric