# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts'
lib6 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts/preprocess_data'
import sys
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
sys.path.append(lib6)
from classifying_cell_types import *
from utilities.general_helpers import *
# ------- SERVER EXTENSIONS ---------
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.manifold import TSNE
from utilities.ML_environment import find_marker_genes_in_cluster
import pickle
# from Bio.Cluster import kcluster
# from Bio.Cluster import kcluster
# import pyclustering
from utilities.ML_environment import find_marker_genes_in_cluster
import os
import numpy as np
import yaml
from os.path import join
import os
import pandas
from collections import Counter
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import sys
import seaborn as sns
# import statsmodels as sm
import scipy.stats as stats
from scipy.stats import rankdata
from sklearn.manifold import TSNE
import pickle
import numpy as np
from utilities.general_helpers import flatten_list
# from utilities.ML_environment import find_marker_genes_in_cluster
from shutil import copyfile
import matplotlib.pyplot as plt
from utilities.clustering_tools import find_marker_genes_in_cluster, find_markers_in_clusters
from utilities.general_helpers import are_the_lists_identical
from utilities.clustering_tools import get_clusters_indices, find_satisfying_list_of_markers_in_clusters



print(f'8.8.21')
COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/immune_cells_26.6.21_4k_genes.pkl'
OUTPUT = r'/storage/md_keren/shitay/garbage/D.pkl'
immune_cells = pickle.load(open(COHORT_PATH, 'rb'))


from scipy.spatial.distance import cdist
D = cdist(immune_cells.counts, immune_cells.counts, 'correlation')

pickle.dump((D), open(OUTPUT, 'wb'), protocol=4)
