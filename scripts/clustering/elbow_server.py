"""
Elbow - estimating the variance explained by a given solution
( Cell, 2018, Figure S8B)

computing the sum of pair-wise distances between all cells in different clusters Dis_b=Σl=1-k(Σi∈cl,j∉cl D(i,j))
and the total distance Dis_t= Σi,j D(i,j). The ratio between these two measures V = Dis_b/Dis_t was used to estimate the
variance explained by a given solution (Figure S8B),
such that in the extreme case where all cells are clustered together or the case where each cell is a single cluster,
this ratio would be 0 and 1, respectively. Exploring this ratio,
we then select the solutions that are near plateau


In order to get elbow plot you can use the output of this script as follows:
plt.pyplot.plot(list(All_Dist_b.keys()), [(ii/2)/Dis_t for ii in All_Dist_b.values()], '-or');
plt.pyplot.ylim((0,1))

Note: Dist_b Needed to be divided by 2.
"""




from os.path import join
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
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
from utilities.general_helpers import create_folder
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
from scipy.stats import pearsonr
import random
from scipy.spatial.distance import cdist






FILTERED_CELLS_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
save_path = r'/storage/md_keren/shitay/outputs/clustering/elbow/elbow.pkl'


D = cdist(filtered_cells.counts, filtered_cells.counts, 'correlation')


Dis_t = 0
for i in range(D.shape[0]):
    for j in range(i, D.shape[0]):
        Dis_t += D[i,j]

All_Dist_b = {}
for K in range(2, 16):
    print(f'K = {K}, ', end='')
    kmeans_clusters_path = fr'/storage/md_keren/shitay/outputs/clustering/kmeans/row_kmeans/kmeans_immune_cells_var0.315_k_{K}.pkl'
    clustering_out = pickle.load(open(kmeans_clusters_path, 'rb'))
    #     clustering_analysis = pickle.load(open(fr'/storage/md_keren/shitay/outputs/clustering/cluster_analysis/cluster_analysis_21.3.21/cluster_analysis_k_{K}.pkl', 'rb'))

    Dist_b = 0
    for idx in range(len(clustering_out['clusters'])):
        cluster_indices = clustering_out['clusters'][idx]
        other_clusters_indices = [ii for ii in range(filtered_cells.number_of_cells) if not ii in cluster_indices]

        for ii in cluster_indices:
            for jj in other_clusters_indices:
                Dist_b += D[ii, jj]

    All_Dist_b[K] = Dist_b


pickle.dump(({'Ks_Dist_b': All_Dist_b, 'Dis_t': Dis_t}), open(save_path, 'wb'))