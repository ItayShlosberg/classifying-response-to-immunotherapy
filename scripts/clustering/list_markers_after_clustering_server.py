from os.path import join
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
import seaborn as sns
# import pyclustering
import matplotlib.pylab as plt
import seaborn as sb
from shutil import copyfile
import matplotlib as plt
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
from utilities.general_helpers import create_folder






OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/markers/21.3.21'
FILTERED_CELLS_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
KMEANS_DIR_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/row_kmeans'
CLUSTERING_ANALYSIS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/cluster_analysis/cluster_analysis_21.3.21'

if __name__ == '__main__':


    NUM_OF_MARKER_GENES = 30

    create_folder(OUTPUT_PATH)
    for K in range(2, 16):
        print(f"working on {K}")
        cluster_dir = join(OUTPUT_PATH, f'markers_cluster_{K}')
        create_folder(cluster_dir)
        print(f"Current K = {K}")
        kmeans_clusters_path = join(KMEANS_DIR_PATH, fr'kmeans_immune_cells_var0.315_k_{K}.pkl')
        clusters_indices = pickle.load(open(kmeans_clusters_path, 'rb'))['clusters']
        # clusters_indices
        clustering_analysis = pickle.load(open(join(CLUSTERING_ANALYSIS_PATH, f'cluster_analysis_k_{K}.pkl'), 'rb'))


        print('num clusters: ', len(clustering_analysis))
        for cls_idx, cluster in enumerate(clustering_analysis):
            pval = cluster['pval']
            logratio = cluster['log ratios']
            gene_ids = np.array(cluster['gene ids'])
            gene_names = np.array(cluster['gene names'])

            logratio_indices = np.flip(np.argsort(logratio))

            gene_names = gene_names[logratio_indices].tolist()
            gene_ids = gene_ids[logratio_indices].tolist()
            logratio = logratio[logratio_indices].tolist()

            aa = np.array([gene_names, pval, logratio]).T
            df = pd.DataFrame(aa, columns=['gene', 'pval adj', 'logFC'])
            df.to_csv(join(cluster_dir, f'markers_cluster_{cls_idx+1}.csv'), index=False)

