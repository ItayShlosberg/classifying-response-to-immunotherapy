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




OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/heatmap/19.3.21_colorbar_1'
FILTERED_CELLS_PATH = fr'/storage/md_keren/shitay/outputs/variance_filtered/immune_cells_var0.315.pkl'
KMEANS_DIR_PATH = r'/storage/md_keren/shitay/outputs/kmeans/row_kmeans'
CLUSTERING_ANALYSIS_PATH = fr'/storage/md_keren/shitay/outputs/cluster_analysis_17.3.21'

if __name__ == '__main__':

    plt.rcParams['figure.dpi'] = 400
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
    NUM_OF_MARKER_GENES = 30
    print_color_heatmap = True
    print_binary_heatmap = False

    for K in range(2, 16):
        print(f"Current K = {K}")
        kmeans_clusters_path = join(KMEANS_DIR_PATH, fr'kmeans_immune_cells_var0.315_k_{K}.pkl')
        clusters_indices = pickle.load(open(kmeans_clusters_path, 'rb'))['clusters']
        # clusters_indices
        clustering_analysis = pickle.load(open(join(CLUSTERING_ANALYSIS_PATH, f'cluster_analysis_k_{K}.pkl', 'rb')))


        features = filtered_cells.features

        print(clustering_analysis[0].keys())
        gene_indices = []
        gene_vertival_lines = [0]  # for heatmap
        cell_horz_lines = [0]  # for heatmap

        print('num markers for each cluster: ', NUM_OF_MARKER_GENES)
        print('num markers for all clusters: ', NUM_OF_MARKER_GENES * len(clustering_analysis))
        print('num clusters: ', len(clustering_analysis))
        for cluster in clustering_analysis:
            gene_indices += [features.index(ii) for ii in cluster['gene ids'][:NUM_OF_MARKER_GENES]]
            gene_vertival_lines += [gene_vertival_lines[-1] + NUM_OF_MARKER_GENES]
        print('Num of repetition in markers:', NUM_OF_MARKER_GENES * K - len(set(gene_indices)), end="\n\n")

        cells_indices = []
        for cls_idx, cluster_indices in enumerate(clusters_indices):
            cells_indices += cluster_indices
            print(f'{cls_idx}. num cells: {len(cluster_indices)}')
            cell_horz_lines += [cell_horz_lines[-1] + len(cluster_indices)]

        arr_heatmap = filtered_cells.counts[cells_indices][:, gene_indices]
        heatmap = np.zeros_like(arr_heatmap)
        heatmap[arr_heatmap > 1] = 1

        ## 19.3.21 using zscore
        arr_heatmap = scipy.stats.zscore(arr_heatmap, axis=0, ddof=1)

        if print_binary_heatmap:

            fig, ax = plt.pyplot.subplots(1)
            fig.set_size_inches(10, 5)

            sb.heatmap(heatmap)

            for ver_line in gene_vertival_lines[1:-1]:
                ax.axvline(ver_line, color='yellow')

            for hor_line in cell_horz_lines[1:-1]:
                ax.axhline(hor_line, color='red')

            ax.set_title(f"Kmeans k={K}")
            ax.set_ylabel('Clusters')
            ax.set_xlabel('Markers')

            ##### save first heatmap
            fig.savefig(join(OUTPUT_PATH, f'binary_heatmap_{K}.png'))

        if print_color_heatmap:
            fig, ax = plt.pyplot.subplots(1)
            fig.set_size_inches(10, 5)

            # sb_out = sb.heatmap(arr_heatmap)
            sb_out = sb.heatmap(arr_heatmap, vmin=-1, vmax=1)

            for ver_line in gene_vertival_lines[1:-1]:
                ax.axvline(ver_line, color='yellow')

            for hor_line in cell_horz_lines[1:-1]:
                ax.axhline(hor_line, color='yellow')

            ax.set_title(f"Kmeans k={K}")
            ax.set_ylabel('Clusters')
            ax.set_xlabel('Markers')

            ##### save second heatmap
            fig.savefig(join(OUTPUT_PATH, f'density_heatmap_{K}.png'))