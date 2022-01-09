"""

7.1.22 - NOT IN USE!!!! print_heatmaps_server is the one we use.


Print Heatmap of gene markers, after finding markers of clusters.
Done for all Ks.
Takes the first NUM_OF_MARKER_GENES (default is 30) gene markers with the highest log_FC.
if there are >= NUM_OF_MARKER_GENES gene markers with pval=0 it will consider only them.

"""

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
import matplotlib
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



OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/heatmap/26.6.21'
FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/cohort_normalized_26.6.21.pkl'
KMEANS_ROW_CLUSTERS_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/26.6.21/row_kmeans'
KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'
CLUSTERING_ANALYSIS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/cluster_analysis/cluster_analysis_26.6.21'


def create_cmap():
    """
    Get a colorful cmap, colors for heatmap in seaborn.
    :return:
    """
    def mult_scalar(a, s):
        size = len(a)
        return tuple([a[i]*s for i in range(size)])

    def subtract(a, b):
        size = len(a)
        return tuple([a[i]-b[i] for i in range(size-1)]+[1])

    newcolors = matplotlib.get_cmap('viridis',100).colors
    for i in range(50):
        newcolors[i, :] = subtract(plt.colors.to_rgba('darkblue'), mult_scalar(matplotlib.colors.to_rgba('darkblue'), (i/50)))
    for i in range(50, 100):
        newcolors[i, :] = subtract(matplotlib.colors.to_rgba('yellow'), mult_scalar(matplotlib.colors.to_rgba('yellow'), ((50-(i-50))/50)))
    cmap = matplotlib.colors.ListedColormap(newcolors)
    return cmap


if __name__ == '__main__':

    plt.rcParams['figure.dpi'] = 1000
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
    filtered_cells = filtered_cells.filter_cells_by_property('is_immune', True)
    cmap = pickle.load(open(r'/storage/md_keren/shitay/outputs/clustering/heatmap/colorbar.pkl', 'rb'))
    NUM_OF_MARKER_GENES = 30
    print_color_heatmap = True
    print_binary_heatmap = False

    create_folder(OUTPUT_PATH)
    for K in range(2, 16):
        print(f"Current K = {K}")

        kmeans_file_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME + f'_k_{K}.pkl')
        print(f'Loading kmeans file:\n{kmeans_file_path}')
        clusters_indices = pickle.load(open(kmeans_file_path, 'rb'))['clusters']
        # markers
        clustering_analysis = pickle.load(open(join(CLUSTERING_ANALYSIS_PATH, f'cluster_analysis_k_{K}.pkl'), 'rb'))


        features = filtered_cells.features

        print(clustering_analysis[0].keys())
        gene_indices = []
        gene_horiz_lines = [0]  # for heatmap
        cell_vertical_lines = [0]  # for heatmap
        top_gene_markers = {}

        print('num markers for each cluster: ', NUM_OF_MARKER_GENES)
        print('num markers for all clusters: ', NUM_OF_MARKER_GENES * len(clustering_analysis))
        print('num clusters: ', len(clustering_analysis))
        for cluster_dic in clustering_analysis:
            cluster_idx = cluster_dic['cluster id']
            cluster = cluster_dic['markers']
            print(f'cluster idx {cluster_idx}')

            #region
            cluster = cluster.sort_values(by=['log_FC'], ascending=False)
            gene_ids = cluster['features'].tolist()
            top_gene_markers[cluster_idx] = cluster['gene names'][:5]
            gene_indices += [features.index(ii) for ii in gene_ids[:NUM_OF_MARKER_GENES]]
            gene_horiz_lines += [gene_horiz_lines[-1] + len(gene_ids[:NUM_OF_MARKER_GENES])]
            #endregion

            #
            # print(cluster)
            # pval = cluster['pval']
            # logratio = cluster['log ratios']
            # gene_ids = np.array(cluster['gene ids'])
            # gene_names = np.array(cluster['gene names'])
            #
            # if sum(pval == 0) >= NUM_OF_MARKER_GENES:
            #     logratio = logratio[pval == 0]
            #     gene_ids = gene_ids[pval == 0]
            #     gene_names = gene_names[pval == 0]
            #
            # logratio_indices = np.flip(np.argsort(logratio))
            #
            # gene_names = gene_names[logratio_indices].tolist()
            # gene_ids = gene_ids[logratio_indices].tolist()
            # logratio = logratio[logratio_indices]
            # top_gene_markers[cluster_idx] = gene_names[:5]
            #
            # gene_indices += [features.index(ii) for ii in gene_ids[:NUM_OF_MARKER_GENES]]
            # gene_horiz_lines += [gene_horiz_lines[-1] + NUM_OF_MARKER_GENES]
        print('Num of repetition in markers:', NUM_OF_MARKER_GENES * K - len(set(gene_indices)), end="\n\n")

        cells_indices = []
        for cls_idx, cluster_indices in enumerate(clusters_indices):
            cells_indices += cluster_indices
            print(f'{cls_idx}. num cells: {len(cluster_indices)}')
            cell_vertical_lines += [cell_vertical_lines[-1] + len(cluster_indices)]

        arr_heatmap = filtered_cells.counts[cells_indices][:, gene_indices]
        heatmap = np.zeros_like(arr_heatmap)
        heatmap[arr_heatmap > 1] = 1

        ## 19.3.21 using zscore
        arr_heatmap = scipy.stats.zscore(arr_heatmap, axis=0, ddof=1)

        if print_binary_heatmap:

            fig, ax = plt.pyplot.subplots(1)
            fig.set_size_inches(10, 5)

            sb.heatmap(heatmap)

            for ver_line in gene_horiz_lines[1:-1]:
                ax.axvline(ver_line, color='yellow')

            for hor_line in cell_vertical_lines[1:-1]:
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
            # sb_out = sb.heatmap(arr_heatmap.T, vmin=-1, vmax=1, cmap=create_cmap())
            # sb_out = sb.heatmap(arr_heatmap.T, vmin=-1, vmax=1, cmap=cmap,
            #                     cbar_kws=dict(shrink=0.7, location="left", use_gridspec=False))
            sb_out = sb.heatmap(arr_heatmap.T, vmin=-1, vmax=1, cmap=cmap)
            cbar = sb_out.collections[0].colorbar
            cbar.set_ticks([-1, 1])
            cbar.set_ticklabels([-1, 1])
            sb_out.set(xticklabels=[])
            sb_out.set(yticklabels=[])
            sb_out.tick_params(bottom=False, left=False)
            for hor_line in gene_horiz_lines[1:-1]:
                ax.axhline(hor_line, color='white')

            for ver_line in cell_vertical_lines[1:-1]:
                ax.axvline(ver_line, color='white')

            ax.set_title(f"Kmeans k={K}")
            ax.set_xlabel('Clusters')
            ax.set_ylabel('Markers')

            ##### save second heatmap
            fig.savefig(join(OUTPUT_PATH, f'density_heatmap_{K}.png'))

            with open(join(OUTPUT_PATH, f'top_5_markers_K_{K}.txt'), 'w') as text_file:
                for key in top_gene_markers.keys():
                    for gene in top_gene_markers[key]:
                        text_file.write(' ' + gene + '\n')
                    text_file.write('\n')

