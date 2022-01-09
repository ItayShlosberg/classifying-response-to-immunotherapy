"""
Print Heatmap of gene markers, after finding markers of clusters.
Done for all Ks.
Takes the first NUM_OF_MARKER_GENES (default is 30) gene markers with the highest log_FC.
if there are >= NUM_OF_MARKER_GENES gene markers with pval=0 it will consider only them.

"""

from os.path import join
import numpy as np
import pandas as pd
import pickle
import os
import numpy as np
import yaml
import os
import pandas
from collections import Counter
from sklearn.cluster import KMeans
import scipy
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib
import seaborn as sb
from shutil import copyfile
import matplotlib as plt

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
from utilities.general_helpers import create_folder
plt.rcParams['figure.dpi'] = 1000
from utilities.droplet_dataset import *


######################################################## Params #######################################################


OUTPUT_PATH = fr'/storage/md_keren/shitay/outputs/clustering/myeloid/summary/subcohort_1.1.22_run_1.1.22'

# FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/CD8_normalized_26.6.21.pkl'
FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/subcohort/normalized/1.1.22/subcohort_immune_cells_normalized_1.1.22_protein_coding_genes.pkl'
SUBSET = 'MYELOIDS'    # None - all cells, MYELOIDS/CYTOTOXIC_T_CELLS

KMEANS_ROW_CLUSTERS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/subcohort_1.1.22_run_1.1.22/row_kmeans'
KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'

# if you used finalize_clusters.ipynb and you use barcode_mapping.csv instead of KMEANS_ROW file it should be the barcode_mapping path else None
IS_USING_BARCODE_MAPPING = r'/storage/md_keren/shitay/outputs/clustering/myeloid/summary/subcohort_1.1.22_run_1.1.22/subcohort_myeloid_1.1.22_clusters_mapping.csv'
CLUSTERING_ANALYSIS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/myeloid/summary/subcohort_1.1.22_run_1.1.22'
# IS_USING_BARCODE_MAPPING = r'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/subcohort_1.1.22_run_1.1.22/barcode_mapping/kmeans_immune_cells_4k_genes_k_10_clusters_mapping.csv'
# CLUSTERING_ANALYSIS_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/cluster_analysis/subcohort_1.1.22_run_1.1.22'

SPECIFIC_K = 9 # None if you want all. (if build with finalize_clusters.ipynb use the K have finished with)
NUM_OF_MARKER_GENES = 30
COLOR_MAP_PATH = r'/storage/md_keren/shitay/outputs/clustering/immune/heatmap/colorbar.pkl'


#######################################################################################################################





def data_loader():
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))
    filtered_cells = filtered_cells.filter_cells_by_property('is_immune', True)
    if SUBSET:
        filtered_cells = get_requested_subset(filtered_cells, SUBSET)
    return filtered_cells


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


def draws_heatmap(arr_heatmap, gene_horiz_lines, cell_vertical_lines, K, cmap):
    # Draws heatmap
    fig, ax = plt.subplots(1)
    fig.set_size_inches(10, 5)
    sb_out = sb.heatmap(arr_heatmap.T, vmin=-1, vmax=1, cmap=cmap);
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
    ax.set_title(f"Kmeans k={K}");
    ax.set_xlabel('Clusters');
    ax.set_ylabel('Markers');
    return fig


def extract_heatmap_values_from_cohort(filtered_cells, cells_indices, gene_indices):
    # Gets heatmap values by giving indices of cells and indices of gene markers.
    arr_heatmap = filtered_cells.counts[cells_indices][:, gene_indices]
    heatmap = np.zeros_like(arr_heatmap)
    heatmap[arr_heatmap > 1] = 1

    ## 19.3.21 using zscore
    arr_heatmap = scipy.stats.zscore(arr_heatmap, axis=0, ddof=1)

    return arr_heatmap


def arrange_heatmap_cells_indices(filtered_cells, clusters_indices, mapping_clusters_to_marker_indices):
    # Loops over clusters indexes and builds a list containing the indexes in the right order
    cell_vertical_lines = [0]  # for heatmap
    cells_indices = []
    for cls_idx, cluster_indices in enumerate(clusters_indices):
        print(f'{cls_idx}. num cells: {len(cluster_indices)}')
        top_markers_indices = mapping_clusters_to_marker_indices[cls_idx]
        mean_values = filtered_cells.counts[cluster_indices][:, top_markers_indices].mean(axis=1)
        cluster_indices = np.array(cluster_indices)[np.flip(np.argsort(mean_values))].tolist()
        cells_indices += cluster_indices
        cell_vertical_lines += [cell_vertical_lines[-1] + len(cluster_indices)]

    return cell_vertical_lines, cells_indices


def arrange_heatmap_genes_indices(filtered_cells, clustering_analysis):
    # Loops over markers DFs and takes the indexes of marker genes from cohort
    features = filtered_cells.features
    print(clustering_analysis[0].keys())
    gene_indices = []
    gene_horiz_lines = [0]  # for heatmap
    print('num markers for each cluster: ', NUM_OF_MARKER_GENES)
    print('num markers for all clusters: ', NUM_OF_MARKER_GENES * len(clustering_analysis))
    print('num clusters: ', len(clustering_analysis))
    mapping_clusters_to_marker_indices = {}
    for idx, cluster_dic in enumerate(clustering_analysis):
        # cluster_idx = cluster_dic['cluster id']
        cluster = cluster_dic['markers']

        cluster = cluster.sort_values(by=['log_FC'], ascending=False)
        gene_ids = cluster['features'].tolist()
        top_markers_indices = [features.index(ii) for ii in gene_ids[:NUM_OF_MARKER_GENES]]
        mapping_clusters_to_marker_indices[idx] = top_markers_indices
        gene_indices += top_markers_indices
        gene_horiz_lines += [gene_horiz_lines[-1] + len(gene_ids[:NUM_OF_MARKER_GENES])]
    print('Num of repetition in markers:', NUM_OF_MARKER_GENES * K - len(set(gene_indices)), end="\n\n")
    return gene_horiz_lines, gene_indices, mapping_clusters_to_marker_indices


def get_clusters_indices(df):
    mapping = list(zip(filtered_cells.samples, filtered_cells.barcodes))
    clusters = sorted(set(df['Cluster']))
    clusters_indexes = []
    for cluster_idx in clusters:
        barcodes_list = df[df['Cluster'] == cluster_idx]['Barcode'].tolist()
        sample_list = df[df['Cluster'] == cluster_idx]['Sample'].tolist()
        cell_idxs = [mapping.index(pair_identifier) for pair_identifier in zip(sample_list, barcodes_list)]
        clusters_indexes.append(cell_idxs)
    return clusters_indexes


if __name__ == '__main__':
    # Loads cohort (and cmap a file for barcolor of heatmap)
    print(f'Params:\nOUTPUT_PATH - {OUTPUT_PATH}\n'
          f'FILTERED_CELLS_PATH - {FILTERED_CELLS_PATH}\n'
          f'KMEANS_ROW_CLUSTERS_PATH - {KMEANS_ROW_CLUSTERS_PATH}\n'
          f'CLUSTERING_ANALYSIS_PATH - {CLUSTERING_ANALYSIS_PATH}')

    print('Loads cohort')
    filtered_cells = data_loader()
    cmap = pickle.load(open(COLOR_MAP_PATH, 'rb'))
    create_folder(OUTPUT_PATH)

    Ks = [SPECIFIC_K] if SPECIFIC_K else list(range(2, 16))
    for K in Ks:
        print(f'Drawing heatmap of K={K}')
        # Loads marker file
        clustering_analysis = pickle.load(open(join(CLUSTERING_ANALYSIS_PATH, f'cluster_analysis_k_{K}.pkl'), 'rb'))

        # loads kmeans_row data
        if IS_USING_BARCODE_MAPPING:
            df = pd.read_csv(IS_USING_BARCODE_MAPPING)
            clusters_indices = get_clusters_indices(df)
        else:
            kmeans_file_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME + f'_k_{K}.pkl')
            print(f'Loading kmeans file:\n{kmeans_file_path}')
            clusters_indices = pickle.load(open(kmeans_file_path, 'rb'))['clusters']

        # Loops over markers DFs and arrange heatmap's  gene marker indices from cohort
        gene_horiz_lines, gene_indices, mapping_clusters_to_marker_indices = arrange_heatmap_genes_indices(filtered_cells, clustering_analysis)

        # Loops over clusters indexes and builds a list containing the indexes in the right order
        cell_vertical_lines, cells_indices = arrange_heatmap_cells_indices(filtered_cells, clusters_indices, mapping_clusters_to_marker_indices)

        # Gets heatmap values by giving indices of cells and indices of gene markers.
        arr_heatmap = extract_heatmap_values_from_cohort(filtered_cells, cells_indices, gene_indices)

        fig = draws_heatmap(arr_heatmap, gene_horiz_lines, cell_vertical_lines, K, cmap)
        # OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/summaries/26.6.21/heatmap.png'

        fig.savefig(join(OUTPUT_PATH, f'heatmap_k_{K}.png'))