"""
After kmeans was done, you can run this script to find markers.
after finding markers run list_markers_after_clustering_server.py on the output of this script to get more convenient
format for the output.

the script get as an input the index of k (the kmeans k) and save the markers gene names the pvals and logFC values.
"""

from os.path import join
import sys
# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/scripts'
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
from utilities.droplet_dataset import *
from utilities.ML_environment import find_markers_in_clusters
from utilities.general_helpers import create_folder


OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/cluster_analysis/cluster_analysis_10.5.21'
FILTERED_CELLS_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/5.21/immune_cells_10.5.21.pkl'
KMEANS_ROW_CLUSTERS_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/10.5.21/row_kmeans'
KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'

if __name__ == '__main__':

    K = int(sys.argv[1])
    print(f'Loading cohort file:\n{FILTERED_CELLS_PATH}')
    filtered_cells = pickle.load(open(FILTERED_CELLS_PATH, 'rb'))

    kmeans_file_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME + f'_k_{K}.pkl')
    print(f'Loading kmeans file:\n{kmeans_file_path}')
    clusters = pickle.load(open(kmeans_file_path, 'rb'))['clusters']

    markers = find_markers_in_clusters(filtered_cells, clusters)
    create_folder(OUTPUT_PATH)
    pickle.dump((markers), open(os.path.join(OUTPUT_PATH, f'cluster_analysis_k_{K}.pkl'), 'wb'))