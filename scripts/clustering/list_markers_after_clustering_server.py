"""
After running Analyze_kmeans_clusters.py tou have the marker gene list.
so to create csv run this script on the output of running Analyze_kmeans_clusters.py.
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

OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/markers/11.8.21'
CLUSTERING_ANALYSIS_PATH = fr'/storage/md_keren/shitay/outputs/clustering/myeloid/cluster_analysis/cluster_analysis_11.8.21' # immune/cytotoxic_t_cells/myeloid

if __name__ == '__main__':


    NUM_OF_MARKER_GENES = 30

    create_folder(OUTPUT_PATH)
    # For each K creates csv file
    for K in range(2, 16):
        print(f"working on {K}")
        cluster_dir = join(OUTPUT_PATH, f'markers_cluster_{K}')
        create_folder(cluster_dir)

        print(f"Current K = {K}")
        clustering_analysis = pickle.load(open(join(CLUSTERING_ANALYSIS_PATH, f'cluster_analysis_k_{K}.pkl'), 'rb'))

        print('num clusters: ', len(clustering_analysis))
        for cluster in clustering_analysis:
            cls_idx = cluster['cluster id']
            markers = cluster['markers']
            markers.to_csv(join(cluster_dir, f'markers_cluster_{cls_idx+1}.csv'), index=False)



    #  5.4.21; create new format of fd using only one csv file and saves it in the main OUTPUT directory
    # sheet (tab) for every K solution and column for each cluster
    writer_path = join(OUTPUT_PATH, 'markers_all_solutions.xlsx')
    writer = pd.ExcelWriter(writer_path, engine='xlsxwriter')
    for K_dir in sorted(os.listdir(OUTPUT_PATH)):
        working_dir = join(OUTPUT_PATH, K_dir)
        K = K_dir.split('_')[-1]
        print(K)

        new_table = {}
        for cluster_marker_list in sorted(os.listdir(working_dir)):
            cluster_idx = cluster_marker_list.split('_')[-1].replace('.csv', '')
            print(cluster_marker_list)
            df = pd.read_csv(join(working_dir, cluster_marker_list))
            new_table[int(cluster_idx)] = list(df['gene names'])[:100]
            if len(new_table[int(cluster_idx)]) < 100:
                new_table[int(cluster_idx)] = new_table[int(cluster_idx)] + [None] * (
                            100 - len(new_table[int(cluster_idx)]))
        print([len(vv) for vv in list(new_table.values())])
        new_sheet = pd.DataFrame(new_table)
        new_sheet = new_sheet.reindex(sorted(new_sheet.columns), axis=1)
        new_sheet.to_excel(writer, sheet_name=K, index=False)
    writer.save()