from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import sys
lib = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\utilities\droplet_dataset'
lib2 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\utilities'
lib3 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\data_analysis'
lib4 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy'
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *
import math
import matplotlib as plt
from sklearn.manifold import TSNE
import statsmodels as sm
import scipy.stats as stats
from Bio.Cluster import kcluster
from sklearn.manifold import TSNE
import pickle
from scipy.stats import pearsonr
import statsmodels as sm
from utilities.ML_environment import multipletests_fdr




######## LOADING DATA ########
sample_id = 'M104'
ROW_SAMPLES_PATH = fr'D:\Technion studies\Keren Laboratory\Data\droplet_seq\ROW_DATA'
SAMPLES_INFORMATION_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\4.3.21'


rna_sample = loading_sample(row_data_path=join(ROW_SAMPLES_PATH, f'{sample_id}.pkl'),
                                cells_information_path=join(SAMPLES_INFORMATION_PATH, f'{sample_id}.pkl'))


######## CLEAN & NORMALIZE DATA ########
VARIANCE = 3
rna_sample = rna_sample.filter_cells_by_property("should_be_removed", False)
rna_sample.normalize_data()
sample = rna_sample.filter_genes_by_variance(VARIANCE, in_place=True).counts
print(f'Sample size after filtering variance: {sample.shape}')


######## KMEANS ########

# user_function = lambda point1, point2: abs(point1[0] - point2[0]) # Manhattan
# user_function = lambda point1, point2: pearsonr(point1, point2) # pearson
def user_function(point1, point2):
#     return abs(point1[0] - point2[0])
    return 1 - pearsonr(point1, point2)[0]


metric = distance_metric(type_metric.USER_DEFINED, func=user_function)

# create K-Means algorithm with specific distance metric
initial_centers = kmeans_plusplus_initializer(sample, 4).initialize()

kmeans_instance = kmeans(sample, initial_centers, metric=metric)

# run cluster analysis and obtain results
kmeans_instance.process()
clusters = kmeans_instance.get_clusters()


from utilities.ML_environment import  find_markers_in_clusters

find_markers_in_clusters(rna_sample, clusters)

