"""
Elbow - estimating the variance explained by a given solution
Used with the following Notebook 'Elbow_plot_&_robustness analysis clustering', run first this script and then
link the output to the notebook.
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
import sys
# ------- SERVER EXTENSIONS ---------
lib = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
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
from scipy.spatial.distance import cdist
from utilities.general_helpers import create_folder
import os.path
from  os.path import join
import os



# COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/6.21/myeloid_normalized_26.6.21_4k_genes.pkl'
COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/subcohort/normalized/1.1.22/subcohort_immune_cells_normalized_1.1.22_4k_genes.pkl'

OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/elbow/subcohort_1.1.22_run_1.1.22'
OUTPUT_FILE_NAME = r'elbow.pkl'
KMEANS_ROW_CLUSTERS_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/kmeans/subcohort_1.1.22_run_1.1.22/row_kmeans'

KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'
SAVE_DISTANCE_MATRIX = True
SUBSET = 'MYELOIDS'    # None - all cells, MYELOIDS/CYTOTOXIC_T_CELLS


def data_loader():
    filtered_cells = pickle.load(open(COHORT_PATH, 'rb'))
    filtered_cells = filtered_cells.filter_cells_by_property('is_immune', True)
    if SUBSET:
        filtered_cells = get_requested_subset(filtered_cells, SUBSET)
    return filtered_cells


if __name__ == '__main__':

    # Loads cohort
    print('Running elbow script')
    create_folder(OUTPUT_PATH)
    print(f'Loading cohort from:\n{COHORT_PATH}')
    cohort = data_loader()

    if SAVE_DISTANCE_MATRIX:
        CDIST_PATH = join(OUTPUT_PATH, r'cdist.pkl')
        if os.path.isfile(CDIST_PATH):
            D = pickle.load(open(CDIST_PATH, 'rb'))
        else:
            D = cdist(cohort.counts, cohort.counts, 'correlation')
            pickle.dump((D), open(CDIST_PATH, 'wb'), protocol=4)

    print('Calculating Dist_t')
    Dis_t = 0
    for i in range(D.shape[0]):
        for j in range(i, D.shape[0]):
            Dis_t += D[i, j]

    All_Dist_b = {}
    for K in range(2, 16):
        print(f'K = {K}, ', end='')
        kmeans_clusters_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME+f'_k_{K}.pkl')
        print(f'Taking file {kmeans_clusters_path}')
        solution_k_clusters = pickle.load(open(kmeans_clusters_path, 'rb'))
        Dist_b = 0
        for idx in range(len(solution_k_clusters['clusters'])):
            cluster_indices = solution_k_clusters['clusters'][idx]
            other_clusters_indices = [ii for ii in range(cohort.number_of_cells) if not ii in cluster_indices]

            for ii in cluster_indices:
                for jj in other_clusters_indices:
                    Dist_b += D[ii, jj]

        All_Dist_b[K] = Dist_b

    pickle.dump(({'Ks_Dist_b': All_Dist_b, 'Dis_t': Dis_t}), open(join(OUTPUT_PATH, OUTPUT_FILE_NAME), 'wb'))
