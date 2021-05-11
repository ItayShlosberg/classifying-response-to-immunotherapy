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


COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/5.21/immune_cells_4k_genes.pkl'
OUTPUT_PATH = r'/storage/md_keren/shitay/outputs/clustering/elbow/10.5.21_cohort/elbow.pkl'
KMEANS_ROW_CLUSTERS_PATH = r'/storage/md_keren/shitay/outputs/clustering/kmeans/10.5.21/row_kmeans'
KMEANS_FILE_NAME = r'kmeans_immune_cells_4k_genes'  # excluding the suffix: '_k_num.pkl'

if __name__ == '__main__':

    # Loads cohort
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    D = cdist(cohort.counts, cohort.counts, 'correlation')

    Dis_t = 0
    for i in range(D.shape[0]):
        for j in range(i, D.shape[0]):
            Dis_t += D[i, j]

    All_Dist_b = {}
    for K in range(2, 16):
        print(f'K = {K}, ', end='')
        kmeans_clusters_path = join(KMEANS_ROW_CLUSTERS_PATH, KMEANS_FILE_NAME+f'_k_{K}.pkl')
        solution_k_clusters = pickle.load(open(kmeans_clusters_path, 'rb'))
        Dist_b = 0
        for idx in range(len(solution_k_clusters['clusters'])):
            cluster_indices = solution_k_clusters['clusters'][idx]
            other_clusters_indices = [ii for ii in range(cohort.number_of_cells) if not ii in cluster_indices]

            for ii in cluster_indices:
                for jj in other_clusters_indices:
                    Dist_b += D[ii, jj]

        All_Dist_b[K] = Dist_b

    pickle.dump(({'Ks_Dist_b': All_Dist_b, 'Dis_t': Dis_t}), open(OUTPUT_PATH, 'wb'))