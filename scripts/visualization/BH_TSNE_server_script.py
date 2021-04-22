"""
Running BH_TSNE.
Conclusion was to use variance 0.315 (~4K genes)
"""

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
import numpy as np
import pickle
from os.path import join
from sklearn.decomposition import PCA
from bhtsne import tsne




OUTPUT_DIR = r'/storage/md_keren/shitay/outputs/TSNE'
FILE_NAME = r'all_cells_bhtsne_21.4.21.pkl'

### cohort should be variance filtered
COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/normalized/cohort_var0.312.pkl'



if __name__ == '__main__':

    print("Running TSNE")
    print(f'File {COHORT_PATH}')
    # print(f'ARG {sys.argv[1]}')

    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    # cohort = cohort.filter_cells_by_property('is_immune', True)
    print(f"Counts shape {cohort.counts.shape}")

    PCs = PCA(n_components=10).fit_transform(cohort.counts)

    print(f"PCs shape {PCs.shape}")

    # cells_embedded = TSNE(n_components=2, perplexity=perplexity, random_state=21).fit_transform()
    cells_embedded = tsne(PCs)

    print(f"TSNE output size {cells_embedded.shape}")
    pickle.dump((cells_embedded), open(join(OUTPUT_DIR, FILE_NAME), 'wb'))


