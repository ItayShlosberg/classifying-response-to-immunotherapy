# # import math
# # import numpy as np
# # import Bio
# # from utilities.general_helpers import *
# # from DL.data_loading import *
# # from utilities.dataset import *
# from matplotlib import pyplot
# import numpy as np
# import scipy
# import pickle
# import matplotlib
# import scrublet as scr
# from DL.data_loading import extract_data_from_pickle as ex
# import pickle
# from DL.data_creation import *
# import matplotlib.pyplot as plt
# from DL.matlab_env.txt_to_python_structures import *
# import random
#
# from scipy.stats import pearsonr
#
# def plot2():
#     hist, bin_edges = np.histogram(cells.sum(axis=1), density=True)
#
#
#
#     hist
#     hist.sum()
#     np.sum(hist * np.diff(bin_edges))
#
#
#
#
#     rng = np.random.RandomState(10)  # deterministic random data
#
#     a = np.hstack((rng.normal(size=1000),
#
#                    rng.normal(loc=5, scale=2, size=1000)))
#
#     # _ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
#     _ = plt.hist(cells.sum(axis=1), bins='auto')  # arguments are passed to np.histogram
#
#     plt.title("Histogram with 'auto' bins")
#
#     plt.show()
#
#
#     breakpoint = 0
#
#
# def distribution(counts, gene_names):
#     colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'lime', 'lavender', 'darkred', 'olive']
#     # for cluster in np.unique(clusters_labels):
#     #     Xi = [X[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
#     #     Yi = [Y[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
#     #     color = colors[cluster]
#     #     plt.plot(Xi, Yi, 'ro', color=color)
#
#     # sums = np.sort(counts.sum(axis=1))
#
#     indexes =
#     random.choices(range(len(counts)), k=50)
#
#
#     sums = counts.sum(axis=1)
#     mito = counts[:, [s.startswith('MT-') for s in gene_names]].sum(axis=1)
#     ones = np.arange(len(sums))
#     print(ones)
#     print(sums[:10])
#     print(sums.shape)
#     print(ones.shape)
#     # plt.plot([1, 2], [1, 4], 'ro')
#     plt.plot(mito[:], sums[:], 'ro')
#     # plt.plot([3, 4], [9, 16], 'ro', color='green')
#     plt.title('Distribution')
#     plt.show()
#
#
# if __name__ == '__main__':
#     sample = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples\M100\enrich_RNA_sample.pkl'
#     data_path = r'..\DATA\rna_seq200k\dataset.pkl'
#
#     data = pickle.load(open(sample, 'rb'))
#     genes = data.gene_names
#     counts = data.counts.T  # Now: Cells (rows) X Genes (Columns)
#     distribution(counts, genes)