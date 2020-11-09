from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib
import scrublet as scr
from DL.data_loading import extract_data_from_pickle as ex
import pickle
from DL.data_creation import *
from DL.matlab_env.txt_to_python_structures import *


sample = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples\M114\enrich_RNA_sample.pkl'
data_path = r'..\DATA\RNAseq_DATA.p'


cells = pickle.load(open(sample, 'rb')).cells.T
# cells, gene_names, patients_information = ex(data_path)
# filter_genes_by_variance(cells, gene_names, required_variance=6)

scrub = scr.Scrublet(cells)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
print(doublet_scores)
print(predicted_doublets)
print(len(doublet_scores))
print(len(predicted_doublets))
print(sum(predicted_doublets))


scrub.plot_histogram()


# Get 2 - D embedding to visualize the results¶
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
# print('Running tSNE...')
# scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
# print('Running ForceAtlas2...')
# scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))

scrub.plot_embedding('UMAP', order_points=True)

pyplot.show()