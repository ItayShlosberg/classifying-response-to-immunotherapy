
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib
import scrublet as scr
import pickle
from DL.data_creation import *
from DL.data_conversions.txt_to_python_structures import *
from os.path import join


SAMPLES = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples'
SAMPLE_ID = 'M114'


def run_scrub_over_all_samples():
    data_path = join(SAMPLES, SAMPLE_ID, 'RNA_sample.pkl')
    rna_sample = pickle.load(open(data_path, 'rb'))
    # cells, gene_names, patients_information = ex(data_path)
    # filter_genes_by_variance(cells, gene_names, required_variance=6)

    scrub = scr.Scrublet(rna_sample.counts)
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


if __name__ == '__main__':
    run_scrub_over_all_samples()