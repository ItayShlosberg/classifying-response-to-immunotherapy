import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *


SAMPLES = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'
OUTPUT = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\Cohort\cohort_all_samples_3.2.21.pkl'

cohort = build_cohort(SAMPLES)#, gene_path =r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\gene_cohort.pkl')

# pickle.dump((cohort), open(OUTPUT, "wb"))
pickle.dump((cohort), open(OUTPUT, 'wb'), protocol=4)#

cc = pickle.load(open(OUTPUT, 'rb'))

_breakpoint = 0
# samples = [subfolder for subfolder in os.listdir(SAMPLES) if not 'csv' in subfolder]
# rna_sample = extract_droplet_data_from_pickle(join(SAMPLES, samples[3]))
# rna_sample = rna_sample[[not aa for aa in rna_sample.cells_information.getattr('should_be_removed')]]
# rna_sample.normalize_data()
#
#
# sample_cell_idx = 50
# cell_barcode = rna_sample.barcodes[sample_cell_idx]
# cohort_cell_index = cohort.barcodes.index(cell_barcode)
#
#
# sample_gene_idx = rna_sample[sample_cell_idx][0].argmax()
# feature_name = rna_sample.features[sample_gene_idx]
# cohort_gene_index = cohort.features.index(feature_name)
#
# check = cohort.counts[cohort_cell_index, cohort_gene_index] == rna_sample.counts[sample_cell_idx, sample_gene_idx]
#
# _breakpoint = 0
