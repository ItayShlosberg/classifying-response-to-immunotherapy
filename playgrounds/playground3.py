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

import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.general_helpers import *
from utilities.droplet_dataset import *
from utilities.droplet_dataset import loading_sample




# SAMPLES = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'

# OUTPUT = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\Cohort\cohort_all_samples_3.2.21.pkl'
# ss = r'C:\Users\itay\Desktop\New folder'
#
#
# # cohort = build_cohort(ss)#, gene_path =r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\gene_cohort.pkl')
#
#
# # pickle.dump((cohort), open(OUTPUT, "wb"))
# # pickle.dump((cohort), open(OUTPUT, 'wb'), protocol=4)#
# # tt = cohort + cohort
# # cc = pickle.load(open(OUTPUT, 'rb'))
#
# _breakpoint = 0
# # samples = [subfolder for subfolder in os.listdir(SAMPLES) if not 'csv' in subfolder]
# # rna_sample = extract_droplet_data_from_pickle(join(SAMPLES, samples[3]))
# # rna_sample = rna_sample[[not aa for aa in rna_sample.cells_information.getattr('should_be_removed')]]
# # rna_sample.normalize_data()
# #
# #
# # sample_cell_idx = 50
# # cell_barcode = rna_sample.barcodes[sample_cell_idx]
# # cohort_cell_index = cohort.barcodes.index(cell_barcode)
# #
# #
# # sample_gene_idx = rna_sample[sample_cell_idx][0].argmax()
# # feature_name = rna_sample.features[sample_gene_idx]
# # cohort_gene_index = cohort.features.index(feature_name)
# #
# # check = cohort.counts[cohort_cell_index, cohort_gene_index] == rna_sample.counts[sample_cell_idx, sample_gene_idx]
# #
# # _breakpoint = 0

# ------- SERVER EXTENSIONS ---------








#
# ROW_SAMPLES_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/ROW_DATA/'
# SAMPLES_INFORMATION_PATH = r'/storage/md_keren/shitay/Data/inferCNV_data/update_runs/4.3.21/'
#
# sample = 'M106.pkl'
#
# rna_sample = loading_sample(row_data_path=join(ROW_SAMPLES_PATH, f'{sample}'),
#                             cells_information_path=join(SAMPLES_INFORMATION_PATH, f'{sample}'))
#
#
# print(rna_sample.number_of_cells)
# # ------- SERVER EXTENSIONS ---------




CELL_INFORMATION_PATH = r'C:\Users\itay\Desktop\test\cells_information'
ROW_DATA_PATH = r'C:\Users\itay\Desktop\test\row_data'
OUTPUT = r'C:\Users\itay\Desktop\test\OUTPUT\cohort.pkl'

gene_ids = build_cohort_gene_list(CELL_INFORMATION_PATH)

cohort = build_cohort(ROW_DATA_PATH, CELL_INFORMATION_PATH, gene_ids)
pickle.dump((cohort), open(OUTPUT, 'wb'), protocol=4)#


cc = pickle.load(open(OUTPUT, 'rb'))
print(cc.number_of_cells)
_breakoint = 0