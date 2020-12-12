from utilities.droplet_dataset import *
from utilities.general_helpers import *
from os.path import join
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib
import scrublet as scr
import pickle
from DL.data_creation import *
from DL.matlab_env.txt_to_python_structures import *
from os.path import join
from utilities.droplet_dataset import *
from utilities import *
import numpy as np
import pickle
import pickle
import pandas as pd
from DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from utilities.general_helpers import *

SOURCE_PATH = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples_5.12.20'
INTERMEDIATE_PATH = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples_10.12.20'
# TARGET_PATH = r'D:\Technion studies\Keren Laboratory\Data\garbage\target'

def first_direction():
    samples = [subfolder for subfolder in os.listdir(SOURCE_PATH)]
    create_folder(INTERMEDIATE_PATH)

    for sample_id in [s for s in samples if not 'csv' in s]:
        print(sample_id)
        data_path = join(SOURCE_PATH, sample_id, f'{sample_id}.pkl')
        rna_sample = extract_droplet_data_from_pickle(data_path)
        sample_v2 = RNAseq_Sample(rna_sample.counts, rna_sample.gene_names, rna_sample.barcodes, rna_sample.features)
        # cell_inf_list = Cell_Inf_List(sample_v2.number_of_cells)
        # cell_inf_list.cells_information_list = rna_sample.cells_information
        # sample_v2.cells_information = cell_inf_list


        create_folder(join(INTERMEDIATE_PATH, sample_id))
        pickle.dump((sample_v2), open(join(INTERMEDIATE_PATH, sample_id, f'{sample_id}.pkl'), 'wb'))
        _breakpoint = 0


def back_direction():
    samples = [subfolder for subfolder in os.listdir(INTERMEDIATE_PATH)]
    create_folder(TARGET_PATH)

    for sample_id in [s for s in samples]:
        print(sample_id)
        data_path = join(INTERMEDIATE_PATH, sample_id, f'{sample_id}.pkl')
        rna_sample = extract_droplet_data_from_pickle(data_path)

        original_sample = RNAseq_Sample(rna_sample.counts, rna_sample.gene_names, rna_sample.samples_list, rna_sample.ens_id_list)
        original_sample.cells_information = rna_sample.cells_information

        create_folder(join(TARGET_PATH, sample_id))
        pickle.dump((original_sample), open(join(TARGET_PATH, sample_id, f'{sample_id}.pkl'), 'wb'))
        _breakpoint = 0

if __name__ == '__main__':
    first_direction()
    # back_direction()