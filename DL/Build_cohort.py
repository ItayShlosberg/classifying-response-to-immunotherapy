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



ROW_DATA_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/ROW_DATA'
CELL_INFORMATION_PATH = r'/storage/md_keren/shitay/Data/inferCNV_data/update_runs/26.6.21'
OUTPUT = r'/storage/md_keren/shitay/Data/droplet_seq/cohort/non_normalized/6.21'
FILE_NAME = r'cohort_non_normalized_26.6.21.pkl'
TO_NORMALIZE = False

if __name__ == '__main__':


    print(f'Run Build_cohort.py')
    print('Parameters:')
    print(f'ROW_DATA_PATH {ROW_DATA_PATH}')
    print(f'CELL_INFORMATION_PATH {CELL_INFORMATION_PATH}')
    print(f'OUTPUT {OUTPUT}')
    print(f'FILE_NAME {FILE_NAME}')
    print(f'TO_NORMALIZE {TO_NORMALIZE}')

    gene_ids = build_cohort_gene_list(CELL_INFORMATION_PATH)

    cohort = build_cohort(ROW_DATA_PATH, CELL_INFORMATION_PATH, gene_ids, to_normalize=TO_NORMALIZE)

    create_folder(OUTPUT)
    pickle.dump((cohort), open(join(OUTPUT, FILE_NAME), 'wb'), protocol=4)
    print(f'OUTPUT has been saved in {join(OUTPUT, FILE_NAME)}')
