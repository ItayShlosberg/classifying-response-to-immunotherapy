""""
Milestone num. 5 issue num. 4 - changing the method of data storing.
We will go over all samples and store the row data.
"""
import sys
# ------- SERVER EXTENSIONS ---------
lib =  r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities/droplet_dataset'
lib2 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/utilities'
lib3 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/data_analysis'
lib4 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy'
lib5 = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/DL'
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
import pandas as pd
import DL
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
import os
from os.path import join

SAMPLES_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\PYTHON_OBJECTS'
OUTPUT_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\ROW_DATA'



samples = [subfolder[:-4] for subfolder in os.listdir(SAMPLES_PATH)]

for sample in samples:
    print(f"Working on {sample}")
    # Extracts one of the samples from PC
    rna_sample = extract_droplet_data_from_pickle(join(SAMPLES_PATH, f'{sample}.pkl'))

    rna_sample.save_row_data(join(OUTPUT_PATH, f'{sample}.pkl'))