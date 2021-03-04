""""
Milestone num. 5 issue num. 4 - changing the method of data storing.
We will go over all samples and store the row data.
"""

import pandas as pd
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
import os
from os.path import join

SAMPLES_PATH = fr'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples_10.12.20'
OUTPUT_PATH = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\ROW_DATA'


samples = [subfolder for subfolder in os.listdir(SAMPLES_PATH)]

for sample in samples:
    print(f"Working on {sample}")
    # Extracts one of the samples from PC
    rna_sample = extract_droplet_data_from_pickle(join(SAMPLES_PATH, sample, f'{sample}.pkl'))

    rna_sample.save_row_data(join(OUTPUT_PATH, f'{sample}.pkl'))