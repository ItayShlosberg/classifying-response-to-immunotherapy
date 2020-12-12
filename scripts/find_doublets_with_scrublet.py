"""
check for doublets over all samples using Scrublet.
Save a summary and  updated PKLs with is_doublet information for each cell.
"""

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


SAMPLES = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\apoptosis\10.12.20'
OUTPUT_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\10.12.20'


def extract_sample(sample_id):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    data_path = join(SAMPLES, sample_id, f'{sample_id}.pkl')
    rna_sample = extract_droplet_data_from_pickle(data_path)
    print(f'sample id {sample_id}')
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')
    return rna_sample


def run_scrub_over_all_samples():
    """
    Main function, run scrub over all samples, for each sample finds doublets and update rna_seq python object's
    meta data. save the updated object in output path.
    Additionally, save a summary Excel conclusion (DF) of all sample.
    """
    summary_df = pd.DataFrame(columns=['sample name',
                                       'number of cells',
                                       'number of scrub doublets',
                                       'p_doublets'])
    samples = [subfolder for subfolder in os.listdir(SAMPLES)]
    create_folder(OUTPUT_PATH)
    for sample in [s for s in samples if not 'csv' in s]:
        # Extracts one of the samples from PC
        rna_sample = extract_sample(sample)
        scrub = scr.Scrublet(rna_sample.counts)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        print(f'num of doublets {sum(predicted_doublets)}')
        print(f'percentage of doublets {sum(predicted_doublets)/rna_sample.number_of_cells}')
        rna_sample.cells_information.setattr('is_doublet', np.arange((rna_sample.number_of_cells)), False)
        rna_sample.cells_information.setattr('is_doublet', np.where(predicted_doublets)[0], True)


        # Append to DF
        summary_df = summary_df.append(pd.DataFrame([[sample,
                                                      rna_sample.number_of_cells,
                                                      sum(predicted_doublets),
                                                      sum(predicted_doublets)/rna_sample.number_of_cells]],
                                       columns=summary_df.columns))

        create_folder(join(OUTPUT_PATH, sample))
        pickle.dump((rna_sample), open(join(OUTPUT_PATH, sample, f'{sample}.pkl'), 'wb'))
    summary_df.to_excel(join(OUTPUT_PATH, r'scrublet_summary.xlsx'), index=False)


if __name__ == '__main__':
    run_scrub_over_all_samples()

