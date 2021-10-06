"""
check for doublets over all samples using Scrublet.
Save a summary and  updated PKLs with is_doublet information for each cell.
"""

import scrublet as scr
import pickle
import pandas as pd
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from utilities.general_helpers import *
from termcolor import colored
from utilities.droplet_dataset import loading_sample

ROW_SAMPLES_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\droplet_seq\new_data_3.10.21\ROW_DATA'
SAMPLES_INFORMATION_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\apoptosis\6.10.21'
OUTPUT_PATH = fr'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\outputs\new_data_3.10.21_outputs\scrublet\6.10.21'

# Union summaries - None if you don't want to combine with older summary,
# otherwise specify the previous summary
UNION_SUMMARY_PATH = None #r'D:\Technion studies\Keren Laboratory\python_playground\outputs\apoptosis\16.12.20_empty_removed\apoptosis_summary.csv'


def extract_sample(sample_id):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    data_path = join(SAMPLES_INFORMATION_PATH, sample_id, f'{sample_id}.pkl')
    rna_sample = extract_droplet_data_from_pickle(data_path)
    print(colored(f'sample id {sample_id}', 'blue'))
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
                                       'n_cells',
                                       'n_scrub_doublets',
                                       'p_doublets'])
    samples = [subfolder for subfolder in os.listdir(SAMPLES_INFORMATION_PATH)]
    create_folder(OUTPUT_PATH)
    for sample in [s for s in samples if not 'csv' in s]:
        print(f"Working on {sample}")
        # Extracts one of the samples from PC
        rna_sample = loading_sample(row_data_path=join(ROW_SAMPLES_PATH, f'{sample}'),
                                    cells_information_path=join(SAMPLES_INFORMATION_PATH, f'{sample}'))
        rna_sample.normalize_data()
        try:
            scrub = scr.Scrublet(rna_sample.counts)
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
        except Exception as ex:
            print(colored(f"Exception has occurred in scrub.scrub_doublets: {ex}. \nTrying with 20 components", 'red'))
            scrub = scr.Scrublet(rna_sample.counts)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=18)
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

        # create_folder(join(OUTPUT_PATH, sample))
        # pickle.dump((rna_sample), open(join(OUTPUT_PATH, sample, f'{sample}.pkl'), 'wb'))
        rna_sample.save_cells_information(join(OUTPUT_PATH, f'{sample}'))

    return summary_df



def union_summaries(summary_df):
    existing_summary_df = pd.read_csv(UNION_SUMMARY_PATH)
    merged_summary = pd.merge(existing_summary_df, summary_df)
    return merged_summary


if __name__ == '__main__':
    summary_df = run_scrub_over_all_samples()
    if UNION_SUMMARY_PATH:
        summary_df = union_summaries(summary_df)
    summary_df.to_excel(join(OUTPUT_PATH, r'scrublet_summary.xlsx'), index=False)