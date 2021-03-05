"""

Summarizing the Quality control. After drawing conclusion with inferCNV and classifying tumor cells and
immune cells using the CNV maps run this script to update the data and save a version of PKL of the
data with the conclusion. The QC includes:

1. classifying cell-types using immune & tumor markers.
2. Finding dying cells (in apoptosis).
3. finding doublets with Scrublet.
4. Using CNV maps of sample to Cluster similar cells together, that way we finished classifying all cells and remove
garbage.

Wiped out (excluded): CellBender.

-------------------------------------

the ouput of the script is Excel file containing :
1. The initial number of droplets.
2. The number of cells after removing dying cells or cells with cancer & immune conflicts (part of it was concluded during the CNV process).
3. Number & percentage of immune, tumor, stromal cells. Each cell has exactly one tag.


"""

import pandas as pd
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from utilities.general_helpers import *
from termcolor import colored
from utilities.droplet_dataset import loading_sample

# the path of the samples we want to summarize:
ROW_SAMPLES_PATH = fr'D:\Technion studies\Keren Laboratory\Data\droplet_seq\ROW_DATA'
SAMPLES_INFORMATION_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\4.3.21'

# where wh excel output will be saved:
OUTPUT_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\QC\QC_summary_5.3.21'


# def extract_sample(sample_id):
#     """
#     Extracts one of the samples from PC
#     :param sample_id: id of rna sample (Mi)
#     :return: rna_sample
#     """
#     data_path = join(ROW_SAMPLES_PATH, sample_id)
#     rna_sample = extract_droplet_data_from_pickle(data_path)
#     print(colored(f'sample id {sample_id}', 'blue'))
#     print(f'count shape {rna_sample.counts.shape}')
#     print(f'number of cells {rna_sample.number_of_cells}')
#     print(f'number of genes {rna_sample.number_of_genes}')
#     return rna_sample


if __name__ == '__main__':
    samples = [subfolder.replace(".pkl", "") for subfolder in os.listdir(ROW_SAMPLES_PATH) if (not 'csv' in subfolder and not 'xlsx' in subfolder)]
    create_folder(OUTPUT_PATH)

    summary_df = pd.DataFrame(columns=['sample',
                                       'n_droplets',
                                       'n_cells',
                                       'p_cells',
                                       'n_immune',
                                       'p_immune',
                                       'n_tumor',
                                       'p_tumor',
                                       'n_stromal',
                                       'p_stromal'])
    for sample_id in samples:
        print(f'Working on {sample_id}')
        # Extracts one of the samples from PC
        rna_sample = loading_sample(row_data_path=join(ROW_SAMPLES_PATH, f'{sample_id}.pkl'),
                                    cells_information_path=join(SAMPLES_INFORMATION_PATH, f'{sample_id}.pkl'))
        n_droplets = rna_sample.number_of_cells

        # remove dying cells and garbage.
        rna_sample = rna_sample.filter_cells_by_property('should_be_removed', False)
        n_cells = rna_sample.number_of_cells
        n_cancer = sum(rna_sample.cells_information.getattr('is_cancer'))
        n_immune = sum(rna_sample.cells_information.getattr('is_immune'))
        n_stromal = sum(rna_sample.cells_information.getattr('is_stromal'))

        # Append to DF
        summary_df = summary_df.append(pd.DataFrame([[sample_id,
                                                      n_droplets,
                                                      n_cells,
                                                      n_cells / n_droplets,
                                                      n_immune,
                                                      n_immune / n_cells,
                                                      n_cancer,
                                                      n_cancer / n_cells,
                                                      n_stromal,
                                                      n_stromal / n_cells]],
                                                    columns=summary_df.columns))

        summary_df.to_excel(join(OUTPUT_PATH, r'QC_summary.xlsx'), index=False)