"""
For droplet_seq csRNAseq 2020 data.
Label Apoptosis cells over all samples. Then save the updated rna sample at the same place (overwrite).

The label process is defined as follow:

1. Take all samples from SAMPLES_PATH. Should be one main folder contains folder with pkl for each sample.

2. Load Excel file contains row for each 'special' sample
(sample that removing apoptosis with some threshold does bad job).

3. For each 'normal' sample (according to privies section) we define
Apoptosis cells as cells with mitochondria content > 0.2.

4. For each 'special' sample we will use Apoptosis clusters using the Excel table.
    all cells from those clusters having mitochondria content > 0.15 will be determined Apoptosis.
    Additionally, the table contains a 'confidence' threshold for apoptosis cells that 'escaped' from the cluster process.

5. We are going back on with cells determined apoptosis which also classified as neutrophils
(they are exceptional and rare and we'd like to use neutrophils as much as we can)


"""

from utilities.droplet_dataset import *
import numpy as np
import pickle
import pandas as pd
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from termcolor import colored

SAMPLES_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\classifying_cell_types\17.2.21'
OUTPUT_PATH = fr'D:\Technion studies\Keren Laboratory\python_playground\outputs\apoptosis\17.2.21'

# Use 10X clusters
CLUSTER_FOLDER_PATH = r'D:\Technion studies\Keren Laboratory\Data\Melanoma\clusters'
# Use APOPTOSIS_CLUSTERS table to see which samples should be treated with clustering
APOPTOSIS_CLUSTERS_PATH = r'D:\Technion studies\Keren Laboratory\Data\tables\apoptosis_remove_using_clusters.xlsx'

# Use the classifying cell-types script for apoptosis summary
# None - if you're not interested in saving a new updated summary.
SUMMARY_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\classifying_cell_types\10.12.20\summary.csv'
MIN_CLUSTER_THRESHOLD = 0.15
NORMAL_SAMPLES_THRESHOLD = 0.2


def extract_sample(sample_id, samples_path=SAMPLES_PATH):
    """
    Extracts one of the samples from PC
    :param sample_id: id of rna sample (Mi)
    :return: rna_sample
    """
    sample_path = join(samples_path, sample_id, f'{sample_id}.pkl')
    rna_sample = extract_droplet_data_from_pickle(sample_path)
    print(f'sample id {sample_id}')
    print(f'count shape {rna_sample.counts.shape}')
    print(f'number of cells {rna_sample.number_of_cells}')
    print(f'number of genes {rna_sample.number_of_genes}')
    return rna_sample


def extract_apoptosis_cluster_indexes(sample_id, intersting_cluster, rna_sample):
    """
    1. Extract 10X cluster mapping table from PC, each cell mapped into its corresponding cluster.
    2. Align RNA sample barcodes to clustering table barcodes.
    3. Check which cells belong to interesting cluster (apoptosis-wise cluster).
    :param sample_id:
    :param intersting_cluster:
    :param rna_sample:
    :return: mapping (boolean list) each cell (in RNAseq sample cells order) belongs to interesting cluster.
    """

    # Extract 10X cluster mapping table from PC, each cell mapped into its corresponding cluster.
    cluster_table_path = join(CLUSTER_FOLDER_PATH, f'Graph_based_{sample_id}.csv')
    cluster_df = pd.read_csv(cluster_table_path)

    # Align RNA sample barcodes to clustering table barcodes.
    clusters_barcodes = cluster_df["Barcode"].tolist()
    rna_sample_alignment_idx_to_clusters = [clusters_barcodes.index(_b) for _b in rna_sample.barcodes]
    cluster_indexes = cluster_df['Graph-based'].to_numpy()[rna_sample_alignment_idx_to_clusters]

    # Convert clustering values to int.
    cluster_df['Graph-based'] = [int(s.split(' ')[1]) for s in cluster_df['Graph-based']]

    # Check which cells belong to interesting cluster (apoptosis-wise cluster).
    cluster_indexes = np.isin(cluster_indexes, intersting_cluster)
    return cluster_indexes


def determine_apoptosis_in_normal_sample():
    """
    At normal samples, said samples that automate removal (only by threshold) doesn't remove more than 10% of the cells.
    we define Apoptosis cells as cells with mitochondria content > 0.2.
    we take all samples that aren't shown in CLUSTER.xlsx file. says they are normal.
    we will save the pkl result for each of the updated samples.
    """
    print(colored('apoptosis_in_normal_samples', 'yellow'))
    apoptosis_clustes_df = pd.read_excel(pd.ExcelFile(APOPTOSIS_CLUSTERS_PATH))

    all_samples = os.listdir(SAMPLES_PATH)
    normal_samples = [sm for sm in all_samples if
                      not sm in apoptosis_clustes_df['sample_id'].tolist() and not '.csv' in sm]

    for sample_id in normal_samples:
        rna_sample = extract_sample(sample_id)

        # Extract mitochondria content.
        counting_reads = rna_sample.counts.sum(axis=1).astype(np.float64)
        mitochondria_content = rna_sample.counts[:, [s.startswith('MT-') for s in rna_sample.gene_names]].sum(
            axis=1).astype(np.float64)
        percetage_of_mitochondria = np.divide(mitochondria_content, counting_reads, out=np.zeros_like(counting_reads),
                                              where=counting_reads != 0)

        # Take all apoptosis cluster cells which their mitochondria content > NORMAL_SAMPLES_THRESHOLD
        apoptosis_indexes = percetage_of_mitochondria > NORMAL_SAMPLES_THRESHOLD

        # Update sample which cells are apoptosis. Exclude Neutrophils.
        for apop_idx in np.where(apoptosis_indexes)[0].tolist():
            if not 'Neutrophils' in rna_sample.cells_information[apop_idx].cell_type_list:
                rna_sample.cells_information[apop_idx].is_apoptosis = True

        # Save rna_sample local. if same path is used it'll overwrite the sample.
        output_folder = join(OUTPUT_PATH, sample_id)
        create_folder(output_folder)
        pickle.dump((rna_sample), open(join(output_folder, f"{sample_id}.pkl"), "wb"))


def determine_apoptosis_in_special_samples():
    """
    Use Apoptosis clusters using the Excel table.
    all cells from those clusters having mitochondria content > 0.15 will be determined Apoptosis.
    Additionally, the table contains a 'confidence' threshold for apoptosis cells that 'escaped' from the
    cluster process.
    At the end of the process the updated samples containing the apoptosis information will be saved.
    :param rna_sample:
    """
    print(colored('apoptosis_in_special_samples', 'blue'))
    apoptosis_clustes_df = pd.read_excel(pd.ExcelFile(APOPTOSIS_CLUSTERS_PATH))
    apoptosis_clustes_df['clusters'] = apoptosis_clustes_df['clusters'].astype(str)
    for row_idx, row in apoptosis_clustes_df.iterrows():
        sample_id = row['sample_id']
        confidence_threshold = row['confidence_threshold']
        clusters = [int(ii) for ii in row['clusters'].split(';')]

        rna_sample = extract_sample(sample_id)
        cluster_indexes = extract_apoptosis_cluster_indexes(sample_id, clusters, rna_sample)

        # Extract mitochondria content.
        counting_reads = rna_sample.counts.sum(axis=1).astype(np.float64)
        mitochondria_content = rna_sample.counts[:, [s.startswith('MT-') for s in rna_sample.gene_names]].sum(
            axis=1).astype(np.float64)
        percetage_of_mitochondria = np.divide(mitochondria_content, counting_reads, out=np.zeros_like(counting_reads),
                                              where=counting_reads != 0)

        # Take all apoptosis cluster cells which their mitochondria content > 0.15
        apoptosis_indexes = np.logical_and(percetage_of_mitochondria > MIN_CLUSTER_THRESHOLD, cluster_indexes)
        # Add cells with mitochondria content > confidence threshold
        apoptosis_indexes = np.logical_or(apoptosis_indexes, percetage_of_mitochondria > confidence_threshold)

        # Update sample which cells are apoptosis. Exclude Neutrophils.
        for apop_idx in np.where(apoptosis_indexes)[0].tolist():
            if not 'Neutrophils' in rna_sample.cells_information[apop_idx].cell_type_list:
                rna_sample.cells_information[apop_idx].is_apoptosis = True

        # Save rna_sample local. if same path is used it'll overwrite the sample.
        output_folder = join(OUTPUT_PATH, sample_id)
        create_folder(output_folder)
        pickle.dump((rna_sample), open(join(output_folder, f"{sample_id}.pkl"), "wb"))


def add_apoptosis_summary():
    """
    Add apoptosis summary to existing summary of classifying_cell_types.py script
    for each sample calculates the number of apoptosis cells and add it to df.
    save a new CSV in OUTPUT path.
    :return:
    """
    print(colored('apoptosis_summary', 'red'))
    existing_summary_df = pd.read_csv(SUMMARY_PATH)
    existing_summary_df['number of apoptosis cells'] = np.zeros(len(existing_summary_df)).astype(np.int)
    updated_apoptosis_df = pd.DataFrame(columns=['sample name',
                                                 'n_cells',
                                                 'p_cancer_or_immune_cells',
                                                 'p_apoptosis_cells',
                                                 'p_cells_with_classification_cancer_immune_apoptosis)',
                                                 'n_immune_cells',
                                                 'n_cancer_cells',
                                                 'n_apoptosis_cells',
                                                 'cancer_immune_conflict'])
    for row_index, row in existing_summary_df.iterrows():
        sample_id = row['sample name']
        # Extract samples from OUTPUT path because that's the path of the UPDATED samples (with the apoptosis)
        rna_sample = extract_sample(sample_id, OUTPUT_PATH)
        number_of_cells = rna_sample.number_of_cells

        number_of_apoptosis_cells = sum([c_inf.is_apoptosis for c_inf in rna_sample.cells_information])
        number_of_immune_cells = sum([c_inf.is_immune for c_inf in rna_sample.cells_information
                                          if not c_inf.is_apoptosis])
        number_of_cancer_cells = sum([c_inf.is_cancer for c_inf in rna_sample.cells_information
                                      if not c_inf.is_apoptosis])
        number_of_cancer_immune_conflict = sum([c_inf.cancer_immune_conflict for c_inf in rna_sample.cells_information
                                                if not c_inf.is_apoptosis])
        percentage_of_cells_with_clss = (number_of_cancer_cells + number_of_immune_cells)/number_of_cells
        percentage_of_apoptosis = number_of_apoptosis_cells/number_of_cells
        percentage_of_controlled_cells = percentage_of_apoptosis + percentage_of_cells_with_clss
        updated_apoptosis_df = updated_apoptosis_df.append(pd.DataFrame([[sample_id,
                                                                          number_of_cells,
                                                                          percentage_of_cells_with_clss,
                                                                          percentage_of_apoptosis,
                                                                          percentage_of_controlled_cells,
                                                                          number_of_immune_cells,
                                                                          number_of_cancer_cells,
                                                                          number_of_apoptosis_cells,
                                                                          number_of_cancer_immune_conflict]],
                                                                          columns=updated_apoptosis_df.columns))


    # for sample_id in [sm for sm in all_samples if not 'csv' in sm]:
    #     rna_sample = extract_sample(sample_id, OUTPUT_PATH)
    #     number_of_apoptosis_cells = sum([c_inf.is_apoptosis for c_inf in rna_sample.cells_information])
    #
    #     # update summary table
    #     row_index = existing_summary_df[existing_summary_df['sample name'] == sample_id].index[0]
    #     existing_summary_df.at[row_index, 'number of apoptosis cells'] = number_of_apoptosis_cells

    updated_apoptosis_df.to_csv(join(OUTPUT_PATH, 'apoptosis_summary.csv'), index= False)


if __name__ == '__main__':
    if not os.path.isdir(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)
    determine_apoptosis_in_special_samples()
    determine_apoptosis_in_normal_sample()
    if SUMMARY_PATH:
        add_apoptosis_summary()
