"""
Going through inferCNV outputs over all samples and updating RNAseq python objects:
For each sample we have infercnv.observations.txt - for cancer cells and cells we marked as cancer due to lack of information.
and infercnv.references.txt  - for immune cells.
Each of those files is a matrix of the modified expressions - rows are genes with CNVs and columns are the cells (barcodes),
so that Vi,j is the CNV value of cell number j for gene i - a value around 1 is normal, V > ~1.1 - Insertion, V < ~0.9 Deletion.

Here, we're gonna update cells information for each of the samples with the CNV information deduced from InferCNV.
"""

from utilities.droplet_dataset import *
import numpy as np
import pickle
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join

SAMPLES_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\10.12.20'
INFERCNV_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\executions\all_data_31.12.20'
OUTPUT_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\31.12.20'


def search_and_update_cnvs_in_file(file_path, updated_rna_sample):
    """
    stromal
    :param file_path: path of infercnv.references\observation.txt file.
    :return:
    """

    # FILE is too big, therefore we're doing not ideal loop.
    # First loop only to extract the size (number of genes) of the matrix.
    number_of_genes = updated_rna_sample.number_of_genes
    matrix = None
    genes = []
    loop_itr = 0
    with open(file_path, 'r') as f:
        for line in f:
            if loop_itr == 0:
                barcodes = [ii.replace('\"', '') for ii in line.split(' ')]
                barcodes[-1] = barcodes[-1][:-1]
            else:
                genes.append(line.split(' ')[0].replace('\"', ''))
                gene_cnv_values = np.array([float(aa) for aa in line.split(' ')[1:]])
                if matrix is None:
                    matrix = np.expand_dims(gene_cnv_values, axis=1)
                else:
                    matrix = np.concatenate([matrix, np.expand_dims(gene_cnv_values, axis=1)], axis=1)
            loop_itr +=1

    cnv_gene_indexes = [updated_rna_sample.gene_names.index(ii) for ii in genes]
    for inferCNV_cell_idx, barcode in enumerate(barcodes):
        barcode_idx_in_sample = updated_rna_sample.barcodes.index(barcode)
        CNV = np.ones([number_of_genes])
        CNV[cnv_gene_indexes] = matrix[inferCNV_cell_idx, :]

        # update properties
        updated_rna_sample.cells_information[barcode_idx_in_sample].inferCNV = CNV
        updated_rna_sample.cells_information[barcode_idx_in_sample].count_insertions = sum(CNV > 1.02)
        updated_rna_sample.cells_information[barcode_idx_in_sample].count_deletions = sum(CNV < 0.98)
        if updated_rna_sample.cells_information[barcode_idx_in_sample].is_cancer and (sum(CNV < 0.98) or sum(CNV > 1.02)):
            updated_rna_sample.cells_information[barcode_idx_in_sample].is_stormal = True

def get_Next_Generation_RNAseq(pr_rna_sample):
    up_rna_sample = RNAseq_Sample(pr_rna_sample.counts,
                                  pr_rna_sample.gene_names,
                                  pr_rna_sample.barcodes,
                                  pr_rna_sample.features)

    up_cells_inf = Cell_Inf_List(pr_rna_sample.number_of_cells)
    for idx, cell_inf in enumerate(pr_rna_sample.cells_information):
        for key, val in cell_inf.__dict__.items():
            up_cells_inf[idx].__setattr__(key, val)
    up_rna_sample.cells_information = up_cells_inf
    return up_rna_sample


def update_sample(sample_id):
    # Extracts one of the samples from PC
    data_path = join(SAMPLES_PATH, sample_id, f'{sample_id}.pkl')
    infer_cnv_cancer_path = join(INFERCNV_PATH, sample_id, 'infercnv.observations.txt')
    infer_cnv_immune_path = join(INFERCNV_PATH, sample_id, 'infercnv.references.txt')

    prev_rna_sample = extract_droplet_data_from_pickle(data_path)
    updated_rna_sample = get_Next_Generation_RNAseq(prev_rna_sample)

    search_and_update_cnvs_in_file(infer_cnv_cancer_path, updated_rna_sample)
    search_and_update_cnvs_in_file(infer_cnv_immune_path, updated_rna_sample)

    create_folder(join(OUTPUT_PATH, sample_id))
    pickle.dump((updated_rna_sample), open(join(OUTPUT_PATH, sample_id, f'{sample_id}.pkl'), 'wb'))


def over_all_samples():
    samples = [subfolder for subfolder in os.listdir(SAMPLES_PATH) if subfolder.startswith('M')]
    create_folder(OUTPUT_PATH)
    for sample_id in samples[20:]:
        print(sample_id)
        update_sample(sample_id)


if __name__ == '__main__':
    over_all_samples()