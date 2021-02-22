import pickle
import numpy as np
import os
from utilities.general_helpers import *
from utilities.droplet_dataset import *

"""
After extracting the data that was collected from MGH (origin Matlab structures) and organizing it in folders.
This script goes through all the folders and builds python structures that will hold them.
"""


PROTEIN_CODING_FILE = r'..\..\Data\gene_ens_map.xlsx'
ROOT_PATH = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples'
PRODUCT_PATH = r'..\..\Data\rna_seq200k'


def str_list_to_float(li):
    return [int(y) for y in li]


def keep_protein_coding_only():
    df = pandas.read_excel(PROTEIN_CODING_FILE)
    kept_genes_names = [gene[0] for gene in df[df.lincRNA == 'protein_coding'][['MIR1302-11']].values]
    # indices_of_protein_coding = [i for i in range(len(gene_names)) if gene_names[i] in kept_genes_names]
    # cells = cd45_cells[:, indices_of_protein_coding]
    # gene_names = operator.itemgetter(*indices_of_protein_coding)(gene_names)
    # return cells, gene_names
    return kept_genes_names


def extract_all_gene_names():
    folders = [(os.path.join(ROOT_PATH, subfolder, 'GeneName.txt'), subfolder) for subfolder in os.listdir(ROOT_PATH)]
    all_genes = []
    for idx, (genes_path, folder) in enumerate(folders):
        # print(f'number {idx+1} folder {folder}')
        with open(genes_path, 'r') as read_file:
            gene_list = read_file.readlines()
            gene_list = [g[:-2] for g in gene_list]
        all_genes.append(gene_list)
    return all_genes


def convert_txt_to_python_structure():
    ens_id_suffix = 'ensID.txt'
    counts_suffix = 'counts.txt'
    gene_name_suffix = 'GeneName.txt'
    sample_name_suffix = 'sample_name.txt'
    folders = [(os.path.join(ROOT_PATH, subfolder), subfolder) for subfolder in os.listdir(ROOT_PATH)]
   # gene_lists = extract_all_gene_names()
   # gene_list = Cohort_RNAseq.uniform_gene_list(gene_lists)
    cohort_RNA_samples = []#Cohort_RNAseq(gene_list)
    for idx, (folder_path, folder) in enumerate(folders):
        if folder != 'M139':
            continue
        print(f'number {idx+1} folder {folder}')
        # geneName
        gene_name_path = os.path.join(folder_path, gene_name_suffix)
        ens_id_path = os.path.join(folder_path, ens_id_suffix)
        counts_path = os.path.join(folder_path, counts_suffix)
        sample_name_path = os.path.join(folder_path, sample_name_suffix)

        with open(gene_name_path, 'r') as read_file:
            gene_lists = read_file.readlines()
            gene_lists = [g[:-2] for g in gene_lists]

        with open(ens_id_path, 'r') as read_file:
            ens_id_list = read_file.readlines()
            ens_id_list = [ens[:-2] for ens in ens_id_list]

        with open(sample_name_path, 'r') as read_file:
            sample_list = read_file.readlines()
            sample_list = [sample[:-2] for sample in sample_list]

        with open(counts_path, 'r') as read_file:
            counts = read_file.readlines()
            counts = np.array([str_list_to_float(counts.split(' ')[:-1]) for counts in counts]).astype(np.uint16).T

        RNA_sample = RNAseq_Sample(counts, gene_lists, sample_list, ens_id_list)
        # cohort_RNA_samples.add_RNAseq(RNA_sample)
        pickle.dump((RNA_sample), open(os.path.join(folder_path, "RNA_sample.pkl"), "wb"))
        # cohort_RNA_samples.append(RNA_sample)
        #cohort_RNA_samples.append(ens_id_list)
    _breakpoint = 0
    # pickle.dump((cohort_RNA_samples), open(os.path.join(PRODUCT_PATH, "dropletRNAseq_dataset.pkl"), "wb"))

    # onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    print(folders)
    return cohort_RNA_samples


if __name__ == '__main__':
    b = convert_txt_to_python_structure()
    print(b)

