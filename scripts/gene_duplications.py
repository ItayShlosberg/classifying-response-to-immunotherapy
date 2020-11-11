import pandas
from collections import Counter
import os
import numpy as np

ROOT_PATH = r'..\Data\rna_seq200k\all_samples'

"""
Finds in each samples gene names which appear multiple times.
"""


def str_list_to_float(li):
    return [int(y) for y in li]


def finds_duplicated_genes():

    counts_suffix = 'counts.txt'
    gene_name_suffix = 'GeneName.txt'
    folders = [(os.path.join(ROOT_PATH, subfolder), subfolder) for subfolder in os.listdir(ROOT_PATH)]
    df = pandas.DataFrame(columns=['File', 'gene', 'number of variations', 'number of cells having a different value in these variations'])
    for idx, (folder_path, folder) in enumerate(folders):

        print(f'number {idx+1} folder {folder}')
        # geneName
        if 'M143'==folder:
            continue
        gene_name_path = os.path.join(folder_path, gene_name_suffix)
        counts_path = os.path.join(folder_path, counts_suffix)


        with open(gene_name_path, 'r') as read_file:
            gene_lists = read_file.readlines()
            gene_lists = [g[:-2] for g in gene_lists]


        with open(counts_path, 'r') as read_file:
            counts = read_file.readlines()
            counts = np.array([str_list_to_float(counts.split(' ')[:-1]) for counts in counts]).astype(np.uint16).T


        gene_names = gene_lists
        duplicates = [(k, v) for k,v in dict(Counter(gene_names)).items() if v>1]
        for gene, n_appearances in duplicates:
            dup_indexes = [(idx,v) for idx, v in enumerate(gene_names) if v==gene]
            difference = sum(counts[:, dup_indexes[0][0]] != counts[:, dup_indexes[1][0]])
            _breakpoint = 0
            df = df.append(pandas.DataFrame([[folder, gene, n_appearances, difference]], columns=list(df.columns)))
        _breakpoint = 0
    df.reset_index(drop=True).to_csv(r'C:\Users\itay\Desktop\out.csv')


if __name__ == '__main__':
    finds_duplicated_genes()