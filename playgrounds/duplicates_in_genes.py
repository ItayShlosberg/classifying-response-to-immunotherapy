# # import math
# # import numpy as np
# # import Bio
# # from utilities.general_helpers import *
# # from DL.data_loading import *
# # from utilities.dataset import *
# from matplotlib import pyplot
# import numpy as np
# import scipy
# import pickle
# import matplotlib
# import scrublet as scr
# from DL.data_loading import extract_data_from_pickle as ex
# import pickle
# from DL.data_creation import *
# from DL.matlab_env.txt_to_python_structures import *
#
# """
# use: pycco file.py
# for documentations
# """
#
#
# # from general_helpers import *
# #
# #
# # # em = Experiments_manager("exp1", "Data")
# # # em.activate_prints_to_file()
# # #     # (r'Data\out.txt')
# # # for i in range(10):
# # #         print(f"print number {i}")
# # #
# # # em.finish_run()
# # def ff():
# #     return 'exp1'
# # E = ff()
# # D = 'Data'
# # @experiment_manager(E, D)
# # def main():
# #     for i in range(25):
# #         print(f"11 decorator print number {i}")
# #
# # # em.finish_run()
# #
# #
# # # print("print number 1")
# #
# # if __name__ == '__main__':
# #     main()
# #     E = 'exp2'
# #     D = 'Data'
# #     main()
#
# # import scipy.io
# # mat_structure_path = r'..\Data\rnaseq2\M97.mat'
# # # mat = scipy.io.loadmat(path)
# # import matlab.engine
# # eng = matlab.engine.start_matlab()
# # Obj = eng.load(mat_structure_path)
# # l = Obj['lib']
# # Names = eng.getfield(Obj, 'name')
# # Strs = eng.getfield(Obj, 'structures')
# # StrLists = eng.getfield(Obj, 'structlist')
# # _breakpoint = 0
#
# # file_path = r'..\Data\rnaseq2\outputs\old\file2 - Copy.txt'
# # genes_path = r'..\Data\rna_seq200k\outputs\old\gene1.txt'
#
# #
# # with open(genes_path, 'r') as f:
# #     s = f.readlines()
# # CONFIG_PATH = r'..\cfg\dummy.yaml'
# # EXPERIMENT_NAME, EXPERIMENTS_FOLDER, config = load_yml(CONFIG_PATH)
# # dataset_config = config['DATASET']
# #
# # data_path = r"..\\"+dataset_config['data_path']
# # # split_data_path = dataset_config['split_data_path']
# # # save_division_path = dataset_config['save_division_path']
# # # test_percent = dataset_config['test_percent']
# # # patients_type = dataset_config['patients_type']
# # # variance = dataset_config['variance']
# # cells, gene_names, patients_information = extract_data_from_pickle(data_path)
# # whole_rna_seq_dataset = RNAseq_Dataset(cells, patients_information, gene_names)
# # _breakpoint = 0
# # _breakpoint = 0
#
#
# # from memory_profiler import profile
# #
# # @profile
# # def allocate():
# #     a = np.zeros((3,4))
# #     b = np.random.rand(20000, 40000)
# #     _breakpoint = 0
# #
# #
# # # allocate()
# # a = np.array([65536]).astype(np.uint16)
# # print(a.dtype)
# # print(a)
#
# sample = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\Data\rna_seq200k\all_samples\M114\enrich_RNA_sample.pkl'
# data_path = r'..\DATA\RNAseq_DATA.p'
#
#
# cells = pickle.load(open(sample, 'rb')).cells.T
#
#
#
#
# hist, bin_edges = np.histogram(cells.sum(axis=1), density=True)
#
#
#
# hist
# hist.sum()
# np.sum(hist * np.diff(bin_edges))
#
#
#
#
# import matplotlib.pyplot as plt
#
# rng = np.random.RandomState(10)  # deterministic random data
#
# a = np.hstack((rng.normal(size=1000),
#
#                rng.normal(loc=5, scale=2, size=1000)))
#
# # _ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
# _ = plt.hist(cells.sum(axis=1), bins='auto')  # arguments are passed to np.histogram
#
# plt.title("Histogram with 'auto' bins")
#
# plt.show()
#
#
# breakpoint = 0
import pandas
from collections import Counter
import os
PROTEIN_CODING_FILE = r'..\Data\gene_ens_map.xlsx'
ROOT_PATH = r'..\Data\rna_seq200k\all_samples'
import numpy as np



def str_list_to_float(li):
    return [int(y) for y in li]

def convert_txt_to_python_structure():

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

convert_txt_to_python_structure()