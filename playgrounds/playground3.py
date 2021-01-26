
import time
import os
from os.path import join
from utilities.droplet_dataset import *
from utilities import *
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import random
from scipy.stats import pearsonr
from matplotlib.pyplot import figure
import numpy as np
import pandas as pd
from os.path import join


print("PlayGround")
path = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\executions\all_data_31.12.20\M110'
# path = r'C:\Users\itay\Desktop\out'
file = r'infercnv.observations.txt'
count = 0
genes_vals =[]
gene_names = []
barcodes = None
s_time = time.time()
with open(join(path, file), 'r') as f:
    for line in f:
        if not count % 1000:
            current_time = time.time() - s_time
            print(f'line num: {count + 1} current time {round(current_time, 2)} sec')
        if count == 0:
            barcodes_length = len(line[:-1].split(' '))
            print(f'number of barcodes {barcodes_length}')
            barcodes = [ii.replace('\"', '') for ii in line.split(' ')]
        else:
            gene = [aa for aa in line[:-1].split(' ')]
            gene_names.append(gene[0].replace('\"', ''))
            genes_vals.append([float(val) for val in gene[1:]])
        if count > 100:
            break
        count +=1

print(f'number of genes {count}')

genes_details = list()
genome_path = r'D:\Technion studies\Keren Laboratory\Data\inferCNV_data\gencode_v19_gene_pos.txt'
with open(genome_path, 'r') as f:
    for line in f:
        genes_details.append(line[:-1].split('\t'))




genome_genes = [gg[0] for gg in genes_details]
gene_names_chr = [(gg, int(genes_details[genome_genes.index(gg)][1].replace('chr', ''))) for gg in gene_names]
chromosomes = list(set([gg[1] for gg in gene_names_chr]))


def retrieve_genes_indices_by_chro(chro_num):
    return [idx for idx, gg in enumerate(gene_names_chr) if gg[1]==chro_num]

_breakpoint = 0
