import os
from os.path import join
import sklearn
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



def check_cross_cell_types(sample_id):
    sample_path = join(SAMPLES_PATH, sample_id, f'{sample_id}.pkl')
    types = ['None', 'Immune_general', 'T cells', 'CD4 helper T cells', 'CD8 Cytotoxic T cells', 'Regulatory T cells',
             'Regulatory CD4 T cells', 'Regulatory CD8 T cells', 'Regulatory CD4_CD8 T cells', 'NKT cells', 'NK cells',
             'B cells', 'Activated T cells', 'Senescence T cells', 'Terminal effector', 'Exhausted T cells',
             'Stem_like T cells', 'Memory T cells', 'Memory CD4 T cells', 'Memory CD8 T cells',
             'Memory CD4_CD8 T cells', 'Macrophage_immature', 'Macrophage_mature', 'Monocyte_immature',
             'Monocyte_mature', 'cDCs_dendritic_cells', 'pDCs', 'myeloid cells_general_immature',
             'myeloid cells_general_mature', 'Neutrophils', 'Granolocytes']
    table = {k1: {k2: 0 for k2 in types} for k1 in types}
    data = pickle.load(open(sample_path, 'rb'))
    cls = [list(set(v)) for v in data.cells_information.getattr('cell_type_list')]
    cls = [v if len(v) else ['None'] for v in cls]

    for cell in cls:
        for t1 in cell:
            for t2 in cell:
                table[t1][t2] += 1
    df = pd.DataFrame(table.values(), index=table.keys())[table.keys()]
    return df