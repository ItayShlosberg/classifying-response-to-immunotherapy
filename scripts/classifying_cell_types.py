
"""
For droplet_seq csRNAseq 2020 data. Classifying cell_types by markers table in ImmuneCellsMarkersUpdated.xlsx.

There are Positive markers so that the cells that are suspected to be classified into them should have value greater
than 1 in order to satisfy the conditions and there are negative markers so that the cells that are suspected to be
classified into them should have value 1.
"""

import sklearn
from utilities.droplet_dataset import *
from utilities import *
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib.pyplot as plt
import pickle
import random
from scipy.stats import pearsonr
from matplotlib.pyplot import figure
import pandas as pd
import os.path as path

PROJECT_PATH = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy'


def dict_append(dic, plc, val):
    dic[plc] = dic[plc] + [val]


def find_indexes_of_markers_in_sample(_sample_genes, markers):
    return [idx1 for idx1, g1 in enumerate(_sample_genes) for idx2, g2 in enumerate(markers) if g1 == g2]


def get_all_possible_combinations_of_markers(cell_type_name, celltype_markers):
    def recursive_comnibation_building(leftover, all_combinations=[[]]):
        def meiosis_of_combination(and_markers, or_markers):
            new_combinations = []
            for or_marker in or_markers:
                for and_marker in and_markers:
                    new_combinations.append(and_marker + [or_marker])
            return new_combinations

        if len(leftover) == 0:
            return all_combinations
        all_currrent_markers = leftover[0]
        del leftover[0]
        all_combination = meiosis_of_combination(all_combinations, all_currrent_markers)
        return recursive_comnibation_building(leftover, all_combination)

    initial_marker_structure = [s.split(';') for s in celltype_markers]
    combinations = recursive_comnibation_building(initial_marker_structure)
    return combinations


