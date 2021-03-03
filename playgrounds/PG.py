import numpy as np
import matplotlib
from utilities.droplet_dataset import  build_cohort
import os
from os.path import join
from DL.Mars_seq_DL.data_loading import *
from utilities.droplet_dataset import *


SAMPLES = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\update_runs\21.2.21'

OUTPUT = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\Cohort\cohort_all_samples_3.2.21.pkl'
# ss = r'C:\Users\itay\Desktop\New folder'

samples = [subfolder for subfolder in os.listdir(SAMPLES) if not 'csv' in subfolder]
rna_sample = extract_droplet_data_from_pickle(join(SAMPLES, samples[0]))

_breakpoint = 0