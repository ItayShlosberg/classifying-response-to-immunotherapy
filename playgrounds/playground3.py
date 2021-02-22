
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
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle


SAMPLE = 'M133'
SAMPLE_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\scrublet\10.12.20'

sample = extract_droplet_data_from_pickle(join(SAMPLE_PATH, SAMPLE, f'{SAMPLE}.pkl'))

ss = sample[[10, 20, 22]]

_breakpoint = 0