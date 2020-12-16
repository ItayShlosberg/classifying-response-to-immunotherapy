from utilities.droplet_dataset import *
from utilities.general_helpers import *
from os.path import join
from matplotlib import pyplot
import numpy as np
import scipy
import pickle
import matplotlib
import scrublet as scr
import pickle
from DL.data_creation import *
from DL.data_conversions.txt_to_python_structures import *
from os.path import join
from utilities.droplet_dataset import *
from utilities import *
import numpy as np
import pickle
import pickle
import pandas as pd
from DL.data_loading import extract_droplet_data_from_pickle
from os.path import join
from utilities.general_helpers import *

SOURCE_PATH = r'D:\PycharmProjects\CellBender\conversions'
OUTPUT_PATH = r'D:\PycharmProjects\CellBender\13.12.20'




samples = [subfolder for subfolder in os.listdir(SOURCE_PATH) if not os.path.isfile(join(SOURCE_PATH,subfolder))]


create_folder(OUTPUT_PATH)
for sample_id in samples:
    path = join(OUTPUT_PATH, sample_id)
    print(path)
    create_folder(path)