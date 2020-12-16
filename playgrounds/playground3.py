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



class A:
    def __init__(self):
        self.x1 = False
        self.x2 = False
        self.x3 = False
        self.x4 = False


a = A()
print(a.__dict__)
a.__setattr__('x1', True)
print(a.__dict__)