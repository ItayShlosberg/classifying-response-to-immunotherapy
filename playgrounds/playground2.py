import pandas as pd
import pandas
from collections import Counter
from pycco import process
import os
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import os
from sklearn.metrics import confusion_matrix, plot_confusion_matrix, ConfusionMatrixDisplay



path = r'C:\Users\itay\Documents\R\win-library\4.0\infercnv\extdata\oligodendroglioma_expression_downsampled.counts.matrix'
path2 = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\M136\matrix.matrix'

path3 = r'C:\Users\itay\Documents\R\win-library\4.0\infercnv\extdata\oligodendroglioma_annotations_downsampled.txt'
path4 = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\inferCNV\M136\annotation.txt'

with open(path, 'rb') as f:
    m1 = f.readlines()

with open(path2, 'rb') as f:
    m2 = f.readlines()


with open(path3, 'rb') as f:
    a1 = f.readlines()

with open(path4, 'rb') as f:
    a2 = f.readlines()

_breakpoint = 0






