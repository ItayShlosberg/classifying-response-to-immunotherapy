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
from playsound import playsound
import time
import random

baby_shark_path = r'C:\Users\itay\Downloads\baby_shark.mp3'
delay = 2 * 60
count = 1

while(True):
    print(count)
    waiting_factor = random.uniform(1, 5) * 60
    print(f'waiting_factor: {waiting_factor/60}')
    time.sleep(waiting_factor)

    playsound(baby_shark_path)

    time.sleep(delay)
    count += 1





