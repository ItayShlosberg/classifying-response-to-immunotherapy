from os.path import join
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score, accuracy_score, recall_score, precision_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib
from sklearn import metrics
from utilities.general_helpers import *
import pickle
from MachineLearning.ML_utilities.dataloder import *
import numpy as np
import pandas as pd


class Metrics:
    """

    """
    def __init__(self, y_true, y_pred):
        self.confusion_matrix = confusion_matrix(y_true, y_pred)
        self.TN, self.FP, self.FN, self.TP = self.confusion_matrix.ravel()

        self.specificity = self.TN / (self.TN + self.FP)
        self.accuracy = accuracy_score(y_true, y_pred)
        self.sensitivity = recall_score(y_true, y_pred)
        self.precision = precision_score(y_true, y_pred)

    def print_scores(self, title=""):
        print('------------------------')
        print(f'{title} scores:')
        print(f'Accuracy: {self.accuracy}')
        print(f'sensitivity: {self.sensitivity}')
        print(f'specificity: {self.specificity}')
        print(f'precision: {self.precision}')
        print(f'Confusion matrix:\n{self.confusion_matrix}')
        print('------------------------')

    def visualization_confusion_matrix(self, title=None, save_path=None, display_labels=None):

        if display_labels:
            disp = ConfusionMatrixDisplay(confusion_matrix=self.confusion_matrix,
                                          display_labels=display_labels)
        else:
            disp = ConfusionMatrixDisplay(confusion_matrix=self.confusion_matrix,
                                          display_labels=['non-response', 'response'])
        _, ax = plt.subplots()
        if title:
            ax.set_title(title)
        disp.plot(include_values=True,
                  cmap='viridis', ax=ax, xticks_rotation='horizontal',
                  values_format=None)
        if save_path:
            plt.savefig(os.path.join(save_path + ".png"))
        plt.show()

        return self.confusion_matrix