import numpy as np
import yaml
import os
import pandas
from collections import Counter
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import sys


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def plot_counter_list(l):
    letter_counts = Counter(l)
    df = pandas.DataFrame.from_dict(letter_counts, orient='index')
    df.plot(kind='bar')


def search_in_list(count_list, key):
    d = {v[0]: v[1] for v in count_list}
    return d.get(key, 0)


def is_there_overlap_in_lists(l1 , l2):
    return len([f for f in l1 if f in l2])!=0


class Experiments_manager:
    def __init__(self, experiment_name, experiment_folder):
        exp_folder_path = os.path.join(experiment_folder, experiment_name)
        if not os.path.isdir(exp_folder_path):
            os.mkdir(exp_folder_path)
        self.print_file_path = os.path.join(exp_folder_path, 'prints.txt')
        self.out_file = open(self.print_file_path, 'w')
        self.orig_stdout = None

    def activate_prints_to_file(self):
        self.orig_stdout = sys.stdout
        sys.stdout = self.out_file
        return self

    def finish_run(self):
        sys.stdout = self.orig_stdout
        self.out_file.close()
        with open(self.print_file_path, 'r+') as f:
            lines = f.readlines()
            for line in lines:
                print(line, end='')

    def print(self, txt, end='\n'):
        print(txt, end=end)
        self.out_file.write(txt+end)


def experiment_manager(experiment_name, experiment_folder):
    def experiment_manager_wrapper(func):
        def inner(*args, **kwargs):
            print(f"Experiment \'{experiment_name}\' has started, prints will be saved in \'{experiment_folder}\'")
            em = Experiments_manager(experiment_name, experiment_folder).activate_prints_to_file()
            output = func(*args, **kwargs)
            em.finish_run()
            return output
        return inner
    return experiment_manager_wrapper


def load_yml(yml_path):
    with open(yml_path, 'r') as f:
        config = yaml.safe_load(f)
    return config['EXPERIMENT']['experiment_name'], config['EXPERIMENT']['experiments_folder'], config


def visualization_confusion_matrix(labels, predictions):
    cm = confusion_matrix(labels, predictions)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=['non-response', 'response'])
    disp.plot(include_values=True,
              cmap='viridis', ax=None, xticks_rotation='horizontal',
              values_format=None)
    plt.show()
    return cm