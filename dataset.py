import numpy as np
from sklearn.metrics import auc
from dataset import *
import numpy as np
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pickle
import random
import pandas as pd
from Bio.Cluster import kcluster
import operator
import seaborn as sns
import time
import sys
from data import filter_genes_by_variance


class RNAseq_Dataset:
    def __init__(self, cells, patient_structure, gene_names):
        self.cells = cells
        self.patients = Patients(patient_structure)
        self.gene_names = gene_names
        self.length = len(self.patients)

    def __len__(self):
        return len(self.patients)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.cells[item], self.patients[item]
        if isinstance(item, slice):
            return RNAseq_Dataset(self.cells[item], self.patients[item], self.gene_names)
        if isinstance(item, list):
            if sum([(ii==0 or ii==1) for ii in item]) == len(item):
                item = [i for i in range(len(self)) if item[i]]
            return RNAseq_Dataset(self.cells[item, :], self.patients[item], self.gene_names)

    def expression_of_genes(self):
        """
        :param cells: PKL format
        :param gene_names: list of gene names in the order shown in the cells.
        :return: sorted avg gene expressions and gene names sorted by avg expression.
        """
        average_val_genes = np.mean(self.cells, axis=0)
        indexes_of_sorted_expression_of_genes = np.argsort(average_val_genes)
        gene_sorted_by_expression = operator.itemgetter(*indexes_of_sorted_expression_of_genes)(self.gene_names)
        sorted_average_expression_genes = average_val_genes[indexes_of_sorted_expression_of_genes]
        return sorted_average_expression_genes, gene_sorted_by_expression

    def get_high_expressed_genes(self):
        """
        Return (only the indexes) of the highest expression genes. TODO: make the genes amount changeable.
        :param cells: pkl format
        :param gene_names: when given, the function return genes names with highest expression also.
        :return: indexes of highest expressed gene and sub-dataset contains cells with the highest expressed genes.
        """
        activated_genes = sum(np.sum(self.cells, axis=0) != 0)
        average_val_genes = np.mean(self.cells, axis=0)
        expression_histogram = np.histogram(average_val_genes)
        expression_histogram_df = np.concatenate((np.expand_dims(expression_histogram[0][:10], axis=0),
                                                  np.expand_dims(expression_histogram[1][:10], axis=0)))
        amount_high_expressed_genes = min(activated_genes * 0.05, sum(expression_histogram[0][-2:]))
        high_expressed_percentage = amount_high_expressed_genes / activated_genes
        indices_of_high_expressed_genes = np.argsort(average_val_genes)[-3:]
        print(f'High expressed genes percent {high_expressed_percentage}')
        return indices_of_high_expressed_genes, \
                RNAseq_Dataset(self.cells[:, indices_of_high_expressed_genes],
                               self.patients,
                                operator.itemgetter(*indices_of_high_expressed_genes)(self.gene_names))

    def filter_cells_by_supervised_classification(self,  required_cell_type="T cells"):
        """
        patients_information contains 'supervised classification' field, which is the classification of the cell
        defined by gene expression in 80% of the cells, and the remaining 20% were done by manual process.
        the function filters cells by cell-type was done by this classification.
        :param required_cell_type: name of the desired cell-type.
        :return: reduced RNAseq dataset contains cells of the desired cell-type
        """
        indexes_list = self.patients.get_cells_belong_to_cells_type(required_cell_type)
        return RNAseq_Dataset(self.cells[indexes_list, :], self.patients[indexes_list], self.gene_names)


class Patients:
    def __init__(self, patient_structure=None):
        self.patients_list = []
        if patient_structure:
            if isinstance(patient_structure[0], dict):
                if patient_structure:
                    for p_dict in patient_structure:
                        self.patients_list.append(Patient_information_cell(p_dict['sample index'],
                                                                           p_dict['cell id'],
                                                                           p_dict['patient details'],
                                                                           p_dict['response'],
                                                                           p_dict['treatment'],
                                                                           p_dict.get('response label', None),
                                                                           p_dict.get('keren cluster', None),
                                                                           p_dict.get('supervised classification', None)))
            elif isinstance(patient_structure[1], Patient_information_cell):
                self.patients_list = patient_structure
        self.length = len(self.patients_list)

    def __len__(self):
        return len(self.patients_list)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.patients_list[item]
        if isinstance(item, slice):
            # Get the start, stop, and step from the slice
            return Patients([self.patients_list[ii] for ii in range(*item.indices(len(self.patients_list)))])
        if isinstance(item, str):
            return [p.__getattribute__(item) for p in self.patients_list]
        if isinstance(item, list):
            if sum([(ii==0 or ii==1) for ii in item]) == len(item):
                item = [i for i in range(len(self)) if item[i]]
            return Patients([self.patients_list[i] for i in item])

    def __getattr__(self, item):
        return [p.__getattribute__(item) for p in self.patients_list]

    def get_cells_belong_to_cells_type(self, cell_type):
        return [p.belongs_to_cell_type(cell_type) for p in self.patients_list]

    def slice_by_binary_indexes(self, binary_indexes):
        self.patients_list = [self.patients_list[i] for i in range(len(self.patients_list)) if binary_indexes[i]]

    def get_patients_names(self):
        return list(set([p.patient_details for p in self.patients_list]))


class Patient_information_cell:
    def __init__(self, sample_index, cell_id, patient_details, response, treatment, response_label, keren_cluster,
                 supervised_classification):
        self.sample_index = sample_index
        self.cell_id = cell_id
        self.patient_details = patient_details
        self.response = response
        self.treatment = treatment
        self.response_label = response_label
        self.keren_cluster = keren_cluster
        self.supervised = supervised_classification

    def belongs_to_cell_type(self, cell_type):
        return cell_type in self.supervised