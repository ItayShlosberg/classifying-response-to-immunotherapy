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


def filter_by_indexes(indices, cells, patients_information=None):
    """
    return reduced cells associated to the corresponding indices.
    :param indices: indexes list of desirable cells.
    :param cells: cells from pkl file.
    :param patients_information: pkl file format
    :return: reduced cells. and their corresponding patients_information.
    """
    if patients_information:
        return cells[indices, :], [patients_information[i] for i in indices]
    return cells[indices, :]


def filter_by_binary_indices(indices, cells, patients_information=None):
    """
    return reduced cells associated to the corresponding indices.
    :param indices: binary list.
    :param cells: cells from pkl file.
    :param patients_information: pkl file format
    :return: reduced cells. and their corresponding patients_information.
    """
    indexes = [i for i in range(len(patients_information)) if indices[i]]
    if patients_information:
        return cells[indexes, :], [patients_information[i] for i in range(len(patients_information)) if indices[i]]
    return cells[indices, :]


def expression_of_genes(cells, gene_names = None):
    """
    :param cells: PKL format
    :param gene_names: list of gene names in the order shown in the cells.
    :return: sorted avg gene expressions and gene names sorted by avg expression.
    """
    average_val_genes = np.mean(cells, axis=0)
    indexes_of_sorted_expression_of_genes = np.argsort(average_val_genes)
    gene_sorted_by_expression = operator.itemgetter(*indexes_of_sorted_expression_of_genes)(gene_names)
    sorted_average_expression_genes = average_val_genes[indexes_of_sorted_expression_of_genes]
    return sorted_average_expression_genes, gene_sorted_by_expression


def filter_cells_by_supervised_classification(cells, patients, required_cell_type="T cells"):
    """
    patients_information contains 'supervised classification' field, which is the classification of the cell
    defined by gene expression in 80% of the cells, and the remaining 20% were done by manual process.
    the function filters cells by cell-type was done by this classification.
    :param cells: pkl format.
    :param patients_information: pkl format.
    :param required_cell_type: name of the desired cell-type.
    :return: cells of the desired cell-type
    """
    indexes_list = patients.get_cells_belong_to_cells_type(required_cell_type)
    cells, patients_information = filter_by_binary_indices(indexes_list, cells, patients_information)
    return cells, patients_information


def get_high_expressed_genes(self):
    """
    Return (only the indexes) of the highest expression genes. TODO: make the genes amount changeable.
    :param cells: pkl format
    :param gene_names: when given, the function return genes names with highest expression also.
    :return: indexes only. if gene names given, the function return genes names with highest expression also.
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
    if self.gene_names:
        return indices_of_high_expressed_genes, operator.itemgetter(*indices_of_high_expressed_genes)(self.gene_names)
    return indices_of_high_expressed_genes
