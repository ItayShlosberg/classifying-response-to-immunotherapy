from utilities.general_helpers import *
import numpy as np


class Cohort_RNAseq:
    def __init__(self, gene_list):
        """
        :param sample_list:
        :param overlapping_genes:
        False - the gene list contains all genes from all cells.
        and if there some cell lacks genes of other cells we complete it to zero.
        True - the gene list contains only genes appearing in all cells.
        """
        self.gene_list = gene_list
        self.sample_list = []
        self.number_of_samples = 0

    def add_RNAseq(self, sample, makeover_required=True):
        """
        A choice - keep only overlapping genes OR complete genes if need, then arrange them in a uniform order.
        :param in_place: True - if you want to save it in self. False - get a new object only.
        :param fixed_length_mode:
        False - the gene list contains all genes from all cells.
        and if there some cell lacks genes of other cells we complete it to zero.
        True - the gene list contains only genes appearing in all cells.
        :return: The new cohort object.
        """
        if makeover_required:
            sample.makeover_genes(self.gene_list)
        self.sample_list.append(sample)
        self.number_of_samples += 1

    @staticmethod
    def uniform_gene_list(gene_lists):
        gene_list = sorted(list(set(flatten_list(gene_lists))))
        return gene_list


class RNAseq_Sample:

    def __init__(self, counts, gene_names, samples_list, ens_id_list):
        """
        :param counts: 2D-numpy array with rows as genes and columns as cells. The created object'll save the count
        table in a manner which rows are the cells and columns are the genes, for convenience.
        :param gene_names: ordered according to the cell order.
        :param ens_id_list: not relevant, used in alignment preprocess.
        :param samples_list: each cell has unique id.
        """
        if counts is not None:
            if counts.shape[0] != len(samples_list) or counts.shape[1] != len(gene_names):
                raise Exception("unsuitable dimensions")
            if ens_id_list and counts.shape[1] != len(ens_id_list):
                raise Exception("unsuitable dimensions")
        self.counts = counts
        self.gene_names = gene_names
        self.samples_list = samples_list
        self.ens_id_list = ens_id_list
        self.cohort_adjustment = False
        self.number_of_genes = len(gene_names)
        self.number_of_cells = len(samples_list)
        self.cells_information = [Cell_information() for _ in range(self.number_of_cells)]

    def makeover_genes(self, ens_id_list):
        all_ens_indexes = self.map_indexes(ens_id_list)
        fixed_counts = np.zeros((self.number_of_cells, len(ens_id_list)))
        fixed_counts[:, all_ens_indexes] = self.counts
        self.ens_id_list = ens_id_list
        self.number_of_genes = len(all_ens_indexes)
        self.counts = fixed_counts

    def map_indexes(self, all_genes):
        current_genes = self.gene_names
        all_genes_indexes = []
        for gene in self.gene_names:
            idx = binary_search(all_genes, gene)
            if idx == -1:
                raise ValueError
            all_genes_indexes.append(idx)

        return all_genes_indexes

    def __getitem__(self, item):
        # if isinstance(item, int):
        #     return self.counts[item], self.cells_information_list[item]
        # if isinstance(item, slice):
        #     return RNAseq_Dataset(self.cells[item], self.cells_information_list[item], self.gene_names)
        # if isinstance(item, list):
        #     # identify if we are dealing with binary indexes or explicit indexes.
        #     if sum([(ii == 0 or ii == 1) for ii in item]) == len(item):
        #         # converts to explicit indexes.s
        #         item = [i for i in range(len(self)) if item[i]]
        #     return RNAseq_Dataset(self.cells[item, :], self.cells_information_list[item], self.gene_names)
        pass


class Cell_information:
    def __init__(self):
        self.is_apoptosis = False     # True after discovery the cell is in apoptosis stage.
        self.cell_type_list = []
        self.is_classified = False
        self.is_cancer = False