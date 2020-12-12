"""
Python structures for droplet data of 2020.
"""
from utilities.general_helpers import *
import numpy as np
import numpy as np


class Alignment_Cohort_RNAseq:
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


class Cohort_RNAseq:
    def __init__(self, gene_list):
        self.sample_list = []
        self.number_of_samples = 0

    def add_RNAseq(self, sample):
        """
        Add new RNAseq to Cohort.
        :param sample: RNAseq object
        :return:
        """
        self.sample_list.append(sample)
        self.number_of_samples += 1


class RNAseq_Sample:

    def __init__(self, counts, gene_names, barcodes, features):
        """
        :param counts: 2D-numpy array with rows as genes and columns as cells. The created object'll save the count
        table in a manner which rows are the cells and columns are the genes, for convenience.
        :param gene_names: ordered according to the cell order.
        :param ens_id_list: not relevant, used in alignment preprocess.
        :param samples_list: each cell has unique id.
        """
        if counts is not None:
            if counts.shape[0] != len(barcodes) or counts.shape[1] != len(gene_names):
                raise Exception("unsuitable dimensions")
            if features and counts.shape[1] != len(features):
                raise Exception("unsuitable dimensions")
        self.counts = counts
        self.gene_names = gene_names
        self.barcodes = barcodes
        self.features = features
        self.cohort_adjustment = False
        self.number_of_genes = len(gene_names)
        self.number_of_cells = len(barcodes)
        self.cells_information = Cell_Inf_List(self.number_of_cells)

    def makeover_genes(self, ens_id_list):
        """
        TODO: build that function from scratch.
        Change the dimension of genes according to the Cohort whole genes. Add zeros in places where there is no
        apparent gene and rearrange the order of the genes.
        :param ens_id_list:
        :return:
        """
        all_ens_indexes = self.map_indexes(ens_id_list)
        fixed_counts = np.zeros((self.number_of_cells, len(ens_id_list)))
        fixed_counts[:, all_ens_indexes] = self.counts
        self.ens_id_list = ens_id_list
        self.number_of_genes = len(all_ens_indexes)
        self.counts = fixed_counts

    def map_indexes(self, all_genes):
        """
        TODO: Not in use, consider removing it.
        :param all_genes:
        :return:
        """
        current_genes = self.gene_names
        all_genes_indexes = []
        for gene in self.gene_names:
            idx = binary_search(all_genes, gene)
            if idx == -1:
                raise ValueError
            all_genes_indexes.append(idx)

        return all_genes_indexes

    def __getitem__(self, item):
        """
        TODO: currently it's only a pattern to build various index functions. use that pattern to build indexing
        :param item:
        :return:
        """
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


class Cell_Inf_List:

    def __init__(self, n_of_cells=0):
        self.cells_information_list = [Cell_information() for _ in range(n_of_cells)]
        self.number_of_cells = n_of_cells

    def __getitem__(self, item):
        if isinstance(item, int) or isinstance(item, np.int16) or isinstance(item, np.int32) or isinstance(item, np.int64
                                                                                                           ):
            return self.cells_information_list[item]
        if isinstance(item, slice):
            return self.cells_information_list[item]
        if isinstance(item, list):
            # identify if we are dealing with binary indexes or explicit indexes.
            if sum([(ii == 0 or ii == 1) for ii in item]) == len(item):
                # converts to explicit indexes
                item = [i for i in range(len(self)) if item[i]]
            return [c_inf for c_idx, c_inf in enumerate(self.cells_information_list) if c_idx in item]
        if isinstance(item, np.ndarray):
            if item.dtype == bool:
                item = np.where(item)[0]
            item = item.tolist()
            return [c_inf for c_idx, c_inf in enumerate(self.cells_information_list) if c_idx in item]

    def setattr(self, attr_name, idx_list, val):
        """
        :param attr_name:
        :param idx_list: Numeric indexes
        :param val:
        :return:
        """
        for idx in idx_list:
            setattr(self.cells_information_list[idx], attr_name, val)

    def setattr_list(self, attr_name, idx_list, val_list):
        for idx, val in zip(idx_list, val_list):
            setattr(self.cells_information_list[idx], attr_name, val)

    def getattr(self, attr_name, idx_list=None):
        if idx_list:
            return [getattr(self.cells_information_list[idx], attr_name) for idx in idx_list]
        else:
            return [getattr(c_inf, attr_name) for c_inf in self.cells_information_list]


class Cell_information:
    def __init__(self):
        """
        conflict_related_cell_types - names of cell-types which the cell was assigned to those cell-types but there was
        a conflict with negative markers. Therefore can't be that there is a cell-type in conflict_related_cell_types
        that also appears in cell_type_list and vise-versa. But, it's possible and actually this is the case many
        times, that cell that has cell-type in conflict_related_cell_types also assigned to be cancer. and it's
        not defined as cancer_immune_conflict if there is not cell-type in cell_type_list.

        is_apoptosis - defined in separated process after assigning cell to cell-types (cancer and immune), therefore
        can be any other case too. For instance - can contain cell-type in cell_type_list and also be apoptosis.
        Note, you have to run remove_apoptosis_cells.py first in order to see a value in that property.

        is_cancer - Only if there was an immune classification (that wasn't removed due to neg-pos conflict) and also
        cancer classification right after.


        is_classified - TRUE if the final classification is immune (not related to cancer classification).

        is_doublet - TRUE is was found to be doublet by Scrublet
        """

        self.cell_type_list = []
        self.conflict_related_cell_types = []
        self.is_apoptosis = False
        self.is_immune = False
        self.is_cancer = False
        self.cancer_immune_conflict = False
        self.is_doublet = False
        self.is_lymphoid = False
        self.is_myeloid = False


