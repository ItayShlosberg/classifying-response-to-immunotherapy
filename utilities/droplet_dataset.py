"""
Python structures for droplet data of 2020.
"""
from utilities.general_helpers import *
import numpy as np
from os.path import join
import numpy as np
import pickle
from DL.Mars_seq_DL.data_loading import *
import pandas as pd
import time
CELL_TYPE_LIST = ['T cells', 'CD4 helper T cells', 'CD8 Cytotoxic T cells', 'Regulatory T cells', 'Regulatory CD4 T cells', 'Regulatory CD8 T cells', 'Regulatory CD4_CD8 T cells', 'NKT cells', 'NK cells', 'B cells', 'Activated T cells', 'Senescence T cells', 'Terminal effector', 'Exhausted T cells', 'Stem_like T cells', 'Memory T cells', 'Memory CD4 T cells', 'Memory CD8 T cells', 'Memory CD4_CD8 T cells', 'Macrophage_immature', 'Macrophage_mature', 'Monocyte_immature', 'Monocyte_mature', 'cDCs_dendritic_cells', 'pDCs', 'myeloid cells_general_immature', 'myeloid cells_general_mature', 'Neutrophils', 'Granolocytes', 'Immune_general']


def normalize_data(counts):
    """
    Will normalize each ***CELL*** separately:
    sum = sum(Cell) # sum up all reads in cell of all genes.
    scaling_factor = sum / 10,000
    For each gene_expression:
        if expression != 0:
            normalized_gene_expression = Log2(gene_expression / scaling_factor + 1)
        if expression = 0:
        normalized_gene_expression = 0
    :return: 1 - if normalizing succeeded, 0 - otherwise (already normalized)
    """

    sum_cells = counts.sum(axis=1)  # for each cell
    scaling_factors = sum_cells / 10000

    normalized_cells = np.log2((counts / scaling_factors[:, None]) + 1)
    normalized_cells[np.isneginf(normalized_cells)] = 0
    return normalized_cells


def build_cohort_gene_list(samples_information_path, save_path=None):

    samples = [subfolder for subfolder in os.listdir(samples_information_path) if not 'csv' in subfolder]

    gene_ids = []
    # loop over all samples and add each of them into the cohort.
    for sample_id in samples:
        # retrieve sample from PC
        sample_information = pickle.load(open(join(samples_information_path, sample_id), 'rb'))
        features = sample_information[2]
        gene_names = sample_information[3]
        gene_ids += list(zip(features, gene_names))

    gene_ids = sorted(list(set(gene_ids)), key=lambda x: x[0])

    if save_path:
        print(f"Saving gene list in {save_path}")
        pickle.dump(gene_ids, open(save_path, 'wb'))

    return gene_ids


def build_cohort(samples_row_data_path, samples_cells_information_path, gene_id_list, to_normalize=True):

    row_samples = [subfolder.replace(".pkl", "") for subfolder in os.listdir(samples_row_data_path) if not 'csv' in subfolder]
    information_samples = [subfolder.replace(".pkl", "") for subfolder in os.listdir(samples_cells_information_path) if not 'csv' in subfolder]
    # Check compatibility of row data with meta data
    if not are_the_lists_identical(row_samples, information_samples):
        raise Exception("Paths are incompatible")

    accumulative_counting_table = None
    cohort_gene_names = [gg[1] for gg in gene_id_list]
    cohort_gene_features = [gg[0] for gg in gene_id_list]
    cohort_barcodes = []
    cohort_mapping_samples = []
    cohort_cells_information = Cell_Inf_List()
    # loop over all samples and add each of them into the cohort.
    for idx, sample_id in enumerate(sorted(row_samples)):

        print(f"\nWorking on {sample_id}, {idx+1}/{len(row_samples)}")
        start_time = time.time()
        # retrieve sample from PC ###
        rna_sample = loading_sample(row_data_path=join(samples_row_data_path, f'{sample_id}.pkl'),
                                    cells_information_path=join(samples_cells_information_path, f'{sample_id}.pkl'))
        print(f'Loading sample time: {time.strftime("%H:%M:%S", time.gmtime(time.time()-start_time))}')

        ### Remove garbage cells ###
        print('Clean data')
        rna_sample = rna_sample.filter_cells_by_property('should_be_removed', False)

        ### Remove garbage cells ###
        if to_normalize:
            print("Normalize Data")
            start_time = time.time()
            rna_sample.normalize_data()
            print(f'Normalizing sample time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')

        sample_number_of_cells = rna_sample.number_of_cells

        ### fill cohort barcodes and mapping_samples ###
        cohort_barcodes += rna_sample.barcodes
        cohort_mapping_samples += [sample_id]*sample_number_of_cells
        cohort_cells_information = cohort_cells_information + rna_sample.cells_information

        ### gene alignment ###
        print("Align sample")
        start_time = time.time()
        # take the corresponding indices of sample_id gene in cohort_gene_features
        cohort_gene_indices = [cohort_gene_features.index(g_id) for g_id in rna_sample.features]
        # build an array with zeros
        aligned_counting_table = np.zeros((sample_number_of_cells, len(gene_id_list)))
        # fill 'aligned counting table' using cohort_gene_indices
        aligned_counting_table[:, cohort_gene_indices] = rna_sample.counts
        print(f'Alignment sample time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')


        ### add current aligned sample into accumulative counting table ###
        print("Add current aligned sample into accumulative counting table")
        start_time = time.time()
        if accumulative_counting_table is not None:
            accumulative_counting_table = np.concatenate((accumulative_counting_table, aligned_counting_table))
        else:
            accumulative_counting_table = aligned_counting_table
        print(f'Concatenate time: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')


        del aligned_counting_table
        del rna_sample
        del cohort_gene_indices


    cohort = Cohort_RNAseq(counts=accumulative_counting_table,
                          gene_names=cohort_gene_names,
                          barcodes=cohort_barcodes,
                          features=cohort_gene_features,
                          samples=cohort_mapping_samples,
                          cells_information=cohort_cells_information)
    return cohort


def loading_sample(row_data_path, cells_information_path=None):
    """
    Loading sample from path,
    note: the sample is stored as row_data (with barcodes for identification).
    and cells_information separately. that way it's possible to store the row_data only once while
    performing multiple analysis

    There is a sanity-check to check whether the barcodes of the cells' information are compatible
    with the order of the barcodes in counts.
    :param row_data_path: the path of the row data on PC.
    :param cells_information_path: the path of the cells' information on PC.
    :return: The RNAseq object
    """

    counts, barcodes_1, features_1, gene_names_1 = pickle.load(open(row_data_path, 'rb'))
    cells_information = None
    if cells_information_path:
        cells_information, barcodes_2, features_2, gene_names_2 = pickle.load(open(cells_information_path, 'rb'))

        if not are_the_lists_identical(features_1, features_2) or not are_the_lists_identical(gene_names_1, gene_names_2) or not are_the_lists_identical(barcodes_1, barcodes_2):
            raise Exception("Paths of sample are incompatible")

    rna_sample = RNAseq_Sample(counts, gene_names_1, barcodes_1, features_1, cells_information)
    return rna_sample


class Cohort_RNAseq:
    """
    Storing all samples in one object. Similar to "RNAseq_Sample", except of having list of
    sample ids corresponding to cells (in the same order), mapping each cell to sample.
    All cells of all samples share the same numpy.array without separation.
    Genes should be known for all samples. It means that before construction samples need to be aligned,
    and if there is a gene of one sample that another sample doesn't have, we need to add a zero-value gene to that
    sample. Only after alignment we can concatenate all cells of all samples together.
    Also, the counting-table should be normalized. the cohort table values are normalized.
    Finally, we have no garbage cell (apoptosis/empty cells). That guarantees the cohort is after QC, and ready for
    analysis.

    External builder function, for first time construction, use 'build_cohort' with path of samples.

    """
    def __init__(self, counts, gene_names, barcodes, features, samples, cells_information=None):
        """
        :param counts: 2D-numpy array with rows as genes and columns as cells. The created object'll save the count
        table in a manner which rows are the cells and columns are the genes, for convenience.
        :param gene_names: ordered according to the cell order.
        :param ens_id_list: not relevant, used in alignment preprocess.
        :param samples_list: each cell has a unique id.
        """
        if counts is not None:
            if counts.shape[0] != len(barcodes) or counts.shape[1] != len(gene_names)\
                    or counts.shape[0] != len(samples):
                raise Exception("unsuitable dimensions")
            if features and counts.shape[1] != len(features):
                raise Exception("unsuitable dimensions")
        self.counts = counts
        self.gene_names = gene_names
        self.barcodes = barcodes
        self.features = features
        self.samples = samples
        self.number_of_genes = len(gene_names)
        self.number_of_cells = len(barcodes)
        self.cells_information = Cell_Inf_List(self.number_of_cells)
        if cells_information:
            self.cells_information = cells_information


    def get_subset_by_identifiers(self, sample_list, barcodes_list):
        """
        On 8.5.21 we noticed that in sample the barcodes of the cells are unique, cell barcodes when looking at the
        entire cohort the barcodes are not unique and , the sample name should also be taken into account.
        therefore here we give an option to get a subet of the cohort by specifying the wanted barcodes and the
        respective sample ids. the lists must be the same length, sample_list[i] is the sam[le_if of barcodes_list[i]
        :param sample_list: list of samples id
        :param barcodes_list:list of barcodes of cells
        :return: a subset cohort, Cohort_RNAseq object of the cells of the given barcodes.
        """
        mapping = list(zip(self.samples, self.barcodes))
        cell_idxs = [mapping.index(pair_identifier) for pair_identifier in zip(sample_list, barcodes_list)]

        return self[cell_idxs]

    def get_subset_by_barcodes(self, barcode_list):
        print("Warning!! You are using get_subset_by_barcodes function to filter cells in COHORT,\ncell barcodes are "
              "not unique and when looking at the entire cohort, the sample name should also be taken into account.")
        indices = [self.barcodes.index(bb) for bb in barcode_list]
        return self[indices]

    def filter_cells_by_property(self, prop_name, value):
        return self[[aa == value for aa in self.cells_information.getattr(prop_name)]]

    def filter_protein_coding_genes(self, in_place=True, PROTEIN_CODING_FILE = r'/storage/md_keren/shitay/Data/tables/gene_ens_map.csv'):
        """
        Keeps only protein coding genes
        :return:
        """
        df = pd.read_csv(PROTEIN_CODING_FILE, header=None, names=['gene_id', 'gene', '1', '2'])
        protein_coding_genes = df[df['2'] == 'protein_coding']['gene'].tolist()
        protein_coding_gene_indices = [idx for idx, g in enumerate(self.gene_names) if
                                       g in protein_coding_genes]
        filtered_cells = self.counts[:, protein_coding_gene_indices]
        filter_features = [self.features[i] for i in protein_coding_gene_indices]
        filtered_genes = [self.gene_names[i] for i in protein_coding_gene_indices]
        if in_place:
            self.counts = filtered_cells
            self.gene_names = filtered_genes
            self.features = filter_features
            self.number_of_genes = len(self.gene_names)
            print(f"Dataset was cleared from non-protein coding genes")
        return RNAseq_Sample(filtered_cells,
                             filtered_genes,
                             self.barcodes,
                             filter_features,
                             self.cells_information)

    def filter_genes_by_variance(self, variance, in_place=True):
        """
        :param variance: genes with variance bigger than given value will stay, otherwise will be removed.
        :param in_place: True - if you want the current working object itself to change,
         or only returning a new calculated object.
        :return: a new object filtered from genes with variance lower than given.
        """

        big_variance_genes = np.var(self.counts, axis=0) > variance
        filtered_cells = self.counts[:, big_variance_genes]
        filterd_features = [self.features[i] for i in range(len(self.gene_names)) if big_variance_genes[i]]
        filtered_genes = [self.gene_names[i] for i in range(len(self.gene_names)) if big_variance_genes[i]]
        if in_place:
            self.counts = filtered_cells
            self.gene_names = filtered_genes
            self.features = filterd_features
            self.number_of_genes = len(self.gene_names)
            print(f"Dataset was cleared from genes with variance of less than {variance}")

        return Cohort_RNAseq(counts=filtered_cells,
                              gene_names=filtered_genes,
                              barcodes=self.barcodes,
                              features=filterd_features,
                              samples=self.samples,
                              cells_information=self.cells_information)

    def __getitem__(self, item):
        """
        :param item: indices (int,list,slice) -
        int: return a value in place.
        list: can be bool (each place - is taken/not taken) or numerical
        (the numbers corresponding the actual places you're interested)

        :return: int - cell vector of genes, and sub list of cells_information.
        otherwise - sub object corresponding to actual indices.
        """
        if isinstance(item, int):
            return self.counts[item], self.cells_information[item]

        if isinstance(item, slice):
            lst_indices = list(
                range(item.indices(self.number_of_cells + 1)[0], item.indices(self.number_of_cells + 1)[1]))
            barcodes = [self.barcodes[ii] for ii in lst_indices]
            samples = [self.samples[ii] for ii in lst_indices]
            cells = self.counts[item, :]
            cells_information = Cell_Inf_List()
            cells_information.cells_information_list = self.cells_information[item]
            cells_information.number_of_cells = len(cells_information.cells_information_list)

            return Cohort_RNAseq(counts=cells,
                                 gene_names=self.gene_names,
                                 barcodes=barcodes,
                                 samples=samples,
                                 features=self.features,
                                 cells_information=cells_information)

        if isinstance(item, np.ndarray):
            item = item.tolist()
            return self[item]

        if isinstance(item, list):
            # identifies if we are dealing with binary indexes or explicit indexes.
            if sum([(ii == 0 or ii == 1) for ii in item]) == len(item):
                # converts to explicit indexes.
                item = [i for i in range(len(item)) if item[i]]

            barcodes = [self.barcodes[ii] for ii in item]
            samples = [self.samples[ii] for ii in item]
            cells = self.counts[item, :]
            cells_information = Cell_Inf_List()
            cells_information.cells_information_list = self.cells_information[item]
            cells_information.number_of_cells = len(cells_information.cells_information_list)

            return Cohort_RNAseq(counts=cells,
                                 gene_names=self.gene_names,
                                 barcodes=barcodes,
                                 samples=samples,
                                 features=self.features,
                                 cells_information=cells_information)

    def __add__(self, other):
        """

        :param other:
        :return:
        """

        if other.number_of_genes != self.number_of_genes or \
                sum([self.features[i]!=self.features[i] for i in range(self.number_of_genes)]) > 0:
            raise AssertionError('genes are not compatible')


        return Cohort_RNAseq(np.concatenate((self.counts, other.counts)),
                             self.gene_names,
                             self.barcodes + other.barcodes,
                             self.features,
                             self.samples + other.samples,
                             self.cells_information + other.cells_information)

    def get_cancer_immune_stroma_map(self):
        """
        :return: bool and str lists according to the cell information object
        """
        is_stroma = self.cells_information.getattr('is_stromal')
        is_cancer = self.cells_information.getattr('is_cancer')
        is_immune = self.cells_information.getattr('is_immune')

        map = [0] * self.number_of_cells
        map_str = ['Nothing'] * self.number_of_cells
        for i in range(self.number_of_cells):
            if is_immune[i]:
                map[i] = 1
                map_str[i] = 'is_immune'
            if is_cancer[i]:
                map[i] = 2
                map_str[i] = 'is_cancer'
            if is_stroma[i]:
                map[i] = 3
                map_str[i] = 'is_stroma'
        return map, map_str

    def get_myeloid_lymphoid_map(self):
        """
        :return: bool and str lists according to the cell information object
        """
        is_myeloid = self.cells_information.getattr('is_myeloid')
        is_lymphoid = self.cells_information.getattr('is_lymphoid')
        map = [0] * self.number_of_cells
        map_str = ['Nothing'] * self.number_of_cells
        for i in range(self.number_of_cells):
            if is_myeloid[i]:
                map[i] = 1
                map_str[i] = 'myeloid'
            if is_lymphoid[i]:
                map[i] = 2
                map_str[i] = 'lymphoid'
            if is_myeloid[i] and is_lymphoid[i]:
                map[i] = 3
                map_str[i] = 'Both'
        return map, map_str


class RNAseq_Sample:

    def __init__(self, counts, gene_names, barcodes, features, cells_information=None):
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
        # self.cohort_adjustment = False
        self.is_normalized = False
        self.number_of_genes = len(gene_names)
        self.number_of_cells = len(barcodes)
        self.cells_information = Cell_Inf_List(self.number_of_cells)
        if cells_information:
            self.cells_information = cells_information

    def get_subset_by_barcodes(self, barcode_list):
        indices = [self.barcodes.index(bb) for bb in barcode_list]
        return self[indices]

    def filter_protein_coding_genes(self, in_place=True, PROTEIN_CODING_FILE = r'/storage/md_keren/shitay/Data/tables/gene_ens_map.csv'):
        """
        Keeps only protein coding genes
        :return:
        """
        df = pd.read_csv(PROTEIN_CODING_FILE, header=None, names=['gene_id', 'gene', '1', '2'])
        protein_coding_genes = df[df['2'] == 'protein_coding']['gene'].tolist()
        protein_coding_gene_indices = [idx for idx, g in enumerate(self.gene_names) if
                                       g in protein_coding_genes]
        filtered_cells = self.counts[:, protein_coding_gene_indices]
        filter_features = [self.features[i] for i in protein_coding_gene_indices]
        filtered_genes = [self.gene_names[i] for i in protein_coding_gene_indices]
        if in_place:
            self.counts = filtered_cells
            self.gene_names = filtered_genes
            self.features = filter_features
            self.number_of_genes = len(self.gene_names)
            print(f"Dataset was cleared from non-protein coding genes")
        return RNAseq_Sample(filtered_cells,
                             filtered_genes,
                             self.barcodes,
                             filter_features,
                             self.cells_information)

    def filter_cells_by_property(self, prop_name, value):
        return self[[aa == value for aa in self.cells_information.getattr(prop_name)]]

    def filter_genes_by_variance(self, variance, in_place=True):
        big_variance_genes = np.var(self.counts, axis=0) > variance
        filtered_cells = self.counts[:, big_variance_genes]
        filter_features = [self.features[i] for i in range(len(self.gene_names)) if big_variance_genes[i]]
        filtered_genes = [self.gene_names[i] for i in range(len(self.gene_names)) if big_variance_genes[i]]
        if in_place:
            self.counts = filtered_cells
            self.gene_names = filtered_genes
            self.features = filter_features
            self.number_of_genes = len(self.gene_names)
            print(f"Dataset was cleared from genes with variance of less than {variance}")
        return RNAseq_Sample(filtered_cells,
                             filtered_genes,
                             self.barcodes,
                             filter_features,
                             self.cells_information)

    def __getitem__(self, item):
        """
        TODO: currently it's only a pattern to build various index functions. use that pattern to build indexing
        :param item:
        :return:
        """
        if isinstance(item, int):
            return self.counts[item], self.cells_information[item]
        if isinstance(item, slice):
            lst_indices = list(
                range(item.indices(self.number_of_cells + 1)[0], item.indices(self.number_of_cells + 1)[1]))
            barcodes = [self.barcodes[ii] for ii in lst_indices]
            cells = self.counts[item, :]
            cells_information = Cell_Inf_List()
            cells_information.cells_information_list = self.cells_information[item]
            cells_information.number_of_cells = len(cells_information.cells_information_list)
            return RNAseq_Sample(counts=cells,
                                 gene_names=self.gene_names,
                                 barcodes=barcodes,
                                 features=self.features,
                                 cells_information=cells_information)
        if isinstance(item, np.ndarray):
            item = item.tolist()
            return self[item]

        if isinstance(item, list):
            # identify if we are dealing with binary indexes or explicit indexes.
            if sum([(ii == 0 or ii == 1) for ii in item]) == len(item):
                # converts to explicit indexes.
                item = [i for i in range(len(item)) if item[i]]

            barcodes = [self.barcodes[ii] for ii in item]
            cells = self.counts[item, :]
            cells_information = Cell_Inf_List()
            cells_information.cells_information_list = self.cells_information[item]
            cells_information.number_of_cells = len(cells_information.cells_information_list)

            return RNAseq_Sample(counts=cells,
                                 gene_names=self.gene_names,
                                 barcodes=barcodes,
                                 features=self.features,
                                 cells_information=cells_information)

    def normalize_data(self):
        """
        in_place function.
        Will normalize each ***CELL*** separately:
        sum = sum(Cell) # sum up all reads in cell of all genes.
        scaling_factor = sum / 10,000
        For each gene_expression:
            if expression != 0:
                normalized_gene_expression = Log2(gene_expression / scaling_factor + 1)
            if expression = 0:
            normalized_gene_expression = 0
        :return: 1 - if normalizing succeeded, 0 - otherwise (already normalized)
        """

        if self.is_normalized:
            return 0
        sum_cells = self.counts.sum(axis=1)  # for each cell
        scaling_factors = sum_cells / 10000

        normalized_cells = np.log2((self.counts / scaling_factors[:, None]) + 1)
        normalized_cells[np.isneginf(normalized_cells)] = 0
        self.counts = normalized_cells
        self.is_normalized = True
        return 1

    def get_statistics(self):
        """
        :return: statistics DFs, one DF of cancer&immune&stroma and one DF of all immune cell-types.
        """
        n_cells = self.number_of_cells
        n_cancer = sum(self.cells_information.getattr('is_cancer'))
        n_stromal = sum(self.cells_information.getattr('is_stromal'))
        immune_bool_indices = self.cells_information.getattr('is_immune')
        n_immune = sum(immune_bool_indices)

        cell_division_quantity = {'immune': n_immune, 'cancer': n_cancer, 'stromal': n_stromal}
        cell_division_df = pd.DataFrame(
            [list(cell_division_quantity.values()), [ii / n_cells for ii in cell_division_quantity.values()]],
            columns=cell_division_quantity.keys())

        immune_sample = self[immune_bool_indices]
        immune_counter = dict(Counter(flatten_list(immune_sample.cells_information.getattr('cell_type_list'))))
        for cell_type in CELL_TYPE_LIST:
            if not cell_type in immune_counter.keys():
                immune_counter[cell_type] = 0
        immune_portions = {key: round(val / n_immune, 3) for key, val in immune_counter.items()}
        immune_cell_types_df = pd.DataFrame([immune_counter.values(), immune_portions.values()], index=['count', 'portion'],
                                            columns=immune_counter.keys())

        return {"cell division": cell_division_df, "immune cell types": immune_cell_types_df}

    def save_cells_information(self, path):
        """
        Saving pickle of the object. default - only the cells information with barcodes
        :return:
        """
        pickle.dump((self.cells_information,
                     self.barcodes,
                     self.features,
                     self.gene_names), open(path, 'wb'))

    def save_row_data(self, path):
        """
        Saving pickle of the object. default - only the cells information with barcodes
        :return:
        """
        pickle.dump((self.counts,
                     self.barcodes,
                     self.features,
                     self.gene_names), open(path, 'wb'))

    def get_cancer_immune_stroma_map(self):
        """
        :return: bool and str lists according to the cell information object
        """
        is_stroma = self.cells_information.getattr('is_stromal')
        is_cancer = self.cells_information.getattr('is_cancer')
        is_immune = self.cells_information.getattr('is_immune')

        map = [0] * self.number_of_cells
        map_str = ['Nothing'] * self.number_of_cells
        for i in range(self.number_of_cells):
            if is_immune[i]:
                map[i] = 1
                map_str[i] = 'is_immune'
            if is_cancer[i]:
                map[i] = 2
                map_str[i] = 'is_cancer'
            if is_stroma[i]:
                map[i] = 3
                map_str[i] = 'is_stroma'
        return map, map_str

    def get_myeloid_lymphoid_map(self):
        """
        :return: bool and str lists according to the cell information object
        """
        is_myeloid = self.cells_information.getattr('is_myeloid')
        is_lymphoid = self.cells_information.getattr('is_lymphoid')
        map = [0] * self.number_of_cells
        map_str = ['Nothing'] * self.number_of_cells
        for i in range(self.number_of_cells):
            if is_myeloid[i]:
                map[i] = 1
                map_str[i] = 'myeloid'
            if is_lymphoid[i]:
                map[i] = 2
                map_str[i] = 'lymphoid'
            if is_myeloid[i] and is_lymphoid[i]:
                map[i] = 3
                map_str[i] = 'Both'
        return map, map_str


class Cell_Inf_List:

    def __init__(self, n_of_cells=0):
        self.cells_information_list = [Cell_information() for _ in range(n_of_cells)]
        self.number_of_cells = n_of_cells

    def __getitem__(self, item):
        if isinstance(item, int) or isinstance(item, np.int16) or isinstance(item, np.int32) or isinstance(item,
                                                                                                           np.int64
                                                                                                           ):
            return self.cells_information_list[item]
        if isinstance(item, slice):
            return self.cells_information_list[item]
        if isinstance(item, list):
            # identify if we are dealing with binary indexes or explicit indexes.
            if sum([(ii == 0 or ii == 1) for ii in item]) == len(item) and len(item) > 1:
                # converts to explicit indexes
                item = [i for i in range(len(self)) if item[i]]
            return [self.cells_information_list[ii] for ii in item]
        if isinstance(item, np.ndarray):
            if item.dtype == bool:
                item = np.where(item)[0]
            item = item.tolist()
            return [self.cells_information_list[ii] for ii in item]

    def __len__(self):
        return self.number_of_cells

    def setattr(self, attr_name, idx_list, val):
        """
        :param attr_name:
        :param idx_list: Numeric indexes. if None all cells' attr_name will be set.
        :param val:
        :return:
        """
        if not idx_list is None:
            for idx in idx_list:
                setattr(self.cells_information_list[idx], attr_name, val)
        else:
            for idx in range(len(self)):
                setattr(self.cells_information_list[idx], attr_name, val)

    def setattr_list(self, attr_name, idx_list, val_list):
        for idx, val in zip(idx_list, val_list):
            setattr(self.cells_information_list[idx], attr_name, val)

    def getattr(self, attr_name, idx_list=None):
        if idx_list:
            return [getattr(self.cells_information_list[idx], attr_name) for idx in idx_list]
        else:
            return [getattr(c_inf, attr_name) for c_inf in self.cells_information_list]

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise Exception(f'There was an attempt to use addition operator with object of type Cell_Inf_List and {type(other)}')

        res = Cell_Inf_List()
        res.number_of_cells = other.number_of_cells + self.number_of_cells
        res.cells_information_list = self.cells_information_list + other.cells_information_list
        return res


class Cell_information:
    def __init__(self):
        """
        cell_type_list - list of all immune cell-types which the cell is associated with, final decision after
        dealing with all kinds of conflicts.

        conflict_related_cell_types - names of cell-types which the cell was assigned to those cell-types but there was
        a conflict with negative markers. Therefore can't be that there is a cell-type in conflict_related_cell_types
        that also appears in cell_type_list and vise-versa. But, it's possible and actually this is the case many
        times, that cell that has cell-type in conflict_related_cell_types also assigned to be cancer. and it's
        not defined as cancer_immune_conflict if there is not cell-type in cell_type_list.

        is_apoptosis - defined in separated process after assigning cell to cell-types (cancer and immune), therefore
        can be any other case too. For instance - can contain cell-type in cell_type_list and also be apoptosis.
        Note, you have to run remove_apoptosis_cells.py first in order to see a value in that property.

        is_immune - TRUE if the final classification is immune (that wasn't removed due to neg_markers&pos_markers conflict)
        and also there wasn't a cancer classification right after. Note: Might be that it is marked as immune and also as apoptosis.

        is_cancer - TRUE if the classification of the cell is cancer. After running classifying_cell_types.py is might be true
        only if a tumor_marker was found to be expressed and there was no conflicts with immune_markers.
        After running use_inferCNV_clustering_to_update_data.py there might be cells that tumor_makers weren't
        found and even though classified as tumor, or cells that might have tumor_markers and immune_markers that before
        weren't defined as something and now are classified as tumor since their CNV map implied they are tumor.
        Note: Might be that it is marked as cancer and also as apoptosis.

        is_stromal - True if the classification of the cell is stroma. After running classifying_cell_types.py
        the cells that didn't have classification as tumor or immune and havn't been removed due to containing
        immune/cancer conflicts. will be probably classified as stroma.

        is_doublet - TRUE is was found to be doublet by Scrublet

        self.is_CelBender_empty - True if in the output of cellBender the cell has been marked as empty cell.
        Note that CellBender output update isn't part of the pipeline and as so, there's need to open the output
        of cellBender separately for analysis.

        should_be_removed - gets a value after running use_inferCNV_clustering_to_update_data. use this flag in order to
        dismiss unwanted cells you wouldn't use in downstream analysis. (it is kind of combination of dying cells
        and cells with conflicts which in inferCNV we concluded that we have no interest in them anymore.
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
        self.is_CelBender_empty = False
        self.is_stromal = False
        self.is_epithelial = False
        self.should_be_removed = False
        self.comment = None
