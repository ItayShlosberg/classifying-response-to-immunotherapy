lib = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/'
import sys

sys.path.append(lib)
from utilities.package_importing import *
import scipy.stats as stats


class Droplet_data_manager:
    """
    Manages dataset for training:
    Contains a droplet_dataset object for splitting the data into train set and test set.
    """

    def __init__(self, cohort, split_data_path, cluster_number):
        """
        :param dataset: Contains a droplet_dataset object for splitting the data into train set and test set.
        """
        self.cohort = cohort
        self.feature_names = cohort.gene_names
        # melanoma_clinical_data = get_constant_cohort(CLINICAL_TABLE_NUM) #
        melanoma_clinical_data = get_clinical_data(71)
        self.clinical_table = melanoma_clinical_data[melanoma_clinical_data['ICI']]

        train_set_barcode_mapping = pd.read_excel(join(split_data_path, f'cluster_{cluster_number}_train_barcodes.xlsx'), engine='openpyxl')
        test_set_barcode_mapping = pd.read_excel(join(split_data_path, f'cluster_{cluster_number}_test_barcodes.xlsx'), engine='openpyxl')
        validation_set_barcode_mapping = pd.read_excel(
            join(split_data_path, f'cluster_{cluster_number}_validation_barcodes.xlsx'), engine='openpyxl')

        self.train_X, self.train_Y = self.get_X_Y(train_set_barcode_mapping)
        self.test_X, self.test_Y = self.get_X_Y(test_set_barcode_mapping)
        self.validation_X, self.validation_Y = self.get_X_Y(validation_set_barcode_mapping)

    def get_X_Y(self, set_barcode_mapping):
        set_cohort = self.cohort.get_subset_by_identifiers(set_barcode_mapping.Sample,
                                                           set_barcode_mapping.Barcode)
        set_X = set_cohort.counts
        set_Y = [1 if self.clinical_table.set_index('Patient id').loc[s].response == 'R' else 0 for s in set_cohort.samples]
        return set_X, set_Y


class OLD_Droplet_data_manager:
    """
    Manages dataset for training:
    Contains a droplet_dataset object for splitting the data into train set and test set.

    """

    def __init__(self, dataset,
                 CLINICAL_LABELS_PATH=r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\tables\clinical_labels.xlsx',
                 MELANOMA_CLINICAL_DATA_PATH=r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\Data\tables\edited_unprotected_Melanoma_clinical_data.xlsx'):
        """
        :param dataset: Contains a droplet_dataset object for splitting the data into train set and test set.
        """
        self.dataset = dataset
        self.clinical_table = self.load_clinical_inf(CLINICAL_LABELS_PATH, MELANOMA_CLINICAL_DATA_PATH)

        self.response_samples = self.clinical_table[self.clinical_table['response'] == 'R'].index.tolist()
        self.non_response_samples = self.clinical_table[self.clinical_table['response'] == 'NR'].index.tolist()
        self.no_used_samples = self.clinical_table[
            self.clinical_table['response'] == 'not in use for now'].index.tolist()

    def load_clinical_inf(self, CLINICAL_LABELS_PATH, MELANOMA_CLINICAL_DATA_PATH):
        # Loads xlsx files
        melanoma_clinical_data = pd.read_excel(MELANOMA_CLINICAL_DATA_PATH)
        clinical_labels = pd.read_excel(CLINICAL_LABELS_PATH)

        # takes nly first 46 samples, fill Nan and creat dictionary mapping
        melanoma_clinical_data = melanoma_clinical_data.iloc[:46, :][
            ['Patient id', 'Clinical response', 'Melanoma type']]
        melanoma_clinical_data['Clinical response'] = melanoma_clinical_data['Clinical response'].fillna('??')
        labels_mapping = {v[0]: v[1] for v in clinical_labels[['Clinical response', 'binary label']].values}

        # adds CRITERIA response into clinical df
        melanoma_clinical_data['response'] = melanoma_clinical_data['Clinical response'].apply(
            lambda x: labels_mapping[x])
        res_melanoma_clinical_data = melanoma_clinical_data.set_index('Patient id')

        # Melanoma type
        melanoma_type_translate = {'Cutaneous': 'Cutaneous', 'Mucosal ': 'Mucosal ', 'Uveal': 'other',
                                   'UN primary': 'other', 'Acral': 'other', 'Unknown': 'other'}
        res_melanoma_clinical_data['Melanoma type'] = res_melanoma_clinical_data['Melanoma type'].apply(
            lambda rr: melanoma_type_translate[rr])

        return res_melanoma_clinical_data

    def train_test_split(self, p_test_size=0.25, shuffle=True, verbose=True):
        train_patient_set, test_patient_set = self._train_test_patient_split_by_response(self.dataset,
                                                                                         p_test_size=p_test_size,
                                                                                         shuffle=shuffle,
                                                                                         verbose=verbose)
        self.train_patient_set = train_patient_set
        self.test_patient_set = test_patient_set

        self.train_cohort = self.dataset.get_subset_by_sample_list(train_patient_set)
        self.test_cohort = self.dataset.get_subset_by_sample_list(test_patient_set)

        labels_response_translator = dict(self.clinical_table['response'])
        labels_str_translator = {'NR': 0, 'R': 1}

        self.X_train = self.train_cohort.counts
        self.X_train_samples = self.train_cohort.samples
        self.X_test = self.test_cohort.counts
        self.X_test_samples = self.test_cohort.samples
        self.y_train = [labels_str_translator[labels_response_translator[ss]] for ss in self.train_cohort.samples]
        self.y_test = [labels_str_translator[labels_response_translator[ss]] for ss in self.test_cohort.samples]

        return self.X_train, self.X_test, self.y_train, self.y_test

    def shuffle_dataset(self):
        indexes = np.arange(self.dataset.number_of_cells)
        np.random.shuffle(indexes)
        self.dataset = self.dataset[indexes]

    def _train_test_patient_split_by_response(self, dataset, p_test_size=0.25, shuffle=True, verbose=True):

        non_response_cohort = dataset.get_subset_by_sample_list(self.non_response_samples)
        response_cohort = dataset.get_subset_by_sample_list(self.response_samples)

        r_train_patient_set, r_test_patient_set = self._train_test_split_patients(response_cohort.samples,
                                                                                  p_test_size=p_test_size,
                                                                                  shuffle=shuffle,
                                                                                  verbose=verbose)
        nr_train_patient_set, nr_test_patient_set = self._train_test_split_patients(non_response_cohort.samples,
                                                                                    p_test_size=p_test_size,
                                                                                    shuffle=shuffle,
                                                                                    verbose=verbose)

        train_patient_set = r_train_patient_set + nr_train_patient_set
        test_patient_set = r_test_patient_set + nr_test_patient_set
        return train_patient_set, test_patient_set

    def _train_test_split_patients(self, samples, p_test_size=0.25, shuffle=True, verbose=True):
        """
        Given a list of samples that represents the mapping samples of a
         list of cells (patient ID repeats more than once), it returns a separation of patients regardless of response
         to therapy.
        :param samples:
        :param p_test_size:
        :param shuffle:
        :param verbose:
        :return:
        """
        n_cells = len(samples)
        test_size = int(np.ceil(p_test_size * n_cells))
        samples = [list(x) for x in list(dict(Counter(samples)).items())]
        if shuffle:
            random.shuffle(samples)

        for i in range(1, len(samples)):
            samples[i][1] += samples[i - 1][1]

        samples = np.array(samples)
        argmin_idx = np.argmin(np.abs(test_size - samples[:, 1].astype(np.int)))
        test_patient_set = samples[0: argmin_idx + 1][:, 0].tolist()
        train_patient_set = samples[argmin_idx + 1:][:, 0].tolist()

        if verbose:
            print(test_patient_set, end='\n')
            print(train_patient_set, end='\n')
            print(f'Test set size: {samples[argmin_idx, 1]}')
            print(f'Train set size: {n_cells - int(samples[argmin_idx, 1])}')

        return train_patient_set, test_patient_set

    def filter(self, filter_by_cell_type, filter_by_immune_cell_type, supervised_classification,
               is_data_prior_to_therapy):
        """
        :param filter_by_cell_type:
        :param filter_by_immune_cell_type:
        :param supervised_classification:
        :param is_data_prior_to_therapy:
        :return:
        """

        if filter_by_cell_type:
            self.dataset = self.dataset.filter_cells_by_property(filter_by_cell_type, True)

        if filter_by_immune_cell_type:
            self.dataset = self.dataset.filter_cells_by_property('is_immune', True)
            indices = [filter_by_immune_cell_type in l for l in
                       self.dataset.cells_information.getattr('cell_type_list')]
            self.dataset = self.dataset[indices]

        if supervised_classification:
            pass

        # . keeps baseline or post patients..
        if is_data_prior_to_therapy:
            pass

    def save_split(self):
        pass



