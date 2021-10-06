

class Droplet_data_manager:
    """
    Manages dataset for training:
    Contains a droplet_dataset object for splitting the data into train set and test set.

    """
    def __init__(self, dataset):
        """
        :param dataset: Contains a droplet_dataset object for splitting the data into train set and test set.
        """
        self.dataset = dataset

    def train_test_split(self, test_size=0.2, shuffle=True, stratify=True):
        """
        :param test_size: in proportion to the number of patients (not proportion to nubber of cells)
        :param shuffle:
        :param random_state:
        :param stratify: Always, equal ration of responder to non-responder between test-set and train-set.
        :return: x_train, x_test, and corresponding indexes
        """

        responders = list(set([p.patient_details for p in self.cells_information_list if p.response_label]))
        non_responders = list(set([p.patient_details for p in self.cells_information_list if not p.response_label]))
        if shuffle:
            random.shuffle(responders)
            random.shuffle(non_responders)

        train_responders = responders[:int(len(responders) * (1 - test_size))]
        test_responders = responders[int(len(responders) * (1 - test_size)):]
        train_non_responders = non_responders[:int(len(non_responders) * (1 - test_size))]
        test_non_responders = non_responders[int(len(non_responders) * (1 - test_size)):]

        train_patients = train_responders + train_non_responders
        test_patients = test_responders + test_non_responders

        train_idxs = [idx for idx, p in enumerate(self.cells_information_list['patient_details']) if p in train_patients]
        test_idxs = [idx for idx, p in enumerate(self.cells_information_list['patient_details']) if p in test_patients]

        x_train = self[train_idxs]
        x_test = self[test_idxs]

        return x_train, x_test, train_idxs, test_idxs
