import pickle


def extract_data_from_pickle(pickle_path):
    """
    Retrieves data from PC located in PICKLE_PATH.
    :return: cells_form, gene_names, patients_information
    """
    cells_form, gene_names, patients_information = pickle.load(open(pickle_path, "rb"))
    return cells_form, gene_names, patients_information
