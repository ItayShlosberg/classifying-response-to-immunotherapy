from Models.enhanced_xgboost import Enhanced_XGboost
from DL.data_loading import *
from utilities.dataset import *
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score
from sklearn import metrics
from utilities.general_helpers import *

# CONFIG_PATH = r'cfg\dummy.yaml'
CONFIG_PATH = r'cfg\factory_cfg\variance_2_test_percent0.30000000000000004_patients_post_cfg.yaml'
CONFIG_PATH = r'cfg\xgboost_1_cfg.yaml'
# CONFIG_PATH = sys.argv[1] if len(sys.argv)>1 else r'cfg\xgboost_1_cfg.yaml' # for terminal with outer config operation

EXPERIMENT_NAME, EXPERIMENTS_FOLDER, config = load_yml(CONFIG_PATH)


def build_datasets(dataset_config):
    # extracts params.
    data_path = dataset_config['data_path']
    split_data_path = dataset_config['split_data_path']
    save_division_path = dataset_config['save_division_path']
    test_percent = dataset_config['test_percent']
    patients_type = dataset_config['patients_type']
    variance = dataset_config['variance']

    cells, gene_names, patients_information = extract_data_from_pickle(data_path)
    whole_rna_seq_dataset = RNAseq_Dataset(cells, patients_information, gene_names)

    # 1. keeps only genes greater than given value.
    if variance:
        whole_rna_seq_dataset.filter_genes_by_variance(variance)

    # 2. keeps baseline or post patients..
    if patients_type == 'post':
        whole_rna_seq_dataset = whole_rna_seq_dataset.get_post_patients_sub_dataset()
    elif patients_type == 'pre':
        whole_rna_seq_dataset = whole_rna_seq_dataset.get_baseline_patients_sub_dataset()

    # 3. if there is already division of patients to datasets.
    if split_data_path:
        print(f"Taking exited data division from {split_data_path}")
        train_patients_names = pickle.load(open(split_data_path, "rb"))
        train_dataset, test_dataset = whole_rna_seq_dataset.split_by_patient_names(train_patients_names)
    else:   # uses the configuration to divide into sets.
        print(f"Dividing data into test/train sets")
        train_dataset, test_dataset, train_idxs, test_idxs = whole_rna_seq_dataset.train_test_split(test_size=test_percent,
                                                                                                    shuffle=True)
        if save_division_path:  # saves the new division for future use.
            pickle.dump((train_dataset.get_all_patients_names()), open(save_division_path, "wb"))
            print(f"New data sets divisions saved in {save_division_path}")

    return train_dataset, test_dataset


@experiment_manager(EXPERIMENT_NAME, EXPERIMENTS_FOLDER)
def main(dataset_config, xgboost_config, experiment_config):
    # Extracts params.
    num_round = xgboost_config['num_round']
    early_stopping_rounds = xgboost_config['early_stopping_rounds']
    k_folds = xgboost_config['k_folds']

    # Builds datasets
    train_dataset, test_dataset = build_datasets(dataset_config)
    print(f'Train dataset patients: {train_dataset.get_all_patients_names()}')
    print(f'Test dataset patients: {train_dataset.get_all_patients_names()}')

    # Builds enhanced XGBoost model.
    model = Enhanced_XGboost(num_round, early_stopping_rounds, k_folds)

    # Trains.
    model.train(train_dataset, True)

    # Inferences on train set.
    avg_prob_cells_predictions, cells_predictions, patients_predictions, df_groupby_patients = model.inference(
        train_dataset)
    patients_labels = df_groupby_patients.values.T[0]
    cells_labels = np.array([p.response_label for p in train_dataset.patients])
    print(f'   train set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    print(f'train set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    print(df_groupby_patients)
    print(f'Test cells predictions CM: {visualization_confusion_matrix(cells_labels, cells_predictions)}')
    print(f'Test patients predictions CM: {visualization_confusion_matrix(patients_labels, patients_predictions)}')

    # Inferences on test set.
    avg_prob_cells_predictions, cells_predictions, patients_predictions, df_groupby_patients = model.inference(
        test_dataset)
    patients_labels = df_groupby_patients.values.T[0]
    cells_labels = np.array([p.response_label for p in test_dataset.patients])
    print(f'   test set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    print(f'test set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    print(df_groupby_patients)
    print(f'Test cells predictions CM: {visualization_confusion_matrix(cells_labels, cells_predictions)}')
    print(f'Test patients predictions CM: {visualization_confusion_matrix(patients_labels, patients_predictions)}')

    # Save model.
    if experiment_config['save_model']:
        model.save_model_in_pkl(os.path.join(experiment_config['experiments_folder'], experiment_config['experiment_name']))


if __name__ == '__main__':
    print(config)
    dataset_config = config['DATASET']
    xgboost_config = config['XGBOOST']
    experiment_config = config['EXPERIMENT']

    main(dataset_config, xgboost_config, experiment_config)
