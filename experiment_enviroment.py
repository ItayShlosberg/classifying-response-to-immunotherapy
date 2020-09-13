from Models.enhanced_xgboost import Enhanced_XGboost
from DL.data_loading import *
from utilities.dataset import *
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score
from sklearn import metrics
from utilities.general_helpers import *


PICKLE_PATH = r'DATA\1-16291cells_all_protein_conding_genes(withoutFilterByVariance).p'
# CONFIG_PATH = r'cfg\dummy.yaml'
CONFIG_PATH = r'cfg\factory_cfg\variance_2_test_percent0.30000000000000004_patients_post_cfg.yaml'
CONFIG_PATH = r'cfg\xgboost_1_cfg.yaml'
# CONFIG_PATH = sys.argv[1] if len(sys.argv)>1 else r'cfg\xgboost_1_cfg.yaml' # for terminal with outer config operation

EXPERIMENT_NAME, EXPERIMENTS_FOLDER, config = load_yml(CONFIG_PATH)


@experiment_manager(EXPERIMENT_NAME, EXPERIMENTS_FOLDER)
def main(test_percent, patients, num_round, early_stopping_rounds, k_folds, variance):
    cells, gene_names, patients_information = extract_data_from_pickle(PICKLE_PATH)
    origin_dataset = RNAseq_Dataset(cells, patients_information, gene_names)
    if patients == 'post':
        origin_dataset = origin_dataset.get_post_patients_sub_dataset()
    elif patients == 'pre':
        origin_dataset = origin_dataset.get_baseline_patients_sub_dataset()
    if variance:
        origin_dataset.filter_genes_by_variance(variance)
    _, _, _, _, train_idxs, test_idxs = origin_dataset.train_test_split(test_size=test_percent, shuffle=True)

    # Builds datasets and enhanced XGBoost model.
    train_dataset = origin_dataset[train_idxs]
    test_dataset = origin_dataset[test_idxs]
    model = Enhanced_XGboost(num_round, early_stopping_rounds, k_folds)

    # train
    model.train(train_dataset, True)

    # Inference on train set
    avg_prob_cells_predictions, cells_predictions, patients_predictions, df_groupby_patients = model.inference(
        train_dataset)
    patients_labels = df_groupby_patients.values.T[0]
    cells_labels = np.array([p.response_label for p in train_dataset.patients])
    print(f'   train set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    print(f'train set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    visualization_confusion_matrix(cells_labels, cells_predictions)
    visualization_confusion_matrix(patients_labels, patients_predictions)

    # Inference on train set
    avg_prob_cells_predictions, cells_predictions, patients_predictions, df_groupby_patients = model.inference(
        test_dataset)
    patients_labels = df_groupby_patients.values.T[0]
    cells_labels = np.array([p.response_label for p in test_dataset.patients])
    print(f'   test set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    print(f'test set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    visualization_confusion_matrix(cells_labels, cells_predictions)
    visualization_confusion_matrix(patients_labels, patients_predictions)

    # model.save_model_in_pkl()


if __name__ == '__main__':
    print(config)
    test_percent = config['DATASET']['test_percent']
    patients = config['DATASET']['patients']
    variance = config['DATASET']['variance']
    num_round = config['XGBOOST']['num_round']
    early_stopping_rounds = config['XGBOOST']['early_stopping_rounds']
    k_folds = config['XGBOOST']['k_folds']
    main(test_percent, patients, num_round, early_stopping_rounds, k_folds, variance)
