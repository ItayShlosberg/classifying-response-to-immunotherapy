"""
Main script to manage experiments with XGBoost.
"""

import sys
# from Models.feature_explorer import Feature_Explorer
lib = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/'
import sys
sys.path.append(lib)
from utilities.package_importing import *

print("\n\n\n\n########## EXPERIMENT HAS STARTED ##########\n\n\n")
# from MachineLearning.Models.enhanced_xgboost import DROPLETseq_Enhanced_XGboost
# from DL.data_loading import *
# from utilities.smart_seq_dataset import *
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score, accuracy_score, recall_score
from sklearn import metrics
# from utilities.general_helpers import *
# from xgboost import plot_importance
import pickle
from MachineLearning.ML_utilities.dataloder import *
from MachineLearning.ML_utilities.general_utilities import *
from MachineLearning.Models.enhanced_xgboost import *

# CONFIG_PATH = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/cfg/server_cfg.yaml'
# CONFIG_PATH = r'cfg/server_cfg.yaml'
# CONFIG_PATH = r'cfg\factory_cfg\variance_2_test_percent0.30000000000000004_patients_post_cfg.yaml'
# CONFIG_PATH = r'C:\Users\KerenYlab\Desktop\Technion studies\Keren laboratory\python_playground\classifying-response-to-immunotherapy\MachineLearning\cfg\dummy.yaml'
CONFIG_PATH = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/MachineLearning/cfg/server_dummy.yaml'
# CONFIG_PATH = sys.argv[1] if len(sys.argv)>1 else r'cfg\dummy.yaml' # for terminal with outer config operation
K = 20

EXPERIMENT_NAME, EXPERIMENTS_FOLDER, config = load_yml(CONFIG_PATH)


def use_markers_to_filter_cohort(cohort, MARKERS_FOLDER_PATH):
    NON_RESPONSE_MARKER_PATH = join(MARKERS_FOLDER_PATH, r'non_response_immune_markers.xlsx')
    RESPONSE_MARKER_PATH = join(MARKERS_FOLDER_PATH, r'response_immune_markers.xlsx')
    non_response_markers = pd.read_excel(NON_RESPONSE_MARKER_PATH)
    response_markers = pd.read_excel(RESPONSE_MARKER_PATH)
    response_markers_indexes = [cohort.features.index(feature) for feature in response_markers.features if feature in cohort.features]
    non_response_markers_indexes = [cohort.features.index(feature) for feature in non_response_markers.features if feature in cohort.features]
    marker_indexes = response_markers_indexes + non_response_markers_indexes[:66]
    cohort.filter_gene_by_indexes(marker_indexes)

def build_datasets(dataset_config):
    # extracts params.
    COHORT_PATH = dataset_config['data_path']
    split_data_path = dataset_config['split_data_path']
    shuffle_dataset = dataset_config['shuffle']
    supervised_classification = dataset_config['supervised_classification']
    MARKERS_FOLDER_PATH = dataset_config['markers_folder_path']
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    use_markers_to_filter_cohort(cohort, MARKERS_FOLDER_PATH)
    data_loader = Droplet_data_manager(cohort, split_data_path, supervised_classification)
    return data_loader


@experiment_manager(EXPERIMENT_NAME, EXPERIMENTS_FOLDER, CONFIG_PATH)
def main(dataset_config, xgboost_config, experiment_config):
    # Extracts params.
    experiment_path = os.path.join(EXPERIMENTS_FOLDER, EXPERIMENT_NAME)
    num_round = xgboost_config['num_round']
    early_stopping_rounds = xgboost_config['early_stopping_rounds']
    max_depth = xgboost_config['max_depth']
    #k_folds = xgboost_config['k_folds']
    model_path = None #xgboost_config['model_path']

    # Builds datasets
    print(f'Load dataset')
    data_loader = build_datasets(dataset_config)
    x_train, x_test, y_train, y_test = data_loader.train_X, data_loader.test_X, data_loader.train_Y, data_loader.test_Y
    x_val, y_val = data_loader.validation_X, data_loader.validation_Y
    feature_names = data_loader.feature_names

    model = NGModel(feature_names, num_round, early_stopping_rounds, max_depth)

    model.train(x_train, y_train, x_val, y_val, verbose=False)
    # xgb.plot_tree(model.model, num_trees=2)
    # plt.show()
    # model.train(X_train, y_train, verbose=False)
    # Inferences on train set.
    print("Train inference")
    y_pred = model.inference(x_train)
    Metrics(y_train, y_pred).print_scores('Train cell')

    print("val inference")
    y_pred = model.inference(x_val)
    Metrics(y_val, y_pred).print_scores('val cell')

    print("Test inference")
    y_pred = model.inference(x_test)
    Metrics(y_test, y_pred).print_scores('Test cell')

    # # patients_preds, cells_preds, df_groupby_patients = model.inference(train_dataset)
    # patients_labels = df_groupby_patients.values.T[0]
    # cells_labels = np.array([p.response_label for p in train_dataset.cells_information_list])
    # # print(f'   train set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    # # print(f'train set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    # print(df_groupby_patients)
    # cells_acc = round(accuracy_score(cells_labels, cells_preds), 3)
    # patients_acc = round(accuracy_score(patients_labels, patients_preds), 3)
    # print(f'Train cells predictions CM:\n{visualization_confusion_matrix(cells_labels, cells_preds, f" {EXPERIMENT_NAME} inference train set cells, accuracy: {cells_acc}", join(experiment_path, f"CM train cells {cells_acc}"))}')
    # print(f'Train patients predictions CM:\n{visualization_confusion_matrix(patients_labels, patients_preds, f"{EXPERIMENT_NAME} inference train set patients, accuracy: {patients_acc}", join(experiment_path, f"CM train patients {patients_acc}"))}')
    # print(f'Train cells classification accuracy:\n{cells_acc}')
    # print(f'Train patients classification accuracy:\n{patients_acc}')
    #
    #
    # # Inferences on test set.
    # print("----------------------------------------------")
    # print("Test inference")
    # patients_preds, cells_preds, df_groupby_patients = model.inference(test_dataset)
    # patients_labels = df_groupby_patients.values.T[0]
    # cells_labels = np.array([p.response_label for p in test_dataset.cells_information_list])
    # # print(f'   test set cells predictions AUC: {roc_auc_score(cells_labels, avg_prob_cells_predictions)}')
    # # print(f'test set patients predictions AUC: {roc_auc_score(patients_labels, patients_predictions)}')
    # print(df_groupby_patients)
    # cells_acc = round(accuracy_score(cells_labels, cells_preds), 3)
    # patients_acc = round(accuracy_score(patients_labels, patients_preds), 3)
    # print(f'Test cells predictions CM:\n{visualization_confusion_matrix(cells_labels, cells_preds, f" {EXPERIMENT_NAME} inference test set cells, accuracy: {cells_acc}", join(experiment_path, f"CM test cells {cells_acc}"))}')
    # print(f'Test patients predictions CM:\n{visualization_confusion_matrix(patients_labels, patients_preds, f"{EXPERIMENT_NAME} inference test set patients, accuracy: {patients_acc}", join(experiment_path, f"CM test patients {patients_acc}"))}')
    # print(f'Test cells classification accuracy:\n{cells_acc}')
    # print(f'Test patients classification accuracy:\n{patients_acc}')
    #
    # print("----------------------------------------------")
    # print("Feature Importance")
    # explorer = Feature_Explorer(model, K)
    # explorer.k_importance_genes(test_dataset.gene_names)
    # plot_importance(model.model_layers[0][0])


if __name__ == '__main__':
    print(config)
    dataset_config = config['DATASET']
    xgboost_config = config['XGBOOST']
    experiment_config = config['EXPERIMENT']

    main(dataset_config, xgboost_config, experiment_config)
