"""
Main script to manage experiments with XGBoost.
"""

import sys
# from Models.feature_explorer import Feature_Explorer
import sys
import os
lib = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\utilities\droplet_dataset'
lib2 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\utilities'
lib3 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\data_analysis'
lib4 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy'
lib5 = r'D:\Technion studies\Keren Laboratory\python_playground\classifying-response-to-immunotherapy\scripts'
sys.path.append(lib)
sys.path.append(lib2)
sys.path.append(lib3)
sys.path.append(lib4)
sys.path.append(lib5)
from os.path import join
import pandas as pd
print("\n\n\n\n########## EXPERIMENT HAS STARTED ##########\n\n\n")
from Models.enhanced_xgboost import DROPLETseq_Enhanced_XGboost
# from DL.data_loading import *
# from utilities.smart_seq_dataset import *
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score, accuracy_score, recall_score
from sklearn import metrics
from utilities.general_helpers import *
from xgboost import plot_importance
import pickle
from MachineLearning.ML_utilities.dataloder import *
from MachineLearning.ML_utilities.general_utilities import *

# CONFIG_PATH = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/cfg/server_cfg.yaml'
# CONFIG_PATH = r'cfg/server_cfg.yaml'
# CONFIG_PATH = r'cfg\factory_cfg\variance_2_test_percent0.30000000000000004_patients_post_cfg.yaml'
# CONFIG_PATH = r'cfg\dummy.yaml'
CONFIG_PATH = r'..\cfg\dummy.yaml'
# CONFIG_PATH = sys.argv[1] if len(sys.argv)>1 else r'cfg\dummy.yaml' # for terminal with outer config operation
K = 20

EXPERIMENT_NAME, EXPERIMENTS_FOLDER, config = load_yml(CONFIG_PATH)


def build_datasets(dataset_config):
    # extracts params.
    COHORT_PATH = dataset_config['data_path']
    split_data_path = dataset_config['split_data_path']
    save_division_path = dataset_config['save_division_path']
    test_percent = dataset_config['test_percent']
    patients_type = dataset_config['patients_type']
    variance = dataset_config['variance']
    protein_coding_genes = dataset_config['protein_coding_genes']
    is_data_prior_to_therapy = dataset_config['prior_to_therapy']
    shuffle_dataset = dataset_config['shuffle']
    supervised_classification = dataset_config['supervised_classification']
    filter_by_cell_type = dataset_config['specific_cell_type']
    filter_by_immune_cell_type = dataset_config['immune_cell_type']


    cohort = pickle.load(open(COHORT_PATH, 'rb'))

    # 1. keeps only genes greater than given value.
    if variance:
        cohort.filter_genes_by_variance(variance)

    # 2. keeps only protein coding genes.
    if protein_coding_genes:
        cohort.filter_protein_coding_genes(protein_coding_genes)

    data_loader = Droplet_data_manager(cohort)


    # . if the dataset is already split into test\train patients.
    if split_data_path:
        print(f"Taking exited data division from {split_data_path}")
        # train_patients_names = pickle.load(open(split_data_path, "rb"))
        # train_dataset, test_dataset = whole_rna_seq_dataset.split_by_patient_names(train_patients_names)
        pass

    else:
        print(f"Dividing data into test/train sets")
        data_loader.filter(filter_by_cell_type,
                           filter_by_immune_cell_type,
                           supervised_classification,
                           is_data_prior_to_therapy)
        if shuffle_dataset:
            data_loader.shuffle_dataset()
        X_train, X_test, y_train, y_test = data_loader.train_test_split(p_test_size=0.1, shuffle=True, verbose=False)


        if save_division_path:  # saves the new division for future use.
            data_loader.save_split(save_division_path)
            print(f"New data sets divisions saved in {save_division_path}")
            # pickle.dump((train_dataset.get_all_patients_names()), open(save_division_path, "wb"))


    return data_loader


@experiment_manager(EXPERIMENT_NAME, EXPERIMENTS_FOLDER, CONFIG_PATH)
def main(dataset_config, xgboost_config, experiment_config):
    # Extracts params.
    experiment_path = os.path.join(EXPERIMENTS_FOLDER, EXPERIMENT_NAME)
    num_round = xgboost_config['num_round']
    early_stopping_rounds = xgboost_config['early_stopping_rounds']
    k_folds = xgboost_config['k_folds']
    model_path = xgboost_config['model_path']

    # Builds datasets
    data_loader = build_datasets(dataset_config)
    X_train, X_test, y_train, y_test = data_loader.X_train, data_loader.X_test, data_loader.y_train, data_loader.y_test
    print(f'Train dataset patients: {data_loader.train_patient_set}')
    print(f'Test dataset patients: {data_loader.test_patient_set}')

    if model_path:
        model = pickle.load(open(model_path, "rb"))
    else:
        # Builds enhanced XGBoost model.
        model = DROPLETseq_Enhanced_XGboost(num_round, early_stopping_rounds, k_folds)

        # Trains.
        model.train(X_train, y_train, verbose=False)


    # Inferences on train set.
    print("Train inference")
    y_pred = model.inference(X_train)
    Metrics(y_train, y_pred).print_scores('Train cell')

    print("Test inference")
    y_pred = model.inference(X_test)
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
