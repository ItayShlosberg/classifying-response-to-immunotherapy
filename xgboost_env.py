import pandas as pd
import pandas
from collections import Counter
from pycco import process
import os
from main import *
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import os
from sklearn.metrics import confusion_matrix, plot_confusion_matrix, ConfusionMatrixDisplay
from sklearn import metrics
os.environ["PATH"] += os.pathsep + r'C:\Program Files\Graphviz\bin'


def visualization_confusion_matrix(labels, predictions):
    cm = confusion_matrix(labels, predictions)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=['non-response', 'response'])
    disp.plot(include_values=True,
              cmap='viridis', ax=None, xticks_rotation='horizontal',
              values_format=None)

    plt.show()


def combine_cells_classification_to_predict_patient_response(dataset, test_idxs, y_pred):
    test_set = dataset[test_idxs]
    ll = [[p.patient_details, p.response_label, y_pred[idx]] for idx, p in enumerate(test_set.patients)]
    df = pd.DataFrame(ll, columns=["p", "labels", 'pred probabilities']).groupby(['p']).mean()
    labels = df.values.T[0]
    predictions_probs = df.values.T[1]
    print(f'TEST PATIENT CLASSIFICATION - AUC: {roc_auc_score(labels, predictions_probs)}')
    np.argmax(metrics.roc_curve(labels, predictions_probs)[1] - metrics.roc_curve(labels, predictions_probs)[0])
    fpr, tpr, thresholds = metrics.roc_curve(labels, predictions_probs)
    best_threshold = thresholds[np.argmax(tpr-fpr)]
    print(f"Best threshold {best_threshold}")
    predictions = threshold_predict(predictions_probs, best_threshold)
    df['final pred'] = predictions
    visualization_confusion_matrix(labels, predictions)
    print(df)


def threshold_predict(a, threshold=0.5):
    return (a>=threshold).astype(np.int)


def majority_vote(all_predictions):
    all_predictions = np.array(all_predictions)
    maj_vote_predictions = np.mean(all_predictions, axis=0)
    return maj_vote_predictions


if __name__ == '__main__':

    k_folds = 5   # num of folds for k folds
    cells, gene_names, patients_information = extract_data_from_pickle()
    origin_dataset = RNAseq_Dataset(cells, patients_information, gene_names)
    # origin_dataset = origin_dataset.get_post_patients_sub_dataset()
    origin_dataset.filter_genes_by_variance(6)
    _, _, _, _, train_idxs, test_idxs = origin_dataset.train_test_split(test_size=0.1, shuffle=True)
    training_dataset = origin_dataset[train_idxs]
    test_dataset = origin_dataset[test_idxs]


    # K fold cross validation
    print("\n\n\n###########################################")

    k_validation = training_dataset.k_fold_cross_validation(k_folds)
    bsts = []
    for x_train, x_val, y_train, y_val, _, _ in k_validation:
        dtrain = xgb.DMatrix(x_train, label=y_train)
        dval = xgb.DMatrix(x_val, label=y_val)

        param = {'max_depth': 20, 'eta': 1, 'objective': 'binary:logistic'}
        param['nthread'] = 4
        param['eval_metric'] = 'auc'
        # evallist = [(dtest, 'eval'), (dtrain, 'train')]
        evallist = [(dtrain, 'train'), (dval, 'validation')]

        # 2. Train model
        num_round = 10
        bst = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=10)
        bsts.append(bst)
        print("\n\n\n###########################################")

    test_labels = np.array([p.response_label for p in test_dataset.patients])
    dtest = xgb.DMatrix(test_dataset.cells, label=test_labels)
    all_predictions = []

    for bst in bsts:
        ypred = bst.predict(dtest)
        all_predictions.append(ypred)
        print(f'TEST CELLS CLASSIFICATION AUC: {roc_auc_score(test_labels, ypred)}')
        ypred = threshold_predict(ypred)
        visualization_confusion_matrix(test_labels, ypred)
        combine_cells_classification_to_predict_patient_response(origin_dataset, test_idxs, ypred)
        print("\n\n\n--------------------------------------------------")

    maj_vote_predictions = majority_vote(all_predictions)
    print(f'maj vote test AUC: {roc_auc_score(test_labels, maj_vote_predictions)}')
    maj_vote_predictions = threshold_predict(maj_vote_predictions)
    visualization_confusion_matrix(test_labels, maj_vote_predictions)
    combine_cells_classification_to_predict_patient_response(origin_dataset, test_idxs, maj_vote_predictions)