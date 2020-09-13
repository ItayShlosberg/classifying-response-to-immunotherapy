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

os.environ["PATH"] += os.pathsep + r'C:\Program Files\Graphviz\bin'

def visualization_confusion_matrix(labels, predictions):
    cm = confusion_matrix(labels, predictions)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=['non-response', 'response'])
    disp.plot(include_values=True,
              cmap='viridis', ax=None, xticks_rotation='horizontal',
              values_format=None)

    plt.show()


def some_func(dataset, test_idxs, y_pred):
    test_set = dataset[test_idxs]
    ll = [[p.patient_details, p.response_label, y_pred[idx]] for idx, p in enumerate(test_set.patients)]
    df = pd.DataFrame(ll, columns=["p", "label", 'pred']).groupby(['p']).mean()
    labels = df.values.T[0]
    predictions = df.values.T[1]
    print(f'THE REAL THING - AUC: {roc_auc_score(labels, predictions)}')
    print(df)
    _breakpoint = 0


def threshold_predict(a):
    return (a>0.5).astype(np.int)


cells, gene_names, patients_information = extract_data_from_pickle()
origin_dataset = RNAseq_Dataset(cells, patients_information, gene_names)
origin_dataset.filter_genes_by_variance(6)
_, _, _, _, train_idxs, test_idxs = origin_dataset.train_test_split(test_size=0.1, shuffle=True)
training_dataset = origin_dataset[train_idxs]
test_dataset = origin_dataset[test_idxs]




# K fold cross validation
print("\n\n\n###########################################")

k_validation = training_dataset.k_fold_cross_validation(5)
bsts = []
for x_train, x_val, y_train, y_val, train_idxs, test_idxs in k_validation:
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

test_labels = np.array([p.response_label for p in test_dataset.patients])
dtest = xgb.DMatrix(test_dataset.cells, label=test_labels)
for bst in bsts:
    ypred = bst.predict(dtest)
    print(f'test AUC: {roc_auc_score(test_labels, ypred)}')
    ypred = threshold_predict(ypred)
    visualization_confusion_matrix(test_labels, ypred)
    some_func(origin_dataset, test_idxs, test_labels)

# all_preictions = []
#
# def majority_vote(all_predictions):
