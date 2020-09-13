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


cells, gene_names, patients_information = extract_data_from_pickle()
dataset = RNAseq_Dataset(cells, patients_information, gene_names)
dataset.filter_genes_by_variance(6)




# _breakpoint = 0
x_train, x_test, y_train, y_test, train_idxs, test_idxs = dataset.train_test_split()


# data = dataset.cells
# labels = np.array(dataset.patients['response_label'])
# X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2, random_state=42, stratify=labels)


dtrain = xgb.DMatrix(x_train, label=y_train)
dval = xgb.DMatrix(x_test, label=y_test)


param = {'max_depth': 20, 'eta': 1, 'objective': 'binary:logistic'}
param['nthread'] = 4
param['eval_metric'] = 'auc'
# evallist = [(dtest, 'eval'), (dtrain, 'train')]
evallist = [(dtrain, 'train'), (dval, 'validation')]


# 2. Train model
num_round = 10
bst = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=10)


def to_int(a):
    return (a>0.5).astype(np.int)

dtest = xgb.DMatrix(x_test)
ypred = bst.predict(dtest)
print(f'test AUC: {roc_auc_score(y_test, ypred)}')
ypred = to_int(ypred)
visualization_confusion_matrix(y_test, ypred)
some_func(dataset, test_idxs, ypred)

ypred = bst.predict(dtrain)
print(f'train AUC: {roc_auc_score(y_train, ypred)}')
ypred = to_int(ypred)
visualization_confusion_matrix(y_train, ypred)




_breakpoint = 0







