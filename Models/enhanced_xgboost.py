import xgboost as xgb
from sklearn import metrics
from utilities.general_helpers import *


def pick_best_threshold(labels, predictions_probs):
    """
    Uses ROC (FTR, TPR) values, in which we decide what is the best threshold. when 1-TPR = 1-FPR for simplicity.
    :param labels: in order to calculate ROC.
    :param predictions_probs: in order to calculate ROC.
    :return: the best threshold.
    """
    fpr, tpr, thresholds = metrics.roc_curve(labels, predictions_probs)
    best_threshold = thresholds[np.argmax(tpr - fpr)]
    return best_threshold


def make_threshold_prediction(probabilities, threshold):
    """
    Binary decision based on threshold.
    :param probabilities: which we check threshold.
    :return: list of 0/1 (greater or smaller than threshold).
    """
    return (probabilities >= threshold).astype(np.int)


def patients_average_cells_predictions(rna_seq_dataset, pred_prob):
    """
    mark probability of each patient as the average of all their prediction cells probabilities.
    :param rna_seq_dataset: the dataset which has been predicted.
    :param pred_prob: prediction cells probabilities of all rna_seq_dataset cells.
    :return: patients labels, patients predictions probs, df_groupby_patients
    """
    match_patient_pred = [[p.patient_details, p.response_label, pred_prob[idx]] for
                          idx, p in enumerate(rna_seq_dataset.patients)]
    df_groupby_patients = pd.DataFrame(match_patient_pred, columns=["patients",
                                                                    "labels",
                                                                    'predictions probabilities']).\
        groupby(['patients']).mean()
    patients_labels = df_groupby_patients.values.T[0]
    patients_predictions_probs = df_groupby_patients.values.T[1]
    return patients_labels, patients_predictions_probs, df_groupby_patients


class Enhanced_XGboost:
    def __init__(self, num_round, early_stopping_rounds, k_folds):
        self.params = {'max_depth': 20,
                       'eta': 1,
                       'objective': 'binary:logistic',
                       'nthread': 4,
                       'eval_metric': 'auc'}
        self.num_round = num_round
        self.early_stopping_rounds = early_stopping_rounds
        self.k_folds = k_folds
        self.model_layers = []
        self.patient_prediction_threshold = None

    def k_fold_cross_val_training(self, rna_seq_dataset, verbose=False):
        self.patient_prediction_threshold = 0
        k_fold_validation = rna_seq_dataset.k_fold_cross_validation(self.k_folds)
        for idx, x_train, x_val, y_train, y_val, train_idxs, val_idxs in enumerate(k_fold_validation):
            if verbose:
                print(f'Train XGBoost on fold number {idx + 1}|{self.k_folds}')
            # Trains current k-1 folds and validate by the k fold.
            dtrain = xgb.DMatrix(x_train, label=y_train)
            dval = xgb.DMatrix(x_val, label=y_val)
            bst = xgb.train(self.param, dtrain, self.num_round,
                            evals=[(dtrain, 'train'), (dval, 'validation')],
                            early_stopping_rounds=self.early_stopping_rounds,
                            verbose_eval=verbose)

            # Defines threshold - majority vote over all cells predictions of patient and then correlation with labels.
            pred_prob_val = bst.predict(x_val)
            val_set = rna_seq_dataset[val_idxs]
            patients_labels, patients_predictions_probs, _ = patients_average_cells_predictions(val_set, pred_prob_val)
            best_threshold_cell_probs = pick_best_threshold(labels=patients_labels,
                                                            predictions_probs=patients_predictions_probs)

            # Adds layer (XGBoost) to the model.
            self.model_layers.append(bst)
            self.patient_prediction_threshold += best_threshold_cell_probs / self.k_folds
            if verbose:
                print("\n\n\n###########################################")

    def inference(self, rna_seq_dataset):
        if len(self.model_layers) == 0:
            return "model hasn't trained"

        data_labels = np.array([p.response_label for p in rna_seq_dataset.patients])
        dtest = xgb.DMatrix(rna_seq_dataset.cells, label=data_labels)
        avg_prob_cells_predictions = np.zeros(len(data_labels))

        for idx, bst in enumerate(self.model_layers):
            pred_probs = bst.predict(dtest)
            avg_prob_cells_predictions += pred_probs / self.k_folds

        patients_labels, patients_predictions_probs, _ = patients_average_cells_predictions(rna_seq_dataset,
                                                                                            avg_prob_cells_predictions)
        patients_predictions = make_threshold_prediction(patients_predictions_probs,
                                                         self.patient_prediction_threshold)

        return patients_predictions, avg_prob_cells_predictions

        visualization_confusion_matrix(test_labels, maj_vote_predictions)

    def save_model_in_pkl(self, path, filename="Enhanced_XGboost_Model.pkl"):
        pickle.dump(self, open(os.path.join(path, filename), "wb"))
