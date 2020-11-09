import numpy as np
import operator


def filter_by_indexes(indices, cells, patients_information=None):
    """
    return reduced cells associated to the corresponding indices.
    :param indices: indexes list of desirable cells.
    :param cells: cells from pkl file.
    :param patients_information: pkl file format
    :return: reduced cells. and their corresponding patients_information.
    """
    if patients_information:
        return cells[indices, :], [patients_information[i] for i in indices]
    return cells[indices, :]


def filter_by_binary_indices(indices, cells, patients_information=None):
    """
    return reduced cells associated to the corresponding indices.
    :param indices: binary list.
    :param cells: cells from pkl file.
    :param patients_information: pkl file format
    :return: reduced cells. and their corresponding patients_information.
    """
    indexes = [i for i in range(len(patients_information)) if indices[i]]
    if patients_information:
        return cells[indexes, :], [patients_information[i] for i in range(len(patients_information)) if indices[i]]
    return cells[indices, :]


def expression_of_genes(cells, gene_names = None):
    """
    :param cells: PKL format
    :param gene_names: list of gene names in the order shown in the cells.
    :return: sorted avg gene expressions and gene names sorted by avg expression.
    """
    average_val_genes = np.mean(cells, axis=0)
    indexes_of_sorted_expression_of_genes = np.argsort(average_val_genes)
    gene_sorted_by_expression = operator.itemgetter(*indexes_of_sorted_expression_of_genes)(gene_names)
    sorted_average_expression_genes = average_val_genes[indexes_of_sorted_expression_of_genes]
    return sorted_average_expression_genes, gene_sorted_by_expression


def filter_cells_by_supervised_classification(cells, patients, required_cell_type="T cells"):
    """
    patients_information contains 'supervised classification' field, which is the classification of the cell
    defined by gene expression in 80% of the cells, and the remaining 20% were done by manual process.
    the function filters cells by cell-type was done by this classification.
    :param cells: pkl format.
    :param patients_information: pkl format.
    :param required_cell_type: name of the desired cell-type.
    :return: cells of the desired cell-type
    """
    indexes_list = patients.get_cells_belong_to_cells_type(required_cell_type)
    cells, patients_information = filter_by_binary_indices(indexes_list, cells, patients_information)
    return cells, patients_information


def get_high_expressed_genes(self):
    """
    Return (only the indexes) of the highest expression genes. TODO: make the genes amount changeable.
    :param cells: pkl format
    :param gene_names: when given, the function return genes names with highest expression also.
    :return: indexes only. if gene names given, the function return genes names with highest expression also.
    """
    activated_genes = sum(np.sum(self.cells, axis=0) != 0)
    average_val_genes = np.mean(self.cells, axis=0)
    expression_histogram = np.histogram(average_val_genes)
    expression_histogram_df = np.concatenate((np.expand_dims(expression_histogram[0][:10], axis=0),
                                              np.expand_dims(expression_histogram[1][:10], axis=0)))
    amount_high_expressed_genes = min(activated_genes * 0.05, sum(expression_histogram[0][-2:]))
    high_expressed_percentage = amount_high_expressed_genes / activated_genes
    indices_of_high_expressed_genes = np.argsort(average_val_genes)[-3:]
    print(f'High expressed genes percent {high_expressed_percentage}')
    if self.gene_names:
        return indices_of_high_expressed_genes, operator.itemgetter(*indices_of_high_expressed_genes)(self.gene_names)
    return indices_of_high_expressed_genes


def combine_cells_classification_to_predict_patient_response(dataset, test_idxs, y_pred):
    test_set = dataset[test_idxs]
    ll = [[p.patient_details, p.response_label, y_pred[idx]] for idx, p in enumerate(test_set.patients)]
    df = pd.DataFrame(ll, columns=["patients", "labels", 'predictions probabilities']).groupby(['patients']).mean()
    labels = df.values.T[0]
    predictions_probs = df.values.T[1]
    print(f'TEST PATIENT CLASSIFICATION - AUC: {roc_auc_score(labels, predictions_probs)}')
    np.argmax(metrics.roc_curve(labels, predictions_probs)[1] - metrics.roc_curve(labels, predictions_probs)[0])
    best_threshold_cell_probs = pick_best_threshold(labels, predictions_probs)
    print(f"Best threshold {best_threshold_cell_probs}")
    predictions = threshold_predict(predictions_probs, best_threshold_cell_probs)
    df['final predictions'] = predictions
    # visualization_confusion_matrix(labels, predictions)
    df = df[["labels", 'final predictions', 'predictions probabilities']]
    print(df)


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
        self.cells_presictions_threshold = None

    def train(self, rna_seq_dataset, verbose=False):
        """
        important Note: it's a k_fold_cross_val training.
        Splits RNAseq dataset into folds and sets each fold as validation set iteratively. Each iteration trains XGBoost
        model on training-set. Then, predicts binary response to validation fold using trained XGBoost,
        calculate average probabilities per patient. and set threshold of each XGBoost model that maximizes TPR and
        minimizes FTR.
        :param rna_seq_dataset: The set in which the k-folds of train/validation sets.
        :param verbose: to print information during training.
        """
        self.patient_prediction_threshold = 0
        self.cells_presictions_threshold = 0
        k_fold_validation = rna_seq_dataset.k_fold_cross_validation(self.k_folds)
        for idx, (x_train, x_val, y_train, y_val, train_idxs, val_idxs) in enumerate(k_fold_validation):
            if verbose:
                print(f'Train XGBoost on fold number {idx + 1}|{self.k_folds}')
            # Trains current k-1 folds and validate by the k fold.
            dtrain = xgb.DMatrix(x_train, label=y_train)
            dval = xgb.DMatrix(x_val, label=y_val)
            bst = xgb.train(self.params, dtrain, self.num_round,
                            evals=[(dtrain, 'train'), (dval, 'validation')],
                            early_stopping_rounds=self.early_stopping_rounds,
                            verbose_eval=verbose)

            # Defines threshold - majority vote over all cells predictions of patient and then correlation with labels.
            pred_cells_prob_val = bst.predict(dval)
            val_set = rna_seq_dataset[val_idxs]
            patients_labels, patients_predictions_probs, _ = patients_average_cells_predictions(val_set, pred_cells_prob_val)
            best_threshold_patient_probs = pick_best_threshold(labels=patients_labels,
                                                            predictions_probs=patients_predictions_probs)
            best_threshold_cell_probs = pick_best_threshold(labels=y_val,
                                                            predictions_probs=pred_cells_prob_val)

            # Adds layer (XGBoost) to the model.
            self.model_layers.append(bst)
            self.patient_prediction_threshold += best_threshold_patient_probs / self.k_folds
            self.cells_presictions_threshold += best_threshold_cell_probs / self.k_folds
            if verbose:
                print("\n\n\n###########################################")

    def inference(self, rna_seq_dataset):
        """
        Model's Inference. Each XGBoost predicts cells probabilities of response. Then activates threshold on the
        average cells response.
        :param rna_seq_dataset: All the data will be predicted.
        :return: patients predictions, avg_prob_cells_predictions, df_groupby_patients['patients']
        """
        if len(self.model_layers) == 0:
            return "model hasn't trained"

        data_labels = np.array([p.response_label for p in rna_seq_dataset.patients])
        dtata = xgb.DMatrix(rna_seq_dataset.cells, label=data_labels)
        avg_prob_cells_predictions = np.zeros(len(data_labels))

        for idx, bst in enumerate(self.model_layers):
            pred_probs = bst.predict(dtata)
            avg_prob_cells_predictions += pred_probs / self.k_folds

        patients_labels, patients_predictions_probs, df_groupby_patients = patients_average_cells_predictions(
            rna_seq_dataset, avg_prob_cells_predictions)
        patients_predictions = make_threshold_prediction(patients_predictions_probs,
                                                         self.patient_prediction_threshold)
        df_groupby_patients['final predictions'] = patients_predictions
        # visualization_confusion_matrix(labels, predictions)
        df_groupby_patients = df_groupby_patients[["labels", 'final predictions', 'predictions probabilities']]
        return avg_prob_cells_predictions, \
               make_threshold_prediction(avg_prob_cells_predictions, self.cells_presictions_threshold),\
               patients_predictions, df_groupby_patients

    def save_model_in_pkl(self, path, filename="Enhanced_XGboost_Model.pkl"):
        """
        Save model local in PKL file.
        :param path: where the model be saved
        :param filename: name of the new file'll be created.
        """
        pickle.dump(self, open(os.path.join(path, filename), "wb"))