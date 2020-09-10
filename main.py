import numpy as np
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pickle
import random
import pandas as pd
from Bio.Cluster import kcluster
import operator
import seaborn as sns
import time
import sys
PICKLE_PATH = r'DATA\1-16291cells.p'
PICKLE_PATH = r'DATA\1-16291cells_supervised_classification.p'
CHECKPOINT_TSNE_PATH = r'DATA\TSNE_Embedded_1-16291cells_randInt21'  # comes into play as import OR export path.
TSNE_IMPORT_EXPORT = False  # FALSE - Import, TRUE - EXPORT


def extract_data_from_pickle():
    """
    Retrieves data from PC located in PICKLE_PATH.
    :return: cells_form, gene_names, patients_information
    """
    cells_form, gene_names, patients_information = pickle.load(open(PICKLE_PATH, "rb"))
    return cells_form, gene_names, patients_information


def cluster(data, n_clusters):
    """
    Activates sklearn.kmeans.
    :param data: cells in PKL format.
    :param n_clusters: desired number of cluster for kmeans
    :return: kmeans labels - list of cluster for each cell.
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(data)

    # print(kmeans.labels_)
    # print(kmeans.predict([[0, 0], [12, 3]]))
    # print(kmeans.cluster_centers_)
    return kmeans.labels_


def TSNE_embedded(cells):
    """
    Activates t-SNE algorithm to find 2D representation of the cells.
    calculate and export if TSNE_IMPORT_EXPORT else import from CHECKPOINT_TSNE_PATH.
    :param cells: pkl format.
    :return: embedded 2D-representation .
    """
    if TSNE_IMPORT_EXPORT:  # perform TSNE and save embedded vector in CHECKPOINT_TSNE_PATH
        cells_embedded = TSNE(n_components=2, random_state=21).fit_transform(cells)
        cells_embedded = cells_embedded.T.tolist()
        pickle.dump(cells_embedded, open(CHECKPOINT_TSNE_PATH, "wb"))
    else:
        cells_embedded = pickle.load(open(CHECKPOINT_TSNE_PATH, "rb"))
    return cells_embedded


def visualize(cells, clusters_labels, title=None, centroids=None):
    """

    :param cells:
    :param clusters_labels:
    :param title:
    :param centroids:
    :return:
    """
    X = cells[0]
    Y = cells[1]

    # plt.plot(X, Y, 'ro')

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'lime', 'lavender', 'darkred', 'olive']
    for cluster in np.unique(clusters_labels):
        Xi = [X[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
        Yi = [Y[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
        color = colors[cluster]
        plt.plot(Xi, Yi, 'ro', color=color)
    if centroids:
        for centroid in centroids:
            plt.plot(centroid[0], centroid[1], 'ro', color='y')
    # plt.plot([1, 2], [1, 4], 'ro')
    # plt.plot([3, 4], [9, 16], 'ro', color='green')
    plt.title(title)
    plt.show()


def build_confusion_matrix(classification1, classification2):
    """
    Given two different classifications of the same group (identical length of list
    with different classifications), builds confusion matrix.
    :param classification1: list in length of the sample size. each place is the classification of the
     corresponding cell. for example L[i] is the classification of cell number i.
    :param classification2: the second classification list. identified to classification1 structure.
    :return: confusion matrix in DataFrame format.
    """
    match_classifications = list(zip(classification1, classification2))
    confusion_matrix = np.zeros((np.max(classification1)+1, np.max(classification2)+1))
    for appearance in match_classifications:
        confusion_matrix[appearance[0], appearance[1]] += 1

    columns = ['V'+str(i) for i in range(max(classification2)+1)]
    index = ['G'+str(i) for i in range(max(classification1)+1)]
    confusion_matrix = pd.DataFrame(confusion_matrix, columns=columns, index=index)

    return confusion_matrix


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
    Return (only the indexes) of the highest expression genes. TODO: make the genes amount changeable.
    :param cells: pkl format
    :param gene_names: when given, the function return genes names with highest expression also.
    :return: indexes only. if gene names given, the function return genes names with highest expression also.
    """
    activated_genes = sum(np.sum(cells, axis=0) !=0)
    average_val_genes = np.mean(cells, axis=0)
    expression_histogram = np.histogram(average_val_genes)
    expression_histogram_df = np.concatenate((np.expand_dims(expression_histogram[0][:10], axis=0),
                        np.expand_dims(expression_histogram[1][:10], axis=0)))
    amount_high_expressed_genes = min(activated_genes * 0.05, sum(expression_histogram[0][-2:]))
    high_expressed_percentage = amount_high_expressed_genes / activated_genes
    indices_of_high_expressed_genes = np.argsort(average_val_genes)[-3:]
    print(f'High expressed genes percent {high_expressed_percentage}')
    if gene_names:
        return indices_of_high_expressed_genes, operator.itemgetter(*indices_of_high_expressed_genes)(gene_names)
    return indices_of_high_expressed_genes


def filter_cells_by_supervised_classification(cells, patients_information, required_cell_type="T cells"):
    """
    patients_information contains 'supervised classification' field, which is the classification of the cell
    defined by gene expression in 80% of the cells, and the remaining 20% were done by manual process.
    the function filters cells by cell-type was done by this classification.
    :param cells: pkl format.
    :param patients_information: pkl format.
    :param required_cell_type: name of the desired cell-type.
    :return: cells of the desired cell-type
    """
    supervised_classification = [p["supervised classification"] for p in patients_information]

    indices_list = [1 if required_cell_type in cl else 0 for cl in supervised_classification]
    cells, patients_information = filter_by_binary_indices(indices_list, cells, patients_information)
    return cells, patients_information


def some_correlation_function(patients_information, response_labels, clusterid):
    x=4
    ratio_non_responder = []
    threshold = [0.05 + (0.0125)*i for i in range(70)]
    threshold_score = [[0,0,0,0] for i in range(70)]
    ratio_responder = []
    for patient in set([p['patient details'] for p in patients_information]):
        indices = [i for i in range(len(patients_information)) if patients_information[i]['patient details']==patient]
        response = response_labels[indices[0]]
        cluster_cells_of_patient = [clusterid[i] for i in indices]
        cluster_1_amount = sum(cluster_cells_of_patient)
        ratio = cluster_1_amount/len(cluster_cells_of_patient)
        if response:
            ratio_responder.append(ratio)
            for i in range(70):
                if threshold[i]>ratio:
                    threshold_score[i][0]+=1
                else:
                    threshold_score[i][1] += 1
        else:
            ratio_non_responder.append(ratio)
            for i in range(70):
                if threshold[i]>ratio:
                    threshold_score[i][2]+=1
                else:
                    threshold_score[i][3] += 1



    print(f'responder:')
    print(ratio_responder)
    print(sum(ratio_responder)/len(ratio_responder))
    print(f'non-responder:')
    print(ratio_non_responder)
    print(sum(ratio_non_responder) / len(ratio_non_responder))
    print(threshold_score)
    print(threshold)


def heatmap_high_epresssed_gene_of_cluster(cells, clusters):
    """
    Shows Seaborn visualized heat-map of clusters of cells and their high expressed genes.
    Made by expression_of_genes function that returns the high expression genes and is activated
    for each cluster filtered by filter_by_indexes function.
    :param cells: pkl format.
    :param clusters: list in cells amount length, that tells for each cell what is its cluster.
    clusters[i] is the cluster of cell i.
    """
    indices_of_high_expressed_genes = []
    heatmap_order = []
    for cluster in set(clusters):
        cluster_indices = [i for i in range(len(clusters)) if cluster==clusters[i]]
        heatmap_order += cluster_indices
        cluster_cells = filter_by_indexes(cluster_indices, cells)
        # cluster_cells[:, expression_of_genes(cluster_cells, gene_names)]
        indices_of_high_expressed_genes += expression_of_genes(cluster_cells).tolist()
    # indices_of_high_expressed_genes = set(indices_of_high_expressed_genes)
    heatmap_x = cells[heatmap_order, :][:, list(indices_of_high_expressed_genes)]
    ax = sns.heatmap(heatmap_x)
    plt.show()


if __name__ == '__main__':
    # path = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\python_playground\Data\smaller_group_of_cells.txt'
    # save_smaller_file(path, 60)
    cells, gene_names, patients_information = extract_data_from_pickle()
    cells, patients_information = filter_cells_by_supervised_classification(cells, patients_information)
    # expression_of_genes(cells, patients_information, gene_names)
    # It has to be decided the order of actions of tSNE and k-means.
    # embedded_cells = TSNE_embedded(cells)
    # cells = np.array(embedded_cells).T    # in order to perform cluster on embedded points.

    keren_clusters = [p['keren cluster'] for p in patients_information]
    response_labels = [p['response label'] for p in patients_information]
    kmeans_clusters = cluster(cells, 2)

    clusterid, error, nfound = kcluster(cells, nclusters=2, dist='c')   # dist 'c' is pearson correlation distance
    # some_correlation_function(patients_information, response_labels, kmeans_clusters)
    some_correlation_function(patients_information, response_labels, clusterid)
    print(build_confusion_matrix(response_labels, clusterid))
    # heatmap(cells, keren_clusters)
    # kmeans_clusters = cluster(cells, 2)
    # print(build_confusion_matrix(response_labels, kmeans_clusters))
    #
    #
    # some_correlation_function(patients_information, response_labels, clusterid)
    # visualize(embedded_cells, keren_clusters, 'keren_clusters')
    # visualize(embedded_cells, response_labels, 'response_labels')
    # print(build_confusion_matrix(keren_clusters, clusterid))


    # Cluster algorithm
    # kmeans_clusters = cluster(cells, 2)
    # clusterid, error, nfound = kcluster(cells, nclusters=2, dist='c')   # dist 'c' is pearson correlation distance
    # clusterid2, error2, nfound2 = kcluster(cells, nclusters=2, dist='e')
    #
    # visualize(embedded_cells, kmeans_clusters, 'euc')
    # visualize(embedded_cells, clusterid, 'pearson')
    # visualize(embedded_cells, clusterid2, 'euc BIO')
    # visualize(cells, labels)

