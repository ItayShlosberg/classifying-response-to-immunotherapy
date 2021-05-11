"""
That script helped to explore at the first time the smart-seq data of 2018.
We used it to cluster, visualize and reproduce steps from the article.
"""

from sklearn.cluster import KMeans
from DL.Mars_seq_DL.data_loading import *
from utilities.general_helpers  import *
from bhtsne import tsne
from sklearn import preprocessing
from sklearn.decomposition import PCA

# PICKLE_PATH = r'DATA\1-16291cells.p'
PICKLE_PATH = r'D:\Technion studies\Keren Laboratory\Data\smart_seq\SmartSeq_RNAseq_DATA.p'
CHECKPOINT_TSNE_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\TSNE_Embedded_correlation.pkl'  # comes into play as import OR export path.
TSNE_IMPORT_EXPORT = True  # FALSE - Import, TRUE - EXPORT




def kmeans(data, n_clusters):
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


# def TSNE_embedded(cells):
#     """
#     Activates t-SNE algorithm to find 2D representation of the cells.
#     calculate and export if TSNE_IMPORT_EXPORT else import from CHECKPOINT_TSNE_PATH.
#     :param cells: pkl format.
#     :return: embedded 2D-representation .
#     """
#     if TSNE_IMPORT_EXPORT:  # perform TSNE and save embedded vector in CHECKPOINT_TSNE_PATH
#         cells_embedded = TSNE(n_components=2, perplexity=30.0, metric="correlation", random_state=21).fit_transform(cells)
#         cells_embedded = cells_embedded.T.tolist()
#         pickle.dump(cells_embedded, open(CHECKPOINT_TSNE_PATH, "wb"))
#     else:
#         cells_embedded = pickle.load(open(CHECKPOINT_TSNE_PATH, "rb"))
#     return cells_embedded


def visualize(cells, clusters_labels, title=None, centroids=None):
    """
    Visualize 2D representation.
    :param cells: embedded cells
    :param clusters_labels: list in number of cells length indicates each cell its cluster.
    :param title: plot title
    :param centroids: of the cluster algorithm
    :return:
    """
    X = cells[0]
    Y = cells[1]

    # plt.plot(X, Y, 'ro')

    colors = ['r', 'g', 'pink', 'c', 'm', 'y', 'b', 'yellow', 'orange', 'purple',  'lime', 'lavender', 'darkred', 'olive', 'k', 'w']
    for cluster, c in zip(np.unique(clusters_labels), colors):
        Xi = [X[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
        Yi = [Y[i] for i in range(len(clusters_labels)) if clusters_labels[i] == cluster]
        plt.plot(Xi, Yi, 'ro', color=c, label=cluster)
    # if centroids:
    #     for centroid in centroids:
    #         plt.plot(centroid[0], centroid[1], 'ro', color='y')
    # plt.plot([1, 2], [1, 4], 'ro')
    # plt.plot([3, 4], [9, 16], 'ro', color='green')
    plt.title(title)
    plt.legend()
    plt.show()



if __name__ == '__main__':


    step = 3


    PCs_path = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\TSNE_smartseq\pcs_scaled.pkl'
    TSNE_path = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\TSNE_smartseq\tsne_scaled.pkl'
    cells_filtered_path = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\TSNE_smartseq\cells_filtered.pkl'


    if step > 0:
        cells_filtered, clusters = pickle.load(open(cells_filtered_path, 'rb'))
    else:
        cells, gene_names, patients_information = extract_smart_seq_data_from_pickle(PICKLE_PATH)
        clusters = [p['general 11 cluster'] for p in patients_information]
        cells_filtered = cells[:, 6 < np.var(cells, axis=0)]
        pickle.dump((cells_filtered, clusters ), open(cells_filtered_path, 'wb'))

    if step > 1:
        PCs = pickle.load(open(PCs_path, 'rb'))
    else:
        PCs = PCA(n_components=10).fit_transform(preprocessing.scale(cells_filtered))
        # PCs = PCA(n_components=10).fit_transform(cells_filtered)
        pickle.dump((PCs), open(PCs_path, 'wb'))

    if step > 2:
        tsne_results = pickle.load(open(TSNE_path, 'rb'))
    else:
        tsne_results = tsne(PCs)
        pickle.dump((tsne_results), open(TSNE_path, 'wb'))

    visualize(tsne_results.T, clusters, 'correlation')