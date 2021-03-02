import statsmodels as sm
import scipy.stats as stats
from Bio.Cluster import kcluster
from sklearn.manifold import TSNE
import pickle


def TSNE_embedded(cells):
    """
    Activates t-SNE algorithm to find 2D representation of the cells.
    calculate and export if TSNE_IMPORT_EXPORT else import from CHECKPOINT_TSNE_PATH.
    :param cells: pkl format.
    :return: embedded 2D-representation .
    """
    pass
    # if TSNE_IMPORT_EXPORT:  # perform TSNE and save embedded vector in CHECKPOINT_TSNE_PATH
    #     cells_embedded = TSNE(n_components=2, random_state=21).fit_transform(cells)
    #     cells_embedded = cells_embedded.T.tolist()
    #     pickle.dump(cells_embedded, open(CHECKPOINT_TSNE_PATH, "wb"))
    # else:
    #     cells_embedded = pickle.load(open(CHECKPOINT_TSNE_PATH, "rb"))
    # return cells_embedded


def kmeans(Data_RNAseq, numer_of_clusters):
    """
    Perform Kmeans algorithm with pearson correlation as distance metric.
    :return: clusters association.
    """
    # dist 'c' is pearson correlation distance
    clusters, error, nfound = kcluster(Data_RNAseq.cells, nclusters=numer_of_clusters, dist='c')
    return clusters


def fing_marker_genes_in_clusters(Data_RNAseq, clusters_list):
    """
    After the clustering process has been done run this function to find marker genes for each cluster.
    The function conducts a Fisher Exact Test for every gene to check whether that gene constitutes a marker
    of one of the clusters. Obviously, we conduct a test per cluster. that means that if we have k clusters and n
    genes we will have k*n tests.
    We test the which the gene expresses in the cluster differently from the other clusters.
    we define expressed as value > 1 (after normalization). and we test the proportion between the number of cells
    expressing each gene compered all other clusters together.
    because we conduct many statistical test we will correct the pvalues using Benjamini/Hochberg correction
    :param Data_RNAseq: Can be cohort object or singular RNAseq object which all its cells having association to one of the clusters.
    :param clusters_list: cluster in size of number of cells such that each place indicate which cluster corresponds
    to that cell (number of cluster).
    :return:
    """

    percentage_voting_expression_in_clusters = []
    p_values = []  # before correction
    for cluster in list(set(clusters_list)):
        print(cluster)
        cluster_indices = [idx for idx, val in enumerate(clusters_list) if val == cluster]
        other_clusters_indices = [idx for idx, val in enumerate(clusters_list) if val != cluster]

        cluster_data = Data_RNAseq[cluster_indices]
        other_clusters_data = Data_RNAseq[other_clusters_indices]
        # statistical test for each gene
        for gene_idx in range(Data_RNAseq.number_of_genes):
            cluster_gene_expression = cluster_data.counts[:, gene_idx] > 1
            other_clusters_gene_expression = other_clusters_data.counts[:, gene_idx] > 1

            number_of_cells_expressing_gene_in_cluster = sum(cluster_gene_expression)
            number_of_cells_expressing_gene_in_other_clusters = sum(other_clusters_gene_expression)

            n_cells_in_cluster = len(cluster_gene_expression)
            n_cells_in_other_clusters = len(other_clusters_gene_expression)

            oddsratio, pvalue = stats.fisher_exact([[number_of_cells_expressing_gene_in_cluster,
                                                     number_of_cells_expressing_gene_in_other_clusters],
                                                    [n_cells_in_cluster, n_cells_in_other_clusters]])

            p_voting_expression_in_cluster = number_of_cells_expressing_gene_in_cluster / n_cells_in_cluster

            percentage_voting_expression_in_clusters.append(p_voting_expression_in_cluster)
            p_values.append(pvalue)

    sm.stats.multitest.multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)