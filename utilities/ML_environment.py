# import statsmodels as sm
import scipy.stats as stats
from scipy.stats import rankdata
from Bio.Cluster import kcluster
import pickle
import pandas as pd
import numpy as np
from utilities.general_helpers import flatten_list

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


def multipletests_fdr(p_vals, alpha=0.05):
    """
    Test results and p-value correction for multiple tests using Benjamini/Hochberg method.
    :param p_vals: np.array of pvals
    :return: pvals_corrected
    """

    ranked_p_values = rankdata(p_vals)
    pvals_corrected = p_vals * len(p_vals) / ranked_p_values
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected < alpha, pvals_corrected


def kmeans(Data_RNAseq, numer_of_clusters):
    """
    Perform Kmeans algorithm with pearson correlation as distance metric.
    :return: clusters association.
    """
    # dist 'c' is pearson correlation distance
    clusters, error, nfound = kcluster(Data_RNAseq.cells, nclusters=numer_of_clusters, dist='c')
    return clusters


def find_marker_genes_in_cluster(cluster_data, other_clusters_data, log_FC_threshold, pval_threshold):
    """
    After the clustering process has been done run this function to find marker genes for each cluster.
    The function conducts a Fisher Exact Test for every gene to check whether that gene constitutes a marker
    of one of the clusters. Obviously, we conduct a test per cluster. that means that if we have k clusters and n
    genes we will have k*n tests.
    We test whether the gene in the cluster expresses differently of the other clusters.
    we define expressed as value > 1 (after normalization). and we test the proportion between the number of cells
    expressing each gene compered to all other clusters together.
    because we conduct many statistical test we will correct the pvalues using Benjamini/Hochberg correction
    :param Data_RNAseq: Can be cohort object or singular RNAseq object which all its cells having association to one of the clusters.
    :param clusters_list: cluster in size of number of cells such that each place indicate which cluster corresponds
    to that cell (number of cluster).

    Note: Marker defined as gene with pval < pval_thresh and log_FC < log_ratio_threshold
    :return:
    gene names - of all gene markers, ordered by pval_corrected
    gene ids - of all gene markers, ordered by pval_corrected
    pval - of all gene markers, ordered by pval_corrected
    log ratios - of all gene markers, ordered by pval_corrected
    """

    # Part 1
    percentage_voting_expression_in_clusters = []
    p_values = []  # before correction

    # statistical test for each gene
    for gene_idx in range(cluster_data.number_of_genes):
        cluster_gene_expression = cluster_data.counts[:, gene_idx] > 1
        other_clusters_gene_expression = other_clusters_data.counts[:, gene_idx] > 1

        number_of_cells_expressing_gene_in_cluster = sum(cluster_gene_expression)
        number_of_cells_expressing_gene_in_other_clusters = sum(other_clusters_gene_expression)

        n_cells_in_cluster = len(cluster_gene_expression)
        n_cells_in_other_clusters = len(other_clusters_gene_expression)

        oddsratio, pvalue = stats.fisher_exact([[number_of_cells_expressing_gene_in_cluster,
                                                 number_of_cells_expressing_gene_in_other_clusters],
                                                [n_cells_in_cluster - number_of_cells_expressing_gene_in_cluster,
                                                 n_cells_in_other_clusters - number_of_cells_expressing_gene_in_other_clusters]])

        percentage_voting_expression_in_clusters.append([[number_of_cells_expressing_gene_in_cluster, number_of_cells_expressing_gene_in_other_clusters],
                                                         [n_cells_in_cluster, n_cells_in_other_clusters]])
        p_values.append(pvalue)

    reject_arr, pvals_corrected = multipletests_fdr(np.array(p_values), alpha=0.05)

    # Part 2
    percentage_voting_expression_in_clusters = np.array(percentage_voting_expression_in_clusters)
    percentage_voting_expression_in_clusters = percentage_voting_expression_in_clusters[:,
                                               0] / percentage_voting_expression_in_clusters[:, 1]

    df = pd.DataFrame(np.array([p_values, pvals_corrected]).T, columns=['pval', 'corected_pval'])
    df['log_FC'] = np.log2((np.mean(cluster_data.counts, axis=0) + 0.01) / (np.mean(other_clusters_data.counts, axis=0) + 0.01))
    df['cluster mean_expression'] = np.mean(cluster_data.counts, axis=0)
    df['other clusters mean expression'] = np.mean(other_clusters_data.counts, axis=0)
    df['n_expressing_cells__cls1 > n_expressing_cells__cls2'] = percentage_voting_expression_in_clusters[:,0] > percentage_voting_expression_in_clusters[:, 1]
    df['features'] = cluster_data.features
    df['gene names'] = cluster_data.gene_names

    # sort log_FC by size, and adjust features/genes and percentage_voting_expression_in_clusters by the log_FC.
    df = df.sort_values(['log_FC'], ascending=False)

    # take only genes with pval < pval_threshold
    df = df[df['corected_pval'] < pval_threshold]
    # take only genes with log_FC > log_FC_threshold
    df = df[df['log_FC'] > log_FC_threshold]
    df = df.reset_index(drop=True)

    return df



def find_markers_in_clusters(data_rna_seq, clusters_indices, log_ratio_threshold = 0, pval_threshold=0.05):
    """

    :param data_rna_seq:
    :param clusters_indices:
    :return:
    """
    cluster_markers_list = []
    for cluster_idx in range(len(clusters_indices)):
        current_cluster_indices = clusters_indices[cluster_idx]
        other_clusters_indices = [ii for ii in flatten_list(clusters_indices) if not ii in clusters_indices[cluster_idx]]
        cluster_markers = find_marker_genes_in_cluster(data_rna_seq[current_cluster_indices], data_rna_seq[other_clusters_indices], log_ratio_threshold, pval_threshold)
        cluster_markers_list.append({'cluster id': cluster_idx, 'markers': cluster_markers})
    return cluster_markers_list



def shows_statistics_in_clusters(Data_RNAseq, clusters_list):
    cluster_samples = {}
    for cluster_id in set(clusters_list):
        cluster_indices = [idx for idx, ii in enumerate(clusters_list) if ii == cluster_id]
        cluster_samples[cluster_id] = Data_RNAseq[cluster_indices]

