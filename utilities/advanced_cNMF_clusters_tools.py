"""
Data clustering/GEPs manager.
Here you can find all methods needed to extract clusters and GEPs information for downstream analysis.
"""
from utilities.package_importing import *

def get_GEP_dfs(res_melanoma_clinical_data,
                EXEC_DIR,
                RUN_RANGE,
                selected_K,
                number_of_genes,
                n_replicates,
                local_density_threshold,
                FEATURE,
                FEATURE_A,
                FEATURE_B):
    RUN_NAME = f'k{RUN_RANGE}_{number_of_genes}genes_{n_replicates}iter'
    USAGES_CONSENSUS_FILE = f'{RUN_NAME}.usages.k_{selected_K}.dt_{str(local_density_threshold).replace(".", "_")}.consensus.txt'

    usage_matrix = pd.read_csv(join(EXEC_DIR, RUN_NAME, USAGES_CONSENSUS_FILE), sep='\t', index_col=0)
    usage_matrix.columns = np.arange(1, selected_K + 1)
    normalized_usage_matrix = usage_matrix.div(usage_matrix.sum(axis=1), axis=0)
    samples = list(set([uu.split('_')[0] for uu in list(normalized_usage_matrix.index)]))
    df = normalized_usage_matrix.copy()
    df['sample'] = [uu.split('_')[0] for uu in list(df.index)]
    df['barcode'] = [uu.split('_')[1] for uu in list(df.index)]
    # keep only data of current patients in cohort
    GEP_all_samples_df = df.copy()
    res_melanoma_clinical_data = res_melanoma_clinical_data.set_index('Patient id')
    df = df[df['sample'].isin(res_melanoma_clinical_data.index.tolist())]

    df[FEATURE] = df['sample'].apply(lambda x: res_melanoma_clinical_data.loc[x][FEATURE])
    ##### Assign each cell one program based on the maximal usage value.
    high_prog = np.argmax(df[list(range(1, selected_K + 1))].values, axis=1) + 1
    df.loc[:, 'associated program'] = high_prog
    df = df.reset_index()[df.columns[-4:].tolist() + df.columns[:-4].tolist()]

    df_r = df[df[FEATURE] == FEATURE_A].reset_index().iloc[:,1:]
    df_nr = df[df[FEATURE] == FEATURE_B].reset_index().iloc[:,1:]
    return df, df_r, df_nr


def get_GEP_fraction_df(GEP_DF):
    n_cells_in_sample_df = GEP_DF.groupby(['sample']).count().reset_index()[['sample', 'barcode']].rename(columns={'barcode': 'count'}).set_index('sample')
    n_cells_in_GEP_df = GEP_DF.groupby(['sample', 'associated program']).count().reset_index()[['sample', 'associated program', 'barcode']].rename(columns={'associated program': 'program', 'barcode': 'count'})

    n_cells_in_GEP_df['n_cells_in_sample'] = n_cells_in_GEP_df['sample'].apply(lambda x: n_cells_in_sample_df.loc[x]['count'])
    n_cells_in_GEP_df['sample_fraction'] = n_cells_in_GEP_df['count'] / n_cells_in_GEP_df['n_cells_in_sample']
    GEP_fraction_df = n_cells_in_GEP_df[['sample', 'program', 'sample_fraction']]

    #### Add sample & cluster pairs which do not appear in fraction_df
    sample_program_pairs = [[sam, cl] for sam in GEP_DF['sample'].unique().tolist() for cl in GEP_DF['associated program'].unique().tolist()]
    zero_fraction_pairs_df = pd.DataFrame([[pair[0],pair[1],0] for pair in sample_program_pairs if not pair in GEP_fraction_df[['sample', 'program']].values.tolist()], columns=['sample', 'program', 'sample_fraction'])
    GEP_fraction_df = GEP_fraction_df.append(zero_fraction_pairs_df).sort_values(['sample', 'program']).reset_index().drop(columns=('index'))
    return GEP_fraction_df


def get_cluster_fraction_df(clusters_barcodes_mapping_df):
    n_cells_in_cluster_df = clusters_barcodes_mapping_df.groupby(['Sample', 'Cluster']).count().reset_index().rename(columns={'Sample':'sample', 'Cluster':'cluster', 'Barcode': 'n_cells_in_cluster'})
    n_cells_in_sample_df = clusters_barcodes_mapping_df.groupby(['Sample']).count().reset_index().rename(columns={'Sample':'sample', 'Cluster':'cluster', 'Barcode': 'n_cells_in_sample'}).iloc[:, :2].set_index('sample')
    n_cells_in_cluster_df['n_cells_in_sample'] = n_cells_in_cluster_df['sample'].apply(lambda x: n_cells_in_sample_df.loc[x]['n_cells_in_sample'])

    n_cells_in_cluster_df['sample_fraction'] = n_cells_in_cluster_df['n_cells_in_cluster'] / n_cells_in_cluster_df['n_cells_in_sample']
    cluster_fraction_df = n_cells_in_cluster_df[['sample', 'cluster', 'sample_fraction']]

    #### Add sample & cluster pairs which do not appear in fraction_df
    sample_cluster_pairs = [[sam, cl] for sam in clusters_barcodes_mapping_df['Sample'].unique().tolist() for cl in clusters_barcodes_mapping_df['Cluster'].unique().tolist()]
    zero_fraction_pairs_df = pd.DataFrame([[pair[0],pair[1],0]   for pair in sample_cluster_pairs if not pair in cluster_fraction_df[['sample', 'cluster']].values.tolist()], columns=['sample', 'cluster', 'sample_fraction'])
    cluster_fraction_df = cluster_fraction_df.append(zero_fraction_pairs_df).sort_values(['sample', 'cluster']).reset_index().drop(columns=('index'))

    return cluster_fraction_df


def find_contaminated_program_in_cNMF_output(GEP_DF, prog_idx_for_removal, output_path, associated_threshold=0.1):
    """
    Given GEP dataframe - cNMF output, and idx of some program found to be contaminated with immune cells. That's it we
    see gene markers of program are immune markers... we will save a csv of all cells associated to that program.
    :param GEP_DF:
    :param prog_idx_for_removal:
    :param output_path:
    :return:
    """
    reduced_GEP_DF = GEP_DF.iloc[:, [0, 1, 2, 3, 17]]
    print(f'{sum(reduced_GEP_DF[14] >= associated_threshold)} cells will be removed')
    contaminated_cells = reduced_GEP_DF[reduced_GEP_DF[prog_idx_for_removal] >= associated_threshold]
    # contaminated_cells.head()
    contaminated_cells.to_csv(output_path)


def filter_contaminated_cells_out_of_GEP_DF(GEP_DF, contaminated_prog_idx=None, contaminated_cells=None, contaminated_cells_path=None):
    """
    Use this function when you want to use GEP_DF (row cNMF matrix output) with the filter of contaminated program -
    Cells of contaminated program removed and contaminated program column taken out.
    :param GEP_DF: cNMF matrix output
    :param contaminated_prog_idx: program column will be removed from GEP_DF (not required).
    :param contaminated_cells: df containing only contaminated cells identifiers (barcodes and samples ids).
    :param contaminated_cells_path: you can use this attribute instead of contaminated_cells (it is the path of
    contaminated_cells df in pc).
    :return: filtered GEP df.
    NOTE: current csv of contaminated cells. use this version for GEP for subcohort_1.1.22:
    r'/storage/md_keren/shitay/Data/tables/GEP/subcohort_1.1.22/contaminated_cells_list/GEP14_contaminated_cells.csv'
    """
    if contaminated_cells_path:
        contaminated_cells = pd.read_csv(contaminated_cells_path)[['sample', 'barcode']]

    # take identifiers of both dfs:
    contaminated_cells_identifiers = [tuple(aa) for aa in contaminated_cells.values.tolist()]
    GEP_identifiers = [tuple(aa) for aa in GEP_DF.drop(columns=[14])[['sample', 'barcode']].values.tolist()]

    # contaminated_indexes stores indexes of contaminated_cells in GEP_df
    contaminated_indexes = [aa in contaminated_cells_identifiers for aa in GEP_identifiers]

    output_df = GEP_DF[~np.array(contaminated_indexes)]
    if contaminated_prog_idx:
        output_df = output_df.drop(columns=[contaminated_prog_idx])
    return output_df