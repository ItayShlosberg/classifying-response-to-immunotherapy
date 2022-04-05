"""
CLUSTERING and GEPS analysis AND differential gene expression analysis
This scripts replace all NBs of response vs. non-response and cutaneous vs. mucosal analysis of immune/tumor/CD8/myeloid
"""


lib = r'/srv01/technion/shitay/Code/classifying_response_to_immunotherapy/'
import sys
sys.path.append(lib)
from utilities.package_importing import *

"""
######################## Arguments ########################
# Specify clinical table and set output path
# clinical table can be 1|2|3
#  1 - NR mucosal vs NR cutaneous
#  2 - NR mucosal vs R cutaneous
#  3 - NR mucosal vs NR & R cutaneous
"""


ROOT_FOLDER_PATH = r'/storage/md_keren/shitay/outputs/response_analysis/subcohort_1.1.22/5.4.22'
EXPERIMENT_NUM = 3

######################## WHAT YOU WANT TO PERFORM ########################
CLUSTERING_ANALYSIS = False  # ranksum-test
GEP_ANALYSIS = False    # ranksum-test
DIFFERENTIAL_GENE_EXP_ANALYSIS = False   # markers and heatmaps
TSNE = False
COMBINING_PVAL_FILES = False
NATMI_ANALYSIS = True
GROUPS = ['immune', 'CD8', 'tumor', 'myeloid'] # immune/tumor/myeloid/CD8

######################## tSNE paths ########################
TSNE_PATH = fr'/storage/md_keren/shitay/outputs/TSNE/subcohort_1.1.22'

######################## gene expression analysis ########################
# COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/cohort/normalized/4.11.21/cohort_normalized_4.11.21_protein_coding_genes.pkl'
COHORT_PATH = r'/storage/md_keren/shitay/Data/droplet_seq/M97_M173/subcohort/normalized/16.3.22/subcohort_normalized_1.1.22_protein_coding_genes.pkl'
######################## Clustering paths ########################

IMMUNE_CLUSTER_BARCODE_MAPPING_PATH = r'/storage/md_keren/shitay/outputs/clustering/immune/summary/subcohort_1.1.22_run_1.1.22/subcohort_immune_1.1.22_clusters_mapping.csv'
CD8_CLUSTER_BARCODE_MAPPING_PATH = r'/storage/md_keren/shitay/outputs/clustering/CD8/summary/subcohort_1.1.22_run_1.1.22/subcohort_CD8_1.1.22_clusters_mapping.csv'
MYELOID_CLUSTER_BARCODE_MAPPING_PATH = r'/storage/md_keren/shitay/outputs/clustering/myeloid/summary/subcohort_1.1.22_run_1.1.22/subcohort_myeloid_1.1.22_clusters_mapping.csv'

######################## NATMI arguments ########################
NATMI_ROOT_PATH = r'/storage/md_keren/shitay/outputs/NATMI/tumor_CD8_myeloid/executions/21.3.22'
NATMI_FILE_NAME = 'Edges_lrc2p.csv'
AVG_SPECIFICITY_THRESHOLD = 0  # 0.01
TOTAL_SPECIFICITY_THRESHOLD = 0  # 0.02
DETECTION_THRESHOLD = 0.7
AVG_EXPRESSION_THRESHOLD = 1
FILTER_GEP_IN_NATNI_DF = False
NATMI_DF_SORT_BY = 'Edge average expression derived specificity'

######################## GEP Arguments ########################
EXEC_DIR = r'/storage/md_keren/shitay/outputs/cNMF/executions/tumor_runs/subcohort_1.1.22'
selected_K = 20
local_density_threshold = '0.10'
number_of_genes = 2000
n_replicates = 200
RUN_RANGE = '20'
N_PROG = 20
NO_PATHWAYS_GEPS = [1]

ACTIVITY_GEP_CELL_USAGE_THRESHOLD = 0.1
ACTIVITY_GEP_SAMPLE_PORTION_THRESHOLD = 0.1
# SUBCOHORT_ACTIVE_PROGRAMS = [2, 4, 6, 7, 8, 9, 10, 11, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30]
######################## Combining PVAL - annotations paths ########################
FILES_EXTENSION = 'csv'#

# annotation paths
immune_annotation_path = r'/storage/md_keren/shitay/Data/tables/clustering_annotations/subcohort_1.1.22/17.1.22/immune.xlsx'
myeloid_annotation_path = r'/storage/md_keren/shitay/Data/tables/clustering_annotations/subcohort_1.1.22/17.1.22/Myeloid.xlsx'
CD8_annotation_path = r'/storage/md_keren/shitay/Data/tables/clustering_annotations/subcohort_1.1.22/17.1.22/CD8.xlsx'
# GEP pathways path
GEP_pathways_path = r'/storage/md_keren/shitay/Data/tables/GEP/subcohort_1.1.22/k20_2000genes_200iter/GEPS_pathways.xlsx'
#######################################################################

if EXPERIMENT_NUM == 1:
    EXPERIMENT_NAME = r'NR_mucosal_vs_NR_cutaneous'
    FEATURE = 'Melanoma type' # 'response' | 'Melanoma type'
elif EXPERIMENT_NUM == 2:
    EXPERIMENT_NAME = r'NR_mucosal_vs_R_cutaneous'
    FEATURE = 'response'# 'response' | 'Melanoma type'
elif EXPERIMENT_NUM == 3:
    EXPERIMENT_NAME = r'NR_mucosal_vs_R_AND_NR_cutaneous'
    FEATURE = 'Melanoma type' # 'response' | 'Melanoma type'

FOLDER_PATH = join(ROOT_FOLDER_PATH, EXPERIMENT_NAME)
create_folder(ROOT_FOLDER_PATH)
create_folder(FOLDER_PATH)

if FEATURE == 'response':
    FEATURE_A = 'R'
    FEATURE_B = 'NR'
    FEATURE_G = 'not in use for now'
else:
    FEATURE_A = 'Cutaneous'
    FEATURE_B = 'Mucosal'
    FEATURE_G = 'other'


def combining_pvals_various_experiments():
    prefix_file_name = '_'.join(FOLDER_PATH.split('/')[-1].split('_')[:-1])
    OUTPUT_FILE_NAME = fr'{prefix_file_name}_pvals_all_statistical_test_results.xlsx'
    ########################################### file paths ###########################################
    GEP_path = f'statistical_tests/GEP_fraction_ranksum_test_pvals.{FILES_EXTENSION}'
    immune_clustering_path = f'statistical_tests/cluster_fractions_ranksum_test_pvals.{FILES_EXTENSION}'
    myeloid_clustering_path = f'myeloids/statistical_tests/cluster_fractions_ranksum_test_pvals.{FILES_EXTENSION}'
    CD8_clustering_path = f'cytotoxic_t_cells/statistical_tests/cluster_fractions_ranksum_test_pvals.{FILES_EXTENSION}'

    if FILES_EXTENSION == 'csv':
        read_func = pd.read_csv
    else:
        read_func = pd.read_excel
    ########################################### dfs ###########################################
    GEP_df = read_func(join(FOLDER_PATH, GEP_path))
    immune_clustering_df = read_func(join(FOLDER_PATH, immune_clustering_path))
    myeloid_clustering_df = read_func(join(FOLDER_PATH, myeloid_clustering_path))
    CD8_clustering_df = read_func(join(FOLDER_PATH, CD8_clustering_path))

    ########################################### annotation dfs ###########################################
    immune_annotation_df = pd.read_excel(join(FOLDER_PATH, immune_annotation_path), engine='openpyxl')
    myeloid_annotation_df = pd.read_excel(join(FOLDER_PATH, myeloid_annotation_path), engine='openpyxl')
    CD8_annotation_df = pd.read_excel(join(FOLDER_PATH, CD8_annotation_path), engine='openpyxl')
    # GEP_pathways_df = pd.read_excel(GEP_pathways_path)

    ########################################### Add annotation to pval dfs ###########################################
    immune_clustering_df['cluster annotation'] = immune_clustering_df['cluster'].apply(
        lambda x: immune_annotation_df.set_index('CLUSTER_IDX').loc[x]['annotations'])
    myeloid_clustering_df['cluster annotation'] = myeloid_clustering_df['cluster'].apply(
        lambda x: myeloid_annotation_df.set_index('cluster').loc[x]['annotation'])
    CD8_clustering_df['cluster annotation'] = CD8_clustering_df['cluster'].apply(
        lambda x: CD8_annotation_df.set_index('cluster').loc[x]['annotation'])

    ########################################### order pval dfs ###########################################
    immune_clustering_df = immune_clustering_df[
        ['cluster', 'cluster annotation', 'pval', 'corrected_pval', 'Median >']].sort_values('pval')
    myeloid_clustering_df = myeloid_clustering_df[
        ['cluster', 'cluster annotation', 'pval', 'corrected_pval', 'Median >']].sort_values('pval')
    CD8_clustering_df = CD8_clustering_df[['cluster', 'cluster annotation', 'pval', 'corrected_pval', 'Median >']].sort_values(
        'pval')
    GEP_df = GEP_df.sort_values('pval')

    ########################################### round pval ###########################################
    pval_round = 4
    immune_clustering_df.pval = round(immune_clustering_df.pval, pval_round)
    immune_clustering_df.corrected_pval = round(immune_clustering_df.corrected_pval, pval_round)

    myeloid_clustering_df.pval = round(myeloid_clustering_df.pval, pval_round)
    myeloid_clustering_df.corrected_pval = round(myeloid_clustering_df.corrected_pval, pval_round)

    CD8_clustering_df.pval = round(CD8_clustering_df.pval, pval_round)
    CD8_clustering_df.corrected_pval = round(CD8_clustering_df.corrected_pval, pval_round)

    GEP_df.pval = round(GEP_df.pval, pval_round)
    GEP_df.corrected_pval = round(GEP_df.corrected_pval, pval_round)

    ####### builds an Excel that will contain all statistical test results in one sheet:
    GEP_df.name = "GEP fraction ranksum test (sorted by pval)"
    immune_clustering_df.name = "Immune clustering fraction ranksum test (sorted by pval)"
    myeloid_clustering_df.name = "myeloid clustering fraction ranksum test (sorted by pval)"
    CD8_clustering_df.name = "CD8 clustering fraction ranksum test (sorted by pval)"

    ########################################### ExcelWriter ###########################################

    writer = pd.ExcelWriter(join(FOLDER_PATH, OUTPUT_FILE_NAME), engine='xlsxwriter')
    workbook = writer.book
    worksheet = workbook.add_worksheet('Result')
    writer.sheets['Result'] = worksheet

    ########################################### write header ###########################################

    worksheet.write_string(0, 0, immune_clustering_df.name)
    worksheet.write_string(immune_clustering_df.shape[0] + 4, 0, myeloid_clustering_df.name)

    worksheet.write_string(0, immune_clustering_df.shape[1] + 1, CD8_clustering_df.name)
    worksheet.write_string(immune_clustering_df.shape[0] + 4, immune_clustering_df.shape[1] + 1, GEP_df.name)

    ########################################### put dfs in xlsx ###########################################
    immune_clustering_df.to_excel(writer, sheet_name='Result', startrow=1, startcol=0, index=False)
    col_fixed(immune_clustering_df, worksheet)

    myeloid_clustering_df.to_excel(writer, sheet_name='Result', startrow=immune_clustering_df.shape[0] + 5, startcol=0,
                                   index=False)

    CD8_clustering_df.to_excel(writer, sheet_name='Result', startrow=1, startcol=immune_clustering_df.shape[1] + 1,
                               index=False)
    col_fixed(CD8_clustering_df, worksheet, col_skip=immune_clustering_df.shape[1] + 1)

    GEP_df.to_excel(writer, sheet_name='Result', startrow=immune_clustering_df.shape[0] + 5,
                    startcol=immune_clustering_df.shape[1] + 1, index=False)

    ########################################### add GEP pathways ###########################################
    pathways_col_skip = immune_clustering_df.shape[1] + GEP_df.shape[1] + 3
    n_p_in_row = 4
    max_n_pathways = 8
    for i in range(1, N_PROG+1):
        if i in NO_PATHWAYS_GEPS:
            pathways = []
        else:
            row_pathways_df = pd.read_excel(GEP_pathways_path, f'GEP{i}', engine='openpyxl')
            pathways = [v[9:].replace('_', ' ').lower() for v in row_pathways_df['Gene Set Name'].values.tolist()][
                       :max_n_pathways]

        pr_pathways_df = pd.DataFrame(pathways, columns=[f'Pr{i}'])
        col_p = (i - 1) % n_p_in_row
        row_p = int((i - 1) / n_p_in_row)
        pr_pathways_df.to_excel(writer, sheet_name='Result', startrow=1 + row_p * (max_n_pathways + 2),
                                startcol=pathways_col_skip + (col_p * 2) + 1, index=False)

    ########################################### save ###########################################
    writer.save()


def print_pvals(pvals, n_groups):
    plt.figure(figsize=(8, 6), dpi=80);
    min_pval = 0.05

    small_pvals = pvals[pvals[:, 1] < min_pval]
    big_pvals = pvals[pvals[:, 1] >= min_pval]

    plt.scatter(x=big_pvals[:, 0], y=big_pvals[:, 1]);
    plt.scatter(x=small_pvals[:, 0], y=small_pvals[:, 1]);

    ax = plt.plot([0, n_groups + 1], [min_pval, min_pval], color='y');
    plt.xticks(np.arange(1, n_groups + 1));
    plt.yticks([min_pval, 1]);
    plt.xlim((0, n_groups + 1));
    plt.title('p-values');

    for coord in pvals:
        plt.text(coord[0], coord[1], '{}'.format(int(coord[0])));


def col_fixed(df, wsheet, col_skip=0):
    """
    Helper function of combining pvals files
    :param df:
    :param wsheet:
    :param col_skip:
    :return:
    """
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
            len(str(series.name))  # len of column name/header
            )) + 1  # adding a little extra space
        wsheet.set_column(idx+col_skip, idx+col_skip, max_len)  # set column width


# def get_coorected_label(x, p_value_dic):
#     pvalue = p_value_dic[x]
#     return str(x) + '\n' + ''.join(['*' for sig in [0.05, 0.005, 0.0005] if pvalue < sig])


def get_corrected_label(x, corrected_p_value_dic, pvalue_dic):
    return str(x) + f'\n{round(pvalue_dic[x], 3)}\n{round(corrected_p_value_dic[x], 3)}\n' + ''.join(
        ['*' for sig in [0.05, 0.005, 0.0005] if corrected_p_value_dic[x] < sig])


def save_pvals(path, pvals, corrected_pvals, GROUP_NAME, comparison_directions, activity_programs_idxs=None):
    pval_df = pd.DataFrame(np.array([pvals[:, 0], pvals[:, 1], corrected_pvals[:, 1]]).T,
                           columns=[GROUP_NAME, 'pval', 'corrected_pval'])  # .sort_values('pval')
    pval_df[GROUP_NAME] = pval_df[GROUP_NAME].astype(int)
    pval_df['Median >'] = list(comparison_directions)
    if GROUP_NAME=='program':
        pval_df['is_activity'] = pval_df[GROUP_NAME].apply(lambda x: 'V' if x in activity_programs_idxs else '')
    pval_df.to_csv(path, index=False)


def create_fraction_df_in_clusters(clusters_barcodes_mapping_df, response_samples, non_response_samples):
    clusters = clusters_barcodes_mapping_df['Cluster'].unique()
    sample_fractions_df = clusters_barcodes_mapping_df.groupby(
        ['Sample', 'Cluster']).count().reset_index()  # .drop(columns=['index'])
    sample_counts_df = clusters_barcodes_mapping_df.groupby(['Sample']).count().reset_index().drop(
        columns=['Cluster'])  # .drop(columns=['index'])
    sample_fractions_df['counts'] = sample_fractions_df.apply(
        lambda x: sample_counts_df.set_index('Sample').loc[x['Sample']][0], axis=1)
    sample_fractions_df['fraction'] = sample_fractions_df['Barcode'] / sample_fractions_df['counts']
    sample_fractions_df[FEATURE] = sample_fractions_df['Sample'].apply(
        lambda x: res_melanoma_clinical_data.loc[x][FEATURE])

    # Add zero fraction samples to DF:
    sample_fractions_df = sample_fractions_df[["Sample", 'Cluster', 'fraction', FEATURE]]
    full_sample_fractions_df = pd.DataFrame(
        sample_fractions_df[sample_fractions_df[FEATURE] != FEATURE_G])
    samples_in_clusters = full_sample_fractions_df[["Sample", 'Cluster']].values.tolist()
    samples = response_samples + non_response_samples
    dic_response = {ss: FEATURE_A for ss in response_samples}
    dic_response.update({ss: FEATURE_B for ss in non_response_samples})
    all_pairs = [[s, cl] for cl in clusters for s in samples]
    pairs_need_to_add = [pair for pair in all_pairs if not pair in samples_in_clusters]
    for sample, cluster in pairs_need_to_add:
        full_sample_fractions_df = full_sample_fractions_df.append(
            pd.DataFrame([[sample, cluster, 0, dic_response[sample]]], columns=full_sample_fractions_df.columns))
    return full_sample_fractions_df


def conducting_statistical_tests(full_sample_fractions_df):
    ### conduct statistical test
    pvals = []
    for cluster in set(full_sample_fractions_df["Cluster"]):
        cluster_df = full_sample_fractions_df[full_sample_fractions_df["Cluster"] == cluster]
        R_fractions = cluster_df[cluster_df[FEATURE] == FEATURE_A]['fraction'].values
        NR_fractions = cluster_df[cluster_df[FEATURE] == FEATURE_B]['fraction'].values
        res = ranksums(R_fractions, NR_fractions)[1]
        pvals.append([cluster, res])
    pvals = np.array(pvals)

    # Correct pvals
    corrected_pvals = np.array([pvals[:, 0], multipletests_fdr(pvals[:, 1])[1]]).T

    return corrected_pvals, pvals


def seaborn_presentation(full_sample_fractions_df, corrected_pvals, pvals, output_path, group):
    ### seaborn plot
    fig = plt.figure(figsize=(6, 3), dpi=80);
    sns.reset_orig();
    full_sample_fractions_df['cluster'] = full_sample_fractions_df['Cluster'].apply(
        lambda x: get_corrected_label(x, dict(corrected_pvals), dict(pvals)))
    full_sample_fractions_df = full_sample_fractions_df.sort_values('Cluster')
    ax = sns.catplot(kind="box", x='cluster', y='fraction', hue=FEATURE, data=full_sample_fractions_df,
                     palette={FEATURE_B: 'r', FEATURE_A: 'b'},
                     height=6, aspect=2.0).set(title=f"{group} cluster fraction distribution");
    sns.stripplot(x='cluster', y='fraction', hue=FEATURE, data=full_sample_fractions_df,
                  jitter=True, dodge=True, linewidth=0.5, palette={FEATURE_B: 'r', FEATURE_A: 'b'});
    # print_pvals(pvals, len(set(sample_fractions_df["Cluster"])))
    # pvals

    ### save plot
    comparison_directions = get_direction_of_ranksum_test(full_sample_fractions_df, 'Cluster')
    save_pvals(join(output_path, 'cluster_fractions_ranksum_test_pvals.csv'), pvals, corrected_pvals,
               'cluster', comparison_directions)
    ax.savefig(join(output_path, f'{group}_cluster_fractions_ranksum_test.png'), dpi=150)


def load_clinical_table():
    ##### Loads clinical table to get labels - R/NR
    melanoma_clinical_data = get_constant_cohort(EXPERIMENT_NUM) #get_clinical_data(71)#ICI=True, after_biopsy='ICI')
    res_melanoma_clinical_data = melanoma_clinical_data.set_index('Patient id')
    print(melanoma_clinical_data.head(20))
    print(f'Num of samples in table is {len(melanoma_clinical_data)}')
    print(res_melanoma_clinical_data[FEATURE].value_counts())
    return res_melanoma_clinical_data, melanoma_clinical_data


def create_sample_lists(melanoma_clinical_data):
    response_samples = melanoma_clinical_data[melanoma_clinical_data[FEATURE] == FEATURE_A]['Patient id'].tolist()
    non_response_samples = melanoma_clinical_data[melanoma_clinical_data[FEATURE] == FEATURE_B]['Patient id'].tolist()
    no_used_samples = melanoma_clinical_data[melanoma_clinical_data[FEATURE] == FEATURE_G]['Patient id'].tolist()

    print(f'number {FEATURE_A} samples: {len(response_samples)}')
    print(f'number {FEATURE_B} samples: {len(non_response_samples)}')
    print(f'number no label samples: {len(no_used_samples)}')
    return response_samples, non_response_samples, no_used_samples


def load_clustering_df(clustering_barcode_path, melanoma_clinical_data):
    clusters_barcodes_mapping_df = pd.read_csv(clustering_barcode_path)
    clusters_barcodes_mapping_df = clusters_barcodes_mapping_df[
        clusters_barcodes_mapping_df["Sample"].isin(melanoma_clinical_data['Patient id'])]
    return clusters_barcodes_mapping_df


def get_direction_of_ranksum_test(full_sample_fractions_df, col_feature):
    directionality_df = full_sample_fractions_df.groupby([col_feature, FEATURE]).median().reset_index().\
                            pivot_table(index=[col_feature], values='fraction', columns=[FEATURE])[[FEATURE_A, FEATURE_B]]
    directions = np.array([FEATURE_A, FEATURE_B])[directionality_df.values.argmax(1)]

    return directions


def clustering_conduct_fraction_test_and_save_output(clustering_barcode_path, melanoma_clinical_data,
                                                     response_samples, non_response_samples, output_path, group):
    # For each cluster, see the number of cells correspond to response patients vs. number of cells corespond to non-response patients.
    clusters_barcodes_mapping_df = load_clustering_df(clustering_barcode_path, melanoma_clinical_data)
    # patient fraction over clusters - ranksum test
    full_sample_fractions_df = create_fraction_df_in_clusters(clusters_barcodes_mapping_df, response_samples, non_response_samples)

    corrected_pvals, pvals = conducting_statistical_tests(full_sample_fractions_df)

    create_folder(output_path)
    seaborn_presentation(full_sample_fractions_df, corrected_pvals, pvals, output_path, group)


# def get_GEP_dfs(res_melanoma_clinical_data):
#     RUN_NAME = f'k{RUN_RANGE}_{number_of_genes}genes_{n_replicates}iter'
#     USAGES_CONSENSUS_FILE = f'{RUN_NAME}.usages.k_{selected_K}.dt_{str(local_density_threshold).replace(".", "_")}.consensus.txt'
#
#     usage_matrix = pd.read_csv(join(EXEC_DIR, RUN_NAME, USAGES_CONSENSUS_FILE), sep='\t', index_col=0)
#     usage_matrix.columns = np.arange(1, selected_K + 1)
#     normalized_usage_matrix = usage_matrix.div(usage_matrix.sum(axis=1), axis=0)
#     samples = list(set([uu.split('_')[0] for uu in list(normalized_usage_matrix.index)]))
#     df = normalized_usage_matrix.copy()
#     df['sample'] = [uu.split('_')[0] for uu in list(df.index)]
#     df['barcode'] = [uu.split('_')[1] for uu in list(df.index)]
#     # keep only data of current patients in cohort
#     GEP_all_samples_df = df.copy()
#     res_melanoma_clinical_data = res_melanoma_clinical_data.set_index('Patient id')
#     df = df[df['sample'].isin(res_melanoma_clinical_data.index.tolist())]
#
#     df[FEATURE] = df['sample'].apply(lambda x: res_melanoma_clinical_data.loc[x][FEATURE])
#
#
#     df_r = df[df[FEATURE] == FEATURE_A]
#     df_nr = df[df[FEATURE] == FEATURE_B]
#     ##### Assign each cell one program based on the maximal usage value.
#     r_high_prog = np.argmax(df_r[list(range(1, N_PROG + 1))].values, axis=1) + 1
#     nr_high_prog = np.argmax(df_nr[list(range(1, N_PROG + 1))].values, axis=1) + 1
#
#     df_r.loc[:, 'associated program'] = r_high_prog
#     df_nr.loc[:, 'associated program'] = nr_high_prog
#     return df, df_r, df_nr


def GEP_create_fraction_df(df_r, df_nr):
    r_associated_program_count_df = pd.DataFrame(df_r.groupby('sample')['associated program'].value_counts()).rename(
        columns={'associated program': 'count'}).reset_index()
    r_samples = r_associated_program_count_df['sample'].unique()

    samples_dic = {}
    for sample in r_samples:
        sample_associated_program_df = r_associated_program_count_df[r_associated_program_count_df['sample'] == sample]
        sample_n_barcodes = sum(sample_associated_program_df['count'])
        count_vector = np.zeros(selected_K + 1)
        for idx, row in sample_associated_program_df.iterrows():
            count_vector[row['associated program']] = row['count'] / sample_n_barcodes
        samples_dic[sample] = count_vector

    programs_r_patients_usage = {p: [samples_dic[sample][p] for sample in r_samples] for p in range(1, selected_K + 1)}

    nr_associated_program_count_df = pd.DataFrame(df_nr.groupby('sample')['associated program'].value_counts()).rename(
        columns={'associated program': 'count'}).reset_index()
    nr_samples = nr_associated_program_count_df['sample'].unique()

    samples_dic = {}
    for sample in nr_samples:
        #     r_associated_program_count_df[r_associated_program_count_df['sample'==sample]]
        sample_associated_program_df = nr_associated_program_count_df[
            nr_associated_program_count_df['sample'] == sample]
        sample_n_barcodes = sum(sample_associated_program_df['count'])
        count_vector = np.zeros(selected_K + 1)
        for idx, row in sample_associated_program_df.iterrows():
            count_vector[row['associated program']] = row['count'] / sample_n_barcodes
        samples_dic[sample] = count_vector

    programs_nr_patients_usage = {p: [samples_dic[sample][p] for sample in nr_samples] for p in
                                  range(1, selected_K + 1)}

    pvals = np.zeros([selected_K, 2])
    for i in range(1, selected_K + 1):
        pvals[i - 1] = ranksums(programs_nr_patients_usage[i], programs_r_patients_usage[i])[1]
        pvals[i - 1, 0] = i

    # Correct pvals
    print(f'With corrected pvals!!')
    corrected_pvals = np.array([pvals[:, 0], multipletests_fdr(pvals[:, 1])[1]]).T

    programs_nr_patients_df = pd.DataFrame(programs_nr_patients_usage)
    programs_nr_patients_df[FEATURE] = FEATURE_B
    programs_r_patients_df = pd.DataFrame(programs_r_patients_usage)
    programs_r_patients_df[FEATURE] = FEATURE_A

    fraction_df = programs_nr_patients_df.melt(id_vars=[FEATURE], var_name="program", value_name="fraction").append(
        programs_r_patients_df.melt(id_vars=[FEATURE], var_name="program", value_name="fraction"))

    return fraction_df, corrected_pvals, pvals


def GEP_seaborn_dispaly(fraction_df, corrected_pvals, pvals, output_path, activity_programs_idxs):
    sns.reset_orig()

    # def get_pval_asterisks(x, p_value_dic):
    #     return str(x) + '\n' + ''.join(['*' for sig in [0.05, 0.005, 0.0005] if p_value_dic[x] < sig])

    fraction_df['Program'] = fraction_df['program'].apply(lambda x: get_corrected_label(x, dict(corrected_pvals), dict(pvals)))
    ax = sns.catplot(kind="box", x='Program', y='fraction', hue=FEATURE, data=fraction_df,
                     palette={FEATURE_B: 'r', FEATURE_A: 'b'}, height=10, aspect=2.9).set(
        title="Programs fraction distribution")
    sns.stripplot(x='Program', y='fraction', hue=FEATURE, data=fraction_df,
                  jitter=True, dodge=True, linewidth=0.5, palette={FEATURE_B: 'r', FEATURE_A: 'b'});
    # print_pvals(corrected_pvals, N_PROG)

    comparison_directions = get_direction_of_ranksum_test(fraction_df, 'program')

    create_folder(output_path)
    save_pvals(join(output_path, 'GEP_fraction_ranksum_test_pvals.csv'), pvals, corrected_pvals, 'program',
               comparison_directions, activity_programs_idxs)
    ax.savefig(join(output_path, 'GEP_fraction_ranksum_test.png'), dpi=150)


def get_activity_GEP_idxs(GEP_all_samples_df):
    # Shows which programs are activity programs, displays the number of samples associated with each program.
    activity_programs_df = GEP_all_samples_df.copy()

    activity_programs_df.iloc[:, :N_PROG] = activity_programs_df.iloc[:, :N_PROG] > ACTIVITY_GEP_CELL_USAGE_THRESHOLD

    sample_counts_df = activity_programs_df.groupby('sample').sum().reset_index().set_index('sample')
    samples_num_barcodes = pd.DataFrame(activity_programs_df.groupby('sample')['barcode'].count()).rename(
        columns={'barcode': 'count'})  # .agg({'barcode': ['count']}).reset_index()
    sample_fraction_df = sample_counts_df.div(samples_num_barcodes.loc[sample_counts_df.index]['count'], axis=0)

    # sample_fraction_df
    num_of_program_in_sample = sample_fraction_df > ACTIVITY_GEP_SAMPLE_PORTION_THRESHOLD

    activity_programs_idx = np.array(num_of_program_in_sample.columns)[num_of_program_in_sample.sum() > 1].astype(int)

    print(f'Number of activity programs (more than 2 samples are associated): {len(activity_programs_idx)}')
    print(f'Activity programs (more than 2 samples are associated): {activity_programs_idx}')
    
    return activity_programs_idx #[np.isin(activity_programs_idx, SUBCOHORT_ACTIVE_PROGRAMS)]


def get_cohort_sample_cells(cohort, response_samples, non_response_samples):
    response_samples_indices = [s in response_samples for s in cohort.samples]
    non_response_samples_indices = [s in non_response_samples for s in cohort.samples]
    no_used_samples_indices = [s in no_used_samples for s in cohort.samples]

    print(f'number of cells in response_samples_indices: {sum(response_samples_indices)}')
    print(f'number of cells in non_response_samples_indices: {sum(non_response_samples_indices)}')
    print(f'number of cells in no_used_samples_indices: {sum(no_used_samples_indices)}')

    response_samples_cells = cohort[response_samples_indices]
    non_response_samples_cells = cohort[non_response_samples_indices]

    return response_samples_cells, non_response_samples_cells


def get_cohort_markers(response_samples_cells, non_response_samples_cells):
    #### Conduct fisher's exact test to find markers of response in cohort cells:
    response_markers = find_satisfying_list_of_markers_in_cluster(response_samples_cells, non_response_samples_cells,
                                                    log_FC_threshold=0.5, pval_threshold=0.05, min_pct=0.1,
                                                    min_diff_pct=0.1, min_markers=50)
    #### Conduct fisher's exact test to find markers of non-response in cohort cells:
    non_response_markers = find_satisfying_list_of_markers_in_cluster(non_response_samples_cells, response_samples_cells,
                                                        log_FC_threshold=0.5, pval_threshold=0.05, min_pct=0.1,
                                                        min_diff_pct=0.1, min_markers=50)
    return response_markers, non_response_markers


def print_heatmap(cohort, response_markers, non_response_markers,
                response_samples_cells, non_response_samples_cells, sub_folder, group_name):
        number_of_markers = 50
        # arrange cells indices by mean expression
        response_markers_indices = np.array(
            [cohort.features.index(f) for f in response_markers.iloc[:number_of_markers]['features']])
        response_cells_indices = np.flip(
            np.argsort(response_samples_cells.counts[:, response_markers_indices].mean(axis=1)))

        non_response_markers_indices = [cohort.features.index(f) for f in
                                        non_response_markers.iloc[:number_of_markers]['features']]
        non_response_cells_indices = np.flip(
            np.argsort(non_response_samples_cells.counts[:, non_response_markers_indices].mean(axis=1)))

        genes_indices = [cohort.features.index(f) for f in response_markers.iloc[:number_of_markers]['features']] + [
            cohort.features.index(f) for f in non_response_markers.iloc[:number_of_markers]['features']]
        arr_heatmap = np.concatenate([response_samples_cells.counts[response_cells_indices][:, genes_indices],
                                      non_response_samples_cells.counts[non_response_cells_indices][:,
                                      genes_indices]])
        heatmap = np.zeros_like(arr_heatmap)
        heatmap[arr_heatmap > 1] = 1
        arr_heatmap = scipy.stats.zscore(arr_heatmap, axis=0, ddof=1)

        cmap = pickle.load(open(r'/storage/md_keren/shitay/outputs/clustering/immune/heatmap/colorbar.pkl', 'rb'))
        fig, ax = plt.subplots(1)
        # fig.set_size_inches(10, 5)
        fig.set_size_inches(32, 16)

        sb_out = sb.heatmap(arr_heatmap.T, vmin=-1, vmax=1, cmap=cmap);
        cbar = sb_out.collections[0].colorbar
        cbar.set_ticks([-1, 1])
        cbar.set_ticklabels([-1, 1])

        sb_out.set_xticks([int(response_samples_cells.number_of_cells / 2),
                           response_samples_cells.number_of_cells + int(
                               non_response_samples_cells.number_of_cells / 2)])  # <--- set the ticks first
        sb_out.set_xticklabels([FEATURE_A, FEATURE_B], rotation='horizontal')

        # sb_out.set(xticklabels=[])
        sb_out.set(yticklabels=[])
        sb_out.tick_params(bottom=False, left=False)

        ax.axhline(number_of_markers, color='white')
        ax.axvline(response_samples_cells.number_of_cells, color='white')
        ax.set_title(f"{group_name} cells markers - {FEATURE_A}/{FEATURE_B}");
        # ax.set_xlabel('response            non-response');
        ax.set_ylabel('Markers');
        file_name = fr'{group_name}_cells_marker_heatmap.png'

        create_folder(join(FOLDER_PATH, sub_folder))
        fig.savefig(join(FOLDER_PATH, sub_folder, file_name))


def differential_gene_exp_analysis(response_samples, non_response_samples):
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    for group in GROUPS:
        sub_folder = ''
        print(f'differential gene expression analysis on {group} cells')
        if group == 'immune':
            filtered_cohort = cohort.filter_cells_by_property('is_immune', True)
        elif group == 'tumor':
            filtered_cohort = cohort.filter_cells_by_property('is_cancer', True)
        elif group == 'CD8':
            filtered_cohort = get_requested_subset(cohort, 'CYTOTOXIC_T_CELLS')
            sub_folder = 'cytotoxic_t_cells'
        elif group == 'myeloid':
            filtered_cohort = get_requested_subset(cohort, 'MYELOIDS')
            sub_folder = 'myeloids'

        one_group_differential_gene_exp_analysis(filtered_cohort, response_samples, non_response_samples, sub_folder, group)


def one_group_differential_gene_exp_analysis(cohort, response_samples, non_response_samples, sub_folder, group_name):
    # differential expression analysis
    response_samples_cells, non_response_samples_cells = get_cohort_sample_cells(cohort, response_samples, non_response_samples)

    # find and save markers
    response_markers, non_response_markers = get_cohort_markers(response_samples_cells, non_response_samples_cells)
    response_markers.to_csv(join(FOLDER_PATH, sub_folder, fr'{FEATURE_A}_{group_name}_markers.csv'), index=False)
    non_response_markers.to_csv(join(FOLDER_PATH, sub_folder, fr'{FEATURE_B}_{group_name}_markers.csv'), index=False)

    # Builds heatmap:
    print_heatmap(cohort, response_markers, non_response_markers, response_samples_cells,
                  non_response_samples_cells, sub_folder, group_name)


def run_GEPS_analysis():
    print('GEP analysis')
    GEP_df, GEP_df_r, GEP_df_nr = get_GEP_dfs(res_melanoma_clinical_data)
    activity_programs_idxs = get_activity_GEP_idxs(GEP_df)
    GEP_fraction_df, corrected_pvals, pvals = GEP_create_fraction_df(GEP_df_r, GEP_df_nr)
    GEP_seaborn_dispaly(GEP_fraction_df, corrected_pvals, pvals, join(FOLDER_PATH, 'statistical_tests'),
                        activity_programs_idxs)


def run_clustering_analysis():
    print('Clustering analysis')
    clustering_conduct_fraction_test_and_save_output(IMMUNE_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data,
                                                     response_samples, non_response_samples,
                                                     join(FOLDER_PATH, 'statistical_tests'), 'immune')

    clustering_conduct_fraction_test_and_save_output(CD8_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data,
                                                     response_samples, non_response_samples,
                                                     join(FOLDER_PATH, 'cytotoxic_t_cells', 'statistical_tests'),
                                                     'CD8')

    clustering_conduct_fraction_test_and_save_output(MYELOID_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data,
                                                     response_samples, non_response_samples,
                                                     join(FOLDER_PATH, 'myeloids', 'statistical_tests'),
                                                     'myeloid')


def print_tSNEs(response_samples, non_response_samples):

    one_group_print_tSNE(response_samples, non_response_samples, '', 'cohort')
    for group in GROUPS:
        print(f'Create a tSNE plot of {group}')
        sub_folder = ''
        if group == 'CD8':
            sub_folder = 'cytotoxic_t_cells'
        elif group == 'myeloid':
            sub_folder = 'myeloids'

        one_group_print_tSNE(response_samples, non_response_samples, sub_folder, group)


def one_group_print_tSNE(response_samples, non_response_samples, sub_folder, group):

    tsne_path = join(TSNE_PATH, group, f'{group}_bhtsne_1.1.22.csv')
    subcohort_tsne_df = pd.read_csv(tsne_path)

    R_tsne_df = subcohort_tsne_df[subcohort_tsne_df.Sample.isin(response_samples)]
    NR_tsne_df = subcohort_tsne_df[subcohort_tsne_df.Sample.isin(non_response_samples)]
    other_samples_tsne_df = subcohort_tsne_df[~subcohort_tsne_df.Sample.isin(res_melanoma_clinical_data.index.tolist())]


    plt.figure(figsize=(16, 8))
    plt.plot(R_tsne_df[['x', 'y']].values[:, 0], R_tsne_df[['x', 'y']].values[:, 1], 'ro', color='b', label=FEATURE_A,
             markersize=1);
    plt.plot(NR_tsne_df[['x', 'y']].values[:, 0], NR_tsne_df[['x', 'y']].values[:, 1], 'ro', label=FEATURE_B,
             markersize=1);
    plt.plot(other_samples_tsne_df[['x', 'y']].values[:, 0], other_samples_tsne_df[['x', 'y']].values[:, 1], 'ro',
             color='gray', label='-', markersize=1);
    plt.legend();
    plt.title(f"{group} tSNE");
    file_name = f'{group}_tSNE.png'
    plt.savefig(join(FOLDER_PATH, sub_folder, file_name))
    print(f'Saved at {join(FOLDER_PATH, sub_folder, file_name)}')


def extract_all_clusters(melanoma_clinical_data):
    ###### extract clusters & GEP
    immune_clusters_barcodes_mapping_df = load_clustering_df(IMMUNE_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data)
    CD8_clusters_barcodes_mapping_df = load_clustering_df(CD8_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data)
    myeloid_clusters_barcodes_mapping_df = load_clustering_df(MYELOID_CLUSTER_BARCODE_MAPPING_PATH, melanoma_clinical_data)
    GEP_df, R_GEP_df, NR_GEP_df = get_GEP_dfs(melanoma_clinical_data, EXEC_DIR, RUN_RANGE, selected_K, number_of_genes, n_replicates, local_density_threshold, FEATURE, FEATURE_A, FEATURE_B)
    GEP_df = filter_contaminated_cells_out_of_GEP_DF(GEP_df, 14, contaminated_cells_path =r'/storage/md_keren/shitay/Data/tables/GEP/subcohort_1.1.22/contaminated_cells_list/GEP14_contaminated_cells.csv')
    R_GEP_df = filter_contaminated_cells_out_of_GEP_DF(R_GEP_df, 14, contaminated_cells_path =r'/storage/md_keren/shitay/Data/tables/GEP/subcohort_1.1.22/contaminated_cells_list/GEP14_contaminated_cells.csv')
    NR_GEP_df = filter_contaminated_cells_out_of_GEP_DF(NR_GEP_df, 14, contaminated_cells_path =r'/storage/md_keren/shitay/Data/tables/GEP/subcohort_1.1.22/contaminated_cells_list/GEP14_contaminated_cells.csv')
    return immune_clusters_barcodes_mapping_df, CD8_clusters_barcodes_mapping_df, myeloid_clusters_barcodes_mapping_df, GEP_df, R_GEP_df, NR_GEP_df


def get_NATMI_FILTERED_DF():
    # READ CSV
    NATMI_df = pd.read_csv(join(NATMI_ROOT_PATH, NATMI_FILE_NAME))

    # filter
    df_filtered = NATMI_df[(NATMI_df['Receptor detection rate'] > DETECTION_THRESHOLD) &
                     (NATMI_df['Ligand detection rate'] > DETECTION_THRESHOLD)]

    df_filtered = df_filtered[(df_filtered['Receptor average expression value'] > AVG_EXPRESSION_THRESHOLD) &
                              (df_filtered['Ligand average expression value'] > AVG_EXPRESSION_THRESHOLD)]
    # df_filtered = df_filtered[(df_filtered['Receptor total expression value']>avg_expression_threshold) &
    #                           (df_filtered['Ligand total expression value']>avg_expression_threshold)]

    df_filtered = df_filtered[
        (df_filtered['Receptor derived specificity of average expression value'] > AVG_SPECIFICITY_THRESHOLD) &
        (df_filtered['Ligand derived specificity of average expression value'] > AVG_SPECIFICITY_THRESHOLD)]

    df_filtered = df_filtered[
        (df_filtered['Receptor derived specificity of total expression value'] > TOTAL_SPECIFICITY_THRESHOLD) &
        (df_filtered['Ligand derived specificity of total expression value'] > TOTAL_SPECIFICITY_THRESHOLD)]

    if FILTER_GEP_IN_NATNI_DF:
        print(f'GEP siganls are filtered out')
        df_filtered = df_filtered[df_filtered[['Sending cluster', 'Target cluster']].apply(
            lambda x: not 'GEP' in x['Sending cluster'] and not 'GEP' in x['Target cluster'], axis=1)]

    print(f'Number of pairs in filtered df is: {len(df_filtered)}')
    ################################ show ################################

    strong_signals_df = df_filtered.sort_values(NATMI_DF_SORT_BY, ascending=False).reset_index(drop=True)

    ####### add parameters

    strong_signals_df[f'Target cluster ({FEATURE_A}) - number of cells'] = None
    strong_signals_df[f'Target cluster ({FEATURE_B}) - number of cells'] = None
    strong_signals_df[f'Target cluster ({FEATURE_A}) - % cells expressing'] = None
    strong_signals_df[f'Target cluster ({FEATURE_B}) - % cells expressing'] = None
    strong_signals_df['Receptor - pval'] = None

    strong_signals_df[f'Sending cluster ({FEATURE_A}) - number of cells'] = None
    strong_signals_df[f'Sending cluster ({FEATURE_B}) - number of cells'] = None
    strong_signals_df[f'Sending cluster ({FEATURE_A}) - % cells expressing'] = None
    strong_signals_df[f'Sending cluster ({FEATURE_B}) - % cells expressing'] = None
    strong_signals_df['Ligand - pval'] = None

    return strong_signals_df


def get_all_cohorts():
    cohort = pickle.load(open(COHORT_PATH, 'rb'))
    immune_cohort = cohort.filter_cells_by_property('is_immune', True)
    tumor_cohort = cohort.filter_cells_by_property('is_cancer', True)
    CD8_cohort = get_requested_subset(cohort, 'CYTOTOXIC_T_CELLS')
    mye_cohort = get_requested_subset(cohort, 'MYELOIDS')

    return cohort, immune_cohort, tumor_cohort, CD8_cohort, mye_cohort


def create_cohorts_of_clusters_per_feature(cohort, clusters_barcodes_mapping_df, clusters_indexes, response_samples, non_response_samples):
    cluster_cohorts_dic = {c_i: {} for c_i in clusters_indexes}
    R_CD8_clusters_df = clusters_barcodes_mapping_df[clusters_barcodes_mapping_df.Sample.isin(response_samples)]
    NR_CD8_clusters_df = clusters_barcodes_mapping_df[
        clusters_barcodes_mapping_df.Sample.isin(non_response_samples)]

    for cluster_curr_index in clusters_indexes:
        R_cells_ids = R_CD8_clusters_df[R_CD8_clusters_df.Cluster == cluster_curr_index][['Sample', 'Barcode']]
        NR_cells_ids = NR_CD8_clusters_df[NR_CD8_clusters_df.Cluster == cluster_curr_index][['Sample', 'Barcode']]

        R_cohort = cohort.get_subset_by_identifiers(R_cells_ids['Sample'], R_cells_ids['Barcode'])
        NR_cohort = cohort.get_subset_by_identifiers(NR_cells_ids['Sample'], NR_cells_ids['Barcode'])

        cluster_cohorts_dic[cluster_curr_index][FEATURE_A] = R_cohort
        cluster_cohorts_dic[cluster_curr_index][FEATURE_B] = NR_cohort

    return cluster_cohorts_dic


def create_cohorts_of_GEPs_per_feature(tumor_cohort, GEP_df, R_GEP_df, NR_GEP_df):
    GEP_indexes = list(GEP_df.columns[4:])
    GEPS_dic = {g_i: {} for g_i in GEP_indexes}
    for GEP_curr_index in GEP_indexes:
        R_cells_ids = R_GEP_df[R_GEP_df['associated program'] == GEP_curr_index][['sample', 'barcode']]
        NR_cells_ids = NR_GEP_df[NR_GEP_df['associated program'] == GEP_curr_index][['sample', 'barcode']]

        NR_cohort = tumor_cohort.get_subset_by_identifiers(NR_cells_ids['sample'], NR_cells_ids['barcode'])
        R_cohort = tumor_cohort.get_subset_by_identifiers(R_cells_ids['sample'], R_cells_ids['barcode'])

        GEPS_dic[GEP_curr_index][FEATURE_A] = R_cohort
        GEPS_dic[GEP_curr_index][FEATURE_B] = NR_cohort
    return GEPS_dic


def conduct_statistical_test_between_expression(R_expression, NR_expression, gene_expression_threshold):
    R_n_cells = len(R_expression)
    NR_n_cells = len(NR_expression)
    R_n_cells_expressing = sum(R_expression > gene_expression_threshold)
    NR_n_cells_expressing = sum(NR_expression > gene_expression_threshold)
    R_per_expressing = R_n_cells_expressing / R_n_cells if R_n_cells != 0 else 0
    NR_per_expressing = NR_n_cells_expressing / NR_n_cells if NR_n_cells != 0 else 0

    oddsratio, pvalue = stats.fisher_exact([[R_n_cells_expressing,
                                             NR_n_cells_expressing],
                                            [R_n_cells - R_n_cells_expressing,
                                             NR_n_cells - NR_n_cells_expressing]])

    return R_n_cells, NR_n_cells, R_per_expressing, NR_per_expressing, pvalue


def conduct_statistical_tests_on_NATMI_signals(cohort, cohorts_dic, strong_signals_df,
                                               gene_expression_threshold=1):
    """
    conduct statistical tests ligands and receptor genes in strong_signals_df and look for differential cases
    between mucosal & cutaneous (or R vs. NR).
    :param cohort:
    :param cohorts_dic: cohorts of all clusters and GEP in a dictionary structure.
    :param strong_signals_df: NATMI df filtered out from weak signals.
    :param gene_expression_threshold: min expression value for a cell to express a gene
    :return:
    """
    print(f'conduct statistical tests on NATMI signals')
    for row_index, current_signal in strong_signals_df.iterrows():

        if not current_signal['Ligand symbol'] in cohort.gene_names or not current_signal[
                                                                               'Receptor symbol'] in cohort.gene_names:
            continue
        ligand_gene_index = cohort.gene_names.index(current_signal['Ligand symbol'])
        receptor_gene_index = cohort.gene_names.index(current_signal['Receptor symbol'])

        cluster_name = current_signal['Sending cluster'].split('_')[0]
        if cluster_name == 'GEP':
            cluster_idx = int(current_signal['Sending cluster'].split('_')[1][1:])
        else:
            cluster_idx = int(current_signal['Sending cluster'].split('_')[2])

        R_counts = cohorts_dic[cluster_name][cluster_idx][FEATURE_A].counts
        NR_counts = cohorts_dic[cluster_name][cluster_idx][FEATURE_B].counts

        R_n_cells, NR_n_cells, R_per_expressing, NR_per_expressing, pval = conduct_statistical_test_between_expression(
            R_counts[:, ligand_gene_index], NR_counts[:, ligand_gene_index], gene_expression_threshold)

        strong_signals_df.at[row_index, [f'Sending cluster ({FEATURE_A}) - number of cells']] = R_n_cells
        strong_signals_df.at[row_index, [f'Sending cluster ({FEATURE_B}) - number of cells']] = NR_n_cells
        strong_signals_df.at[
            row_index, [f'Sending cluster ({FEATURE_A}) - % cells expressing']] = R_per_expressing
        strong_signals_df.at[
            row_index, [f'Sending cluster ({FEATURE_B}) - % cells expressing']] = NR_per_expressing
        strong_signals_df.at[row_index, ['Ligand - pval']] = pval

        cluster_name = current_signal['Target cluster'].split('_')[0]
        if cluster_name == 'GEP':
            cluster_idx = int(current_signal['Target cluster'].split('_')[1][1:])
        else:
            cluster_idx = int(current_signal['Target cluster'].split('_')[2])

        R_counts = cohorts_dic[cluster_name][cluster_idx][FEATURE_A].counts
        NR_counts = cohorts_dic[cluster_name][cluster_idx][FEATURE_B].counts

        R_n_cells, NR_n_cells, R_per_expressing, NR_per_expressing, pval = conduct_statistical_test_between_expression(
            R_counts[:, receptor_gene_index], NR_counts[:, receptor_gene_index], gene_expression_threshold)
        strong_signals_df.at[row_index, [f'Target cluster ({FEATURE_A}) - number of cells']] = R_n_cells
        strong_signals_df.at[row_index, [f'Target cluster ({FEATURE_B}) - number of cells']] = NR_n_cells
        strong_signals_df.at[
            row_index, [f'Target cluster ({FEATURE_A}) - % cells expressing']] = R_per_expressing
        strong_signals_df.at[
            row_index, [f'Target cluster ({FEATURE_B}) - % cells expressing']] = NR_per_expressing
        strong_signals_df.at[row_index, ['Receptor - pval']] = pval

    # remove null lines:
    strong_signals_df = strong_signals_df[~((strong_signals_df['Receptor - pval'].isnull()) | (
        strong_signals_df['Ligand - pval'].isnull()))]

    # correct pvals
    strong_signals_df['Receptor - qval'] = multipletests_fdr(strong_signals_df['Receptor - pval'])[1]
    strong_signals_df['Ligand - qval'] = multipletests_fdr(strong_signals_df['Ligand - pval'])[1]
    # take rows with qval < 0.05
    strong_signals_df = strong_signals_df[
        (strong_signals_df['Receptor - qval'] < 0.05) | (strong_signals_df['Ligand - qval'] < 0.05)]
    strong_signals_df = strong_signals_df.iloc[:,
                               [0, 3, 1, 2] + list(range(25, 29)) + list(range(20, 24)) + [31, 30] +
                               list(range(4, 19)) + [29, 24]]
    return strong_signals_df

def infer_NATMI_cellular_communication(response_samples, non_response_samples, melanoma_clinical_data):
    """
    run NATMI analysis - create a DF of ligand-receptor pairs in myeloid,CD8 clusters and tumor GEPS.
    then, conduct statistical tests to examine which pairs are significant different between mucosal & cutaneous samples.
    return a df of significant cases sorted by specificity of edges.
    :param response_samples:
    :param non_response_samples:
    :param melanoma_clinical_data:
    :return:
    """
    print(f'infer NATMI cellular communication networks')
    immune_clusters_barcodes_mapping_df, CD8_clusters_barcodes_mapping_df, myeloid_clusters_barcodes_mapping_df, GEP_df, R_GEP_df, NR_GEP_df = extract_all_clusters(melanoma_clinical_data)
    strong_signals_df = get_NATMI_FILTERED_DF()
    cohort, immune_cohort, tumor_cohort, CD8_cohort, mye_cohort = get_all_cohorts()

    # Create cohorts of clusters per feature
    CD8_dic = create_cohorts_of_clusters_per_feature(CD8_cohort, CD8_clusters_barcodes_mapping_df,
                                                    list(set(CD8_clusters_barcodes_mapping_df.Cluster)),
                                                    response_samples, non_response_samples)
    mye_dic = create_cohorts_of_clusters_per_feature(mye_cohort, myeloid_clusters_barcodes_mapping_df,
                                                    list(set(myeloid_clusters_barcodes_mapping_df.Cluster)),
                                                    response_samples, non_response_samples)
    GEPs_dic = create_cohorts_of_GEPs_per_feature(tumor_cohort, GEP_df, R_GEP_df, NR_GEP_df)

    cohorts_dic = {'GEP':GEPs_dic, 'Myeloid':mye_dic, 'CD8':CD8_dic}

    # get siginificant cases
    sig_cases_strong_signals_df = conduct_statistical_tests_on_NATMI_signals(cohort, cohorts_dic, strong_signals_df,
                                               gene_expression_threshold=1)
    sig_cases_strong_signals_df.to_excel(join(FOLDER_PATH, r'NATMI_strong_signal_sig_cases.xlsx'))


if __name__ == '__main__':

    create_folder(FOLDER_PATH)
    create_folder(join(FOLDER_PATH, 'cytotoxic_t_cells'))
    create_folder(join(FOLDER_PATH, 'myeloids'))

    # Loads clinical table to get labels - R/NR
    res_melanoma_clinical_data, melanoma_clinical_data = load_clinical_table()

    # Builds response/non_response sample list
    response_samples, non_response_samples, no_used_samples = create_sample_lists(melanoma_clinical_data)

    if TSNE:
        print_tSNEs(response_samples, non_response_samples)

    # Association with clusters - statistical tests: immune\CD8\myeloid
    if CLUSTERING_ANALYSIS:
        run_clustering_analysis()

    # Association with GEPs - statistical tests
    if GEP_ANALYSIS:
        run_GEPS_analysis()

    # differential gene expression analysis
    if DIFFERENTIAL_GENE_EXP_ANALYSIS:
        differential_gene_exp_analysis(response_samples, non_response_samples)

    # Combine pvals files
    if COMBINING_PVAL_FILES:
        combining_pvals_various_experiments()

    if NATMI_ANALYSIS:
        infer_NATMI_cellular_communication(response_samples, non_response_samples, melanoma_clinical_data)