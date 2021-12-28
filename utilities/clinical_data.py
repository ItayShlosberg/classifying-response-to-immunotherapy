from utilities.general_helpers import *
import numpy as np
from os.path import join
import numpy as np
import pickle
from DL.Mars_seq_DL.data_loading import *
import pandas as pd
import time


def therapy_ICI_translator(x):
    if not x or str(x) == 'nan':
        return 'nothing'
    x = str(x).lower()
    if 'ipi' in x or 'pembro' in x or 'pembo' in x or 'nivo' in x or 'apd' in x or 'pd1' in x:
        return 'ICI'

    elif 'enco' in x or 'd+t' in x:
        return 'other therapy'

    elif 'radiation' in x:
        return 'other therapy'

    elif 'carbo' in x or 'imatinib' in x or 'tmz' in x or 'tvec' in x or 'cdk4/6' in x:
        return 'other immune therapy'

    elif 'surgical ' in x or 'surgery ' in x or 'neodajuvant' in x:
        return 'other therapy'

    else:
        return x


def therapy_IPI_transletror(x):
    if not x or str(x)=='nan':
        return 'nothing'
    x = str(x).lower()
    if 'ipi' in x and 'nivo' in x:
        return 'ICI'
    else: return 'other therapy'


def get_clinical_data(n_samples=71, ICI=None, melanoma_type=None, prior_biopsy=None, after_biopsy=None,
                      only_metastasis_sample=False, response=None, therapy_translator=therapy_ICI_translator):
    #  Loads xlsx files
    CLINICAL_LABELS_PATH = r'/storage/md_keren/shitay/Data/tables/clinical_labels.xlsx'
    MELANOMA_CLINICAL_DATA_PATH = r'/storage/md_keren/shitay/Data/tables/Melanoma_clinical_data_12.21_unportected.xlsx'
    print(f'Using clinical table in path:\n {MELANOMA_CLINICAL_DATA_PATH}\n\nand labels:\n{CLINICAL_LABELS_PATH}')

    melanoma_clinical_data = pd.read_excel(MELANOMA_CLINICAL_DATA_PATH)
    clinical_labels = pd.read_excel(CLINICAL_LABELS_PATH)
    # takes only first 71 samples, fill Nan and creat dictionary mapping
    melanoma_clinical_data = melanoma_clinical_data.iloc[:n_samples][
        ['Patient id', 'Clinical response', 'Melanoma type', 'Therapy(ies) prior to biopsy', 'Therapy after biopsy',
         'Primary=1; Metastasis=0', 'Genotype ']]
    # fill Nans
    melanoma_clinical_data['Melanoma type'] = melanoma_clinical_data['Melanoma type'].fillna('??')
    melanoma_clinical_data['Clinical response'] = melanoma_clinical_data['Clinical response'].fillna('??')
    # Take only samples that are metastasis (if requested)
    if only_metastasis_sample:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['Primary=1; Metastasis=0'] != 1]
    melanoma_clinical_data = melanoma_clinical_data.drop(columns='Primary=1; Metastasis=0')
    # Convert response label mapping using clinical_label table
    labels_mapping = {v[0]: v[1] for v in clinical_labels[['Clinical response', 'binary label']].values}
    labels_mapping['R '] = 'R'
    labels_mapping['PD (NR; for pembro) NR for d+t'] = 'NR'
    melanoma_clinical_data['response'] = melanoma_clinical_data['Clinical response'].apply(lambda x: labels_mapping[x])
    # Convert melanoma type
    melanoma_type_translate = {'Cutaneous': 'Cutaneous', 'Mucosal ': 'Mucosal', 'Uveal': 'other', 'UN primary': 'other',
                               'Acral': 'other', 'Unknown': 'other', '??': 'other'}
    melanoma_clinical_data['Melanoma type'] = melanoma_clinical_data['Melanoma type'].apply(
        lambda rr: melanoma_type_translate[rr])
    # Convert therapy prior/after biopsy
    melanoma_clinical_data['prior to biopsy'] = melanoma_clinical_data['Therapy(ies) prior to biopsy'].apply(
        lambda x: therapy_translator(x))
    melanoma_clinical_data['after biopsy'] = melanoma_clinical_data['Therapy after biopsy'].apply(
        lambda x: therapy_translator(x))
    melanoma_clinical_data['ICI'] = (melanoma_clinical_data['prior to biopsy'] == 'ICI') | (
            melanoma_clinical_data['after biopsy'] == 'ICI')
    melanoma_clinical_data = melanoma_clinical_data.drop(
        columns=['Therapy(ies) prior to biopsy', 'Therapy after biopsy', 'Clinical response'])
    # Convert Genotype - save information only for BRAF: is BRAF mutated or not BRAF mutated
    melanoma_clinical_data['BRAF'] = melanoma_clinical_data['Genotype '].astype(str).apply(lambda x: True if 'BRAF' in x else False)
    melanoma_clinical_data = melanoma_clinical_data.drop(columns=['Genotype '])

    if not ICI is None:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['ICI'] == ICI]
    if not after_biopsy is None:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['after biopsy'] == after_biopsy]
    if not prior_biopsy is None:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['prior to biopsy'] == prior_biopsy]
    if not melanoma_type is None:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['Melanoma type'].isin(melanoma_type)]
    if not response is None:
        melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['response'] == response]

    return melanoma_clinical_data.reset_index()


def get_clinical_subdata(n_R_muc, n_NR_muc, n_R_cut, n_NR_cut, n_samples=71, ICI=None, melanoma_type=None, prior_biopsy=None, after_biopsy=None,
                      only_metastasis_sample=False, response=None, therapy_translator=therapy_ICI_translator):
    melanoma_clinical_data = get_clinical_data(n_samples=n_samples, ICI=ICI, melanoma_type=melanoma_type, prior_biopsy=prior_biopsy, after_biopsy=after_biopsy,
                      only_metastasis_sample=only_metastasis_sample, response=response, therapy_translator=therapy_translator)

    try:
        R_Cut_df = melanoma_clinical_data[(melanoma_clinical_data['response'] == 'R').values &
                                          (melanoma_clinical_data['Melanoma type'] == 'Cutaneous').values].sample(
            n_R_cut)

        NR_Cut_df = melanoma_clinical_data[(melanoma_clinical_data['response'] == 'NR').values &
                                           (melanoma_clinical_data['Melanoma type'] == 'Cutaneous').values].sample(
            n_NR_cut)

        R_Muc_df = melanoma_clinical_data[(melanoma_clinical_data['response'] == 'R').values &
                                          (melanoma_clinical_data['Melanoma type'] == 'Mucosal').values].sample(n_R_muc)

        NR_Muc_df = melanoma_clinical_data[(melanoma_clinical_data['response'] == 'NR').values &
                                           (melanoma_clinical_data['Melanoma type'] == 'Mucosal').values].sample(
            n_NR_muc)

        return pd.concat([R_Muc_df, NR_Muc_df, R_Cut_df, NR_Cut_df]).reset_index()
    except Exception:
        print('There are no enough samples as you requested! valid the number that you ask!')

def get_clinical_table_size(n_samples=71, ICI=None, melanoma_type=None, prior_biopsy=None, after_biopsy=None,
                      only_metastasis_sample=False, response=None, therapy_translator=therapy_ICI_translator):
    melanoma_clinical_data = get_clinical_data(n_samples=n_samples, ICI=ICI, melanoma_type=melanoma_type, prior_biopsy=prior_biopsy, after_biopsy=after_biopsy,
                      only_metastasis_sample=only_metastasis_sample, response=response, therapy_translator=therapy_translator)
    return len(melanoma_clinical_data)

def get_constant_cohort(comparison_type=1):
    """
    Constant cohort for analysis of R vs NR
    Taking from the cutaneous group 9 NR and 5 R for the analysis.
    Taking from the mucosal group 4 R for the analysis.
    :param cutaneous_comparison:
    1. NR_Mucosal VS. NR_Cutaneous (drop R cutaneous samples).
    2. NR_Mucosal VS. R_Cutaneous (drop NR cutaneous samples).
    3. NR_Mucosal VS. NR_Cutaneous and R_Cutaneous (all sub cohort).

    :return:
    """
    samples = ['M128', 'M136', 'M111', 'M99', 'M137', 'M147', 'M131', 'M130', 'M162', 'M145',
               'M153', 'M141', 'M107', 'M146', 'M106', 'M161', 'M120', 'M151']
    melanoma_clinical_data = get_clinical_data(therapy_translator=therapy_IPI_transletror)
    melanoma_clinical_data = melanoma_clinical_data[melanoma_clinical_data['Patient id'].isin(samples)]
    response = melanoma_clinical_data['response'].values
    melanoma_type = melanoma_clinical_data['Melanoma type'].values

    if comparison_type == 1:
        return melanoma_clinical_data[(melanoma_type == 'Mucosal') | ((response == 'NR') & (melanoma_type == 'Cutaneous'))]
    if comparison_type == 2:
        return melanoma_clinical_data[(melanoma_type == 'Mucosal') | ((response == 'R') & (melanoma_type == 'Cutaneous'))]
    if comparison_type == 3:
        return melanoma_clinical_data