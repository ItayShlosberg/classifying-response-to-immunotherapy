import pandas as pd
import numpy as np
import pickle
from collections import Counter

"""
Smart_seq scRNAseqBuilding 2018 Data. classifying cells into cell-types.
"""

TABLE_PATH = r'D:\Technion studies\Keren Laboratory\research\articles\Tables and files article5\additional files\supervised_analysis.xlsx'
ORIGIN_PICKLE_PATH = r'D:\Technion studies\Keren Laboratory\Data\smart_seq\SmartSeq_RNAseq_DATA.p'
ADDED_INFORMATION_PICKLE_PATH = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\smart_seq\SmartSeq_RNAseq_DATA_3.1.2021.p'
BARCODES_ORDER_PATH = r'C:\Users\itay\Downloads\updated files 4.1.20\sample_name_all.xlsx'

def save_to_pickle(cells_form, gene_names, patients_information):
    pickle.dump((cells_form, gene_names, patients_information), open(ADDED_INFORMATION_PICKLE_PATH, "wb"))


def extract_data_from_pickle():
    cells_form, gene_names, patients_information = pickle.load(open(ORIGIN_PICKLE_PATH, "rb"))
    return cells_form, gene_names, patients_information


def get_supervised_cell_types_list():
    xls = pd.ExcelFile(TABLE_PATH)
    df = pd.read_excel(xls)

    # The columns of the table are actually the first rows, therefore we detach it and add it to the end of the table
    columns = list(df.columns)
    new_row = [columns[0]] + [x if type(x) == int else 0 for x in columns[1:]]
    df.loc[29] = new_row
    # Now, changes the columns to integers to avoid confusion.
    new_columns = {columns[i]:i for i in range(len(columns))}
    cell_types_list = df.rename(columns=new_columns).values

    # cell_types_list holds each cell_type and indexes belong to the CD45 cells indexes (1-16291).
    # there are cells that were wiped out during QC also (16,292-~20k)
    def remove_zeros_from_list(l):
        return [i for i in l if i!=0]
    cell_types_list = [{cell_type[0]: remove_zeros_from_list(cell_type[1:])} for cell_type in cell_types_list]
    return cell_types_list

def extract_barcodes_order():
    xls = pd.ExcelFile(BARCODES_ORDER_PATH)
    annotations_df = pd.read_excel(xls)
    ann = annotations_df.columns.tolist() + [a[0] for a in annotations_df.values.tolist()]
    ann_dic = {idx+1: barcode for idx,barcode in enumerate(ann)}
    return ann_dic


def add_cell_type_to_cell_patients(patients_information, cell_types_list):
    # adds the types of every cell to patients_information.
    for patient in patients_information:
        patient['supervised classification'] = []

    ann_order_dic = extract_barcodes_order()
    barcode_to_idx = {pt['cell id']: idx for idx, pt in enumerate(patients_information)}
    for cell_type_dict in cell_types_list:
        for cell_type, cells_indexes in cell_type_dict.items():
            for cell_idx in cells_indexes:
                barcode = ann_order_dic[cell_idx]
                if barcode in barcode_to_idx.keys():
                    corresponding_idx = barcode_to_idx[barcode]
                    patients_information[corresponding_idx]['supervised classification'].append(cell_type)
                else:
                    print(f"Barcode: {barcode} has not been found in patients_information list and has been noted in cell-types list as {cell_type}")

    # for cell in range(1, 16292):
    #     print(f'working on cell number: {cell}')
    #     corresponding_cell_idx = cell-1     # (1-16291) to (0-16290)
    #     current_cell_types = []
    #     for cell_type in cell_types_list:
    #         if cell in cell_type[1]:
    #             current_cell_types.append(cell_type[0])


    return patients_information


def add_supervised_cell_types_to_patients(patients_information):
    cell_types_list = get_supervised_cell_types_list()
    patients_information = add_cell_type_to_cell_patients(patients_information, cell_types_list)
    return patients_information


if __name__ == '__main__':

    cells, gene_names, patients_information = extract_data_from_pickle()

    patients_information = add_supervised_cell_types_to_patients(patients_information)

    # saves the new information added
    save_to_pickle(cells, gene_names, patients_information)
    print(f"saved in {ADDED_INFORMATION_PICKLE_PATH}")
