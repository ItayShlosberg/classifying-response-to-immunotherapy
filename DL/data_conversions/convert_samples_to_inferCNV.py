"""
Use that script to create input to InferCNV.
You should define the input path of the samples in PKL format.
InferCNV annotation: immune -> immune. NOT immune - Cancer (even if there is no tag at all).
"""

from utilities.droplet_dataset import *
from DL.Mars_seq_DL.data_loading import extract_droplet_data_from_pickle
from os.path import join


INPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\17.2.21'
OUTPUT_DIR = r'D:\Technion studies\Keren Laboratory\python_playground\outputs\temporal garbage\converted'

def convert_matrix(rna_sample, sample_id):

    gene_duplications = {g:False for g in rna_sample.gene_names}

    writer_path = join(OUTPUT_DIR, sample_id)
    create_folder(writer_path)
    with open(join(writer_path, 'matrix.matrix'), "wb") as writer:
        header = str.encode('\t'.join(rna_sample.barcodes)+'\n')
        writer.write(header)
        for idx, gene_name in enumerate(rna_sample.gene_names):
            if gene_duplications[gene_name]:
                continue
            gene_duplications[gene_name] = True
            # print(f'{idx}/{len(rna_sample.gene_names)}')
            values = '\t'.join(rna_sample.counts[:, idx].astype(str).tolist())
            line = gene_name+'\t'+values + '\n'
            writer.write(str.encode(line))


    _breakpoint = 0


def create_annotations(rna_sample, sample_id):
    immune_indexes = [idx for idx, ci in enumerate(rna_sample.cells_information) if ci.is_immune]
    cancer_indexes = [idx for idx, ci in enumerate(rna_sample.cells_information) if not ci.is_immune]

    immune_barcodes = [rna_sample.barcodes[i] for i in immune_indexes]
    cancer_barcodes = [rna_sample.barcodes[i] for i in cancer_indexes]

    writer_path = join(OUTPUT_DIR, sample_id)
    create_folder(writer_path)
    with open(join(writer_path, 'annotation.txt'), "wb") as writer:
        for barcode in immune_barcodes:
            st = str.encode(barcode + '\t' + 'immune' + '\n')
            writer.write(st)

        for barcode in cancer_barcodes:
            st = str.encode(barcode + '\t' + 'cancer' + '\n')
            writer.write(st)


def convert_all_samples():
    samples = [subfolder for subfolder in os.listdir(INPUT_DIR)]

    # Extract ImmuneCellsMarkersUpdated Excel file from PC and load it into DataFrame.

    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    for sample_id in samples:
        print(sample_id)
        # Extracts one of the samples from PC
        sample_path = join(INPUT_DIR, sample_id, f'{sample_id}.pkl')
        rna_sample = extract_droplet_data_from_pickle(sample_path)

        convert_matrix(rna_sample, sample_id)
        create_annotations(rna_sample, sample_id)


if __name__ == '__main__':

    convert_all_samples()

