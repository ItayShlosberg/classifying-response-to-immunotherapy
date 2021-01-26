
from utilities.droplet_dataset import *
from os.path import join

SOURCE_PATH = r'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples'
PKL_name = r'enrich_RNA_sample.pkl'

def remove():
    samples = [subfolder for subfolder in os.listdir(SOURCE_PATH)]

    for sample_id in samples:
        file = join(SOURCE_PATH, sample_id, PKL_name)
        if os.path.isfile(file):
            print(file)
            os.remove(file)


if __name__ == '__main__':
    remove()