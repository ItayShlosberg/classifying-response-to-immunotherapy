from os.path import join
from utilities.general_helpers import *

SOURCE_PATH = r'D:\PycharmProjects\CellBender\conversions'
OUTPUT_PATH = r'D:\PycharmProjects\CellBender\13.12.20'




samples = [subfolder for subfolder in os.listdir(SOURCE_PATH) if not os.path.isfile(join(SOURCE_PATH,subfolder))]


create_folder(OUTPUT_PATH)
for sample_id in samples:
    path = join(OUTPUT_PATH, sample_id)
    print(path)
    create_folder(path)