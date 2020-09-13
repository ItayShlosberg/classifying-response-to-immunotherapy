import math
import numpy as np
import Bio
"""
use: pycco main.py
for documentations
"""


from general_helpers import *


# em = Experiments_manager("exp1", "Data")
# em.activate_prints_to_file()
#     # (r'Data\out.txt')
# for i in range(10):
#         print(f"print number {i}")
#
# em.finish_run()
def ff():
    return 'exp1'
E = ff()
D = 'Data'
@experiment_manager(E, D)
def main():
    for i in range(25):
        print(f"11 decorator print number {i}")

# em.finish_run()


# print("print number 1")

if __name__ == '__main__':
    main()
    E = 'exp2'
    D = 'Data'
    main()

