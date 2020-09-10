import math
import numpy as np
import Bio


# for x in range(1, 18):
#     y = x + 1
#     z = math.log(y, 2)
#     print(f" number {x},   log2(TPM+1): {z}")


path = r'DATA\high expressed genes list based on article5.txt'

with open(path, 'r') as f:
    lines = f.readlines()

lines = [l.replace("\n", "") for l in lines if len(l)>1]

print(lines)
