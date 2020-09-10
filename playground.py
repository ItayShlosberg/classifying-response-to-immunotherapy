import math
import numpy as np
import Bio

l = [[1,2], [3,4], [1,2], [3,4]]
data = np.array(l)

# print(a)
# print(a.T)


for x in range(1, 18):
    y = x + 1
    z = math.log(y, 2)
    print(f" number {x},   log2(TPM+1): {z}")


from Bio.Cluster import distancematrix
matrix = distancematrix(data)
print(matrix)
