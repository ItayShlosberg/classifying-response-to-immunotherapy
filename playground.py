import math
import numpy as np
import Bio
"""
use: pycco main.py
for documentations
"""

# for x in range(1, 18):
#     y = x + 1
#     z = math.log(y, 2)
#     print(f" number {x},   log2(TPM+1): {z}")


from sklearn.model_selection import train_test_split
import numpy as np
n_samples, n_features, n_classes = 10, 2, 2
data = np.random.randn(n_samples, n_features)  # 10 training examples
labels = np.random.randint(n_classes, size=n_samples)  # 10 labels
indices = np.arange(n_samples)
x1, x2, y1, y2, idx1, idx2 = train_test_split(
    data, labels, indices, test_size=0.2)


_breakpoint = 0