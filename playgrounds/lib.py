from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import numpy as np

print("kmeans")
X = np.array([[1, 2], [1, 4], [1, 0],
                  [10, 2], [10, 4], [10, 0]])

kmeans = KMeans(n_clusters=2, random_state=0).fit(X)

print(kmeans.labels_)
# array([1, 1, 1, 0, 0, 0], dtype=int32)

print(kmeans.predict([[0, 0], [12, 3]]))
# array([1, 0], dtype=int32)

print(kmeans.cluster_centers_)
# array([[10.,  2.],
#        [ 1.,  2.]])

print()
print()
print("tSNE")
X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
X_embedded = TSNE(n_components=2).fit_transform(X)
print(X_embedded)
print(X_embedded.T.tolist())


print()
print()
print()
import matplotlib.pyplot as plt
plt.plot([1, 2], [1, 4], 'ro')
plt.plot([3, 4], [9, 16], 'ro', color='green')
# plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
# plt.axis([0, 6, 0, 20])
plt.show()