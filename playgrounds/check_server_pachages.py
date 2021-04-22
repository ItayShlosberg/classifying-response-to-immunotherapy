
""""
TODO: remove this script. only to test packages installation.
"""
import numpy as np
import pandas as pd
import os, errno
import datetime
import uuid
import itertools
import yaml
import subprocess
import scipy.sparse as sp


from scipy.spatial.distance import squareform
from sklearn.decomposition import non_negative_factorization
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.utils import sparsefuncs



from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list

import matplotlib.pyplot as plt

import scanpy as sc


from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list

import matplotlib.pyplot as plt

import scanpy as sc

import palettable
from bhtsne import tsne

from sklearn import preprocessing
from sklearn.decomposition import PCA
import scanpy as sc

from IPython.display import Image
