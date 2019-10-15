from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import FastICA
from sklearn.decomposition import LatentDirichletAllocation
from sklearn.decomposition import NMF

from sklearn.manifold import MDS
from sklearn.manifold import TSNE
from sklearn.manifold import Isomap

from hdbscan import HDBSCAN

from granatum_clustering.clustering_utils import NMFClustering

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN

from sklearn.mixture import GaussianMixture

import numpy as np

from granatum_clustering.clustering_utils import _silhouette_analysis
from granatum_clustering.clustering_utils import _silhouette_DBSCAN_analysis
from granatum_clustering.clustering_utils import _mixture_bic_analysis


CLUSTERING = {
    'kmeans': KMeans(n_jobs=-1),
    'DBSCAN': DBSCAN(n_jobs=-1),
    'WARD': AgglomerativeClustering(),
    'HDBSCAN': HDBSCAN(core_dist_n_jobs=-1),
    # 'NMF': NMFClustering(),
    'Agglomerative_correlation': AgglomerativeClustering(
    affinity='correlation', linkage='average'),
    'Gaussian_Mixture': GaussianMixture(covariance_type='spherical'),
}

EMBEDDING = {
    'PCA': PCA(),
    'IncrementalPCA': IncrementalPCA(),
    'FactorAnalysis': FactorAnalysis(),
    'FastICA': FastICA(),
    'NMF': NMF(),
    'LatentDirichletAllocation': LatentDirichletAllocation(
        learning_method='online'),
    'MDS': MDS(n_jobs=-1),
    None:None,
}

PLOTTING_EMBEDDING = {
    'PCA': PCA(),
    'IncrementalPCA': IncrementalPCA(),
    'FactorAnalysis': FactorAnalysis(),
    'FastICA': FastICA(),
    'ISOMAP': Isomap(n_neighbors=10),
    'NMF': NMF(),
    'TSNE': TSNE(),
    'MDS': MDS(n_jobs=-1),
}

FIND_BEST_CLUSTER_RANGE = range(2,30)
FIND_BEST_DBSCAN_CLUSTER_RANGE = np.linspace(0.1, 5.0, 30)

FIND_BEST_CLUSTER = {
    'kmeans': {'method': _silhouette_analysis,
               'range': FIND_BEST_CLUSTER_RANGE},
    'WARD': {'method': _silhouette_analysis,
               'range': FIND_BEST_CLUSTER_RANGE},
    'DBSCAN': {'method': _silhouette_DBSCAN_analysis,
               'range': FIND_BEST_DBSCAN_CLUSTER_RANGE,
               'params': {'scale':'eps'}},
    'HDBSCAN': {'method': _silhouette_DBSCAN_analysis,
               'range': FIND_BEST_DBSCAN_CLUSTER_RANGE,
                'params': {'scale':'alpha'}},
    'Agglomerative_correlation': {'method': _silhouette_analysis,
               'range': FIND_BEST_CLUSTER_RANGE},
    'Gaussian_Mixture': {'method': _mixture_bic_analysis,
               'range': FIND_BEST_CLUSTER_RANGE},
}
