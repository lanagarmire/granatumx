from sklearn.decomposition import NMF

from sklearn.metrics import silhouette_score

import numpy as np


class NMFClustering():
    """
    """
    def __init__(self, n_clusters=2):
        """ """
        self.n_clusters = n_clusters
        self.nmf = NMF(n_components=n_clusters)

    def fit(self, X):
        """
        """
        self.nmf.fit(X)

    def predict(self, X):
        """
        """
        return np.argmax(self.nmf.transform(X), axis=1)

    def fit_predict(self, X):
        """
        """
        self.fit(X)
        return self.predict(X)

    def set_params(self, **kwargs):
        """
        """
        for param in kwargs:
            setattr(self, param, kwargs[param])


def _silhouette_analysis(algo, matrix, clusters, **kwargs):
    """
    """
    max_score = 0
    max_n_clusters = None

    n_samples = len(matrix)

    for n_clusters in clusters:
        algo.set_params(n_clusters=n_clusters)
        clusters = algo.fit_predict(matrix)

        if len(set(clusters)) == 1:
            continue

        score = silhouette_score(matrix, clusters)

        if score >= max_score:
            max_score = score
            max_n_clusters = n_clusters

        if n_samples - 1 <= n_clusters:
            break

    algo.set_params(n_clusters=max_n_clusters)

    return max_n_clusters

def _silhouette_DBSCAN_analysis(algo, matrix, scale_list, scale='eps'):
    """
    """
    max_score = None
    max_eps = None
    number_of_clusters = None
    number_of_samples = len(matrix)

    max_scale = 30

    while True:
        for eps in scale_list:
            params = {scale: eps}
            algo.set_params(**params)
            clusters = algo.fit_predict(matrix)
            outliers = clusters == -1
            number_of_clusters = len(set(clusters))

            if number_of_clusters <= 1:
                continue

            score = silhouette_score(matrix, clusters)
            score = score * float(len(outliers)) / number_of_samples

            if max_score is None:
                max_score = score

            if score >= max_score:
                max_eps = eps

        if max_eps is not None:
            break
        else:
            if scale_list.max() > max_scale:
                break

            scale_list = scale_list + scale_list.max()

    if max_eps is None:
        max_eps = scale_list[-1]

    algo.set_params(**{scale:max_eps})

    return max_eps

def _mixture_bic_analysis(algo, matrix, clusters, **kwargs):
    """
    """
    min_score = None
    max_n_clusters = None
    n_samples = len(matrix)

    assert(algo.bic)

    for n_clusters in clusters:
        algo.set_params(n_components=n_clusters)
        algo.fit(matrix)
        clusters = algo.predict(matrix)

        score = algo.bic(matrix)

        if min_score is None:
            min_score = score

        if score <= min_score:
            min_score = score
            max_n_clusters = n_clusters

        if n_samples - 1 <= n_clusters:
            break

    algo.set_params(n_components=max_n_clusters)


    return max_n_clusters
