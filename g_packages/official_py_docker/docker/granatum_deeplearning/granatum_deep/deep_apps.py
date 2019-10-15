from granatum_clustering.clustering_apps import GranatumClustering
from granatum_clustering.clustering_utils import _silhouette_analysis

from granatum_clustering.granatum_config import FIND_BEST_CLUSTER_RANGE

from granatum_deep.deepmodel_base import DeepBase

from granatum_deep.deep_config import EPOCHS
from granatum_deep.deep_config import LEVEL_DIMS_IN
from granatum_deep.deep_config import LEVEL_DIMS_OUT
from granatum_deep.deep_config import DROPOUT

from granatum_clustering.granatum_config import CLUSTERING
from granatum_clustering.granatum_config import EMBEDDING


class GranatumDeepClustering(GranatumClustering):
    """
    Clustering and embedding of the data using Deep neural-networks

    This package use the original `granatum_clustering` package and all the clustering
    algorithms can be used. Moreover, The `autoencoder` algorithm can be used to cluster
    or embed, prior to clustering, the data.


    CURRENT CLUSTERING ALGORITHM :      kmeans, WARD, DBSCAN, HDBSCAN, autoencoder
                                        Agglomerative_correlation, Gaussian_Mixture
    CURRENT EMBEDDING ALGORITHM :       autoencoder

        Parameters
    ----------
    selected_clustering : str <'kmeans': default>        Selected clustering
    n_components : Int <30: default>                     Hidden layer in the autoencoder network
    n_clusters : Int  <2: default>                       only used if find_best_number_of_cluster
                                                         is set to False
    epochs : Int <10: default>                           Number of epoch used bu the DLLs
    dropout : 0.0 < float < 1.0  <0.5: default>          Number of epoch used bu the DLLs

    find_best_number_of_cluster: Bool <True: default>    automatically find the optimal
                                                         number of clusters
    embedding_with_autoencoder: Bool <True: default>     Transform the data, prior to clustering,
                                                         using autoencoder
    clustering_with_autoencoder: Bool <False: default>   Use autoencoder to find the best cluster
                                                         for each sample
    level_dims_in: list(Int) <default: [50]>             Defines the levels and the dimensions of
                                                         the hidden layers before the bottleneck
    level_dims_out: list(Int) <default: [50]>            Defines the levels and the dimensions of
                                                         the hidden layers after the bottleneck


    Return (when fit is called)
    ----------
    json : {
        'clusters' :        cluster dictionnary
        'clusters_array' :        cluster array
        'n_clusters' :      number of clusters
        'n_components':     number of components used for the embbedding
        'embedding':        embedding algorithm used
        'outliers':         outliers array
        'embedding_for_plotting':    embedding algorithm used for plotting
        'clustering_algorithm':      clustering algorithm used
        'plot_figure':               is figure plotted
        'figure_png':                figure plotted (png) converted into base64
        'figure_html':                figure plotted html
     }

    """
    def __init__(self,
                 n_components=10,
                 n_clusters=3,
                 selected_embedding='autoencoder',
                 level_dims_in=LEVEL_DIMS_IN,
                 level_dims_out=LEVEL_DIMS_OUT,
                 dropout=DROPOUT,
                 epochs=EPOCHS,
                 _clustering=CLUSTERING,
                 _embedding=EMBEDDING,
                 **kwargs):
        """
        """
        self.level_dims_in = level_dims_in
        self.level_dims_out = level_dims_out
        self.dropout = dropout
        self.epochs = epochs

        _embedding['autoencoder'] = DeepBase(
            dropout=self.dropout,
            epochs=self.epochs,
            level_dims_in=self.level_dims_in,
            level_dims_out=self.level_dims_out,
            n_components=n_components)


        _clustering['autoencoder'] = DeepBase(
            dropout=self.dropout,
            epochs=self.epochs,
            level_dims_in=self.level_dims_in,
                level_dims_out=self.level_dims_out,
                n_components=n_clusters)

        GranatumClustering.__init__(self,
                                    selected_embedding=selected_embedding,
                                    n_components=n_components,
                                    n_clusters=n_clusters,
                                    _embedding=_embedding,
                                    _clustering=_clustering,
                                    **kwargs)

        self._plotting_embedding['autoencoder'] = DeepBase(
            dropout=self.dropout,
            epochs=self.epochs,
            level_dims_in=self.level_dims_in,
            level_dims_out=self.level_dims_out,
            n_components=2)

        self._find_best_cluster['autoencoder'] = {
            'method': _silhouette_analysis, 'range': FIND_BEST_CLUSTER_RANGE}

def main():
    """ """
    from sklearn.datasets import make_blobs

    test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)

    metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(test_dataset))}
    plot_figures_path = '/home/opoirion/code/d3visualisation/'

    granatum = GranatumDeepClustering(
        n_clusters=5,
        selected_embedding='PCA',
        selected_clustering='kmeans',
    )

    granatum.fit(test_dataset, metadata=metadata)
    granatum.plot(
        local_path_to_save_plots=plot_figures_path)


if __name__ == '__main__':
    main()
