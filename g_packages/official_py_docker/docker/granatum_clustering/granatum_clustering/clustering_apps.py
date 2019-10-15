from io import BytesIO

from granatum_clustering.granatum_config import CLUSTERING
from granatum_clustering.granatum_config import EMBEDDING
from granatum_clustering.granatum_config import PLOTTING_EMBEDDING
from granatum_clustering.granatum_config import FIND_BEST_CLUSTER

from os import remove

import numpy as np

import json

from collections import defaultdict

import base64


class GranatumClustering():
    """
    CURRENT CLUSTERING ALGORITHM : kmeans, WARD, DBSCAN, HDBSCAN,
                                   Agglomerative_correlation, Gaussian_Mixture
    CURRENT EMBEDDING ALGORITHM  : PCA, IncrementalPCA, FactorAnalysis, FastICA,
                                   NMF, LatentDirichletAllocation, MDS, None



    Parameters
    ----------
    selected_embedding          : str <'PCA' : default>
    selected_clustering         : str <'kmeans' : default>
    n_components                : int <30 : default>
    n_clusters                  : int <2 : default> only used if find_best_number_of_cluster is set to False
    plot_figure                 : Bool <True : default> plot figure in png
    plot_figure_html            : Bool <True : default> plot figure in html
    find_best_number_of_cluster : Bool <True : default> automatically find the optimal number of clusters
    figsize                     : tuple (int, int)  <(12, 12) : default> matplotlig size of the figure

    Return (when fit is called)
    ----------
    json : {
        'clusters'             : cluster dictionnary
        'clusters_array'       : cluster array
        'n_clusters'           : number of clusters
        'n_components'         : number of components used for the embbedding
        'embedding'            : embedding algorithm used
        'outliers'             : outliers array
        'clustering_algorithm' : clustering algorithm used
     }

    """
    def __init__(self,
                 selected_embedding='PCA',
                 selected_clustering='kmeans',
                 find_best_number_of_cluster=True,
                 n_components=10,
                 n_clusters=2,
                 _plotting_embedding=PLOTTING_EMBEDDING,
                 _find_best_cluster=FIND_BEST_CLUSTER,
                 _clustering=CLUSTERING,
                 _embedding=EMBEDDING,
                 **kwargs):
        """
        """
        self._plotting_embedding = _plotting_embedding
        self._find_best_cluster = _find_best_cluster
        self._clustering = _clustering
        self._embedding = _embedding

        self.selected_plotting_embedding = None
        self.local_path_to_save_plots = None

        self.selected_embedding = selected_embedding
        self.selected_clustering = selected_clustering

        self.find_best_number_of_cluster = find_best_number_of_cluster

        assert(selected_embedding in self._embedding)
        assert(selected_clustering in self._clustering)

        self.n_components = n_components
        self.n_clusters = n_clusters

        if selected_embedding is not None:
            self.embedding = self._embedding[selected_embedding].set_params(
                n_components=self.n_components)

        self.clustering = self._clustering[selected_clustering]

        if selected_clustering not in ['DBSCAN', 'HDBSCAN', 'Gaussian_Mixture']:
            self.clustering.set_params(n_clusters=self.n_clusters)

        if selected_clustering == 'Gaussian_Mixture':
            self.clustering.set_params(n_components=self.n_clusters)

        self._init_values()

    def _init_values(self):
        """
        """
        self.clusters = None
        self.metadata = {}
        self.sample_ids = None
        self.matrix = None
        self.figure = None
        self.figure_html = None
        self.figure_png = None

        self.outliers = []

    def fit(self, matrix, sample_ids=None, metadata={}, jsonify=False):
        """
        main function to transform an input matrix.
        return a json object.
        In case of metadata, each sample should have the same metadatas

        input
        ----------
        matrix     : (sample x gene) matrix <float>
        sample_ids : array<str>     OPTIONAL
        jsonify    : bool    OPTIONAL
        metadata   : dict {'sample_ids' : {metadata_key : metadata_value}}    OPTIONAL

        Return
        ----------
        json : {
            'clusters'               : cluster dictionnary
            'clusters_array'         : cluster array
            'n_clusters'             : number of clusters
            'n_components'           : number of components used for the embbedding
            'embedding'              : embedding algorithm used
            'outliers'               : outliers array
            'embedding_for_plotting' : embedding algorithm used for plotting
            'clustering_algorithm'   : clustering algorithm used

         }

        """
        self._init_values()

        self.sample_ids = sample_ids
        self.metadata = metadata

        if not sample_ids:
            self.sample_ids = range(len(matrix))

        if self.selected_embedding is not None:
            if self.selected_embedding in ['NMF', 'LatentDirichletAllocation']:
                matrix = self._shift_data_to_positive(matrix)

            matrix = self.embedding.fit_transform(matrix)

        if self.find_best_number_of_cluster:
            self._compute_best_cluster(matrix)

        self.clusters = self._fit_predict(matrix)
        self.clusters_array = self.clusters[:]

        self._format_found_clusters()

        results = {
            'clusters': self.clusters,
            'clusters_array': self.clusters_array.tolist(),
            'n_clusters': len(set(self.clusters.values())),
            'n_components': self.n_components,
            'embedding': self.selected_embedding,
            'outliers': self.outliers,
            'clustering_algorithm': self.selected_clustering,
        }

        if jsonify:
            results = json.dumps(results, indent=2)

        self.matrix = matrix

        return results

    def plot(self,
             embedding='PCA',
             local_path_to_save_plots=None,
             plot_figure_png=True,
             plot_figure_html=True,
             figsize=(600, 600 ),
             jsonify=False):
        """
        CURRENT PLOTTING ALGORITHM : PCA, IncrementalPCA, FactorAnalysis, FastICA,
                                     NMF, TSNE, MDS, ISOMAP

        input
        ----------
        local_path_to_save_plots    : str  <None default>    Local path where to save the plot
        figsize                     : tuple <int : 600, int : 600 default>      Number of pixels for width / height
        plot_figure_png             : bool OPTIONAL
        plot_figure_html            : bool OPTIONAL
        jsonify                     : bool OPTIONAL
        selected_plotting_embedding : str <'PCA' : default>

        Return
        ----------
        json : {
            'figure_png'             : figure plotted (png) converted into base64
            'figure_html'            : figure plotted in html
            'embedding_for_plotting' : embedding algorithm used for plotting
         }
        """
        assert(embedding in self._plotting_embedding)
        self.local_path_to_save_plots = local_path_to_save_plots

        if embedding is not None:
            self.embedding_for_plotting = self._plotting_embedding[
                embedding].set_params(
                    n_components=2)

        self.selected_plotting_embedding = embedding

        self.figsize = figsize

        if plot_figure_png:
            self._create_fig(self.matrix)

        if plot_figure_html:
            self._create_html_fig(self.matrix)

        results = {
            'plot_figure_html': self.figure_html,
            'plot_figure_png': self.figure_png,
            'embedding_for_plotting': self.selected_plotting_embedding,
        }

        return results

    def _shift_data_to_positive(self, matrix):
        """
        """
        matrix_min = matrix.min()

        if matrix_min < 0:
            matrix = matrix - matrix_min

        return matrix

    def _fit_predict(self, matrix):
        """
        """
        if self.selected_clustering != 'Gaussian_Mixture':
            clusters = self.clustering.fit_predict(matrix)

        else:
            self.clustering.fit(matrix)
            clusters = self.clustering.predict(matrix)

        return clusters

    def _format_found_clusters(self):
        """
        """
        clusters = []

        self.clusters = self.clusters.tolist()

        for sample, cluster in zip(self.sample_ids, self.clusters):
            if cluster == -1 :
                clusters.append('outliers')
                self.outliers.append(sample)

            else:
                clusters.append('cluster_{0}'.format(cluster))

        clusters = {sample: cluster for sample, cluster in
                    zip(self.sample_ids, clusters)}

        self.clusters = clusters
        self.clusters_array = np.array([self.clusters[sample]
                                        for sample in self.sample_ids])

    def _compute_best_cluster(self, matrix):
        """
        """
        find_best = self._find_best_cluster[self.selected_clustering]['method']
        clusters = self._find_best_cluster[self.selected_clustering]['range']

        params = {}

        if 'params' in self._find_best_cluster[self.selected_clustering]:
            params = self._find_best_cluster[self.selected_clustering]['params']

        find_best(self.clustering, matrix, clusters, **params)

    def _create_fig(self, matrix):
        """
        """
        import matplotlib as mpl
        mpl.use('Agg')
        import pylab as plt

        f_io = BytesIO()
        cluster_set = set(self.clusters_array)

        C1, C2, color_dict = self._get_components_and_color_dict(matrix, cluster_set)

        fig, ax = plt.subplots(figsize=(float(self.figsize[0]) / 75,
                                        float(self.figsize[1]) / 75))

        for cluster in cluster_set:
            index = self.clusters_array == cluster
            ax.scatter(C1[index], C2[index], color=color_dict[cluster], label=cluster)

        ax.legend()
        ax.set_xlabel('Component 1')
        ax.set_ylabel('Component 2')

        ax.set_title('Projection: {0}    clustering algorithms: {1}    embbeding: {2}'.format(
            self.selected_plotting_embedding,
            self.selected_clustering,
            self.selected_embedding)
        )

        fig.savefig(f_io, format='png')

        if self.local_path_to_save_plots:
            try:
                fig.savefig('{0}/granatum_scatter_plot.png'.format(
                    self.local_path_to_save_plots))
                print('fig saved at: {0}/granatum_scatter_plot.png'.format(
                    self.local_path_to_save_plots))
            except Exception:
                pass

        f_io.seek(0)
        self.figure_png =base64.b64encode(f_io.read())

    def _create_html_fig(self, matrix):
        """
        """
        from bokeh.plotting import figure
        from bokeh.plotting import save
        from bokeh.resources import CDN
        from bokeh.embed import file_html

        from bokeh.models import ColumnDataSource
        from bokeh.models import HoverTool

        from jinja2 import Environment, PackageLoader, Markup


        f_io_name = '{0}/granatum_scatter_plot.html'.format(
            self.local_path_to_save_plots)

        f_io_name_tmp = './tmp_bokeh_scatter_{0}.html'.format(
            np.random.randint(0, 1000))

        cluster_set = set(self.clusters_array)
        C1, C2, color_dict = self._get_components_and_color_dict(matrix, cluster_set)

        color_list = [color_dict[cluster] for cluster in self.clusters_array]

        metadatas = defaultdict(list)

        if self.metadata:
            meta_ref = self.metadata[self.sample_ids[0]].keys()

            for sample in self.sample_ids:
                assert(meta_ref == self.metadata[sample].keys())

                for meta in self.metadata[sample]:
                    metadatas[meta].append(self.metadata[sample][meta])

        data = {
            'x': C1,
            'y': C2,
            'cluster':self.clusters_array,
            'color': color_list,
            'sample_ids': self.sample_ids
        }

        for meta in metadatas:
            data[meta] = metadatas[meta]

        source = ColumnDataSource(data=data)

        hover = HoverTool()
        hover.tooltips = [("SampleID", "@sample_ids"),
                          ("clusterID", "@cluster")]

        for meta in metadatas:
            hover.tooltips.append((meta, "@{0}".format(meta)))

        title = 'Projection: {0}    clustering algorithms: {1}    embbeding: {2}'.format(
            self.selected_plotting_embedding,
            self.selected_clustering,
            self.selected_embedding)

        fig = figure(plot_width=self.figsize[0],
                     plot_height=self.figsize[1],
                     tools=[hover, 'pan', 'box_zoom', 'reset', 'save'], title=title)

        fig.xaxis.axis_label = 'Component 1'
        fig.yaxis.axis_label = 'Component 2'

        fig.circle(x='x', y='y',
                   source=source,
                   size=7,
                   line_color='black',
                   fill_alpha=0.7,
                   color='color',
                   alpha=0.7)

        _env = Environment(loader=PackageLoader('granatum_clustering', '_templates'))
        template = _env.get_template("template.html")

        self.figure_html = file_html(fig, resources=CDN, template=template)

    def _get_components_and_color_dict(self, matrix, cluster_set):
        """
        """
        from colour import Color

        if self.selected_plotting_embedding in ['NMF']:
            matrix = self._shift_data_to_positive(matrix)

        C1, C2 = self.embedding_for_plotting.fit_transform(matrix).T

        color_list = [color.get_hex_l() for color in
                      Color('blue').range_to('red', len(cluster_set))]

        color_dict = {cluster: color for cluster, color in zip(cluster_set, color_list)}

        if 'outliers' in color_dict:
            color_dict['outliers'] = '#C8C8C8'

        return C1, C2, color_dict


def main():
    """ """
    from sklearn.datasets import make_blobs
    plot_figures_path = '/home/opoirion/code/d3visualisation/'

    test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)

    metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(test_dataset))}

    granatum = GranatumClustering(find_best_number_of_cluster=False, n_clusters=30)
    results = granatum.fit(test_dataset, metadata=metadata)
    print('number of clusters: {0}'.format(results['n_clusters']))

    plot = granatum.plot(
        embedding='NMF',
        local_path_to_save_plots=plot_figures_path)

    assert('embedding_for_plotting' in plot)


if __name__ == '__main__':
    main()
