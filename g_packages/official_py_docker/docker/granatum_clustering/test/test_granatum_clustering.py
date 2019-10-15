from granatum_clustering.clustering_apps import GranatumClustering
from granatum_clustering.clustering_apps import EMBEDDING
from granatum_clustering.clustering_apps import PLOTTING_EMBEDDING
from granatum_clustering.clustering_apps import CLUSTERING


from sklearn.datasets import make_blobs

TEST_DATASET, REF_ARRAY = make_blobs(n_samples=100, n_features=10, centers=8)


def test_0_default():
    """
    test default parameter on the default test dataset and the plotting
    """
    granatum = GranatumClustering(n_components=5, selected_clustering='DBSCAN')
    results = granatum.fit(TEST_DATASET)

    assert('clusters' in results)
    assert('n_clusters' in results)
    assert('n_components' in results)
    assert('embedding' in results)
    assert('clustering_algorithm' in results)

    plots = granatum.plot(plot_figure_png=True,
                          plot_figure_html=True,
                          embedding='PCA',
                          jsonify=False)

    assert('embedding_for_plotting' in plots)
    assert('plot_figure_html' in plots)
    assert('plot_figure_png' in plots)


def test_2_default():
    """
    test all algorithm for embedding
    """
    for embedd in EMBEDDING:
        print('#### embedding tested: {0}'.format(embedd))

        granatum = GranatumClustering(
            n_components=5,
            find_best_number_of_cluster=False,
            selected_embedding=embedd)

        results = granatum.fit(TEST_DATASET)

        assert('clusters' in results)
        assert('n_clusters' in results)
        assert('n_components' in results)
        assert('embedding' in results)
        assert('clustering_algorithm' in results)

def test_3_default():
    """
    test all clustering_algorithm
    """
    for clustering in CLUSTERING:
        print('#### clustering algo tested: {0}'.format(clustering))

        granatum = GranatumClustering(
            n_components=5,
            selected_plotting_embedding='PCA',
            selected_clustering=clustering,
            find_best_number_of_cluster=False)
        results = granatum.fit(TEST_DATASET)

        assert('clusters' in results)
        assert('n_clusters' in results)
        assert('n_components' in results)
        assert('embedding' in results)
        assert('clustering_algorithm' in results)

def test_4_default():
    """
    test all clustering algorithm with optimum cluster selection
    """
    for clustering in CLUSTERING:
        print('#### clustering algo tested: {0}'.format(clustering))

        granatum = GranatumClustering(
            n_components=5,
            selected_clustering=clustering)

        results = granatum.fit(TEST_DATASET)
        print('number of clusters found: {0}'.format(results['n_clusters']))

        assert('clusters' in results)
        assert('n_clusters' in results)
        assert('n_components' in results)
        assert('embedding' in results)
        assert('clustering_algorithm' in results)

def test_5_default():
    """
    test all plotting embbeding algorithm
    """
    granatum = GranatumClustering(
        n_components=5)
    granatum.fit(TEST_DATASET)

    for embbed in PLOTTING_EMBEDDING:
        print('#### plotting embbeding used tested: {0}'.format(embbed))

        plots = granatum.plot(embedding=embbed)

        assert('embedding_for_plotting' in plots)
        assert('plot_figure_html' in plots)
        assert('plot_figure_png' in plots)


if __name__ == '__main__':
    test_0_default()
    test_2_default()
    test_3_default()
    test_4_default()
    test_5_default()
