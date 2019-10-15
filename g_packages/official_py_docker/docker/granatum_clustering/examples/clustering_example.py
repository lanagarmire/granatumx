from granatum_clustering.clustering_apps import GranatumClustering


def main():
    """ """
    from sklearn.datasets import make_blobs

    test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)
    metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(test_dataset))}
    sample_ids = list(range(len(test_dataset)))

    granatum = GranatumClustering(
        selected_embedding='PCA',
        selected_clustering='WARD')

    results = granatum.fit(matrix=test_dataset,
                           metadata=metadata,
                           sample_ids=sample_ids)

    assert('clusters' in results)
    assert('n_clusters' in results)
    assert('n_components' in results)
    assert('embedding' in results)
    assert('clustering_algorithm' in results)

    plots = granatum.plot(
        figsize=(600, 600),
        embedding='PCA',
        plot_figure_png=True,
        plot_figure_html=True,
        jsonify=False)

    assert('plot_figure_html' in plots)
    assert('plot_figure_png' in plots)


if __name__ == '__main__':
    main()
