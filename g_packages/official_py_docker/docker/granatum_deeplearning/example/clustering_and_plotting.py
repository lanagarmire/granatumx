from granatum_deep.deep_apps import GranatumDeepClustering
from sklearn.datasets import make_blobs


def example_1():
    """
    use of granatum_deep to embedd with DNN, clusters with classic algo and
    plot either with PCA or DNN
    """
    test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)

    metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(test_dataset))}
    plot_figures_path = '/home/opoirion/code/d3visualisation/'

    granatum = GranatumDeepClustering(
        selected_clustering='DBSCAN', # same clustering algorithms used by granatum_clustering
        selected_embedding='autoencoder', # here a DNN is used to embed the data
    )

    granatum.fit(test_dataset, metadata=metadata)

    # classic use of embedding
    granatum.plot(
        embedding='PCA',
        local_path_to_save_plots=plot_figures_path)

    granatum.plot(
        embedding='autoencoder', # here we use the DNN embedding to visualise the data
        local_path_to_save_plots=plot_figures_path)

def example_2():
    """
    use of granatum_deep  clusters with DNN with and without DNN prior embedding
    """
    test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)

    granatum = GranatumDeepClustering(
        n_clusters=8,
        selected_embedding='PCA', # here a DNN is used to embed the data
        find_best_number_of_cluster=False,
        selected_clustering='autoencoder', # It is possible to use this option but quite long
    )

    granatum.fit(test_dataset)

    granatum = GranatumDeepClustering(
        n_clusters=8,
        selected_embedding='autoencoder', # here a DNN is used to embed the data
        find_best_number_of_cluster=True, # Prepare to wait!
        selected_clustering='autoencoder',  # It is possible to use this option but quite long
    )

    granatum.fit(test_dataset)


if __name__ == '__main__':
    example_1()
    example_2()
