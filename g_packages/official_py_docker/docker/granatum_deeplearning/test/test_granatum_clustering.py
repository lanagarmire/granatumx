from granatum_deep.deep_apps import GranatumDeepClustering
from sklearn.datasets import make_blobs

TEST_DATASET, REF_ARRAY = make_blobs(n_samples=300, n_features=200, centers=8)


def test_1():
    """
    use of granatum_deep to embedd with DNN, clusters with classic algo and
    plot either with PCA or DNN
    """
    metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(TEST_DATASET))}

    granatum = GranatumDeepClustering(
        selected_clustering='DBSCAN', # same clustering algorithms used by granatum_clustering
        selected_embedding='autoencoder', # here a DNN is used to embed the data
    )

    granatum.fit(TEST_DATASET, metadata=metadata)

    # classic use of embedding
    granatum.plot(embedding='PCA')

    granatum.plot(
        embedding='autoencoder', # here we use the DNN embedding to visualise the data
    )

def test_2():
    """
    use of granatum_deep  clusters with DNN without DNN prior embedding
    """
    granatum = GranatumDeepClustering(
        n_clusters=8,
        selected_embedding='PCA', # here a DNN is used to embed the data
        find_best_number_of_cluster=False,
        selected_clustering='autoencoder', # It is possible to use this option but quite long
    )

    granatum.fit(TEST_DATASET)

def test_3():
    """
    use of granatum_deep  clusters with DNN with and DNN prior embedding
    and automatic nb of cluster selection
    """

    granatum = GranatumDeepClustering(
        n_clusters=8,
        selected_embedding='autoencoder', # here a DNN is used to embed the data
        find_best_number_of_cluster=True, # Prepare to wait!
        selected_clustering='autoencoder',  # It is possible to use this option but quite long
    )
    granatum._find_best_cluster['autoencoder']['range'] = range(2,5)

    granatum.fit(TEST_DATASET)


if __name__ == '__main__':
    test_1()
    test_2()
    test_3()
