# Clustering module for GRANATUM

clustering module for granatum encoding results in json and accepting inputs as json parsed from the terminal

## Requirements
* Linux working environment
* [python 2 (>=2.7)](https://www.python.org/download/releases/2.7.2/)
* Python libraries (automatically installed with the pip install command):
  * [Scikit-learn](http://scikit-learn.org/) (version = 0.19)
  * colour
  * matplotlib

## installation (local)

Tested using python2(.7) but should work also using python3

```bash
git clone https://github.com/lanagarmire/granatum_clustering.git
cd granatum_clustering
pip2 install -r requirements.txt --user # python 2.7.X must be used
```

## test

```bash
pip2 install nose --user
pip3 install nose --user
nosetests-2.7 -v
nosetests-3.X -v
```

## config
* The meta configs (i.e. classifiers used and specific classifier parameters) are defined into granatum_clustering/granatum_config.py

* The parameters that can be provided by the user are defined in granatum_clustering/clustering_apps.py

## parameters
    CURRENT CLUSTERING ALGORITHM :      kmeans, WARD, DBSCAN, HDBSCAN,
                                        Agglomerative_correlation, Gaussian_Mixture
    CURRENT EMBEDDING ALGORITHM :       PCA, IncrementalPCA, FactorAnalysis, FastICA,
                                        NMF, LatentDirichletAllocation, MDS, None



        Parameters
    ----------
    selected_embedding : str <'PCA': default>
    selected_clustering : str <'kmeans': default>
    n_components : int <30: default>
    n_clusters : int  <2: default>    only used if find_best_number_of_cluster
                                      is set to False
    plot_figure: Bool <True: default>                    plot figure in png
    plot_figure_html: Bool <True: default>               plot figure in html
    find_best_number_of_cluster: Bool <True: default>    automatically find the optimal
                                                         number of clusters

    Return (when fit is called)
    ----------
    json : {
        'clusters' :        cluster dictionnary
        'clusters_array' :        cluster array
        'n_clusters' :      number of clusters
        'n_components':     number of components used for the embbedding
        'embedding':        embedding algorithm used
        'outliers':         outliers array
        'clustering_algorithm':      clustering algorithm used
     }

## usage
* Using python

```python
from granatum_clustering.clustering_apps import GranatumClustering
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

```

## Developpers

* main developper: Olivier Poirion (o.poirion@gmail.com)
* frontend/integration developper: Xun Zhu (zhuxun2@gmail.com )