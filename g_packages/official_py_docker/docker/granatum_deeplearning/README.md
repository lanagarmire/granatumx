# Deep-learning module for GRANATUM

Clustering and embedding of the data using Deep neural-networks
This package use the original `granatum_clustering` package and all the clustering
algorithms can be used. Moreover, The `autoencoder` algorithm can be used to cluster
or embed, prior to clustering, the data.

## Requirements
* Linux working environment
* [python 2 (>=2.7)](https://www.python.org/download/releases/2.7.2/)
* Python libraries (automatically installed with the pip install command):
  * [Scikit-learn](http://scikit-learn.org/) (version = 0.19)
  * theano (tested with version 0.8.2)
  * keras  (tested with version 2.1.3)

## installation (local)

Tested using python2(.7)

```bash
git clone https://github.com/lanagarmire/granatum_deeplearning.git
cd granatum_deeplearning
pip2 install -r requirements.txt --user # python 2.7.X must be used
pip2 install -e . --user # python 2.7.X must be used
nano ~/.keras/keras.json # edit keras.json
```
to use theano, keras.json should be configured as bellow:

```json
{
    "image_dim_ordering": "th",
    "epsilon": 1e-07,
    "floatx": "float32",
    "backend": "theano"
}

```

## test

```bash
pip2 install nose --user
nosetests-2.7 -v

```

## config
* The extensive neural network parameters can be found in the deep_config.py file

```python
## Dimensions of the intermediate layers before and after the middle hidden layer
# level on the dims of the hidden layers BEFORE new dim
LEVEL_DIMS_IN = [50]
# Number of nodes in the middle hidden layer
# (i.e. the new dimensions of the transformed data)
N_COMPONENTS = 100
# level on the dims of the hidden layers AFTER new dim
LEVEL_DIMS_OUT = [50]
# Percentage of edges being dropout at each training iteration (None for no dropout)
DROPOUT = 0.5
# L2 Regularization constant on the node activity
ACT_REG = False
# L1 Regularization constant on the weight
W_REG = False
# Fraction of the dataset to be used as test set when building the autoencoder
DATA_SPLIT = None
# activation function
ACTIVATION = 'relu'
# Number of epoch
EPOCHS = 10
# Loss function to minimize
LOSS = 'binary_crossentropy'
# Optimizer (sgd for Stochastic Gradient Descent)
OPTIMIZER = 'adam'
# train/test split
DATA_SPLIT = None
# Random seed
SEED = None
```

* The parameters that can be provided by the user are defined in granatum_clustering/clustering_apps.py

## parameters
    CURRENT CLUSTERING ALGORITHM :      kmeans, WARD, DBSCAN, HDBSCAN, autoencoder
                                        Agglomerative_correlation, Gaussian_Mixture
    CURRENT EMBEDDING ALGORITHM :       autoencoder, PCA, IncrementalPCA, FactorAnalysis,
                                        NMF, LatentDirichletAllocation, MDS, FastICA, None

        Parameters
    ----------
    selected_clustering : str <'kmeans': default>        Selected clustering
    n_components : int <30: default>                     Hidden layer in the autoencoder network
    n_clusters : int  <2: default>                       only used if find_best_number_of_cluster
                                                         is set to False
    epochs : int <10: default>                           Number of epoch used bu the DLLs
    dropout : 0.0 < float < 1.0  <0.5: default>          Number of epoch used bu the DLLs

    find_best_number_of_cluster: Bool <True: default>    automatically find the optimal
                                                         number of clusters
    level_dims_in: list(Int) <default: [50]>             Defines the levels and the dimensions of
                                                         the hidden layers before the bottleneck
    level_dims_out: list(Int) <default: [50]>             Defines the levels and the dimensions of
                                                         the hidden layers after the bottleneck

## usage
* Using python

```python
from granatum_deep.deep_apps import GranatumDeepClustering
from sklearn.datasets import make_blobs

test_dataset, ref_array = make_blobs(n_samples=300, n_features=200, centers=8)

metadata = {i: {'dummy':'dummy:{0}'.format(i)} for i in range(len(test_dataset))}
plot_figures_path = '/home/opoirion/code/d3visualisation'

granatum = GranatumDeepClustering(
    selected_clustering='DBSCAN', # same clustering algorithms used by granatum_clustering
    selected_embedding='autoencoder', # here a DNN is used to embed the data
    find_best_number_of_cluster=True,
    clustering_with_autoencoder=False
    )

granatum.fit(test_dataset, metadata=metadata)

# classic use of embedding
granatum.plot(
    embedding='PCA',
    local_path_to_save_plots=plot_figures_path)
```

## Developpers

* main developper: Olivier Poirion (o.poirion@gmail.com)
* frontend/integration developper: Xun Zhu (zhuxun2@gmail.com )