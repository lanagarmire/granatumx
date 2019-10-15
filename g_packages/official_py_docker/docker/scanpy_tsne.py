#!/usr/bin/env python

from itertools import combinations
import random

import scanpy.api as sc
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

# import pandas as pd
# import seaborn as sns


def main():
  gn = Granatum()

  adata = gn.ann_data_from_assay(gn.get_import('assay'))
  random_seed = gn.get_arg('random_seed')

  sc.tl.tsne(adata, 20, random_state=random_seed)

  X_tsne = adata.obsm['X_tsne']

  plt.figure()
  plt.scatter(X_tsne[:, 0], X_tsne[:, 1], 5000 / adata.shape[0])
  plt.xlabel('t-SNE dim. 1')
  plt.ylabel('t-SNE dim. 2')
  plt.tight_layout()
  gn.add_current_figure_to_results('t-SNE plot: each dot represents a cell', dpi=75)

  pca_export = {
    'dimNames': ['t-SNE dim. 1', 't-SNE dim. 2'],
    'coords': {sample_id: X_tsne[i, :].tolist() for i, sample_id in enumerate(adata.obs_names)},
  }
  gn.export_statically(pca_export, 't-SNE coordinates')

  gn.commit()


if __name__ == '__main__':
  main()
