#!/usr/bin/env python

import multiprocessing

from MulticoreTSNE import MulticoreTSNE as TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

# import pandas as pd
# import seaborn as sns

def main():
  gn = Granatum()

  assay = gn.get_import('assay')

  matrix = np.array(assay.get('matrix'))
  sample_ids = assay.get('sampleIds')

  num_samples = matrix.shape[1]

  # ---- PCA --------------------------------------------------------------------

  X = np.transpose(matrix)
  model = PCA(n_components=2)
  Y_pca = model.fit_transform(X)

  pca_export = {
    'dimNames': ['PCA-1', 'PCA-2'],
    'coords':
    {sample_id: Y_pca[i, :].tolist()
     for i, sample_id in enumerate(sample_ids)},
  }
  gn.export_statically(pca_export, 'pca')

  plt.figure()
  plt.scatter(Y_pca[:, 0], Y_pca[:, 1], 5000/num_samples)
  plt.tight_layout()

  gn.add_current_figure_to_results('Principal Component Analysis (PCA) scatter-plot', dpi=75)

  # ---- T-SNE ------------------------------------------------------------------

  X = np.transpose(matrix)
  model = TSNE(n_jobs=multiprocessing.cpu_count())
  Y_tsne = model.fit_transform(X)

  tsne_export = {
    'dimNames': ['tSNE-1', 'tSNE-2'],
    'coords': {
      sample_id: Y_tsne[i, :].tolist()
      for i, sample_id in enumerate(sample_ids)
    },
  }
  gn.export_statically(tsne_export, 'tsne')

  plt.figure()
  plt.scatter(Y_tsne[:, 0], Y_tsne[:, 1], s=5000/num_samples)
  plt.tight_layout()

  gn.add_current_figure_to_results('t-Distributed Stochastic Neighbor Embedding (t-SNE) scatter-plot', dpi=75)

  gn.commit()


if __name__ == '__main__':
  main()
