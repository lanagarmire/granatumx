#!/usr/bin/env python

import scanpy.api as sc

import numpy as np

import matplotlib.pyplot as plt

from granatum_sdk import Granatum


def main():
  gn = Granatum()

  adata = gn.ann_data_from_assay(gn.get_import('assay'))
  sample_coords = gn.get_import('sampleCoords')
  random_seed = gn.get_arg('random_seed')

  sc.pp.neighbors(adata, n_neighbors=20, use_rep='X', method='gauss')
  sc.tl.louvain(adata, random_state=random_seed)

  cluster_assignment = dict(zip(adata.obs_names, ['Cluster {}'.format(int(c) + 1) for c in adata.obs['louvain']]))
  gn.export_statically(cluster_assignment, 'Cluster assignment')

  dim_names = sample_coords.get('dimNames')
  coords_dict = sample_coords.get('coords')

  plt.figure()
  clusters = adata.obs['louvain'].cat.categories
  for c in clusters:
    cell_ids = adata.obs_names[adata.obs['louvain'] == c]
    coords = [coords_dict.get(x) for x in cell_ids]
    coords_x = [x[0] for x in coords]
    coords_y = [x[1] for x in coords]
    plt.scatter(coords_x, coords_y, label='Cluster {}'.format(int(c) + 1))

  plt.xlabel(dim_names[0])
  plt.ylabel(dim_names[1])
  plt.legend()
  plt.tight_layout()

  gn.add_current_figure_to_results(
    'Scatter-plot using imported cell coordinates. Each dot represents a cell. The colors indicate the indentified cell clusters.',
    dpi=75
  )

  gn.commit()


if __name__ == '__main__':
  main()
