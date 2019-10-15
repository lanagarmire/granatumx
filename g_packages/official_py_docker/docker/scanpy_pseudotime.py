#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scanpy.api as sc

from granatum_sdk import Granatum

def main():
  gn = Granatum()

  n_neighbors = gn.get_arg('nNeighbors', 15)
  neighbor_method = gn.get_arg('neighborMethod', 'gauss')

  assay = gn.get_import('assay')

  adata = sc.AnnData(np.array(assay.get('matrix')).transpose())
  adata.var_names = assay.get('geneIds')
  adata.obs_names = assay.get('sampleIds')

  sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X', method=neighbor_method)
  sc.tl.dpt(adata, n_branchings=1)

  gn._pickle(adata, 'adata')

  # dpt_groups

  for spec in [{'col': 'dpt_order', 'caption': 'Cell order'}, {'col': 'dpt_groups', 'caption': 'Cell groups'}]:
    plt.figure()
    sc.pl.diffmap(adata, color=spec['col'])
    gn.add_current_figure_to_results(spec['caption'])
    gn.export_statically(dict(zip(adata.obs_names.tolist(), adata.obs[spec['col']].values.tolist())), spec['col'])

  gn.commit()


if __name__ == '__main__':
  main()
