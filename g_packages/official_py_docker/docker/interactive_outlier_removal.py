#!/usr/bin/env python

import json
from sys import stdin, stdout
import numpy as np
import scanpy.api as sc

from granatum_sdk import Granatum


def remove_from(x, y):
  try:
    y.remove(x)
  except ValueError:
    return


def main():
  gn = Granatum()

  adata = gn.ann_data_from_assay(gn.get_import('assay'))
  outliers = gn.get_arg('outliers')

  num_cells_before = adata.shape[0]

  kept_cell_ids = adata.obs_names.drop(outliers, errors='ignore').values

  adata = adata[kept_cell_ids, :]

  gn.export_statically(gn.assay_from_ann_data(adata), 'Outlier removed assay')
  gn.add_result(
    'You removed {} outliers from {} cells, the result assay has {} cells (and {} genes).'.format(
      len(outliers), num_cells_before, adata.shape[0], adata.shape[1]
    ),
    type='markdown'
  )

  gn.commit()


if __name__ == '__main__':
  main()
