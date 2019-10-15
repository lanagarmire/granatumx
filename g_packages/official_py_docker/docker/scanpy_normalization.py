from itertools import combinations

import multiprocessing

import scanpy.api as sc
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

# import pandas as pd
# import seaborn as sns

nans = np.array([np.nan, np.nan])


def trim_extreme(x, a, b):
    low = np.percentile(x, a)
    high = np.percentile(x, b)
    filtered = x[(x > low) & (x < high)]
    return filtered.copy()


def make_plot(adata, log_trans=False):
    violin_data = []
    for cell in adata.X:
        filtered = trim_extreme(cell, 5, 95)
        if log_trans:
            filtered = np.log1p(filtered)
        if filtered.shape[0] == 0:
            filtered = nans
        violin_data.append(filtered)

    plt.figure()
    plt.boxplot(violin_data)
    plt.xlabel('Cells')
    plt.ylabel('Expression lvl (log transformed)')
    plt.tight_layout()


def quantile_normalization(mat):
    # double argsort for getting the corresponding ranks for
    # each element in the vector
    rank_mat = np.argsort(np.argsort(mat, 1), 1)
    medians = np.median(np.sort(mat, 1), 0)
    normalized = np.zeros_like(mat)

    for i in range(rank_mat.shape[0]):
        normalized[i, :] = medians[rank_mat[i, :]]

    return normalized


def main():
    gn = Granatum()

    adata = gn.ann_data_from_assay(gn.get_import('assay'))
    num_cells_to_sample = gn.get_arg('num_cells_to_sample')
    method = gn.get_arg('method')
    log_trans_when_plot = gn.get_arg('log_trans_when_plot')

    if num_cells_to_sample > adata.shape[0]:
        num_cells_to_sample = adata.shape[0]

    sampled_cells_idxs = np.sort(np.random.choice(adata.shape[0], num_cells_to_sample, replace=False))

    make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
    gn.add_current_figure_to_results(
        'Before normalization: Each bar in the box plot represents one cell. Only expression levels between the 5 and 95 percentiles (exclusive) are plotted.',
        height=350,
        dpi=75 * 40 / max(40, num_cells_to_sample)
    )

    if method == 'quantile':
        adata.X = quantile_normalization(adata.X)
    elif method == 'scanpy':
        sc.pp.normalize_per_cell(adata)
    else:
        raise ValueError()

    make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
    gn.add_current_figure_to_results(
        'After normalization: Each bar in the box plot represents one cell. Only expression levels between the 5 and 95 percentiles (exclusive) are plotted.',
        height=350,
        dpi=75 * 40 / max(40, num_cells_to_sample)
    )

    gn.export_statically(gn.assay_from_ann_data(adata), 'Normalized assay')

    gn.commit()


if __name__ == '__main__':
    main()
