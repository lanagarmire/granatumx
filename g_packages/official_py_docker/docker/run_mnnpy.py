from itertools import combinations

import multiprocessing
import pandas as pd

import scanpy.api as sc
import matplotlib.pyplot as plt
import numpy as np
from mnnpy import mnn_correct

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


def main():
    gn = Granatum()

    adata = gn.ann_data_from_assay(gn.get_import('assay'))
    batch_dict = gn.get_import('batch_dict')
    num_cells_to_sample = gn.get_arg('num_cells_to_sample')
    log_trans_when_plot = gn.get_arg('log_trans_when_plot')
    k = gn.get_arg('k')
    sigma = gn.get_arg('sigma')

    if num_cells_to_sample > adata.shape[0]:
        num_cells_to_sample = adata.shape[0]

    sampled_cells_idxs = np.sort(np.random.choice(adata.shape[0], num_cells_to_sample, replace=False))

    make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
    gn.add_current_figure_to_results(
        'Before Batch-effect removal: Each bar in the box plot represents one cell. Only expression levels between the 5 and 95 percentiles (exclusive) are plotted.',
        height=350,
        dpi=75 * 40 / max(40, num_cells_to_sample)
    )

    # -- main --
    batch_vec = pd.Series(batch_dict).reindex(adata.obs_names)
    batches = [adata[batch_vec == batch_id, :] for batch_id in batch_vec.unique()]
    adata, mnn_list, angle_list = mnn_correct(*batches, k=k, sigma=sigma, index_unique=None)
    # -- main end --

    gn.add_markdown(f"""\
```
    batches: {batch_vec.unique()}
```
""")

    make_plot(adata[sampled_cells_idxs, :], log_trans=log_trans_when_plot)
    gn.add_current_figure_to_results(
        'After Batch-effect removal: Each bar in the box plot represents one cell. Only expression levels between the 5 and 95 percentiles (exclusive) are plotted.',
        height=350,
        dpi=75 * 40 / max(40, num_cells_to_sample)
    )

    gn.export_statically(gn.assay_from_ann_data(adata), 'Batch corrected assay')
    for i, mnn in enumerate(mnn_list[1:]):
        print(mnn)
        gn.export(mnn.to_csv(), f"MNN_{i}.csv", raw=True)

    gn.commit()


if __name__ == '__main__':
    main()
