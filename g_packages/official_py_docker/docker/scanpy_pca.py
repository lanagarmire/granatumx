#!/usr/bin/env python

from itertools import combinations

import multiprocessing

import scanpy.api as sc
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

# import pandas as pd
# import seaborn as sns


def to_percentage(x):
    return "%{:3.2f}".format(x * 100)


def main():
    gn = Granatum()

    adata = gn.ann_data_from_assay(gn.get_import("assay"))
    num_top_comps = gn.get_arg("num_top_comps")

    sc.pp.pca(adata, 20)

    variance_ratios = adata.uns["pca"]["variance_ratio"]
    pc_labels = ["PC{}".format(x + 1) for x in range(len(variance_ratios))]

    plt.figure()
    plt.bar(pc_labels, variance_ratios)
    plt.tight_layout()
    gn.add_current_figure_to_results("Explained variance (ratio) by each Principal Component (PC)", height=350, dpi=75)

    X_pca = adata.obsm["X_pca"]

    for i, j in combinations(range(num_top_comps), 2):
        xlabel = "PC{}".format(i + 1)
        ylabel = "PC{}".format(j + 1)

        plt.figure()
        plt.scatter(X_pca[:, i], X_pca[:, j], s=5000 / adata.shape[0])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        gn.add_current_figure_to_results("PC{} vs. PC{}".format(i + 1, j + 1), dpi=75)

        pca_export = {
            "dimNames": [xlabel, ylabel],
            "coords": {sample_id: X_pca[k, [i, j]].tolist() for k, sample_id in enumerate(adata.obs_names)},
        }
        gn.export(pca_export, "PC{} vs. PC{}".format(i + 1, j + 1), kind="sampleCoords", meta={})

    gn.commit()


if __name__ == "__main__":
    main()
