#!/usr/bin/env python

import math
import random

import scanpy.api as sc
import numpy as np
from granatum_sdk import Granatum


def main():
    gn = Granatum()

    adata = gn.ann_data_from_assay(gn.get_import("assay"))
    num_cells_to_sample = gn.get_arg("num_cells_to_sample")
    random_seed = gn.get_arg("random_seed")

    np.random.seed(random_seed)

    num_cells_before = adata.shape[0]
    num_genes_before = adata.shape[1]

    if num_cells_to_sample > 0 and num_cells_to_sample < 1:
        num_cells_to_sample = round(num_cells_before * num_cells_to_sample)
    else:

        num_cells_to_sample = round(num_cells_to_sample)

    if num_cells_to_sample > num_cells_before:
        num_cells_to_sample = num_cells_before

    if num_cells_to_sample < 1:
        num_cells_to_sample = 1

    sampled_cells_idxs = np.sort(np.random.choice(num_cells_before, num_cells_to_sample, replace=False))

    adata = adata[sampled_cells_idxs, :]

    gn.add_result(
        "\n".join(
            [
                "The assay before down-sampling has **{}** cells and {} genes.".format(
                    num_cells_before, num_genes_before
                ),
                "",
                "The assay after down-sampling has **{}** cells and {} genes.".format(adata.shape[0], adata.shape[1]),
            ]
        ),
        type="markdown",
    )

    gn.export(gn.assay_from_ann_data(adata), "Down-sampled Assay", dynamic=False)

    gn.commit()


if __name__ == "__main__":
    main()
