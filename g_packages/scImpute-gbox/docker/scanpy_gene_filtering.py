#!/usr/bin/env python

import math

import scanpy.api as sc
import numpy as np
from granatum_sdk import Granatum


def main():
    gn = Granatum()

    adata = gn.ann_data_from_assay(gn.get_import("assay"))
    min_cells_expressed = gn.get_arg("min_cells_expressed")
    min_mean = gn.get_arg("min_mean")
    max_mean = gn.get_arg("max_mean")
    min_disp = gn.get_arg("min_disp")
    max_disp = gn.get_arg("max_disp")

    num_genes_before = adata.shape[1]

    sc.pp.filter_genes(adata, min_cells=min_cells_expressed)

    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, min_mean=math.log(min_mean), max_mean=math.log(max_mean), min_disp=min_disp, max_disp=max_disp
    )
    adata = adata[:, filter_result.gene_subset]

    sc.pl.filter_genes_dispersion(filter_result)
    gn.add_current_figure_to_results(
        "Each dot represent a gene. The gray dots are the removed genes. The x-axis is log-transformed.",
        zoom=3,
        dpi=50,
        height=400,
    )

    gn.add_result(
        "\n".join(
            [
                "Number of genes before filtering: **{}**".format(num_genes_before),
                "",
                "Number of genes after filtering: **{}**".format(adata.shape[1]),
            ]
        ),
        type="markdown",
    )

    gn.export(gn.assay_from_ann_data(adata), "Filtered Assay", dynamic=False)

    gn.commit()


if __name__ == "__main__":
    main()
