#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

from granatum_sdk import Granatum

COLORS = ["#3891ea", "#29ad19", "#ac2d58", "#db7580", "#ed2310", "#ca2dc2", "#5f7575", "#7cc1b5", "#c3bd78", "#4ffa24"]


def main():
    gn = Granatum()

    sample_coords = gn.get_import("viz_data")
    df = gn.pandas_from_assay(gn.get_import("assay"))
    gene_id = gn.get_arg("gene_id")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")

    scatter_df = pd.DataFrame(
        {"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": df.loc[gene_id, :]},
        index=coords.keys(),
    )

    plt.scatter(x=scatter_df["x"], y=scatter_df["y"], s=5000 / scatter_df.shape[0], c=scatter_df["value"], cmap="Reds")
    plt.colorbar()

    plt.xlabel(dim_names[0])
    plt.ylabel(dim_names[1])
    plt.tight_layout()

    gn.add_current_figure_to_results("Scatter-plot", dpi=75)

    gn.commit()


if __name__ == "__main__":
    main()
