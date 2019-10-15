#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

from granatum_sdk import Granatum

COLORS = ["#3891ea", "#29ad19", "#ac2d58", "#db7580", "#ed2310", "#ca2dc2", "#5f7575", "#7cc1b5", "#c3bd78", "#4ffa24"]


def main():
    gn = Granatum()
    sample_coords = gn.get_import("viz_data")
    value = gn.get_import("value")
    print(value)
    coloring_type = gn.get_arg("coloring_type")

    coords = sample_coords.get("coords")
    dim_names = sample_coords.get("dimNames")

    df = pd.DataFrame(
        {"x": [a[0] for a in coords.values()], "y": [a[1] for a in coords.values()], "value": pd.Series(value)},
        index=coords.keys()
    )

    print(df)

    if coloring_type == "categorical":
        for i, cat in enumerate(df["value"].unique()):
            dff = df[df["value"] == cat]
            plt.scatter(x=dff["x"], y=dff["y"], s=5000 / df.shape[0], c=COLORS[i % len(COLORS)], label=cat)
        lgd = plt.legend(markerscale=60 / (5000 / df.shape[0]), loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    elif coloring_type == "continuous":
        plt.scatter(x=df["x"], y=df["y"], s=5000 / df.shape[0], c=df["value"], cmap="Reds")
        plt.colorbar()

    plt.xlabel(dim_names[0])
    plt.ylabel(dim_names[1])
    # plt.tight_layout()

    gn.add_current_figure_to_results(
        "Scatter-plot",
        dpi=75,
        width=750,
        height=650,
        savefig_kwargs={'bbox_extra_artists': (lgd,), 'bbox_inches': 'tight'}
    )

    gn.commit()


if __name__ == "__main__":
    main()
