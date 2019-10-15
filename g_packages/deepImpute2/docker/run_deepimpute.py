import granatum_sdk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from deepimpute.deepImpute import deepImpute 


def main():
    gn = granatum_sdk.Granatum()

    assay = gn.get_import("assay")
    data = np.array(assay.get("matrix")).T

    seed = gn.get_arg("seed")
    checkbox = gn.get_arg("use_auto_limit")
    cell_subset = gn.get_arg("cell_subset")

    NN_lim = {False: gn.get_arg("NN_lim"), True: "auto"}.get(checkbox, True)

    #model = MultiNet(n_cores="all", seed=seed)
    #model.fit(data, NN_lim=NN_lim, cell_subset=cell_subset)

    frameddata = pd.DataFrame(data)
    imputed, model = deepImpute(frameddata, NN_lim=NN_lim, cell_subset=cell_subset)

    LABELS_PARAMS = {"fontsize": 14, "fontweight": "bold", "fontname": "serif"}

    vmax = np.percentile(np.log10(1 + data.flatten()), 99)

    print("Generating Heatmap")
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(np.log10(1 + data), aspect="auto", vmax=vmax)
    ax[1].imshow(np.log10(1 + imputed), aspect="auto", vmax=vmax)
    ax[0].set_xlabel("Genes", **LABELS_PARAMS)
    ax[1].set_xlabel("Genes", **LABELS_PARAMS)
    ax[0].set_ylabel("Cells", **LABELS_PARAMS)
    ax[0].set_title("raw (log)", **LABELS_PARAMS)
    ax[1].set_title("imputed (log)", **LABELS_PARAMS)

    gn.add_current_figure_to_results("Heatmaps")
    nb_genes = len(set(model.targets.flatten()))
    #nb_genes = np.sum([len(net.targetGenes) for net in model.targets])

    def calc_dropout(matrix):
        return np.sum((np.array(matrix) == 0)) * 1. / data.size

    r, p = model.score(frameddata)
    rows, cols = frameddata.shape

    message = "\n".join(
        [
            "  - Data frame number of rows: **{0}**",
            "  - Data frame number of columns: **{1}**",
            "  - Number of imputed genes: **{2}**",
            "  - Percentage of dropout entries *before* imputation: **{3:.2f}%**",
            "  - Percentage of dropout entries *after* imputation: **{4:.2f}%**",
            "  - Accuracy (correlation) on masked data: **{5:.2f}**"
        ]
    ).format(
        rows,
        cols,
        nb_genes,
        calc_dropout(data) * 100,
        calc_dropout(imputed.to_numpy()) * 100,
	r
    )

    gn.add_result(message, data_type="markdown")

    assay["matrix"] = imputed.T.to_numpy().tolist()
    gn.export_statically(assay, "Imputed assay")

    gn.commit()


if __name__ == "__main__":
    main()
