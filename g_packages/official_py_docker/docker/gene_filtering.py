import math

import scanpy.api as sc
import numpy as np
from granatum_sdk import Granatum, eprint
import pandas as pd


def main():
    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import("assay"))

    epsilon = gn.get_arg('epsilon')
    min_cells_expressed = gn.get_arg('min_cells_expressed')

    filter_df = pd.DataFrame({'gene': df.index})
    filter_df['sum_expr'] = [sum(df.values[i, :]) for i in range(df.shape[0])]
    filter_df['avg_expr'] = filter_df['sum_expr'] / df.shape[1]
    filter_df['num_expressed_genes'] = [sum([x > epsilon for x in df.values[i, :]]) for i in range(df.shape[0])]
    filter_df['removed'] = filter_df['num_expressed_genes'] < min_cells_expressed

    new_df = df.loc[np.logical_not(filter_df['removed'].values), :]

    gn.add_result(
        "\n".join(
            [
                "Number of genes before filtering: **{}**".format(df.shape[0]),
                "",
                "Number of genes after filtering: **{}**".format(new_df.shape[0]),
            ]
        ),
        type="markdown",
    )

    if filter_df.shape[0] > 0:
        filter_df_deleted = filter_df.loc[filter_df['removed'].values, :].drop('removed', axis=1)
        gn.add_result(
            {
                'title': f"Removed genes ({filter_df_deleted.shape[0]})",
                'orient': 'split',
                'columns': filter_df_deleted.columns.values.tolist(),
                'data': filter_df_deleted.values.tolist(),
            },
            data_type='table',
        )
    else:
        gn.add_result(
            f"No genes were removed. All {df.shape[0]} genes were kept. "
            f"See attachment **gene_selection.csv** for detail.",
            'markdown',
        )

    gn.export(filter_df.to_csv(index=False), 'gene_selection.csv', kind='raw', meta=None, raw=True)
    gn.export(gn.assay_from_pandas(new_df), "Filtered Assay", dynamic=False)

    gn.commit()


if __name__ == "__main__":
    main()
