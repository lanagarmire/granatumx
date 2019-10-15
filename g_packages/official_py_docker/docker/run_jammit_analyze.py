import datetime
import json

import numpy as np
from granatum_sdk import Granatum
from numpy import ones, concatenate, linspace, tile, shape
from numpy.random import standard_normal, choice
import matplotlib.pyplot as plt
import seaborn as sns

# from jammit.efdr_ssvdanyN import efdr_ssvdanyN
from jammit import JAMMIT

import pandas as pd

EPSILON = 1e-7


def main():
    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    alpha = gn.get_arg('alpha')

    jammit = JAMMIT.from_dfs([df])

    res = jammit.run_for_one_alpha(
        alpha,
        verbose=1,
        convergence_threshold=0.000000001,
    )

    u = res['u']
    v = res['v']

    gn.export(dict(zip(df.index, u)), 'Genes loadings', kind='geneMeta')
    gn.export(dict(zip(df.columns, v)), 'Sample scores', kind='sampleMeta')

    gene_df = pd.DataFrame({'id_': df.index, 'abs_loading': abs(u), 'loading': u})
    gene_df = gene_df[['id_', 'abs_loading', 'loading']]
    gene_df = gene_df.loc[gene_df['loading'].abs() > EPSILON]
    gene_df = gene_df.sort_values('abs_loading', ascending=False)

    gn.add_result(
        {
            'title': f"Signal genes ({len(gene_df)})",
            'orient': 'split',
            'columns': gene_df.columns.values.tolist(),
            'data': gene_df.values.tolist(),
        },
        data_type='table',
    )
    gn.export(gene_df.to_csv(index=False), 'signal_genes.csv', kind='raw', meta=None, raw=True)

    sample_df = pd.DataFrame({'id_': df.columns, 'abs_score': abs(v), 'score': v})
    sample_df = sample_df[['id_', 'abs_score', 'score']]
    sample_df = sample_df.loc[sample_df['score'].abs() > EPSILON]
    sample_df = sample_df.sort_values('abs_score', ascending=False)

    gn.add_result(
        {
            'title': f"Signal samples ({len(sample_df)})",
            'orient': 'split',
            'columns': sample_df.columns.values.tolist(),
            'data': sample_df.values.tolist(),
        },
        data_type='table',
    )
    gn.export(sample_df.to_csv(index=False), 'signal_samples.csv', kind='raw', meta=None, raw=True)

    subset_df = df.loc[gene_df['id_'], sample_df['id_']]
    gn.export(gn.assay_from_pandas(subset_df), 'Assay with only signal genes and samples', kind='assay')

    sns.clustermap(subset_df, cmap='RdBu')
    gn.add_current_figure_to_results(
        description='Cluster map of the signal genes and signal samples',
        zoom=2,
        width=750,
        height=850,
        dpi=50,
    )
    plt.close()

    plt.figure()
    plt.scatter(range(len(u)), u, s=2, c='red')
    plt.xlabel('index')
    plt.ylabel('value in u')
    gn.add_current_figure_to_results(
        description='The *u* vector (loadings for genes) plotted as a scatter plot.',
        zoom=2,
        width=750,
        height=450,
        dpi=50,
    )
    plt.close()

    plt.figure()
    plt.plot(range(len(v)), v)
    plt.scatter(range(len(v)), v, s=6, c='red')
    plt.xlabel('index')
    plt.ylabel('value in v')
    gn.add_current_figure_to_results(
        description='The *v* vector (scores for samples) plotted as a line plot.',
        zoom=2,
        width=750,
        height=450,
        dpi=50,
    )
    plt.close()

    # gn.export_current_figure(
    #     'cluster_map.pdf',
    #     zoom=2,
    #     width=750,
    #     height=850,
    #     dpi=50,
    # )

    gn.commit()


if __name__ == '__main__':
    main()
