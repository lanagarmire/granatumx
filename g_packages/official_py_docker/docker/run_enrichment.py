import math

import pandas as pd
import scanpy.api as sc
import numpy as np
from granatum_sdk import Granatum, guess_gene_id_type, convert_gene_ids, get_all_genes
from collections import Counter
import zgsea
import sys


def progress(x, *args, **kwargs):
    print(x, *args, file=sys.stderr, **kwargs)


def main():
    gn = Granatum()

    gene_scores_dict = gn.get_import("gene_scores")
    species = gn.get_arg("species")
    gset_group_id = gn.get_arg("gset_group_id")
    threshold = gn.get_arg("threshold")
    use_abs = gn.get_arg("use_abs")
    background = gn.get_arg("background")

    gene_ids = list(gene_scores_dict.keys())
    gene_scores = list(gene_scores_dict.values())

    gene_id_type = guess_gene_id_type(list(gene_ids)[:5])
    if gene_id_type != 'symbol':
        gene_ids = convert_gene_ids(gene_ids, gene_id_type, 'symbol', species)

    if species == "human":
        pass
    elif species == "mouse":
        gene_ids = zgsea.to_human_homolog(gene_ids, "mouse")
    else:
        raise ValueError()

    if use_abs:
        input_list = np.array(gene_ids)[np.abs(np.array(gene_scores)) >= threshold]
    else:
        input_list = np.array(gene_ids)[np.array(gene_scores) >= threshold]

    gn.add_result(
        f"""\
Number of genes after thresholding: {len(input_list)} (out of original {len(gene_ids)}).

Please see the attachment `list_of_genes.csv` for the list of genes considered in this enrichment analysis.""",
        'markdown',
    )

    gn.export(pd.Series(input_list).to_csv(index=False), 'list_of_genes.csv', kind='raw', meta=None, raw=True)

    if background == 'all':
        background_list = get_all_genes('human')
    elif background == 'from_gene_sets':
        background_list = None
    elif background == 'from_input':
        background_list = gene_ids
    else:
        raise ValueError()

    result_df = zgsea.simple_fisher(input_list, gset_group_id, background_list=background_list)
    result_df = result_df.sort_values('fdr')
    result_df = result_df[[
        'gene_set_name',
        'size',
        'p_val',
        'fdr',
        'odds_ratio',
        'n_overlaps',
        'overlapping_genes',
    ]]
    result_df.columns = [
        'Gene set',
        'Gene set size',
        'p-value',
        'FDR',
        'Odds ratio',
        'Number of overlapping genes',
        'Overlapping genes',
    ]

    gn.add_pandas_df(result_df)
    gn.export(result_df.to_csv(index=False), 'enrichment_results.csv', kind='raw', meta=None, raw=True)

    gn.commit()


if __name__ == "__main__":
    main()
