#!/usr/bin/env python

import scanpy.api as sc

import numpy as np
import pandas as pd
import math

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    assay = gn.get_import('assay')
    sample_ids = assay.get('sampleIds')
    group_dict = gn.get_import('groupVec')
    group_vec = pd.Categorical([group_dict.get(x) for x in sample_ids])
    num_groups = len(group_vec.categories)
    figheight = 400 * (math.floor((num_groups - 1) / 7) + 1)

    adata = sc.AnnData(np.array(assay.get('matrix')).transpose())
    adata.var_names = assay.get('geneIds')
    adata.obs_names = assay.get('sampleIds')
    adata.obs['groupVec'] = group_vec

    sc.pp.neighbors(adata, n_neighbors=20, use_rep='X', method='gauss')
    sc.tl.rank_genes_groups(adata, 'groupVec', n_genes=100000)
    sc.pl.rank_genes_groups(adata, n_genes=20)
    gn.add_current_figure_to_results('One-vs-rest marker genes', dpi=50, height=figheight)

    gn._pickle(adata, 'adata')

    rg_res = adata.uns['rank_genes_groups']

    for group in rg_res['names'].dtype.names:
        genes_names = [str(x[group]) for x in rg_res['names']]
        scores = [float(x[group]) for x in rg_res['scores']]
        gn.export(dict(zip(genes_names, scores)), 'Marker score ({} vs. rest)'.format(group), kind='geneMeta')

    # cluster_assignment = dict(zip(adata.obs_names, adata.obs['louvain'].values.tolist()))
    # gn.export_statically(cluster_assignment, 'cluster_assignment')

    gn.commit()


if __name__ == '__main__':
    main()
