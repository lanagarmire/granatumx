#!/usr/bin/env python

from sys import stdin
from sys import stdout
from os import path
import json

#from granatum_clustering.clustering_apps import GranatumClustering
from granatum_deep.deep_apps import GranatumDeepClustering
import numpy as np

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    assay = gn.get_import('assay')

    args_for_init = {
        'selected_embedding': gn.get_arg('selectedEmbedding'),
        'selected_clustering': gn.get_arg('selectedClustering'),
        'n_components': gn.get_arg('nComponents'),
        'n_clusters': gn.get_arg('nClusters'),
        'find_best_number_of_cluster': gn.get_arg('findBestNumberOfCluster'),
    }

    args_for_fit = {
        'matrix': np.transpose(np.array(assay.get('matrix'))),
        'sample_ids': assay.get('sampleIds'),
    }

    granatum_clustering = GranatumDeepClustering(**args_for_init)
    fit_results = granatum_clustering.fit(**args_for_fit)

    gn.export_statically(fit_results.get('clusters'), 'Cluster assignment')

    md_str = f"""\
## Results

  * Cluster array: `{fit_results.get('clusters_array')}`
  * Cluster array: `{fit_results.get('clusters_array')}`
  * nClusters: {fit_results.get('n_clusters')}
  * Number of components: {fit_results.get('n_components')}
  * Outliers: {fit_results.get('outliers')}"""
    # gn.add_result(md_str, 'markdown')

    gn.add_result(
        {
            'orient': 'split',
            'columns': ['Sample ID', 'Cluster Assignment'],
            'data': [[x, y] for x, y in zip(assay.get('sampleIds'), fit_results.get('clusters_array'))],
        },
        'table',
    )

    gn.commit()


if __name__ == '__main__':
    main()
