import datetime
import json

import numpy as np
from granatum_sdk import Granatum
from numpy import ones, concatenate, linspace, tile, shape
from numpy.random import standard_normal, choice
import matplotlib.pyplot as plt

# from jammit.efdr_ssvdanyN import efdr_ssvdanyN
from jammit import JAMMIT

import pandas as pd


def main():
    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    n_steps = gn.get_arg('n_steps')
    min_theta = gn.get_arg('min_theta')
    max_theta = gn.get_arg('max_theta')

    jammit = JAMMIT.from_dfs([df])

    jammit.scan(
        thetas=np.linspace(min_theta, max_theta, n_steps),
        calculate_fdr=True,
        n_perms=10,
        verbose=1,
        convergence_threshold=0.000000001,
    )

    jammit_result = jammit.format(columns=['theta', 'alpha', 'n_sigs', 'fdr'])
    jammit_result['theta'] = jammit_result['theta'].round(3)
    jammit_result['alpha'] = jammit_result['alpha'].round(3)

    plt.plot(jammit_result['alpha'], jammit_result['fdr'])
    plt.xlabel('alpha')
    plt.ylabel('FDR')
    gn.add_current_figure_to_results('FDR plotted against alpha', height=400)

    gn.add_result(
        {
            'pageSize': n_steps,
            'orient': 'split',
            'columns': [{'name': h, 'type': 'number', 'round': 3} for h in jammit_result.columns],
            'data': jammit_result.values.tolist(),
        },
        data_type='table',
    )

    gn.commit()


if __name__ == '__main__':
    main()
