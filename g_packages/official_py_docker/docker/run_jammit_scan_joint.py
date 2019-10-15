import datetime
import json

import numpy as np
from granatum_sdk import Granatum, eprint
from numpy import ones, concatenate, linspace, tile, shape
from numpy.random import standard_normal, choice
import matplotlib.pyplot as plt

# from jammit.efdr_ssvdanyN import efdr_ssvdanyN
from jammit import JAMMIT

import pandas as pd


def main():
    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    track_dict = gn.get_import('track_dict')
    n_steps = gn.get_arg('n_steps')
    min_theta = gn.get_arg('min_theta')
    max_theta = gn.get_arg('max_theta')

    track_vec = np.array([track_dict.get(k) for k in df.index])
    jammit = JAMMIT.from_df_and_track_vec(df, track_vec)

    jammit.scan(
        thetas=np.linspace(min_theta, max_theta, n_steps),
        calculate_fdr=True,
        n_perms=10,
        verbose=1,
        convergence_threshold=0.000000001,
    )

    jammit_result = jammit.format(show_individual_tracks=True)
    # jammit_result['theta'] = jammit_result['theta'].round(3)
    # jammit_result['alpha'] = jammit_result['alpha'].round(3)
    eprint(f"jammit_result = {jammit_result}")

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

    gn.export(jammit_result.to_csv(index=False), 'jammit_scan_result.csv', kind='raw', meta=None, raw=True)

    gn.commit()


if __name__ == '__main__':
    main()
