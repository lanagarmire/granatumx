#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    assay = gn.get_import('assay')
    matrix = np.array(assay.get('matrix'))

    log_base = gn.get_arg('logBase')
    pseudo_counts = gn.get_arg('pseudoCounts')

    transformed_matrix = np.log(matrix + pseudo_counts) / np.log(log_base)

    non_zero_values_before = matrix.flatten()
    # non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5)) &
    #                                                 (non_zero_values_before < np.percentile(non_zero_values_before, 95))]
    non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]

    non_zero_values_after = transformed_matrix.flatten()
    # non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5)) &
    #                                               (non_zero_values_after < np.percentile(non_zero_values_after, 95))]
    non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

    plt.figure()

    plt.subplot(2, 1, 1)
    plt.title('Before log transformation')
    plt.hist(non_zero_values_before, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.subplot(2, 1, 2)
    plt.title('After log transformation')
    plt.hist(non_zero_values_after, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.tight_layout()

    caption = (
        'The distribution of expression level before and after log transformation. Only the values greater '
        'than the 5 percentile (usually zero in single-cell data) and lower than 95 percentile are considered.'
    )
    gn.add_current_figure_to_results(caption, zoom=2, dpi=50)

    assay['matrix'] = transformed_matrix.tolist()
    gn.export_statically(assay, 'Log transformed assay')

    gn.commit()


if __name__ == '__main__':
    main()
