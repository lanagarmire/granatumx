import numpy as np

from granatum_sdk import Granatum
import matplotlib.pyplot as plt
import numpy as np


def plot_distribution_comparison(m1, m2, gn):
    non_zero_values_before = m1.flatten()
    non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]

    non_zero_values_after = m2.flatten()
    non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

    plt.figure()

    plt.subplot(2, 1, 1)
    plt.title('Before gene centering')
    plt.hist(non_zero_values_before, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.subplot(2, 1, 2)
    plt.title('After gene centering')
    plt.hist(non_zero_values_after, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.tight_layout()

    caption = (
        'The distribution of expression level before and after gene centering. Only the values greater '
        'than the 5 percentile (usually zero in single-cell data) and lower than 95 percentile are considered.'
    )
    gn.add_current_figure_to_results(caption, zoom=2, dpi=50)


def main():
    gn = Granatum()

    assay = gn.get_import('assay')

    matrix = np.array(assay.get('matrix'))

    transformed_matrix = matrix - matrix.mean(axis=1, keepdims=True)
    assay['matrix'] = transformed_matrix.tolist()

    plot_distribution_comparison(matrix, transformed_matrix, gn)

    gn.export_statically(assay, 'Gene centered assay')

    gn.commit()


if __name__ == '__main__':
    main()
