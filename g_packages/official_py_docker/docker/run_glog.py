#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    assay = gn.get_import('assay')
    x = np.array(assay.get('matrix'))
    log_base = gn.get_arg('log_base')
    n_top = gn.get_arg('n_top')
    n_bottom = gn.get_arg('n_bottom')
    which_mid = gn.get_arg('which_mid')

    gene_df = pd.DataFrame(
        {
            'row_num': range(x.shape[0]),
            'gene_id': assay.get('geneIds'),
            'exp_mean': np.mean(x, axis=1),
            'exp_std': np.std(x, axis=1),
        }
    )
    gene_df = gene_df.sort_values('exp_mean', ascending=False)
    top_gene_row = gene_df.head(n_top).sort_values('exp_std', ascending=False).iloc[0]
    bottom_gene_row = gene_df.tail(n_bottom).sort_values('exp_std').iloc[0]

    hk_gene = x[top_gene_row['row_num'], :]
    neg_gene = x[bottom_gene_row['row_num'], :]

    if which_mid == 'mean':
        alphabk = np.mean(neg_gene[:])
    elif which_mid == 'median':
        alphabk = np.median(neg_gene[:])
    else:
        raise ValueError()

    loghkdatabk = np.log((hk_gene - alphabk).astype(complex)) / np.log(log_base)

    c = (np.std(neg_gene[:], ddof=1) / np.std(loghkdatabk, ddof=1))**2

    xbk = x - alphabk
    transformed_matrix = np.log((xbk + np.sqrt(xbk**2 + c)) / 2) / np.log(log_base)

    gn.add_result(
        '\n'.join(
            [
                f"Selected benchmarking genes:",
                f"  * housekeeping gene: **{top_gene_row['gene_id']}** "
                f"(mean: {top_gene_row['exp_mean']}, std: {top_gene_row['exp_std']}) ",
                f"  * negative control gene: **{bottom_gene_row['gene_id']}**"
                f"(mean: {bottom_gene_row['exp_mean']}, std: {bottom_gene_row['exp_std']})",
                f"",
                f"Final formula is `y = log{log_base}((z + sqrt(z^2 + c))/2)`, where `z = x - {alphabk}` and `c = {c}`.",
            ]
        ), 'markdown'
    )

    non_zero_values_before = x.flatten()
    non_zero_values_before = non_zero_values_before[(non_zero_values_before > np.percentile(non_zero_values_before, 5))]

    non_zero_values_after = transformed_matrix.flatten()
    non_zero_values_after = non_zero_values_after[(non_zero_values_after > np.percentile(non_zero_values_after, 5))]

    plt.figure()

    plt.subplot(2, 1, 1)
    plt.title('Before glog transformation')
    plt.hist(non_zero_values_before, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.subplot(2, 1, 2)
    plt.title('After glog transformation')
    plt.hist(non_zero_values_after, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Expression level')

    plt.tight_layout()

    caption = (
        'The distribution of expression level before and after glog transformation. Only the values greater '
        'than the 5 percentile (usually zero in single-cell data) and lower than 95 percentile are considered.'
    )
    gn.add_current_figure_to_results(caption, zoom=2, dpi=50)

    assay['matrix'] = transformed_matrix.tolist()
    gn.export_statically(assay, 'GLog transformed assay')

    gn.commit()


if __name__ == '__main__':
    main()
