#!/usr/bin/env python

from itertools import combinations
import random

import scanpy.api as sc
import matplotlib.pyplot as plt
import numpy as np

from granatum_sdk import Granatum

import pandas as pd
import seaborn as sns


def main():
    gn = Granatum()

    tb1 = gn.pandas_from_assay(gn.get_import('assay1'))
    tb2 = gn.pandas_from_assay(gn.get_import('assay2'))
    label1 = gn.get_arg('label1')
    label2 = gn.get_arg('label2')
    direction = gn.get_arg('direction')
    normalization = gn.get_arg('normalization')

    if direction == 'samples':
        tb1 = tb1.T
        tb2 = tb2.T

    overlapped_index = set(tb1.index) & set(tb2.index)
    tb1.index = [f"{label1}_{x}" if x in overlapped_index else x for x in tb1.index]
    tb2.index = [f"{label2}_{x}" if x in overlapped_index else x for x in tb2.index]

    if normalization == 'none':
        tb = pd.concat([tb1, tb2], axis=0)
    elif normalization == 'frobenius':
        ntb1 = np.linalg.norm(tb1)
        ntb2 = np.linalg.norm(tb2)
        ntb = np.mean([ntb1, ntb2])
        fct1 = ntb / ntb1
        fct2 = ntb / ntb2
        tb = pd.concat([tb1 * fct1, tb2 * fct2], axis=0)
        gn.add_markdown(
            f"""\

Normalization info:

  - Assay **{label1}** is multiplied by {fct1}
  - Assay **{label2}** is multiplied by {fct2}
"""
        )
    elif normalization == 'mean':
        ntb1 = np.mean(tb1)
        ntb2 = np.mean(tb2)
        ntb = np.mean([ntb1, ntb2])
        fct1 = ntb / ntb1
        fct2 = ntb / ntb2
        tb = pd.concat([tb1 * fct1, tb2 * fct2], axis=0)

        gn.add_markdown(
            f"""\

Normalization info:",

  - Assay **{label1}** is multiplied by {fct1}
  - Assay **{label2}** is multiplied by {fct2}
"""
        )
    else:
        raise ValueError()

    if direction == 'samples':
        tb = tb.T

    gn.add_markdown(
        f"""\
You combined the following assays:

  - Assay 1 (with {tb1.shape[0]} genes and {tb1.shape[1]} cells)
  - Assay 2 (with {tb2.shape[0]} genes and {tb2.shape[1]} cells)

into:

  - Combined Assay (with {tb.shape[0]} genes and {tb.shape[1]} cells)
"""
    )

    gn.export_statically(gn.assay_from_pandas(tb), 'Combined assay')

    if direction == 'samples':
        meta_type = 'sampleMeta'
    elif direction == 'genes':
        meta_type = 'geneMeta'
    else:
        raise ValueError()

    gn.export({**{x: label1 for x in tb1.index}, **{x: label2 for x in tb2.index}}, 'Assay label', meta_type)

    gn.commit()


if __name__ == '__main__':
    main()
