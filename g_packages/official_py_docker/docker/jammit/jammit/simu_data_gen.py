import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.random import choice, randn

from .utils import center_rows


def simu_data_gen(
    n_features=20,
    n_samples=10,
    n_sig_features=0.25,
    n_sig_support_samples=0.5,
    sig_intensity=10,
    noise_intensity=1,
    features_as_rows=True,
    randomize_sig_idxs=False,
    randomize_sig_support_idxs=False,
    centered=True,
):
    """This function generates a simulated matrix with a given size and signal
    """
    assert n_sig_features <= n_features
    assert n_sig_support_samples <= n_samples

    if n_sig_features < 1:
        n_sig_features = math.floor(n_features * n_sig_features)

    if n_sig_support_samples < 1:
        n_sig_support_samples = math.floor(n_sig_support_samples * n_samples)

    if randomize_sig_idxs:
        sig_feature_idxs = choice(range(n_features), n_sig_features, replace=False)
    else:
        sig_feature_idxs = np.array(range(n_sig_features))

    if randomize_sig_support_idxs:
        sig_support_sample_idxs = choice(
            range(n_samples), n_sig_support_samples, replace=False
        )
    else:
        sig_support_sample_idxs = np.array(range(n_sig_support_samples))

    mat = noise_intensity * randn(n_features, n_samples)

    true_u = np.full(n_features, 0)
    true_v = np.full(n_samples, 0)

    for i in sig_feature_idxs:
        true_u[i] = 1

    for j in sig_support_sample_idxs:
        true_v[j] = 1

    for i in sig_feature_idxs:
        for j in sig_support_sample_idxs:
            mat[i, j] += sig_intensity

    if centered:
        mat = center_rows(mat)

    if features_as_rows:
        return mat, true_u, true_v
    else:
        return mat.T, true_u, true_v


def main():
    mat, _, _ = simu_data_gen(
        randomize_sig_idxs=True, randomize_sig_support_idxs=True, centered=False
    )
    plt.imshow(mat)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    main()
