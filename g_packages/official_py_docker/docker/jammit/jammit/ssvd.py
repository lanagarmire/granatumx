import sys
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy

from .simu_data_gen import simu_data_gen
from .utils import all_equal, eprint, normalize


def ssvd(xx, tu=None, alpha=0, max_n_iter=5000, verbose=1, convergence_threshold=0.0001):
    # estimate the u0 from classical SVD
    if tu is None:
        tu, _, _ = np.linalg.svd(xx, full_matrices=False)

    # calcuate sparse, rank-1 SVD
    u0 = tu[:, 0]
    ud = 1
    iter_cnt = 0
    while ud > convergence_threshold:
        iter_cnt += 1
        v = normalize(xx.T @ u0)
        u = xx @ v
        u = np.sign(u) * np.maximum(np.abs(u) - alpha, 0)
        u = normalize(u)
        ud = np.linalg.norm(u0 - u)
        u0 = u
        if iter_cnt > max_n_iter:
            if verbose >= 1:
                eprint(f"Failed to converge in {max_n_iter} iterations.")
            break

    # signing - make sure that v is positively correlated with colmean(xx)
    if scipy.stats.pearsonr(v, np.mean(xx, axis=0))[0] < 0:
        u = -u
        v = -v

    return (u, v)
