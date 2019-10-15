import math
import sys

import numpy as np
from colors import color

BLOCKS = " ▁▂▃▄▅▆▇█"


def all_equal(v):
    return all([x == v[0] for x in v[1:]])


def number_to_block(x, b):
    if x > 0:
        c = "red"
    else:
        c = "blue"

    if math.isnan(x):
        return "n"
    else:
        return color(BLOCKS[int(round(abs(x) / b * 8))], fg=c)


def array_visualizer_abs(xs):
    return "".join([number_to_block(x, np.max(np.abs(xs))) for x in xs])


def array_visualizer(xs):
    return "".join(["n" if np.isnan(x) else "." if x == 0 else "+" if x > 0 else "-" for x in xs])


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def center_rows(xx):
    return xx - np.mean(xx, axis=1, keepdims=True)


def normalize(v, epsilon=1e-10):
    return v / (np.linalg.norm(v) + epsilon)


def split_by_lens(xs, lens):
    return np.split(xs, np.cumsum(lens))


def permute_each_row(xx, replace=False):
    """Permute each row independently for a given matrix and return the
    result as a new matrix

    Args:
        xx (np.array with dim=2):
            The matrix to be permuted

    Kwargs:
        replace (bool):
            Whether to have replacement

    Returns:
        The permuted matrix, with the same shape as xx
    """

    if replace:
        return np.array([np.random.choice(x, len(x)) for x in xx])
    else:
        return np.array([np.random.permutation(x) for x in xx])
