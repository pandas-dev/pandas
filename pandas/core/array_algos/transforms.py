"""
transforms.py is for shape-preserving functions.
"""

import numpy as np

from pandas.core.dtypes.common import ensure_platform_int


def shift(values: np.ndarray, periods: int, axis: int, fill_value) -> np.ndarray:
    new_values = values

    # make sure array sent to np.roll is c_contiguous
    f_ordered = values.flags.f_contiguous
    if f_ordered:
        new_values = new_values.T
        axis = new_values.ndim - axis - 1

    if np.prod(new_values.shape):
        new_values = np.roll(new_values, ensure_platform_int(periods), axis=axis)

    axis_indexer = [slice(None)] * values.ndim
    if periods > 0:
        axis_indexer[axis] = slice(None, periods)
    else:
        axis_indexer[axis] = slice(periods, None)
    new_values[tuple(axis_indexer)] = fill_value

    # restore original order
    if f_ordered:
        new_values = new_values.T

    return new_values
