"""
Implementations of high-level numpy functions that are ExtensionArray-compatible.
"""
import numpy as np

from pandas._typing import ArrayLike


def tile(arr: ArrayLike, shape) -> ArrayLike:
    raise NotImplementedError


def broadcast_to(arr: ArrayLike, shape, orient=None) -> ArrayLike:
    if isinstance(arr, np.ndarray):
        values = arr
    else:
        # ExtensionArray
        values = arr._values_for_factorize()[0]

    # TODO: if we are ndim==size==1 it shouldnt matter whether rowlike/columnlike?
    if values.ndim == 1 and orient is not None:
        # SUpport treating a 1-dimensional array as either a row or column
        assert orient in ["rowlike", "columnlike"]
        if orient == "rowlike":
            values = values.reshape(1, -1)
        else:
            values = values.reshape(-1, 1)

    btvalues = np.broadcast_to(values, shape)
    if isinstance(arr, np.ndarray):
        result = btvalues
    else:
        result = type(arr)._from_factorized(btvalues, arr)
    return result
