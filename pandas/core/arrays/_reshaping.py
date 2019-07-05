"""
Utilities for implementing 2D compatibility for 1D ExtensionArrays.
"""
from functools import wraps
from typing import Tuple

import numpy as np

from pandas._libs.lib import is_integer


def implement_2d(cls):
    """
    A decorator to take a 1-dimension-only ExtensionArray subclass and make
    it support limited 2-dimensional operations.
    """
    from pandas.core.arrays import ExtensionArray

    # For backwards-compatibility, if an EA author implemented __len__
    #  but not size, we use that __len__ method to get an array's size.
    has_size = cls.size is not ExtensionArray.size
    has_shape = cls.shape is not ExtensionArray.shape
    has_len = cls.__len__ is not ExtensionArray.__len__

    if not has_size and has_len:
        cls.size = property(cls.__len__)
        cls.__len__ = ExtensionArray.__len__

    elif not has_size and has_shape:

        @property
        def size(self) -> int:
            return np.prod(self.shape)

        cls.size = size

    orig_copy = cls.copy

    @wraps(orig_copy)
    def copy(self):
        result = orig_copy(self)
        result._shape = self._shape
        return result

    cls.copy = copy

    return cls


def can_safe_ravel(shape: Tuple[int, ...]) -> bool:
    """
    Check if an array with the given shape can be ravelled unambiguously
    regardless of column/row order.

    Parameters
    ----------
    shape : tuple[int]

    Returns
    -------
    bool
    """
    if len(shape) == 1:
        return True
    if len(shape) > 2:
        raise NotImplementedError(shape)
    if shape[0] == 1 or shape[1] == 1:
        # column-like or row-like
        return True
    return False


def tuplify_shape(size: int, shape, restrict=True) -> Tuple[int, ...]:
    """
    Convert a passed shape into a valid tuple.
    Following ndarray.reshape, we accept either `reshape(a, b)` or
    `reshape((a, b))`, the latter being canonical.

    Parameters
    ----------
    size : int
    shape : tuple
    restrict : bool, default True
        Whether to restrict to shapes (N), (1,N), and (N,1)

    Returns
    -------
    tuple[int, ...]
    """
    if len(shape) == 0:
        raise ValueError("shape must be a non-empty tuple of integers", shape)

    if len(shape) == 1:
        if is_integer(shape[0]):
            pass
        else:
            shape = shape[0]
            if not isinstance(shape, tuple):
                raise ValueError("shape must be a non-empty tuple of integers", shape)

    if not all(is_integer(x) for x in shape):
        raise ValueError("shape must be a non-empty tuple of integers", shape)

    if any(x < -1 for x in shape):
        raise ValueError("Invalid shape {shape}".format(shape=shape))

    if -1 in shape:
        if shape.count(-1) != 1:
            raise ValueError("Invalid shape {shape}".format(shape=shape))
        idx = shape.index(-1)
        others = [n for n in shape if n != -1]
        prod = np.prod(others)
        dim = size // prod
        shape = shape[:idx] + (dim,) + shape[idx + 1 :]

    if np.prod(shape) != size:
        raise ValueError(
            "Product of shape ({shape}) must match "
            "size ({size})".format(shape=shape, size=size)
        )

    num_gt1 = len([x for x in shape if x > 1])
    if num_gt1 > 1 and restrict:
        raise ValueError(
            "The default `reshape` implementation is limited to "
            "shapes (N,), (N,1), and (1,N), not {shape}".format(shape=shape)
        )
    return shape
