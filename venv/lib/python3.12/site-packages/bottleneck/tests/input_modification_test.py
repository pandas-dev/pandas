"""Test functions."""

import warnings

import numpy as np
from numpy.testing import assert_equal
import bottleneck as bn
from .util import DTYPES
import pytest


def arrays(dtypes):
    """Iterator that yield arrays to use for unit testing."""
    ss = {}
    ss[1] = {"size": 4, "shapes": [(4,)]}
    ss[2] = {"size": 6, "shapes": [(2, 3)]}
    ss[3] = {"size": 6, "shapes": [(1, 2, 3)]}
    rs = np.random.RandomState([1, 2, 3])
    for ndim in ss:
        size = ss[ndim]["size"]
        shapes = ss[ndim]["shapes"]
        for dtype in dtypes:
            a = np.arange(size, dtype=dtype)
            if issubclass(a.dtype.type, np.inexact):
                idx = rs.rand(*a.shape) < 0.2
                a[idx] = np.inf
                idx = rs.rand(*a.shape) < 0.2
                a[idx] = np.nan
                idx = rs.rand(*a.shape) < 0.2
                a[idx] *= -1
            for shape in shapes:
                a = a.reshape(shape)
                yield a


@pytest.mark.parametrize("func", bn.get_functions("all"), ids=lambda x: x.__name__)
def test_modification(func):
    """Test that bn.xxx gives the same output as np.xxx."""
    name = func.__name__
    if name == "replace":
        return
    msg = "\nInput array modified by %s.\n\n"
    msg += "input array before:\n%s\nafter:\n%s\n"
    for i, a in enumerate(arrays(DTYPES)):
        axes = list(range(-a.ndim, a.ndim))
        if all(x not in name for x in ["push", "move", "sort", "partition"]):
            axes += [None]

        second_arg = 1
        if "partition" in name:
            second_arg = 0

        for axis in axes:
            with np.errstate(invalid="ignore"):
                a1 = a.copy()
                a2 = a.copy()
                if any(x in name for x in ["move", "sort", "partition"]):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        func(a1, second_arg, axis=axis)
                else:
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            func(a1, axis=axis)
                    except ValueError as e:
                        if name.startswith(
                            "nanarg"
                        ) and "All-NaN slice encountered" in str(e):
                            continue
                assert_equal(a1, a2, msg % (name, a1, a2))
