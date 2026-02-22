"""Check that functions can handle list input"""

import warnings

import numpy as np
from numpy.testing import assert_array_almost_equal
import bottleneck as bn
from .util import DTYPES
import pytest


def lists(dtypes=DTYPES):
    """Iterator that yields lists to use for unit testing."""
    ss = {}
    ss[1] = {"size": 4, "shapes": [(4,)]}
    ss[2] = {"size": 6, "shapes": [(1, 6), (2, 3)]}
    ss[3] = {"size": 6, "shapes": [(1, 2, 3)]}
    ss[4] = {"size": 24, "shapes": [(1, 2, 3, 4)]}
    for ndim in ss:
        size = ss[ndim]["size"]
        shapes = ss[ndim]["shapes"]
        a = np.arange(size)
        for shape in shapes:
            a = a.reshape(shape)
            for dtype in dtypes:
                yield a.astype(dtype).tolist()


@pytest.mark.parametrize("func", bn.get_functions("all"), ids=lambda x: x.__name__)
def test_list_input(func):
    """Test that bn.xxx gives the same output as bn.slow.xxx for list input."""
    msg = "\nfunc %s | input %s (%s) | shape %s\n"
    msg += "\nInput array:\n%s\n"
    name = func.__name__
    if name == "replace":
        return
    func0 = eval("bn.slow.%s" % name)
    for i, a in enumerate(lists()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                actual = func(a)
                desired = func0(a)
            except TypeError:
                actual = func(a, 2)
                desired = func0(a, 2)
        a = np.array(a)
        tup = (name, "a" + str(i), str(a.dtype), str(a.shape), a)
        err_msg = msg % tup
        assert_array_almost_equal(actual, desired, err_msg=err_msg)
