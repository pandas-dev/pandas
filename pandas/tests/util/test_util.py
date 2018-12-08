# -*- coding: utf-8 -*-
import pytest

from pandas.util._decorators import deprecate_kwarg, make_signature
from pandas.util._validators import validate_kwargs

import pandas.util.testing as tm


def test_rands():
    r = tm.rands(10)
    assert(len(r) == 10)


def test_rands_array_1d():
    arr = tm.rands_array(5, size=10)
    assert(arr.shape == (10,))
    assert(len(arr[0]) == 5)


def test_rands_array_2d():
    arr = tm.rands_array(7, size=(10, 10))
    assert(arr.shape == (10, 10))
    assert(len(arr[1, 1]) == 7)


def test_numpy_err_state_is_default():
    # The defaults since numpy 1.6.0
    expected = {"over": "warn", "divide": "warn",
                "invalid": "warn", "under": "ignore"}
    import numpy as np

    # The error state should be unchanged after that import.
    assert np.geterr() == expected


@pytest.mark.parametrize("func,expected", [
    # Case where the func does not have default kwargs.
    (validate_kwargs, (["fname", "kwargs", "compat_args"],
                       ["fname", "kwargs", "compat_args"])),

    # Case where the func does have default kwargs.
    (deprecate_kwarg, (["old_arg_name", "new_arg_name",
                        "mapping=None", "stacklevel=2"],
                       ["old_arg_name", "new_arg_name",
                        "mapping", "stacklevel"]))
])
def test_make_signature(func, expected):
    # see gh-17608
    assert make_signature(func) == expected
