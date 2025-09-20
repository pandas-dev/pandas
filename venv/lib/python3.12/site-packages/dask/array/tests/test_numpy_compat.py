from __future__ import annotations

import numpy as np
import pytest

import dask.array as da
from dask.array.utils import assert_eq


@pytest.fixture(
    params=[
        [("A", ("f4", (3, 2))), ("B", ("f4", 3)), ("C", ("f8", 3))],
        [("A", ("i4", (3, 2))), ("B", ("f4", 3)), ("C", ("S4", 3))],
    ]
)
def dtype(request):
    return np.dtype(request.param)


@pytest.fixture(params=[["A"], ["A", "B"], ["A", "B", "C"]])
def index(request):
    return request.param


def test_basic():
    # sanity check
    dtype = [("a", "f8"), ("b", "f8"), ("c", "f8")]
    x = np.ones((5, 3), dtype=dtype)
    dx = da.ones((5, 3), dtype=dtype, chunks=3)
    result = dx[["a", "b"]]
    expected = x[["a", "b"]]
    assert_eq(result, expected)


def test_min_max_round_funcs():
    # Regression test for gh-5031
    image = da.from_array(np.array([[0, 1], [1, 2]]), chunks=(1, 2))
    # These use __array_function__ (and min/max/round are aliased,
    # to amin/amax/round_ in numpy)
    assert int(np.min(image)) == 0
    assert int(np.max(image)) == 2
    assert np.round(image)[1, 1] == 2
