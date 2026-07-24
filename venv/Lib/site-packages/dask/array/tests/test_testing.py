from __future__ import annotations

import sys

import numpy as np
import pytest

import dask.array as da
from dask.array.utils import assert_eq


@pytest.mark.skipif(bool(sys.flags.optimize), reason="Assertions disabled.")
def test_assert_eq_checks_scalars():
    # https://github.com/dask/dask/issues/2680
    with pytest.raises(AssertionError):
        assert_eq(np.array(0), np.array(1))

    a = da.from_array(np.array([0]), 1)[0]
    b = np.array([1])[0]
    with pytest.raises(AssertionError):
        assert_eq(a, b)


def test_assert_eq_array():
    x = np.arange(6).reshape((2, 3))
    d = da.from_array(x, chunks=(1, 3))
    assert_eq(d, x)
    assert_eq(d, d)
    with pytest.raises(AssertionError):
        assert_eq(d, x + 1)


def test_assert_eq_checks_shape():
    a = da.ones(3, chunks=2)
    assert_eq(a, a)
    with pytest.raises(AssertionError):
        assert_eq(da.ones(2), da.ones(3))


def test_assert_eq_check_dtype():
    a = da.ones(3, chunks=2, dtype="i4")
    b = da.ones(3, chunks=2, dtype="f8")
    with pytest.raises(AssertionError):
        assert_eq(a, b)
    # values are equal, so disabling the dtype check passes
    assert_eq(a, b, check_dtype=False)


def test_assert_eq_equal_nan():
    a = da.from_array(np.array([1.0, np.nan]), chunks=2)
    b = da.from_array(np.array([1.0, np.nan]), chunks=2)
    # NaNs compare equal by default
    assert_eq(a, b)
    with pytest.raises(AssertionError):
        assert_eq(a, b, equal_nan=False)
