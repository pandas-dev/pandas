"""Test moving window functions."""

import numpy as np
from numpy.testing import assert_equal, assert_array_almost_equal, assert_raises
import bottleneck as bn
from .util import arrays, array_order
import pytest


@pytest.mark.parametrize("func", bn.get_functions("move"), ids=lambda x: x.__name__)
def test_move(func):
    """Test that bn.xxx gives the same output as a reference function."""
    fmt = (
        "\nfunc %s | window %d | min_count %s | input %s (%s) | shape %s | "
        "axis %s | order %s\n"
    )
    fmt += "\nInput array:\n%s\n"
    aaae = assert_array_almost_equal
    func_name = func.__name__
    func0 = eval("bn.slow.%s" % func_name)
    if func_name == "move_var":
        decimal = 3
    else:
        decimal = 5
    for i, a in enumerate(arrays(func_name)):
        axes = range(-1, a.ndim)
        for axis in axes:
            windows = range(1, a.shape[axis])
            for window in windows:
                min_counts = list(range(1, window + 1)) + [None]
                for min_count in min_counts:
                    actual = func(a, window, min_count, axis=axis)
                    desired = func0(a, window, min_count, axis=axis)
                    tup = (
                        func_name,
                        window,
                        str(min_count),
                        "a" + str(i),
                        str(a.dtype),
                        str(a.shape),
                        str(axis),
                        array_order(a),
                        a,
                    )
                    err_msg = fmt % tup
                    aaae(actual, desired, decimal, err_msg)
                    err_msg += "\n dtype mismatch %s %s"
                    da = actual.dtype
                    dd = desired.dtype
                    assert_equal(da, dd, err_msg % (da, dd))


# ---------------------------------------------------------------------------
# Test argument parsing


@pytest.mark.parametrize("func", bn.get_functions("move"), ids=lambda x: x.__name__)
def test_arg_parsing(func, decimal=5):
    """test argument parsing."""

    name = func.__name__
    func0 = eval("bn.slow.%s" % name)

    a = np.array([1.0, 2, 3])

    fmt = "\n%s" % func
    fmt += "%s\n"
    fmt += "\nInput array:\n%s\n" % a

    actual = func(a, 2)
    desired = func0(a, 2)
    err_msg = fmt % "(a, 2)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, 2, 1)
    desired = func0(a, 2, 1)
    err_msg = fmt % "(a, 2, 1)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, window=2)
    desired = func0(a, window=2)
    err_msg = fmt % "(a, window=2)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, window=2, min_count=1)
    desired = func0(a, window=2, min_count=1)
    err_msg = fmt % "(a, window=2, min_count=1)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, window=2, min_count=1, axis=0)
    desired = func0(a, window=2, min_count=1, axis=0)
    err_msg = fmt % "(a, window=2, min_count=1, axis=0)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, min_count=1, window=2, axis=0)
    desired = func0(a, min_count=1, window=2, axis=0)
    err_msg = fmt % "(a, min_count=1, window=2, axis=0)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, axis=-1, min_count=None, window=2)
    desired = func0(a, axis=-1, min_count=None, window=2)
    err_msg = fmt % "(a, axis=-1, min_count=None, window=2)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a=a, axis=-1, min_count=None, window=2)
    desired = func0(a=a, axis=-1, min_count=None, window=2)
    err_msg = fmt % "(a=a, axis=-1, min_count=None, window=2)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    if name in ("move_std", "move_var"):
        actual = func(a, 2, 1, -1, ddof=1)
        desired = func0(a, 2, 1, -1, ddof=1)
        err_msg = fmt % "(a, 2, 1, -1, ddof=1)"
        assert_array_almost_equal(actual, desired, decimal, err_msg)

    # regression test: make sure len(kwargs) == 0 doesn't raise
    args = (a, 1, 1, -1)
    kwargs = {}
    func(*args, **kwargs)


@pytest.mark.parametrize("func", bn.get_functions("move"), ids=lambda x: x.__name__)
def test_arg_parse_raises(func):
    """test argument parsing raises in move"""
    a = np.array([1.0, 2, 3])
    assert_raises(TypeError, func)
    assert_raises(TypeError, func, axis=a)
    assert_raises(TypeError, func, a, 2, axis=0, extra=0)
    assert_raises(TypeError, func, a, 2, axis=0, a=a)
    assert_raises(TypeError, func, a, 2, 2, 0, 0, 0)
    assert_raises(TypeError, func, a, 2, axis="0")
    assert_raises(TypeError, func, a, 1, min_count="1")
    if func.__name__ not in ("move_std", "move_var"):
        assert_raises(TypeError, func, a, 2, ddof=0)


# ---------------------------------------------------------------------------
# move_median.c is complicated. Let's do some more testing.
#
# If you make changes to move_median.c then do lots of tests by increasing
# range(100) in the two functions below to range(10000). And for extra credit
# increase size to 30. With those two changes the unit tests will take a
# LONG time to run.


def test_move_median_with_nans():
    """test move_median.c with nans"""
    fmt = "\nfunc %s | window %d | min_count %s\n\nInput array:\n%s\n"
    aaae = assert_array_almost_equal
    min_count = 1
    size = 10
    func = bn.move_median
    func0 = bn.slow.move_median
    rs = np.random.RandomState([1, 2, 3])
    for i in range(100):
        a = np.arange(size, dtype=np.float64)
        idx = rs.rand(*a.shape) < 0.1
        a[idx] = np.inf
        idx = rs.rand(*a.shape) < 0.2
        a[idx] = np.nan
        rs.shuffle(a)
        for window in range(2, size + 1):
            actual = func(a, window=window, min_count=min_count)
            desired = func0(a, window=window, min_count=min_count)
            err_msg = fmt % (func.__name__, window, min_count, a)
            aaae(actual, desired, decimal=5, err_msg=err_msg)


def test_move_median_without_nans():
    """test move_median.c without nans"""
    fmt = "\nfunc %s | window %d | min_count %s\n\nInput array:\n%s\n"
    aaae = assert_array_almost_equal
    min_count = 1
    size = 10
    func = bn.move_median
    func0 = bn.slow.move_median
    rs = np.random.RandomState([1, 2, 3])
    for i in range(100):
        a = np.arange(size, dtype=np.int64)
        rs.shuffle(a)
        for window in range(2, size + 1):
            actual = func(a, window=window, min_count=min_count)
            desired = func0(a, window=window, min_count=min_count)
            err_msg = fmt % (func.__name__, window, min_count, a)
            aaae(actual, desired, decimal=5, err_msg=err_msg)


# ----------------------------------------------------------------------------
# Regression test for square roots of negative numbers


def test_move_std_sqrt():
    """Test move_std for neg sqrt."""

    a = [
        0.0011448196318903589,
        0.00028718669878572767,
        0.00028718669878572767,
        0.00028718669878572767,
        0.00028718669878572767,
    ]
    err_msg = "Square root of negative number. ndim = %d"
    b = bn.move_std(a, window=3)
    assert np.isfinite(b[2:]).all(), err_msg % 1

    a2 = np.array([a, a])
    b = bn.move_std(a2, window=3, axis=1)
    assert np.isfinite(b[:, 2:]).all(), err_msg % 2

    a3 = np.array([[a, a], [a, a]])
    b = bn.move_std(a3, window=3, axis=2)
    assert np.isfinite(b[:, :, 2:]).all(), err_msg % 3
