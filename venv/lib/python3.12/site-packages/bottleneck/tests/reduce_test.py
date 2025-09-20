"""Test reduce functions."""

import warnings
import traceback

import numpy as np
from numpy.testing import assert_equal, assert_raises, assert_array_almost_equal

import bottleneck as bn
from .util import arrays, array_order, DTYPES
import pytest


@pytest.mark.parametrize("func", bn.get_functions("reduce"), ids=lambda x: x.__name__)
def test_reduce(func):
    """test reduce functions"""
    return unit_maker(func)


def unit_maker(func, decimal=5, skip_dtype=("nansum", "ss")):
    """Test that bn.xxx gives the same output as bn.slow.xxx."""
    fmt = "\nfunc %s | input %s (%s) | shape %s | axis %s | order %s\n"
    fmt += "\nInput array:\n%s\n"
    name = func.__name__
    func0 = eval("bn.slow.%s" % name)
    for i, a in enumerate(arrays(name)):
        if a.ndim == 0:
            axes = [None]  # numpy can't handle e.g. np.nanmean(9, axis=-1)
        else:
            axes = list(range(-1, a.ndim)) + [None]
        for axis in axes:
            actual = "Crashed"
            desired = "Crashed"
            actualraised = False
            try:
                # do not use a.copy() here because it will C order the array
                actual = func(a, axis=axis)
            except:  # noqa
                actualraised = True
            desiredraised = False
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    desired = func0(a, axis=axis)
            except:  # noqa
                desiredraised = True
            if actualraised and desiredraised:
                pass
            else:
                tup = (
                    name,
                    "a" + str(i),
                    str(a.dtype),
                    str(a.shape),
                    str(axis),
                    array_order(a),
                    a,
                )
                err_msg = fmt % tup
                if actualraised != desiredraised:
                    if actualraised:
                        fmt2 = "\nbn.%s raised\nbn.slow.%s ran\n\n%s"
                    else:
                        fmt2 = "\nbn.%s ran\nbn.slow.%s raised\n\n%s"
                    msg = fmt2 % (name, name, traceback.format_exc())
                    err_msg += msg
                    assert False, err_msg
                assert_array_almost_equal(actual, desired, decimal, err_msg)
                err_msg += "\n dtype mismatch %s %s"
                if name not in skip_dtype:
                    if hasattr(actual, "dtype") and hasattr(desired, "dtype"):
                        da = actual.dtype
                        dd = desired.dtype
                        assert_equal(da, dd, err_msg % (da, dd))


# ---------------------------------------------------------------------------
# Test argument parsing


@pytest.mark.parametrize("func", bn.get_functions("reduce"), ids=lambda x: x.__name__)
def test_arg_parsing(func):
    """test argument parsing"""
    return unit_maker_argparse(func)


def unit_maker_argparse(func, decimal=5):
    """test argument parsing."""

    name = func.__name__
    func0 = eval("bn.slow.%s" % name)

    a = np.array([1.0, 2, 3])

    fmt = "\n%s" % func
    fmt += "%s\n"
    fmt += "\nInput array:\n%s\n" % a

    actual = func(a)
    desired = func0(a)
    err_msg = fmt % "(a)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, 0)
    desired = func0(a, 0)
    err_msg = fmt % "(a, 0)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, None)
    desired = func0(a, None)
    err_msg = fmt % "(a, None)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, axis=0)
    desired = func0(a, axis=0)
    err_msg = fmt % "(a, axis=0)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, axis=None)
    desired = func0(a, axis=None)
    err_msg = fmt % "(a, axis=None)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a=a)
    desired = func0(a=a)
    err_msg = fmt % "(a)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    # regression test: make sure len(kwargs) == 0 doesn't raise
    args = (a, 0)
    kwargs = {}
    func(*args, **kwargs)


@pytest.mark.parametrize("func", bn.get_functions("reduce"))
def test_arg_parse_raises(func):
    """test argument parsing raises in reduce"""
    return unit_maker_argparse_raises(func)


def unit_maker_argparse_raises(func):
    """test argument parsing raises in reduce"""
    a = np.array([1.0, 2, 3])
    assert_raises(TypeError, func)
    assert_raises(TypeError, func, axis=a)
    assert_raises(TypeError, func, a, axis=0, extra=0)
    assert_raises(TypeError, func, a, axis=0, a=a)
    assert_raises(TypeError, func, a, 0, 0, 0, 0, 0)
    assert_raises(TypeError, func, a, axis="0")
    if func.__name__ not in ("nanstd", "nanvar"):
        assert_raises(TypeError, func, a, ddof=0)
    assert_raises(TypeError, func, a, a)
    # assert_raises(TypeError, func, None) results vary


# ---------------------------------------------------------------------------
# Check that exceptions are raised


def test_nanmax_size_zero(dtypes=DTYPES):
    """Test nanmax for size zero input arrays."""
    shapes = [(0,), (2, 0), (1, 2, 0)]
    for shape in shapes:
        for dtype in dtypes:
            a = np.zeros(shape, dtype=dtype)
            assert_raises(ValueError, bn.nanmax, a)
            assert_raises(ValueError, bn.slow.nanmax, a)


def test_nanmin_size_zero(dtypes=DTYPES):
    """Test nanmin for size zero input arrays."""
    shapes = [(0,), (2, 0), (1, 2, 0)]
    for shape in shapes:
        for dtype in dtypes:
            a = np.zeros(shape, dtype=dtype)
            assert_raises(ValueError, bn.nanmin, a)
            assert_raises(ValueError, bn.slow.nanmin, a)


# ---------------------------------------------------------------------------
# nanstd and nanvar regression test (issue #60)


def test_nanstd_issue60():
    """nanstd regression test (issue #60)"""

    f = bn.nanstd([1.0], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanstd([1.0], ddof=1)
    assert_equal(f, s, err_msg="bn.nanstd([1.0], ddof=1) wrong")

    f = bn.nanstd([1], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanstd([1], ddof=1)
    assert_equal(f, s, err_msg="bn.nanstd([1], ddof=1) wrong")

    f = bn.nanstd([1, np.nan], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanstd([1, np.nan], ddof=1)
    assert_equal(f, s, err_msg="bn.nanstd([1, nan], ddof=1) wrong")

    f = bn.nanstd([[1, np.nan], [np.nan, 1]], axis=0, ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanstd([[1, np.nan], [np.nan, 1]], axis=0, ddof=1)
    assert_equal(f, s, err_msg="issue #60 regression")


def test_nanvar_issue60():
    """nanvar regression test (issue #60)"""

    f = bn.nanvar([1.0], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanvar([1.0], ddof=1)
    assert_equal(f, s, err_msg="bn.nanvar([1.0], ddof=1) wrong")

    f = bn.nanvar([1], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanvar([1], ddof=1)
    assert_equal(f, s, err_msg="bn.nanvar([1], ddof=1) wrong")

    f = bn.nanvar([1, np.nan], ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanvar([1, np.nan], ddof=1)
    assert_equal(f, s, err_msg="bn.nanvar([1, nan], ddof=1) wrong")

    f = bn.nanvar([[1, np.nan], [np.nan, 1]], axis=0, ddof=1)
    with np.errstate(invalid="ignore"):
        s = bn.slow.nanvar([[1, np.nan], [np.nan, 1]], axis=0, ddof=1)
    assert_equal(f, s, err_msg="issue #60 regression")


@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize("func", (bn.nanstd, bn.nanvar), ids=lambda x: x.__name__)
def test_ddof_nans(func, dtype):
    array = np.ones((1, 1), dtype=dtype)
    for axis in [None, 0, 1, -1]:
        result = func(array, axis=axis, ddof=3)
        assert np.isnan(result)
