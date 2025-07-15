import numpy as np
from numpy.testing import (
    assert_equal,
    assert_array_equal,
    assert_array_almost_equal,
    assert_raises,
)

import bottleneck as bn
from .reduce_test import (
    unit_maker as reduce_unit_maker,
    unit_maker_argparse as unit_maker_parse_rankdata,
)
from .util import arrays, array_order, DTYPES
import pytest

# ---------------------------------------------------------------------------
# partition, argpartition


@pytest.mark.parametrize(
    "func", (bn.partition, bn.argpartition), ids=lambda x: x.__name__
)
def test_partition_and_argpartition(func):
    """test partition or argpartition"""

    msg = "\nfunc %s | input %s (%s) | shape %s | n %d | axis %s | order %s\n"
    msg += "\nInput array:\n%s\n"

    name = func.__name__
    func0 = eval("bn.slow.%s" % name)

    rs = np.random.RandomState([1, 2, 3])
    for i, a in enumerate(arrays(name)):
        if a.ndim == 0 or a.size == 0 or a.ndim > 3:
            continue
        for axis in list(range(-1, a.ndim)) + [None]:
            if axis is None:
                nmax = a.size - 1
            else:
                nmax = a.shape[axis] - 1
            if nmax < 1:
                continue
            n = rs.randint(nmax)
            s0 = func0(a, n, axis)
            s1 = func(a, n, axis)
            if name == "argpartition":
                s0 = complete_the_argpartition(s0, a, n, axis)
                s1 = complete_the_argpartition(s1, a, n, axis)
            else:
                s0 = complete_the_partition(s0, n, axis)
                s1 = complete_the_partition(s1, n, axis)
            tup = (
                name,
                "a" + str(i),
                str(a.dtype),
                str(a.shape),
                n,
                str(axis),
                array_order(a),
                a,
            )
            err_msg = msg % tup
            assert_array_equal(s1, s0, err_msg)


def complete_the_partition(a, n, axis):
    def func1d(a, n):
        a[:n] = np.sort(a[:n])
        a[n + 1 :] = np.sort(a[n + 1 :])
        return a

    a = a.copy()
    ndim = a.ndim
    if axis is None:
        if ndim != 1:
            raise ValueError("`a` must be 1d when axis is None")
        axis = 0
    elif axis < 0:
        axis += ndim
        if axis < 0:
            raise ValueError("`axis` out of range")
    a = np.apply_along_axis(func1d, axis, a, n)
    return a


def complete_the_argpartition(index, a, n, axis):
    a = a.copy()
    ndim = a.ndim
    if axis is None:
        if index.ndim != 1:
            raise ValueError("`index` must be 1d when axis is None")
        axis = 0
        ndim = 1
        a = a.reshape(-1)
    elif axis < 0:
        axis += ndim
        if axis < 0:
            raise ValueError("`axis` out of range")
    if ndim == 1:
        a = a[index]
    elif ndim == 2:
        if axis == 0:
            for i in range(a.shape[1]):
                a[:, i] = a[index[:, i], i]
        elif axis == 1:
            for i in range(a.shape[0]):
                a[i] = a[i, index[i]]
        else:
            raise ValueError("`axis` out of range")
    elif ndim == 3:
        if axis == 0:
            for i in range(a.shape[1]):
                for j in range(a.shape[2]):
                    a[:, i, j] = a[index[:, i, j], i, j]
        elif axis == 1:
            for i in range(a.shape[0]):
                for j in range(a.shape[2]):
                    a[i, :, j] = a[i, index[i, :, j], j]
        elif axis == 2:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    a[i, j, :] = a[i, j, index[i, j, :]]
        else:
            raise ValueError("`axis` out of range")
    else:
        raise ValueError("`a.ndim` must be 1, 2, or 3")
    a = complete_the_partition(a, n, axis)
    return a


def test_transpose():
    """partition transpose test"""
    a = np.arange(12).reshape(4, 3)
    actual = bn.partition(a.T, 2, -1).T
    desired = bn.slow.partition(a.T, 2, -1).T
    assert_equal(actual, desired, "partition transpose test")


# ---------------------------------------------------------------------------
# rankdata, nanrankdata, push


@pytest.mark.parametrize(
    "func", (bn.rankdata, bn.nanrankdata, bn.push), ids=lambda x: x.__name__
)
def test_nonreduce_axis(func):
    """Test nonreduce axis functions"""
    return reduce_unit_maker(func)


def test_push():
    """Test push"""
    ns = (0, 1, 2, 3, 4, 5, None)
    a = np.array([np.nan, 1, 2, np.nan, np.nan, np.nan, np.nan, 3, np.nan])
    for n in ns:
        actual = bn.push(a.copy(), n=n)
        desired = bn.slow.push(a.copy(), n=n)
        assert_array_equal(actual, desired, "failed on n=%s" % str(n))


# ---------------------------------------------------------------------------
# Test argument parsing


@pytest.mark.parametrize(
    "func", bn.get_functions("nonreduce_axis"), ids=lambda x: x.__name__
)
def test_arg_parsing(func):
    """test argument parsing in nonreduce_axis"""
    name = func.__name__
    if name in ("partition", "argpartition"):
        return unit_maker_parse(func)
    elif name in ("push"):
        return unit_maker_parse(func)
    elif name in ("rankdata", "nanrankdata"):
        return unit_maker_parse_rankdata(func)
    else:
        fmt = "``%s` is an unknown nonreduce_axis function"
        raise ValueError(fmt % name)
    return unit_maker_raises(func)


def unit_maker_parse(func, decimal=5):
    """test argument parsing."""

    name = func.__name__
    func0 = eval("bn.slow.%s" % name)

    a = np.array([1.0, 2, 3])

    fmt = "\n%s" % func
    fmt += "%s\n"
    fmt += "\nInput array:\n%s\n" % a

    actual = func(a, 1)
    desired = func0(a, 1)
    err_msg = fmt % "(a, 1)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    actual = func(a, 1, axis=0)
    desired = func0(a, 1, axis=0)
    err_msg = fmt % "(a, 1, axis=0)"
    assert_array_almost_equal(actual, desired, decimal, err_msg)

    if name != "push":

        actual = func(a, 2, None)
        desired = func0(a, 2, None)
        err_msg = fmt % "(a, 2, None)"
        assert_array_almost_equal(actual, desired, decimal, err_msg)

        actual = func(a, 1, axis=None)
        desired = func0(a, 1, axis=None)
        err_msg = fmt % "(a, 1, axis=None)"
        assert_array_almost_equal(actual, desired, decimal, err_msg)

        # regression test: make sure len(kwargs) == 0 doesn't raise
        args = (a, 1, -1)
        kwargs = {}
        func(*args, **kwargs)

    else:

        # regression test: make sure len(kwargs) == 0 doesn't raise
        args = (a, 1)
        kwargs = {}
        func(*args, **kwargs)


def unit_maker_raises(func):
    """test argument parsing raises in nonreduce_axis"""
    a = np.array([1.0, 2, 3])
    assert_raises(TypeError, func)
    assert_raises(TypeError, func, axis=a)
    assert_raises(TypeError, func, a, axis=0, extra=0)
    assert_raises(TypeError, func, a, axis=0, a=a)
    if func.__name__ in ("partition", "argpartition"):
        assert_raises(TypeError, func, a, 0, 0, 0, 0, 0)
        assert_raises(TypeError, func, a, axis="0")


@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize(
    "func", (bn.partition, bn.argpartition), ids=lambda x: x.__name__
)
def test_out_of_bounds_raises(func, dtype):
    array = np.ones((10, 10), dtype=dtype)
    for axis in [None, 0, 1, -1]:
        with pytest.raises(ValueError, match="must be between"):
            func(array, 1000, axis=axis)

        with pytest.raises(ValueError, match="must be between"):
            func(array, -1, axis=axis)
