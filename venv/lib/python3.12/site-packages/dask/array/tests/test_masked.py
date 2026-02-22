from __future__ import annotations

import random
import sys
from copy import deepcopy
from itertools import product

import numpy as np
import pytest

import dask.array as da
from dask.array.numpy_compat import ComplexWarning
from dask.array.utils import assert_eq
from dask.base import tokenize
from dask.utils import typename

pytest.importorskip("dask.array.ma")


def test_tokenize_masked_array():
    m = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=10)
    m2 = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=0)
    m3 = np.ma.masked_array([1, 2, 3], mask=False, fill_value=10)
    assert tokenize(m) == tokenize(m)
    assert tokenize(m2) == tokenize(m2)
    assert tokenize(m3) == tokenize(m3)
    assert tokenize(m) != tokenize(m2)
    assert tokenize(m) != tokenize(m3)


def test_from_array_masked_array():
    m = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=10)
    dm = da.from_array(m, chunks=(2,), asarray=False)
    assert_eq(dm, m)


def test_copy_deepcopy():
    t = np.ma.masked_array([1, 2], mask=[0, 1])
    x = da.from_array(t, chunks=t.shape, asarray=False)
    # x = da.arange(5, chunks=(2,))
    y = x.copy()
    memo = {}
    y2 = deepcopy(x, memo=memo)

    xx = da.ma.masked_where([False, True], [1, 2])
    assert_eq(x, xx)

    assert_eq(y, t)
    assert isinstance(y.compute(), np.ma.masked_array)
    assert_eq(y2, t)
    assert isinstance(y2.compute(), np.ma.masked_array)


functions = [
    lambda x: x,
    lambda x: da.expm1(x),
    lambda x: 2 * x,
    lambda x: x / 2,
    lambda x: x**2,
    lambda x: x + x,
    lambda x: x * x,
    lambda x: x[0],
    lambda x: x[:, 1],
    lambda x: x[:1, None, 1:3],
    lambda x: x.T,
    lambda x: da.transpose(x, (1, 2, 0)),
    lambda x: x.sum(),
    lambda x: x.dot(np.arange(x.shape[-1])),
    lambda x: x.dot(np.eye(x.shape[-1])),
    lambda x: da.tensordot(x, np.ones(x.shape[:2]), axes=[(0, 1), (0, 1)]),
    lambda x: x.sum(axis=0),
    lambda x: x.max(axis=0),
    lambda x: x.sum(axis=(1, 2)),
    lambda x: x.astype(np.complex128),
    lambda x: x.map_blocks(lambda x: x * 2),
    lambda x: x.round(1),
    lambda x: x.reshape((x.shape[0] * x.shape[1], x.shape[2])),
    lambda x: abs(x),
    lambda x: x > 0.5,
    lambda x: x.rechunk((4, 4, 4)),
    lambda x: x.rechunk((2, 2, 1)),
]


@pytest.mark.parametrize("func", functions)
def test_basic(func):
    x = da.random.default_rng().random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.4] = 0

    y = da.ma.masked_equal(x, 0)

    xx = func(x)
    yy = func(y)

    assert_eq(xx, da.ma.filled(yy, 0))

    if yy.shape:
        zz = yy.compute()
        assert isinstance(zz, np.ma.masked_array)


def test_tensordot():
    rng = da.random.default_rng()
    x = rng.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.4] = 0
    y = rng.random((4, 3, 2), chunks=(2, 2, 1))
    y[y < 0.4] = 0

    xx = da.ma.masked_equal(x, 0)
    yy = da.ma.masked_equal(y, 0)

    assert_eq(
        da.tensordot(x, y, axes=(2, 0)),
        da.ma.filled(da.tensordot(xx, yy, axes=(2, 0)), 0),
    )
    assert_eq(
        da.tensordot(x, y, axes=(1, 1)),
        da.ma.filled(da.tensordot(xx, yy, axes=(1, 1)), 0),
    )
    assert_eq(
        da.tensordot(x, y, axes=((1, 2), (1, 0))),
        da.ma.filled(da.tensordot(xx, yy, axes=((1, 2), (1, 0))), 0),
    )


@pytest.mark.parametrize("func", functions)
@pytest.mark.filterwarnings(f"ignore::{typename(ComplexWarning)}")  # abs() in assert_eq
def test_mixed_concatenate(func):
    rng = da.random.default_rng()
    x = rng.random((2, 3, 4), chunks=(1, 2, 2))
    y = rng.random((2, 3, 4), chunks=(1, 2, 2))

    y[y < 0.4] = 0
    yy = da.ma.masked_equal(y, 0)

    d = da.concatenate([x, y], axis=0)
    s = da.concatenate([x, yy], axis=0)

    dd = func(d)
    ss = func(s)
    assert_eq(dd, ss, check_meta=False, check_type=False)


@pytest.mark.parametrize("func", functions)
@pytest.mark.filterwarnings(f"ignore::{typename(ComplexWarning)}")  # abs() in assert_eq
def test_mixed_random(func):
    d = da.random.default_rng().random((4, 3, 4), chunks=(1, 2, 2))
    d[d < 0.4] = 0

    fn = lambda x: np.ma.masked_equal(x, 0) if random.random() < 0.5 else x
    s = d.map_blocks(fn)

    dd = func(d)
    ss = func(s)

    assert_eq(dd, ss, check_meta=False, check_type=False)


def test_mixed_output_type():
    y = da.random.default_rng().random((10, 10), chunks=(5, 5))
    y[y < 0.4] = 0

    y = da.ma.masked_equal(y, 0)
    x = da.zeros((10, 1), chunks=(5, 1))

    z = da.concatenate([x, y], axis=1)
    assert z.shape == (10, 11)
    zz = z.compute()
    assert isinstance(zz, np.ma.masked_array)


def test_creation_functions():
    x = np.array([-2, -1, 0, 1, 2] * 20).reshape((10, 10))
    y = np.array([-2, 0, 1, 1, 0] * 2)
    dx = da.from_array(x, chunks=5)
    dy = da.from_array(y, chunks=4)

    sol = np.ma.masked_greater(x, y)
    for a, b in product([dx, x], [dy, y]):
        assert_eq(da.ma.masked_greater(a, b), sol)

    # These are all the same as masked_greater, just check for correct op
    assert_eq(da.ma.masked_greater(dx, 0), np.ma.masked_greater(x, 0))
    assert_eq(da.ma.masked_greater_equal(dx, 0), np.ma.masked_greater_equal(x, 0))
    assert_eq(da.ma.masked_less(dx, 0), np.ma.masked_less(x, 0))
    assert_eq(da.ma.masked_less_equal(dx, 0), np.ma.masked_less_equal(x, 0))
    assert_eq(da.ma.masked_equal(dx, 0), np.ma.masked_equal(x, 0))
    assert_eq(da.ma.masked_not_equal(dx, 0), np.ma.masked_not_equal(x, 0))

    # masked_where
    assert_eq(da.ma.masked_where(False, dx), np.ma.masked_where(False, x))
    assert_eq(da.ma.masked_where(dx > 2, dx), np.ma.masked_where(x > 2, x))

    with pytest.raises(IndexError):
        da.ma.masked_where((dx > 2)[:, 0], dx)

    assert_eq(da.ma.masked_inside(dx, -1, 1), np.ma.masked_inside(x, -1, 1))
    assert_eq(da.ma.masked_outside(dx, -1, 1), np.ma.masked_outside(x, -1, 1))
    assert_eq(da.ma.masked_values(dx, -1), np.ma.masked_values(x, -1))

    # masked_equal and masked_values in numpy sets the fill_value to `value`,
    # which can sometimes be an array. This is hard to support in dask, so we
    # forbid it. Check that this isn't supported:
    with pytest.raises(ValueError):
        da.ma.masked_equal(dx, dy)

    with pytest.raises(ValueError):
        da.ma.masked_values(dx, dy)

    y = x.astype("f8")
    y[0, 0] = y[7, 5] = np.nan
    dy = da.from_array(y, chunks=5)

    assert_eq(da.ma.masked_invalid(dy), np.ma.masked_invalid(y))

    my = np.ma.masked_greater(y, 0)
    dmy = da.ma.masked_greater(dy, 0)

    assert_eq(da.ma.fix_invalid(dmy, fill_value=0), np.ma.fix_invalid(my, fill_value=0))


def test_filled():
    x = np.array([-2, -1, 0, 1, 2] * 20).reshape((10, 10))
    dx = da.from_array(x, chunks=5)

    mx = np.ma.masked_equal(x, 0)
    mdx = da.ma.masked_equal(dx, 0)

    assert_eq(da.ma.filled(mdx), np.ma.filled(mx))
    assert_eq(da.ma.filled(mdx, -5), np.ma.filled(mx, -5))


def assert_eq_ma(a, b):
    res = a.compute()
    if res is np.ma.masked:
        assert res is b
    else:
        assert type(res) == type(b)
        if hasattr(res, "mask"):
            np.testing.assert_equal(res.mask, b.mask)
            a = da.ma.filled(a)
            b = np.ma.filled(b)
        assert_eq(a, b, equal_nan=True)


@pytest.mark.parametrize("dtype", ("i8", "f8"))
@pytest.mark.parametrize(
    "reduction", ["sum", "prod", "mean", "var", "std", "min", "max", "any", "all"]
)
def test_reductions(dtype, reduction):
    x = (np.random.default_rng(42).random((11, 11)) * 10).astype(dtype)
    dx = da.from_array(x, chunks=(4, 4))
    mx = np.ma.masked_greater(x, 5)
    mdx = da.ma.masked_greater(dx, 5)

    dfunc = getattr(da, reduction)
    func = getattr(np, reduction)

    assert_eq_ma(dfunc(mdx), func(mx))
    assert_eq_ma(dfunc(mdx, axis=0), func(mx, axis=0))
    assert_eq_ma(dfunc(mdx, keepdims=True, split_every=4), func(mx, keepdims=True))
    assert_eq_ma(dfunc(mdx, axis=0, split_every=2), func(mx, axis=0))
    assert_eq_ma(
        dfunc(mdx, axis=0, keepdims=True, split_every=2),
        func(mx, axis=0, keepdims=True),
    )
    assert_eq_ma(dfunc(mdx, axis=1, split_every=2), func(mx, axis=1))
    assert_eq_ma(
        dfunc(mdx, axis=1, keepdims=True, split_every=2),
        func(mx, axis=1, keepdims=True),
    )


@pytest.mark.parametrize("dtype", ("i8", "f8"))
@pytest.mark.parametrize(
    "reduction", ["sum", "prod", "mean", "var", "std", "min", "max", "any", "all"]
)
def test_reductions_allmasked(dtype, reduction):
    x = np.ma.masked_array([1, 2], dtype=dtype, mask=True)
    dx = da.from_array(x, asarray=False)

    dfunc = getattr(da, reduction)
    func = getattr(np, reduction)

    assert_eq_ma(dfunc(dx), func(x))


@pytest.mark.parametrize("reduction", ["argmin", "argmax"])
def test_arg_reductions(reduction):
    x = np.random.default_rng().random((10, 10, 10))
    dx = da.from_array(x, chunks=(3, 4, 5))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    dfunc = getattr(da, reduction)
    func = getattr(np, reduction)

    assert_eq_ma(dfunc(dmx), func(mx))
    assert_eq_ma(dfunc(dmx, 0), func(mx, 0))
    assert_eq_ma(dfunc(dmx, 1), func(mx, 1))
    assert_eq_ma(dfunc(dmx, 2), func(mx, 2))


def test_cumulative():
    x = np.random.default_rng(0).random((20, 24, 13))
    dx = da.from_array(x, chunks=(6, 5, 4))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    for axis in [0, 1, 2]:
        assert_eq_ma(dmx.cumsum(axis=axis), mx.cumsum(axis=axis))
        assert_eq_ma(dmx.cumprod(axis=axis), mx.cumprod(axis=axis))


def test_accessors():
    x = np.random.default_rng().random((10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    assert_eq(da.ma.getmaskarray(dmx), np.ma.getmaskarray(mx))
    assert_eq(da.ma.getmaskarray(dx), np.ma.getmaskarray(x))
    assert_eq(da.ma.getdata(dmx), np.ma.getdata(mx))
    assert_eq(da.ma.getdata(dx), np.ma.getdata(x))


def test_masked_array():
    x = np.random.default_rng().random((10, 10)).astype("f4")
    dx = da.from_array(x, chunks=(3, 4))
    f1 = da.from_array(np.array(1), chunks=())

    fill_values = [(None, None), (0.5, 0.5), (1, f1)]
    for data, (df, f) in product([x, dx], fill_values):
        assert_eq(
            da.ma.masked_array(data, fill_value=df), np.ma.masked_array(x, fill_value=f)
        )
        assert_eq(
            da.ma.masked_array(data, mask=data > 0.4, fill_value=df),
            np.ma.masked_array(x, mask=x > 0.4, fill_value=f),
        )
        assert_eq(
            da.ma.masked_array(data, mask=data > 0.4, fill_value=df),
            np.ma.masked_array(x, mask=x > 0.4, fill_value=f),
        )
        assert_eq(
            da.ma.masked_array(data, fill_value=df, dtype="f8"),
            np.ma.masked_array(x, fill_value=f, dtype="f8"),
        )

    with pytest.raises(ValueError):
        da.ma.masked_array(dx, fill_value=dx)

    with pytest.raises(np.ma.MaskError):
        da.ma.masked_array(dx, mask=dx[:3, :3])


def test_set_fill_value():
    x = np.random.default_rng().integers(0, 10, (10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    mx = np.ma.masked_greater(x, 3)
    dmx = da.ma.masked_greater(dx, 3)

    da.ma.set_fill_value(dmx, -10)
    np.ma.set_fill_value(mx, -10)
    assert_eq_ma(dmx, mx)

    da.ma.set_fill_value(dx, -10)
    np.ma.set_fill_value(x, -10)
    assert_eq_ma(dx, x)

    with pytest.raises(TypeError):
        da.ma.set_fill_value(dmx, 1e20)

    with pytest.raises(ValueError):
        da.ma.set_fill_value(dmx, dx)


@pytest.mark.parametrize("keepdims", [False, True])
def test_average_weights_with_masked_array(keepdims):
    mask = np.array([[True, False], [True, True], [False, True]])
    data = np.arange(6).reshape((3, 2))
    a = np.ma.array(data, mask=mask)
    d_a = da.ma.masked_array(data=data, mask=mask, chunks=2)

    weights = np.array([0.25, 0.75])
    d_weights = da.from_array(weights, chunks=2)

    da_avg = da.ma.average(d_a, weights=d_weights, axis=1, keepdims=keepdims)
    assert_eq(da_avg, np.ma.average(a, weights=weights, axis=1, keepdims=keepdims))


def test_arithmetic_results_in_masked():
    mask = np.array([[True, False], [True, True], [False, True]])
    x = np.arange(6).reshape((3, 2))
    masked = np.ma.array(x, mask=mask)
    dx = da.from_array(x, chunks=(2, 2))

    res = dx + masked
    sol = x + masked
    assert_eq(res, sol)
    assert isinstance(res.compute(), np.ma.masked_array)


def test_count():
    data = np.arange(120).reshape((12, 10))
    mask = (data % 3 == 0) | (data % 4 == 0)
    x = np.ma.masked_where(mask, data)
    dx = da.from_array(x, chunks=(2, 3))

    for axis in (None, 0, 1):
        res = da.ma.count(dx, axis=axis)
        sol = np.ma.count(x, axis=axis)
        assert_eq(res, sol)

    res = da.ma.count(dx, keepdims=True)
    sol = np.ma.count(x, keepdims=True)
    assert_eq(res, sol)

    # Test all masked
    x = np.ma.masked_all((12, 10))
    dx = da.from_array(x, chunks=(2, 3))
    assert_eq(da.ma.count(dx), np.ma.count(x))

    # Test on non-masked array
    x = np.arange(120).reshape((12, 10))
    dx = da.from_array(data, chunks=(2, 3))
    for axis in (None, 0, 1):
        res = da.ma.count(dx, axis=axis)
        sol = np.ma.count(x, axis=axis)
        assert_eq(res, sol, check_dtype=sys.platform != "win32")


@pytest.mark.parametrize("funcname", ["ones_like", "zeros_like", "empty_like"])
def test_like_funcs(funcname):
    mask = np.array([[True, False], [True, True], [False, True]])
    data = np.arange(6).reshape((3, 2))
    a = np.ma.array(data, mask=mask)
    d_a = da.ma.masked_array(data=data, mask=mask, chunks=2)

    da_func = getattr(da.ma, funcname)
    np_func = getattr(np.ma.core, funcname)

    res = da_func(d_a)
    sol = np_func(a)

    if "empty" in funcname:
        assert_eq(da.ma.getmaskarray(res), np.ma.getmaskarray(sol))
    else:
        assert_eq(res, sol)


def test_nonzero():
    data = np.arange(9).reshape((3, 3))
    mask = np.array([[True, False, False], [True, True, False], [True, False, True]])
    a = np.ma.array(data, mask=mask)
    d_a = da.ma.masked_array(data=data, mask=mask, chunks=2)

    for c1, c2 in [
        (a > 4, d_a > 4),
        (a, d_a),
        (a <= -2, d_a <= -2),
        (a == 0, d_a == 0),
    ]:
        sol = np.ma.nonzero(c1)
        res = da.ma.nonzero(c2)

        assert isinstance(res, type(sol))
        assert len(res) == len(sol)

        for i in range(len(sol)):
            assert_eq(res[i], sol[i])


def test_where():
    rng = np.random.default_rng()
    # Copied and adapted from the da.where test.
    x = rng.integers(10, size=(15, 14))
    mask = rng.choice(a=[False, True], size=(15, 14), p=[0.5, 0.5])
    x[5, 5] = x[4, 4] = 0  # Ensure some false elements
    d = da.ma.masked_array(x, mask=mask, chunks=(4, 5))
    x = np.ma.array(x, mask=mask)
    y = rng.integers(10, size=15).astype(np.uint8)
    e = da.from_array(y, chunks=(4,))

    # Nonzero test
    sol = np.ma.where(x)
    res = da.ma.where(d)
    for i in range(len(sol)):
        assert_eq(res[i], sol[i])

    for c1, c2 in [
        (d > 5, x > 5),
        (d, x),
        (1, 1),
        (5, 5),
        (True, True),
        (np.True_, np.True_),
        (0, 0),
        (False, False),
        (np.False_, np.False_),
    ]:
        for b1, b2 in [(0, 0), (-e[:, None], -y[:, None]), (e[:14], y[:14])]:
            w1 = da.ma.where(c1, d, b1)
            w2 = np.ma.where(c2, x, b2)
            assert_eq(w1, w2)
