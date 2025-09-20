from __future__ import annotations

import os
import warnings
from contextlib import nullcontext as does_not_warn
from itertools import permutations, zip_longest

import pytest

np = pytest.importorskip("numpy")

import itertools

import dask.array as da
import dask.config as config
from dask.array.numpy_compat import ComplexWarning
from dask.array.utils import assert_eq, same_keys
from dask.core import get_deps


@pytest.mark.parametrize("dtype", ["f4", "i4"])
@pytest.mark.parametrize("keepdims", [True, False])
@pytest.mark.parametrize("nan", [True, False])
def test_numel(dtype, keepdims, nan):
    x = np.ones((2, 3, 4))
    if nan:
        y = np.random.default_rng().uniform(-1, 1, size=(2, 3, 4))
        x[y < 0] = np.nan
        numel = da.reductions.nannumel

        def _sum(arr, **kwargs):
            n = np.sum(np.ma.masked_where(np.isnan(arr), arr), **kwargs)
            return n.filled(0) if isinstance(n, np.ma.MaskedArray) else n

    else:
        numel = da.reductions.numel
        _sum = np.sum

    assert_eq(
        numel(x, axis=(), keepdims=keepdims, dtype=dtype),
        _sum(x, axis=(), keepdims=keepdims, dtype=dtype),
    )
    assert_eq(
        numel(x, axis=0, keepdims=keepdims, dtype=dtype),
        _sum(x, axis=0, keepdims=keepdims, dtype=dtype),
    )

    for length in range(x.ndim):
        for sub in itertools.combinations([d for d in range(x.ndim)], length):
            assert_eq(
                numel(x, axis=sub, keepdims=keepdims, dtype=dtype),
                _sum(x, axis=sub, keepdims=keepdims, dtype=dtype),
            )

    for length in range(x.ndim):
        for sub in itertools.combinations([d for d in range(x.ndim)], length):
            ssub = np.random.default_rng().shuffle(list(sub))
            assert_eq(
                numel(x, axis=ssub, keepdims=keepdims, dtype=dtype),
                _sum(x, axis=ssub, keepdims=keepdims, dtype=dtype),
            )


def reduction_0d_test(da_func, darr, np_func, narr):
    expected = np_func(narr)
    actual = da_func(darr)

    assert_eq(actual, expected)
    assert_eq(da_func(narr), expected)  # Ensure Dask reductions work with NumPy arrays
    assert actual.size == 1


def test_reductions_0D():
    x = np.int_(3)  # np.int_ has a dtype attribute, np.int does not.
    a = da.from_array(x, chunks=(1,))

    reduction_0d_test(da.sum, a, np.sum, x)
    reduction_0d_test(da.prod, a, np.prod, x)
    reduction_0d_test(da.mean, a, np.mean, x)
    reduction_0d_test(da.var, a, np.var, x)
    reduction_0d_test(da.std, a, np.std, x)
    reduction_0d_test(da.min, a, np.min, x)
    reduction_0d_test(da.max, a, np.max, x)
    reduction_0d_test(da.any, a, np.any, x)
    reduction_0d_test(da.all, a, np.all, x)

    reduction_0d_test(da.nansum, a, np.nansum, x)
    reduction_0d_test(da.nanprod, a, np.nanprod, x)
    reduction_0d_test(da.nanmean, a, np.mean, x)
    reduction_0d_test(da.nanvar, a, np.var, x)
    reduction_0d_test(da.nanstd, a, np.std, x)
    reduction_0d_test(da.nanmin, a, np.nanmin, x)
    reduction_0d_test(da.nanmax, a, np.nanmax, x)


def reduction_1d_test(da_func, darr, np_func, narr, use_dtype=True, split_every=True):
    assert_eq(da_func(darr), np_func(narr))
    assert_eq(
        da_func(narr), np_func(narr)
    )  # Ensure Dask reductions work with NumPy arrays
    assert_eq(da_func(darr, keepdims=True), np_func(narr, keepdims=True))
    assert_eq(da_func(darr, axis=()), np_func(narr, axis=()))
    assert same_keys(da_func(darr), da_func(darr))
    assert same_keys(da_func(darr, keepdims=True), da_func(darr, keepdims=True))
    if use_dtype:
        with pytest.warns(ComplexWarning) if np.iscomplexobj(narr) else does_not_warn():
            assert_eq(da_func(darr, dtype="f8"), np_func(narr, dtype="f8"))
            assert_eq(da_func(darr, dtype="i8"), np_func(narr, dtype="i8"))
            assert same_keys(da_func(darr, dtype="i8"), da_func(darr, dtype="i8"))
    if split_every:
        a1 = da_func(darr, split_every=2)
        a2 = da_func(darr, split_every={0: 2})
        assert same_keys(a1, a2)
        assert_eq(a1, np_func(narr))
        assert_eq(a2, np_func(narr))
        assert_eq(
            da_func(darr, keepdims=True, split_every=2), np_func(narr, keepdims=True)
        )


@pytest.mark.parametrize("dtype", ["f4", "i4", "c8"])
def test_reductions_1D(dtype):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ComplexWarning)
        x = (np.arange(5) + 1j * np.arange(5)).astype(dtype)
        a = da.from_array(x, chunks=(2,))

    reduction_1d_test(da.sum, a, np.sum, x)
    reduction_1d_test(da.prod, a, np.prod, x)
    reduction_1d_test(da.mean, a, np.mean, x)
    reduction_1d_test(da.var, a, np.var, x)
    reduction_1d_test(da.std, a, np.std, x)
    reduction_1d_test(da.min, a, np.min, x, False)
    reduction_1d_test(da.max, a, np.max, x, False)
    reduction_1d_test(da.any, a, np.any, x, False)
    reduction_1d_test(da.all, a, np.all, x, False)

    reduction_1d_test(da.nansum, a, np.nansum, x)
    reduction_1d_test(da.nanprod, a, np.nanprod, x)
    reduction_1d_test(da.nanmean, a, np.mean, x)
    reduction_1d_test(da.nanvar, a, np.var, x)
    reduction_1d_test(da.nanstd, a, np.std, x)
    reduction_1d_test(da.nanmin, a, np.nanmin, x, False)
    reduction_1d_test(da.nanmax, a, np.nanmax, x, False)


def reduction_2d_test(da_func, darr, np_func, narr, use_dtype=True, split_every=True):
    assert_eq(da_func(darr), np_func(narr))
    assert_eq(da_func(darr, keepdims=True), np_func(narr, keepdims=True))
    assert_eq(da_func(darr, axis=()), np_func(narr, axis=()))
    assert_eq(da_func(darr, axis=0), np_func(narr, axis=0))
    assert_eq(da_func(darr, axis=1), np_func(narr, axis=1))
    assert_eq(da_func(darr, axis=-1), np_func(narr, axis=-1))
    assert_eq(da_func(darr, axis=-2), np_func(narr, axis=-2))
    assert_eq(
        da_func(darr, axis=1, keepdims=True), np_func(narr, axis=1, keepdims=True)
    )
    assert_eq(
        da_func(darr, axis=(), keepdims=True), np_func(narr, axis=(), keepdims=True)
    )
    assert_eq(da_func(darr, axis=(1, 0)), np_func(narr, axis=(1, 0)))

    assert same_keys(da_func(darr, axis=()), da_func(darr, axis=()))
    assert same_keys(da_func(darr, axis=1), da_func(darr, axis=1))
    assert same_keys(da_func(darr, axis=(1, 0)), da_func(darr, axis=(1, 0)))

    if use_dtype:
        with pytest.warns(ComplexWarning) if np.iscomplexobj(narr) else does_not_warn():
            assert_eq(da_func(darr, dtype="f8"), np_func(narr, dtype="f8"))
            assert_eq(da_func(darr, dtype="i8"), np_func(narr, dtype="i8"))

    if split_every:
        a1 = da_func(darr, split_every=4)
        a2 = da_func(darr, split_every={0: 2, 1: 2})
        assert same_keys(a1, a2)
        assert_eq(a1, np_func(narr))
        assert_eq(a2, np_func(narr))
        assert_eq(
            da_func(darr, keepdims=True, split_every=4),
            np_func(narr, keepdims=True),
        )
        assert_eq(da_func(darr, axis=(), split_every=2), np_func(narr, axis=()))
        assert_eq(da_func(darr, axis=0, split_every=2), np_func(narr, axis=0))
        assert_eq(
            da_func(darr, axis=(), keepdims=True, split_every=2),
            np_func(narr, axis=(), keepdims=True),
        )
        assert_eq(
            da_func(darr, axis=0, keepdims=True, split_every=2),
            np_func(narr, axis=0, keepdims=True),
        )
        assert_eq(da_func(darr, axis=1, split_every=2), np_func(narr, axis=1))
        assert_eq(
            da_func(darr, axis=1, keepdims=True, split_every=2),
            np_func(narr, axis=1, keepdims=True),
        )


def test_reduction_errors():
    x = da.ones((5, 5), chunks=(3, 3))
    with pytest.raises(ValueError):
        x.sum(axis=2)
    with pytest.raises(ValueError):
        x.sum(axis=-3)


@pytest.mark.slow
@pytest.mark.parametrize("dtype", ["f4", "i4", "c8"])
def test_reductions_2D(dtype):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ComplexWarning)
        x = (np.arange(1, 122) + 1j * np.arange(1, 122)).reshape((11, 11)).astype(dtype)
    a = da.from_array(x, chunks=(4, 4))

    b = a.sum(keepdims=True)
    assert b.__dask_keys__() == [[(b.name, 0, 0)]]

    reduction_2d_test(da.sum, a, np.sum, x)
    reduction_2d_test(da.mean, a, np.mean, x)
    reduction_2d_test(da.var, a, np.var, x, False)  # Difference in dtype algo
    reduction_2d_test(da.std, a, np.std, x, False)  # Difference in dtype algo
    reduction_2d_test(da.min, a, np.min, x, False)
    reduction_2d_test(da.max, a, np.max, x, False)
    reduction_2d_test(da.any, a, np.any, x, False)
    reduction_2d_test(da.all, a, np.all, x, False)

    reduction_2d_test(da.nansum, a, np.nansum, x)
    reduction_2d_test(da.nanmean, a, np.mean, x)
    reduction_2d_test(da.nanvar, a, np.nanvar, x, False)  # Difference in dtype algo
    reduction_2d_test(da.nanstd, a, np.nanstd, x, False)  # Difference in dtype algo
    reduction_2d_test(da.nanmin, a, np.nanmin, x, False)
    reduction_2d_test(da.nanmax, a, np.nanmax, x, False)

    # prod/nanprod overflow for data at this size, leading to warnings about
    # overflow/invalid values.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        reduction_2d_test(da.prod, a, np.prod, x)
        reduction_2d_test(da.nanprod, a, np.nanprod, x)


@pytest.mark.parametrize(
    ["dfunc", "func"],
    [
        (da.argmin, np.argmin),
        (da.argmax, np.argmax),
        (da.nanargmin, np.nanargmin),
        (da.nanargmax, np.nanargmax),
    ],
)
def test_arg_reductions(dfunc, func):
    x = np.random.default_rng().random((10, 10, 10))
    a = da.from_array(x, chunks=(3, 4, 5))

    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))
    assert_eq(dfunc(a, 1), func(x, 1))
    assert_eq(dfunc(a, 2), func(x, 2))
    with config.set(split_every=2):
        assert_eq(dfunc(a), func(x))
        assert_eq(dfunc(a, 0), func(x, 0))
        assert_eq(dfunc(a, 1), func(x, 1))
        assert_eq(dfunc(a, 2), func(x, 2))
    assert_eq(dfunc(a, keepdims=True), func(x, keepdims=True))

    pytest.raises(ValueError, lambda: dfunc(a, 3))
    pytest.raises(TypeError, lambda: dfunc(a, (0, 1)))

    x2 = np.arange(10)
    a2 = da.from_array(x2, chunks=3)
    assert_eq(dfunc(a2), func(x2))
    assert_eq(dfunc(a2, 0), func(x2, 0))
    assert_eq(dfunc(a2, 0, split_every=2), func(x2, 0))

    x3 = np.array(1)
    a3 = da.from_array(x3)
    assert_eq(dfunc(a3), func(x3))


@pytest.mark.parametrize(
    ["dfunc", "func"], [(da.nanmin, np.nanmin), (da.nanmax, np.nanmax)]
)
def test_nan_reduction_warnings(dfunc, func):
    x = np.random.default_rng().random((10, 10, 10))
    x[5] = np.nan
    a = da.from_array(x, chunks=(3, 4, 5))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # All-NaN slice encountered
        expected = func(x, 1)
    assert_eq(dfunc(a, 1), expected)


@pytest.mark.parametrize(
    ["dfunc", "func"], [(da.nanargmin, np.nanargmin), (da.nanargmax, np.nanargmax)]
)
def test_nanarg_reductions(dfunc, func):
    x = np.random.default_rng().random((10, 10, 10))
    x[5] = np.nan
    a = da.from_array(x, chunks=(3, 4, 5))
    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # All-NaN slice encountered
        with pytest.raises(ValueError):
            dfunc(a, 1).compute()

        with pytest.raises(ValueError):
            dfunc(a, 2).compute()

        x[:] = np.nan
        a = da.from_array(x, chunks=(3, 4, 5))
        with pytest.raises(ValueError):
            dfunc(a).compute()


@pytest.mark.parametrize(["dfunc", "func"], [(da.min, np.min), (da.max, np.max)])
def test_min_max_empty_chunks(dfunc, func):
    x1 = np.arange(10)
    a1 = da.from_array(x1, chunks=1)
    assert_eq(dfunc(a1[a1 < 2]), func(x1[x1 < 2]))

    x2 = np.arange(10)
    a2 = da.from_array(x2, chunks=((5, 0, 5),))
    assert_eq(dfunc(a2), func(x2))

    x3 = np.array([[1, 1, 2, 3], [1, 1, 4, 0]])
    a3 = da.from_array(x3, chunks=1)
    assert_eq(dfunc(a3[a3 >= 2]), func(x3[x3 >= 2]))

    a4 = da.arange(10)
    with pytest.raises(
        ValueError
    ):  # Checking it mimics numpy behavior when all chunks are empty
        dfunc(a4[a4 < 0]).compute()


@pytest.mark.parametrize("func", ["argmax", "nanargmax"])
def test_arg_reductions_unknown_chunksize(func):
    x = da.arange(10, chunks=5)
    x = x[x > 1]

    with pytest.raises(ValueError) as info:
        getattr(da, func)(x)

    assert "unknown chunksize" in str(info.value)


@pytest.mark.parametrize("func", ["argmax", "nanargmax"])
def test_arg_reductions_unknown_chunksize_2d(func):
    x = da.ones((10, 10), chunks=(5, 5))
    x = x[x[0, :] > 0, :]  # unknown chunks in first dimension only

    with pytest.raises(ValueError):
        getattr(da, func)(x, axis=0)

    getattr(da, func)(x, axis=1).compute()


@pytest.mark.parametrize("func", ["argmax", "nanargmax"])
def test_arg_reductions_unknown_single_chunksize(func):
    x = da.ones((10, 10), chunks=(10, 10))
    x = x[x[0, :] > 0, :]  # unknown chunks in first dimension only

    getattr(da, func)(x, axis=0).compute()
    getattr(da, func)(x, axis=1).compute()


def test_reductions_2D_nans():
    # chunks are a mix of some/all/no NaNs
    x = np.full((4, 4), np.nan)
    x[:2, :2] = np.array([[1, 2], [3, 4]])
    x[2, 2] = 5
    x[3, 3] = 6
    a = da.from_array(x, chunks=(2, 2))

    reduction_2d_test(da.sum, a, np.sum, x, False, False)
    reduction_2d_test(da.prod, a, np.prod, x, False, False)
    reduction_2d_test(da.mean, a, np.mean, x, False, False)
    reduction_2d_test(da.var, a, np.var, x, False, False)
    reduction_2d_test(da.std, a, np.std, x, False, False)
    reduction_2d_test(da.min, a, np.min, x, False, False)
    reduction_2d_test(da.max, a, np.max, x, False, False)
    reduction_2d_test(da.any, a, np.any, x, False, False)
    reduction_2d_test(da.all, a, np.all, x, False, False)

    reduction_2d_test(da.nansum, a, np.nansum, x, False, False)
    reduction_2d_test(da.nanprod, a, np.nanprod, x, False, False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        reduction_2d_test(da.nanmean, a, np.nanmean, x, False, False)
        reduction_2d_test(da.nanvar, a, np.nanvar, x, False, False)
        reduction_2d_test(da.nanstd, a, np.nanstd, x, False, False)
        reduction_2d_test(da.nanmin, a, np.nanmin, x, False, False)
        reduction_2d_test(da.nanmax, a, np.nanmax, x, False, False)

        assert_eq(da.argmax(a), np.argmax(x))
        assert_eq(da.argmin(a), np.argmin(x))
        assert_eq(da.nanargmax(a), np.nanargmax(x))
        assert_eq(da.nanargmin(a), np.nanargmin(x))

        assert_eq(da.argmax(a, axis=0), np.argmax(x, axis=0))
        assert_eq(da.argmin(a, axis=0), np.argmin(x, axis=0))
        assert_eq(da.nanargmax(a, axis=0), np.nanargmax(x, axis=0))
        assert_eq(da.nanargmin(a, axis=0), np.nanargmin(x, axis=0))

        assert_eq(da.argmax(a, axis=1), np.argmax(x, axis=1))
        assert_eq(da.argmin(a, axis=1), np.argmin(x, axis=1))
        assert_eq(da.nanargmax(a, axis=1), np.nanargmax(x, axis=1))
        assert_eq(da.nanargmin(a, axis=1), np.nanargmin(x, axis=1))


def test_moment():
    def moment(x, n, axis=None):
        return ((x - x.mean(axis=axis, keepdims=True)) ** n).sum(
            axis=axis
        ) / np.ones_like(x).sum(axis=axis)

    # Poorly conditioned
    x = np.array([1.0, 2.0, 3.0] * 10).reshape((3, 10)) + 1e8
    a = da.from_array(x, chunks=5)
    assert_eq(a.moment(2), moment(x, 2))
    assert_eq(a.moment(3), moment(x, 3))
    assert_eq(a.moment(4), moment(x, 4))

    x = np.arange(1, 122).reshape((11, 11)).astype("f8")
    a = da.from_array(x, chunks=(4, 4))
    assert_eq(a.moment(4, axis=1), moment(x, 4, axis=1))
    assert_eq(a.moment(4, axis=(1, 0)), moment(x, 4, axis=(1, 0)))

    # Tree reduction
    assert_eq(a.moment(order=4, split_every=4), moment(x, 4))
    assert_eq(a.moment(order=4, axis=0, split_every=4), moment(x, 4, axis=0))
    assert_eq(a.moment(order=4, axis=1, split_every=4), moment(x, 4, axis=1))


def test_reductions_with_negative_axes():
    x = np.random.default_rng().random((4, 4, 4))
    a = da.from_array(x, chunks=2)

    assert_eq(a.argmin(axis=-1), x.argmin(axis=-1))
    assert_eq(a.argmin(axis=-1, split_every=2), x.argmin(axis=-1))

    assert_eq(a.sum(axis=-1), x.sum(axis=-1))
    assert_eq(a.sum(axis=(0, -1)), x.sum(axis=(0, -1)))


def test_nan():
    x = np.array([[1, np.nan, 3, 4], [5, 6, 7, np.nan], [9, 10, 11, 12]])
    d = da.from_array(x, chunks=(2, 2))

    assert_eq(np.nansum(x), da.nansum(d))
    assert_eq(np.nansum(x, axis=0), da.nansum(d, axis=0))
    assert_eq(np.nanmean(x, axis=1), da.nanmean(d, axis=1))
    assert_eq(np.nanmin(x, axis=1), da.nanmin(d, axis=1))
    assert_eq(np.nanmax(x, axis=(0, 1)), da.nanmax(d, axis=(0, 1)))
    assert_eq(np.nanvar(x), da.nanvar(d))
    assert_eq(np.nanstd(x, axis=0), da.nanstd(d, axis=0))
    assert_eq(np.nanargmin(x, axis=0), da.nanargmin(d, axis=0))
    assert_eq(np.nanargmax(x, axis=0), da.nanargmax(d, axis=0))
    assert_eq(np.nanprod(x), da.nanprod(d))


@pytest.mark.parametrize("func", ["nansum", "sum", "nanmin", "min", "nanmax", "max"])
def test_nan_object(func):
    with warnings.catch_warnings():
        if os.name == "nt" and func in {"min", "max"}:
            # RuntimeWarning: invalid value encountered in reduce in wrapreduction
            # from NumPy.
            warnings.simplefilter("ignore", RuntimeWarning)

        x = np.array([[1, np.nan, 3, 4], [5, 6, 7, np.nan], [9, 10, 11, 12]]).astype(
            object
        )
        d = da.from_array(x, chunks=(2, 2))

        if func in {"nanmin", "nanmax"}:
            warnings.simplefilter("ignore", RuntimeWarning)

        assert_eq(getattr(np, func)(x, axis=()), getattr(da, func)(d, axis=()))

        if func in {"nanmin", "nanmax"}:
            warnings.simplefilter("default", RuntimeWarning)

        if func in {"min", "max"}:
            warnings.simplefilter("ignore", RuntimeWarning)
        assert_eq(getattr(np, func)(x, axis=0), getattr(da, func)(d, axis=0))
        if os.name != "nt" and func in {"min", "max"}:
            warnings.simplefilter("default", RuntimeWarning)

        assert_eq(getattr(np, func)(x, axis=1), getattr(da, func)(d, axis=1))
        # wrap the scalar in a numpy array since the dask version cannot know dtype
        assert_eq(np.array(getattr(np, func)(x)).astype(object), getattr(da, func)(d))


def test_0d_array():
    x = da.mean(da.ones(4, chunks=4), axis=()).compute()
    x = da.mean(da.ones(4, chunks=4), axis=0).compute()
    y = np.mean(np.ones(4))
    assert type(x) == type(y)

    x = da.sum(da.zeros(4, chunks=1)).compute()
    y = np.sum(np.zeros(4))
    assert type(x) == type(y)


def test_reduction_on_scalar():
    x = da.from_array(np.array(1.0), chunks=())
    assert (x == x).all()


def test_reductions_with_empty_array():
    dx1 = da.ones((10, 0, 5), chunks=4)
    x1 = dx1.compute()
    dx2 = da.ones((0, 0, 0), chunks=4)
    x2 = dx2.compute()

    for dx, x in [(dx1, x1), (dx2, x2)]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # Mean of empty slice
            assert_eq(dx.mean(), x.mean())
            assert_eq(dx.mean(axis=()), x.mean(axis=()))
            assert_eq(dx.mean(axis=0), x.mean(axis=0))
            assert_eq(dx.mean(axis=1), x.mean(axis=1))
            assert_eq(dx.mean(axis=2), x.mean(axis=2))


def assert_max_deps(x, n, eq=True):
    dependencies, dependents = get_deps(x.dask)
    if eq:
        assert max(map(len, dependencies.values())) == n
    else:
        assert max(map(len, dependencies.values())) <= n


def test_tree_reduce_depth():
    # 2D
    x = da.from_array(np.arange(242).reshape((11, 22)), chunks=(3, 4))
    thresh = {0: 2, 1: 3}
    assert_max_deps(x.sum(split_every=thresh), 2 * 3)
    assert_max_deps(x.sum(axis=(), split_every=thresh), 1)
    assert_max_deps(x.sum(axis=0, split_every=thresh), 2)
    assert_max_deps(x.sum(axis=1, split_every=thresh), 3)
    assert_max_deps(x.sum(split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(), split_every=20), 1)
    assert_max_deps(x.sum(axis=0, split_every=20), 4)
    assert_max_deps(x.sum(axis=1, split_every=20), 6)

    # 3D
    x = da.from_array(np.arange(11 * 22 * 29).reshape((11, 22, 29)), chunks=(3, 4, 5))
    thresh = {0: 2, 1: 3, 2: 4}
    assert_max_deps(x.sum(split_every=thresh), 2 * 3 * 4)
    assert_max_deps(x.sum(axis=(), split_every=thresh), 1)
    assert_max_deps(x.sum(axis=0, split_every=thresh), 2)
    assert_max_deps(x.sum(axis=1, split_every=thresh), 3)
    assert_max_deps(x.sum(axis=2, split_every=thresh), 4)
    assert_max_deps(x.sum(axis=(0, 1), split_every=thresh), 2 * 3)
    assert_max_deps(x.sum(axis=(0, 2), split_every=thresh), 2 * 4)
    assert_max_deps(x.sum(axis=(1, 2), split_every=thresh), 3 * 4)
    assert_max_deps(x.sum(split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(), split_every=20), 1)
    assert_max_deps(x.sum(axis=0, split_every=20), 4)
    assert_max_deps(x.sum(axis=1, split_every=20), 6)
    assert_max_deps(x.sum(axis=2, split_every=20), 6)
    assert_max_deps(x.sum(axis=(0, 1), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(0, 2), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(1, 2), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(0, 1), split_every=40), 4 * 6)
    assert_max_deps(x.sum(axis=(0, 2), split_every=40), 4 * 6)
    assert_max_deps(x.sum(axis=(1, 2), split_every=40), 6 * 6)


def test_tree_reduce_set_options():
    x = da.from_array(np.arange(242).reshape((11, 22)), chunks=(3, 4))
    with config.set(split_every={0: 2, 1: 3}):
        assert_max_deps(x.sum(), 2 * 3)
        assert_max_deps(x.sum(axis=()), 1)
        assert_max_deps(x.sum(axis=0), 2)


def test_reduction_names():
    x = da.ones(5, chunks=(2,))
    assert x.sum().name.startswith("sum")
    assert "max" in x.max().name.split("-")[0]
    assert x.var().name.startswith("var")
    assert x.all().name.startswith("all")
    assert any(k[0].startswith("nansum") for k in da.nansum(x).dask)
    assert x.mean().name.startswith("mean")


def test_general_reduction_names():
    dtype = int
    a = da.reduction(
        da.ones(10, dtype, chunks=2), np.sum, np.sum, dtype=dtype, name="foo"
    )
    names, tokens = list(zip_longest(*[key[0].rsplit("-", 1) for key in a.dask]))
    assert set(names) == {"ones_like", "foo", "foo-partial", "foo-aggregate"}
    assert all(tokens)


@pytest.mark.parametrize("func", [np.sum, np.argmax])
def test_array_reduction_out(func):
    x = da.arange(10, chunks=(5,))
    y = da.ones((10, 10), chunks=(4, 4))
    func(y, axis=0, out=x)
    assert_eq(x, func(np.ones((10, 10)), axis=0))


@pytest.mark.parametrize("func", ["cumsum", "cumprod", "nancumsum", "nancumprod"])
@pytest.mark.parametrize("use_nan", [False, True])
@pytest.mark.parametrize("axis", [None, 0, 1, -1])
@pytest.mark.parametrize("method", ["sequential", "blelloch"])
def test_array_cumreduction_axis(func, use_nan, axis, method):
    np_func = getattr(np, func)
    da_func = getattr(da, func)

    s = (10, 11, 12)
    a = np.arange(np.prod(s), dtype=float).reshape(s)
    if use_nan:
        a[1] = np.nan
    d = da.from_array(a, chunks=(4, 5, 6))
    if func in ["cumprod", "nancumprod"] and method == "blelloch" and axis is None:
        with pytest.warns(RuntimeWarning):
            da_func(d, axis=axis, method=method).compute()
            return

    a_r = np_func(a, axis=axis)
    d_r = da_func(d, axis=axis, method=method)

    assert_eq(a_r, d_r)


@pytest.mark.parametrize("func", [np.cumsum, np.cumprod])
def test_array_cumreduction_out(func):
    x = da.ones((10, 10), chunks=(4, 4))
    func(x, axis=0, out=x)
    assert_eq(x, func(np.ones((10, 10)), axis=0))


@pytest.mark.parametrize(
    "npfunc,daskfunc", [(np.sort, da.topk), (np.argsort, da.argtopk)]
)
@pytest.mark.parametrize("split_every", [None, 2, 4, 8])
def test_topk_argtopk1(npfunc, daskfunc, split_every):
    # Test data
    k = 5
    # Test at least 3 levels of aggregation when split_every=2
    # to stress the different chunk, combine, aggregate kernels
    rng = np.random.default_rng()
    npa = rng.random(800)
    npb = rng.random((10, 20, 30))

    a = da.from_array(npa, chunks=((120, 80, 100, 200, 300),))
    b = da.from_array(npb, chunks=(4, 8, 8))

    # 1-dimensional arrays
    # top 5 elements, sorted descending
    assert_eq(npfunc(npa)[-k:][::-1], daskfunc(a, k, split_every=split_every))
    # bottom 5 elements, sorted ascending
    assert_eq(npfunc(npa)[:k], daskfunc(a, -k, split_every=split_every))

    # n-dimensional arrays
    # also testing when k > chunk
    # top 5 elements, sorted descending
    assert_eq(
        npfunc(npb, axis=0)[-k:, :, :][::-1, :, :],
        daskfunc(b, k, axis=0, split_every=split_every),
    )
    assert_eq(
        npfunc(npb, axis=1)[:, -k:, :][:, ::-1, :],
        daskfunc(b, k, axis=1, split_every=split_every),
    )
    assert_eq(
        npfunc(npb, axis=-1)[:, :, -k:][:, :, ::-1],
        daskfunc(b, k, axis=-1, split_every=split_every),
    )
    with pytest.raises(ValueError):
        daskfunc(b, k, axis=3, split_every=split_every)

    # bottom 5 elements, sorted ascending
    assert_eq(
        npfunc(npb, axis=0)[:k, :, :], daskfunc(b, -k, axis=0, split_every=split_every)
    )
    assert_eq(
        npfunc(npb, axis=1)[:, :k, :], daskfunc(b, -k, axis=1, split_every=split_every)
    )
    assert_eq(
        npfunc(npb, axis=-1)[:, :, :k],
        daskfunc(b, -k, axis=-1, split_every=split_every),
    )
    with pytest.raises(ValueError):
        daskfunc(b, -k, axis=3, split_every=split_every)


@pytest.mark.parametrize(
    "npfunc,daskfunc", [(np.sort, da.topk), (np.argsort, da.argtopk)]
)
@pytest.mark.parametrize("split_every", [None, 2, 3, 4])
@pytest.mark.parametrize("chunksize", [1, 2, 3, 4, 5, 10])
def test_topk_argtopk2(npfunc, daskfunc, split_every, chunksize):
    """Fine test use cases when k is larger than chunk size"""
    npa = np.random.default_rng().random((10,))
    a = da.from_array(npa, chunks=chunksize)
    k = 5

    # top 5 elements, sorted descending
    assert_eq(npfunc(npa)[-k:][::-1], daskfunc(a, k, split_every=split_every))
    # bottom 5 elements, sorted ascending
    assert_eq(npfunc(npa)[:k], daskfunc(a, -k, split_every=split_every))


def test_topk_argtopk3():
    a = da.random.default_rng().random((10, 20, 30), chunks=(4, 8, 8))

    # As Array methods
    assert_eq(a.topk(5, axis=1, split_every=2), da.topk(a, 5, axis=1, split_every=2))
    assert_eq(
        a.argtopk(5, axis=1, split_every=2), da.argtopk(a, 5, axis=1, split_every=2)
    )


@pytest.mark.parametrize(
    "func",
    [da.cumsum, da.cumprod, da.argmin, da.argmax, da.min, da.max, da.nansum, da.nanmax],
)
@pytest.mark.parametrize("method", ["sequential", "blelloch"])
def test_regres_3940(func, method):
    if func in {da.cumsum, da.cumprod}:
        kwargs = {"method": method}
    else:
        kwargs = {}
    a = da.ones((5, 2), chunks=(2, 2))
    assert func(a, **kwargs).name != func(a + 1, **kwargs).name
    assert func(a, axis=0, **kwargs).name != func(a, **kwargs).name
    assert func(a, axis=0, **kwargs).name != func(a, axis=1, **kwargs).name
    if func not in {da.cumsum, da.cumprod, da.argmin, da.argmax}:
        assert func(a, axis=()).name != func(a).name
        assert func(a, axis=()).name != func(a, axis=0).name


def test_trace():
    def _assert(a, b, *args, **kwargs):
        return assert_eq(a.trace(*args, **kwargs), b.trace(*args, **kwargs))

    b = np.arange(12).reshape((3, 4))
    a = da.from_array(b, 1)
    _assert(a, b)
    _assert(a, b, 0)
    _assert(a, b, 1)
    _assert(a, b, -1)

    b = np.arange(8).reshape((2, 2, 2))
    a = da.from_array(b, 2)
    _assert(a, b)
    _assert(a, b, 0)
    _assert(a, b, 1)
    _assert(a, b, -1)
    _assert(a, b, 0, 0, 1)
    _assert(a, b, 0, 0, 2)
    _assert(a, b, 0, 1, 2, int)
    _assert(a, b, 0, 1, 2, float)
    _assert(a, b, offset=1, axis1=0, axis2=2, dtype=int)
    _assert(a, b, offset=1, axis1=0, axis2=2, dtype=float)


@pytest.mark.parametrize("func", ["median", "nanmedian"])
@pytest.mark.parametrize("axis", [0, [0, 1], 1, -1])
@pytest.mark.parametrize("keepdims", [True, False])
def test_median(axis, keepdims, func):
    x = np.arange(100).reshape((2, 5, 10))
    d = da.from_array(x, chunks=2)
    assert_eq(
        getattr(da, func)(d, axis=axis, keepdims=keepdims),
        getattr(np, func)(x, axis=axis, keepdims=keepdims),
    )


@pytest.mark.parametrize("func", ["median", "nanmedian"])
@pytest.mark.parametrize("axis", [0, [0, 2], 1])
def test_median_does_not_rechunk_if_whole_axis_in_one_chunk(axis, func):
    x = np.arange(100).reshape((2, 5, 10))
    d = da.from_array(x, chunks=(2, 1, 10))

    actual = getattr(da, func)(d, axis=axis)
    expected = getattr(np, func)(x, axis=axis)
    assert_eq(actual, expected)
    does_rechunk = "rechunk" in str(dict(actual.__dask_graph__()))
    if axis == 1:
        assert does_rechunk
    else:
        assert not does_rechunk


@pytest.mark.parametrize("method", ["sum", "mean", "prod"])
def test_object_reduction(method):
    arr = da.ones(1).astype(object)
    result = getattr(arr, method)().compute()
    assert result == 1


@pytest.mark.parametrize("func", ["nanmin", "nanmax"])
def test_empty_chunk_nanmin_nanmax(func):
    # see https://github.com/dask/dask/issues/8352
    x = np.arange(10).reshape(2, 5)
    d = da.from_array(x, chunks=2)
    x = x[x > 4]
    d = d[d > 4]
    block_lens = np.array([len(x.compute()) for x in d.blocks])
    assert 0 in block_lens
    with pytest.raises(ValueError) as err:
        getattr(da, func)(d)
    assert "Arrays chunk sizes are unknown" in str(err)
    d = d.compute_chunk_sizes()
    assert_eq(getattr(da, func)(d), getattr(np, func)(x))


@pytest.mark.parametrize("func", ["nanmin", "nanmax"])
def test_empty_chunk_nanmin_nanmax_raise(func):
    # see https://github.com/dask/dask/issues/8352
    x = np.arange(10).reshape(2, 5)
    d = da.from_array(x, chunks=2)
    d = d[d > 9]
    x = x[x > 9]
    d = d.compute_chunk_sizes()
    with pytest.raises(ValueError) as err_np:
        getattr(np, func)(x)
    with pytest.raises(ValueError) as err_da:
        d = getattr(da, func)(d)
        d.compute()
    assert str(err_np.value) == str(err_da.value)


def test_mean_func_does_not_warn():
    # non-regression test for https://github.com/pydata/xarray/issues/5151
    xr = pytest.importorskip("xarray")
    a = xr.DataArray(da.from_array(np.full((10, 10), np.nan)))

    with warnings.catch_warnings(record=True) as rec:
        a.mean().compute()
    assert not rec  # did not warn


@pytest.mark.parametrize("func", ["nanvar", "nanstd"])
def test_nan_func_does_not_warn(func):
    # non-regression test for #6105
    x = np.ones((10,)) * np.nan
    x[0] = 1
    x[1] = 2
    d = da.from_array(x, chunks=2)
    with warnings.catch_warnings(record=True) as rec:
        getattr(da, func)(d).compute()
    assert not rec  # did not warn


@pytest.mark.parametrize("chunks", list(permutations(((2, 1) * 8, (3,) * 8, (6,) * 4))))
@pytest.mark.parametrize("split_every", [2, 4])
@pytest.mark.parametrize(
    "axes", list(permutations((0, 1, 2), 2)) + list(permutations((0, 1, 2)))
)
def test_chunk_structure_independence(axes, split_every, chunks):
    # Reducing an array should not depend on its chunk-structure!!!
    # See Issue #8541: https://github.com/dask/dask/issues/8541
    shape = tuple(np.sum(s) for s in chunks)
    np_array = np.arange(np.prod(shape)).reshape(*shape)
    x = da.from_array(np_array, chunks=chunks)
    reduced_x = da.reduction(
        x,
        lambda x, axis, keepdims: x,
        lambda x, axis, keepdims: x,
        keepdims=True,
        axis=axes,
        split_every=split_every,
        dtype=x.dtype,
        meta=x._meta,
    )
    assert_eq(reduced_x, np_array, check_chunks=False, check_shape=False)


def test_weighted_reduction():
    # Weighted reduction
    def w_sum(x, weights=None, dtype=None, computing_meta=False, **kwargs):
        """`chunk` callable for (weighted) sum"""
        if computing_meta:
            return x
        if weights is not None:
            x = x * weights
        return np.sum(x, dtype=dtype, **kwargs)

    # Arrays
    a = 1 + np.ma.arange(60).reshape(6, 10)
    a[2, 2] = np.ma.masked
    dx = da.from_array(a, chunks=(4, 5))
    # Weights
    w = np.linspace(1, 2, 6).reshape(6, 1)

    # No weights (i.e. normal sum)
    x = da.reduction(dx, w_sum, np.sum, dtype=dx.dtype)
    assert_eq(x, np.sum(a), check_shape=True)

    # Weighted sum
    x = da.reduction(dx, w_sum, np.sum, dtype="f8", weights=w)
    assert_eq(x, np.sum(a * w), check_shape=True)

    # Non-broadcastable weights (short axis)
    with pytest.raises(ValueError):
        da.reduction(dx, w_sum, np.sum, weights=[1, 2, 3])

    # Non-broadcastable weights (too many dims)
    with pytest.raises(ValueError):
        da.reduction(dx, w_sum, np.sum, weights=[[[2]]])


def test_cumreduction_no_rechunk_on_1d_array():
    x = da.ones((5,))
    y = da.cumsum(x)
    no_rechunk = "rechunk" not in str(dict(y.__dask_graph__()))
    assert no_rechunk


@pytest.mark.parametrize("axis", [3, 0, [1, 3]])
@pytest.mark.parametrize("q", [0.75, [0.75], [0.75, 0.4]])
@pytest.mark.parametrize("rechunk", [True, False])
def test_nanquantile(rechunk, q, axis):
    shape = 7, 10, 7, 10
    arr = np.random.randn(*shape)
    indexer = np.random.randint(0, 10, size=shape)
    arr[indexer >= 8] = np.nan
    arr[:, :, :, 1] = 1
    arr[1, :, :, :] = 1

    darr = da.from_array(arr, chunks=(2, 3, 4, (5 if rechunk else -1)))
    assert_eq(da.nanquantile(darr, q, axis=axis), np.nanquantile(arr, q, axis=axis))
    assert_eq(
        da.nanquantile(darr, q, axis=axis, keepdims=True),
        np.nanquantile(arr, q, axis=axis, keepdims=True),
    )
    assert_eq(
        da.nanpercentile(darr, q * 100, axis=axis),
        np.nanpercentile(arr, q * 100, axis=axis),
    )
    assert_eq(
        da.nanpercentile(darr, q * 100, axis=axis, keepdims=True),
        np.nanpercentile(arr, q * 100, axis=axis, keepdims=True),
    )


@pytest.mark.parametrize("axis", [3, [1, 3]])
@pytest.mark.parametrize("q", [0.75, [0.75]])
@pytest.mark.parametrize("rechunk", [True, False])
def test_quantile(rechunk, q, axis):
    shape = 10, 15, 20, 15
    arr = np.random.randn(*shape)
    indexer = np.random.randint(0, 10, size=shape)
    arr[indexer >= 8] = np.nan

    darr = da.from_array(arr, chunks=(2, 3, 4, (5 if rechunk else -1)))
    assert_eq(da.quantile(darr, q, axis=axis), np.quantile(arr, q, axis=axis))
    assert_eq(
        da.quantile(darr, q, axis=axis, keepdims=True),
        np.quantile(arr, q, axis=axis, keepdims=True),
    )
    assert_eq(da.percentile(darr, q, axis=axis), np.percentile(arr, q, axis=axis))
    assert_eq(
        da.percentile(darr, q, axis=axis, keepdims=True),
        np.percentile(arr, q, axis=axis, keepdims=True),
    )


@pytest.mark.parametrize("func", [da.quantile, da.nanquantile, da.nanpercentile])
def test_quantile_func_family_with_axis_none(func):

    # Check that these functions raise a NotImplementedError
    # when axis=None and more than one chunk is present
    # along at least one dimension
    darr = da.ones((3, 3), chunks=(2, 2))
    with pytest.raises(
        NotImplementedError, match="The full algorithm is difficult to do in parallel"
    ):
        func(darr, 0.5, axis=None)

    # Check that the functions behave as expected
    # when axis=None and the array is a single chunk
    darr = da.from_array([-1, 0, 1])
    assert_eq(func(darr, 0.0, axis=None), -1.0)


def test_nanquantile_all_nan():
    shape = 10, 15, 20, 15
    arr = np.random.randn(*shape)
    arr[:] = np.nan
    darr = da.from_array(arr, chunks=(2, 3, 4, -1))
    da.nanquantile(darr, 0.75, axis=-1).compute()
    with pytest.raises(RuntimeWarning):
        assert_eq(
            da.nanquantile(darr, 0.75, axis=-1), np.nanquantile(arr, 0.75, axis=-1)
        )
        assert_eq(da.percentile(darr, 0.75, axis=-1), np.percentile(arr, 0.75, axis=-1))


def test_nanquantile_method():
    shape = 10, 15, 20, 15
    arr = np.random.randn(*shape)
    indexer = np.random.randint(0, 10, size=shape)
    arr[indexer >= 8] = np.nan
    darr = da.from_array(arr, chunks=(2, 3, 4, -1))
    assert_eq(
        da.nanquantile(darr, 0.75, axis=-1, method="weibull"),
        np.nanquantile(arr, 0.75, axis=-1, method="weibull"),
    )
    assert_eq(
        da.nanpercentile(darr, 0.75, axis=-1, method="weibull"),
        np.nanpercentile(arr, 0.75, axis=-1, method="weibull"),
    )


def test_nanquantile_one_dim():
    arr = np.random.randn(10)
    darr = da.from_array(arr, chunks=(2,))
    assert_eq(da.nanquantile(darr, 0.75, axis=-1), np.nanquantile(arr, 0.75, axis=-1))


def test_nanquantile_two_dims():
    arr = np.random.randn(10, 10)
    darr = da.from_array(arr, chunks=(2, -1))
    assert_eq(da.nanquantile(darr, 0.75, axis=-1), np.nanquantile(arr, 0.75, axis=-1))
    assert_eq(
        da.nanpercentile(darr, 0.75, axis=-1), np.nanpercentile(arr, 0.75, axis=-1)
    )
