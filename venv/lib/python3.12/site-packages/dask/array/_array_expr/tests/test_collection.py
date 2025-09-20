from __future__ import annotations

import operator

import numpy as np
import pytest

import dask.array as da
from dask import is_dask_collection
from dask.array import Array, assert_eq
from dask.array._array_expr._rechunk import Rechunk


@pytest.fixture()
def arr():
    return da.random.random((10, 10), chunks=(5, 6))


@pytest.mark.array_expr
@pytest.mark.parametrize(
    "op",
    [
        "__add__",
        "__sub__",
        "__mul__",
        "__truediv__",
        "__floordiv__",
        "__pow__",
        "__radd__",
        "__rsub__",
        "__rmul__",
        "__rtruediv__",
        "__rfloordiv__",
        "__rpow__",
    ],
)
def test_arithmetic_ops(arr, op):
    result = getattr(arr, op)(2)
    expected = getattr(arr.compute(), op)(2)
    assert_eq(result, expected)


@pytest.mark.array_expr
def test_rechunk(arr):
    result = arr.rechunk((7, 3))
    expected = arr.compute()
    assert_eq(result, expected)


@pytest.mark.array_expr
def test_blockwise():
    x = da.random.random((10, 10), chunks=(5, 5))
    z = da.blockwise(operator.add, "ij", x, "ij", 100, None, dtype=x.dtype)
    assert_eq(z, x.compute() + 100)

    x = da.random.random((10, 10), chunks=(5, 5))
    z = da.blockwise(operator.add, "ij", x, "ij", x, "ij", dtype=x.dtype)
    expr = z.expr.optimize()
    assert len(list(expr.find_operations(Rechunk))) == 0
    assert_eq(z, x.compute() * 2)

    # align
    x = da.random.random((10, 10), chunks=(5, 5))
    y = da.random.random((10, 10), chunks=(7, 3))
    z = da.blockwise(operator.add, "ij", x, "ij", y, "ij", dtype=x.dtype)
    expr = z.expr.optimize()
    assert len(list(expr.find_operations(Rechunk))) > 0
    assert_eq(z, x.compute() + y.compute())


@pytest.mark.parametrize("func", ["min", "max", "sum", "prod", "mean", "any", "all"])
def test_reductions(arr, func):
    # var and std need __array_function__
    result = getattr(arr, func)(axis=0)
    expected = getattr(arr.compute(), func)(axis=0)
    assert_eq(result, expected)


@pytest.mark.parametrize(
    "func",
    [
        "sum",
        "mean",
        "any",
        "all",
        "max",
        "min",
        "nanmin",
        "nanmax",
        "nanmean",
        "nansum",
        "nanprod",
    ],
)
def test_reductions_toplevel(arr, func):
    # var and std need __array_function__
    result = getattr(da, func)(arr, axis=0)
    expected = getattr(np, func)(arr.compute(), axis=0)
    assert_eq(result, expected)


def test_from_array():
    x = np.random.random((10, 10))
    d = da.from_array(x, chunks=(5, 5))
    assert_eq(d, x)
    assert d.chunks == ((5, 5), (5, 5))


@pytest.mark.array_expr
def test_is_dask_collection_doesnt_materialize():
    class ArrayTest(Array):
        def __dask_graph__(self):
            raise NotImplementedError

    arr = ArrayTest(da.random.random((10, 10), chunks=(5, 5)).expr)
    assert is_dask_collection(arr)
    with pytest.raises(NotImplementedError):
        arr.__dask_graph__()


def test_astype():
    x = da.random.randint(1, 100, (10, 10), chunks=(5, 5))
    result = x.astype(np.float64)
    expected = x.compute().astype(np.float64)
    assert_eq(result, expected)


def test_stack_promote_type():
    i = np.arange(10, dtype="i4")
    f = np.arange(10, dtype="f4")
    di = da.from_array(i, chunks=5)
    df = da.from_array(f, chunks=5)
    res = da.stack([di, df])
    assert_eq(res, np.stack([i, f]))
