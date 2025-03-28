from __future__ import annotations

import pytest

import dask.array._array_expr as da
from dask.array import assert_eq


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
