from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.utils import assert_eq, same_keys

cupy = pytest.importorskip("cupy")


def test_percentile():
    d = da.from_array(cupy.ones((16,)), chunks=(4,))
    qs = np.array([0, 50, 100])

    result = da.percentile(d, qs, method="midpoint")
    assert_eq(result, np.array([1, 1, 1], dtype=d.dtype), check_type=False)

    x = cupy.array([0, 0, 5, 5, 5, 5, 20, 20])
    d = da.from_array(x, chunks=(3,))

    result = da.percentile(d, qs, method="midpoint")
    assert_eq(result, np.array([0, 5, 20], dtype=result.dtype), check_type=False)

    assert not same_keys(
        da.percentile(d, qs, "midpoint"),
        da.percentile(d, [0, 50], "midpoint"),
    )


def test_percentile_tokenize():
    d = da.from_array(cupy.ones((16,)), chunks=(4,))
    qs = np.array([0, 50, 100])
    assert same_keys(da.percentile(d, qs), da.percentile(d, qs))


def test_percentiles_with_empty_arrays():
    x = da.from_array(cupy.ones(10), chunks=((5, 0, 5),))
    result = da.percentile(x, [10, 50, 90], method="midpoint")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, np.array([1, 1, 1], dtype=x.dtype), check_type=False)


def test_percentiles_with_empty_q():
    x = da.from_array(cupy.ones(10), chunks=((5, 0, 5),))
    result = da.percentile(x, [], method="midpoint")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, np.array([], dtype=x.dtype), check_type=False)


@pytest.mark.parametrize("q", [5, 5.0, np.int64(5), np.float64(5)])
def test_percentiles_with_scaler_percentile(q):
    # Regression test to ensure da.percentile works with scalar percentiles
    # See #3020
    d = da.from_array(cupy.ones((16,)), chunks=(4,))
    result = da.percentile(d, q, method="midpoint")
    assert type(result._meta) == cupy.ndarray
    assert_eq(result, result)  # Check that _meta and computed arrays match types
    assert_eq(result, np.array([1], dtype=d.dtype), check_type=False)


def test_percentiles_with_unknown_chunk_sizes():
    rng = da.random.default_rng(cupy.random.default_rng())
    x = rng.random(1000, chunks=(100,))
    x._chunks = ((np.nan,) * 10,)

    result = da.percentile(x, 50, method="midpoint").compute()
    assert type(result) == cupy.ndarray
    assert 0.1 < result < 0.9

    a, b = da.percentile(x, [40, 60], method="midpoint").compute()
    assert type(a) == cupy.ndarray
    assert type(b) == cupy.ndarray
    assert 0.1 < a < 0.9
    assert 0.1 < b < 0.9
    assert a < b
