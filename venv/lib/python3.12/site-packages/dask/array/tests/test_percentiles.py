from __future__ import annotations

import pytest

pytest.importorskip("numpy")

import numpy as np

import dask.array as da
from dask.array.utils import assert_eq, same_keys

try:
    import crick
except ImportError:
    crick = None


percentile_internal_methods = pytest.mark.parametrize(
    "internal_method",
    [
        pytest.param(
            "tdigest", marks=pytest.mark.skipif(not crick, reason="Requires crick")
        ),
        "dask",
    ],
)


@percentile_internal_methods
def test_percentile(internal_method):
    d = da.ones((16,), chunks=(4,))
    qs = [0, 50, 100]

    assert_eq(
        da.percentile(d, qs, internal_method=internal_method),
        np.array([1, 1, 1], dtype=d.dtype),
    )

    x = np.array([0, 0, 5, 5, 5, 5, 20, 20])
    d = da.from_array(x, chunks=(3,))

    result = da.percentile(d, qs, internal_method=internal_method)
    assert_eq(result, np.array([0, 5, 20], dtype=result.dtype))

    assert same_keys(
        da.percentile(d, qs, internal_method=internal_method),
        da.percentile(d, qs, internal_method=internal_method),
    )
    assert not same_keys(
        da.percentile(d, qs, internal_method=internal_method),
        da.percentile(d, [0, 50], internal_method=internal_method),
    )

    if internal_method != "tdigest":
        x = np.array(["a", "a", "d", "d", "d", "e"])
        d = da.from_array(x, chunks=(3,))
        assert_eq(
            da.percentile(d, [0, 50, 100]), np.array(["a", "d", "e"], dtype=x.dtype)
        )


@pytest.mark.skip
def test_percentile_with_categoricals():
    try:
        import pandas as pd
    except ImportError:
        return
    x0 = pd.Categorical(["Alice", "Bob", "Charlie", "Dennis", "Alice", "Alice"])
    x1 = pd.Categorical(["Alice", "Bob", "Charlie", "Dennis", "Alice", "Alice"])

    dsk = {("x", 0): x0, ("x", 1): x1}

    x = da.Array(dsk, "x", chunks=((6, 6),))

    p = da.percentile(x, [50])
    assert (p.compute().categories == x0.categories).all()
    assert (p.compute().codes == [0]).all()
    assert same_keys(da.percentile(x, [50]), da.percentile(x, [50]))


@percentile_internal_methods
def test_percentiles_with_empty_arrays(internal_method):
    x = da.ones(10, chunks=((5, 0, 5),))
    assert_eq(
        da.percentile(x, [10, 50, 90], internal_method=internal_method),
        np.array([1, 1, 1], dtype=x.dtype),
    )


@percentile_internal_methods
def test_percentiles_with_empty_q(internal_method):
    x = da.ones(10, chunks=((5, 0, 5),))
    assert_eq(
        da.percentile(x, [], internal_method=internal_method),
        np.array([], dtype=x.dtype),
    )


@percentile_internal_methods
@pytest.mark.parametrize("q", [5, 5.0, np.int64(5), np.float64(5)])
def test_percentiles_with_scaler_percentile(internal_method, q):
    # Regression test to ensure da.percentile works with scalar percentiles
    # See #3020
    d = da.ones((16,), chunks=(4,))
    assert_eq(
        da.percentile(d, q, internal_method=internal_method),
        np.array([1], dtype=d.dtype),
    )


@percentile_internal_methods
def test_unknown_chunk_sizes(internal_method):
    x = da.random.default_rng().random(1000, chunks=(100,))
    x._chunks = ((np.nan,) * 10,)

    result = da.percentile(x, 50, internal_method=internal_method).compute()
    assert 0.1 < result < 0.9

    a, b = da.percentile(x, [40, 60], internal_method=internal_method).compute()
    assert 0.1 < a < 0.9
    assert 0.1 < b < 0.9
    assert a < b
