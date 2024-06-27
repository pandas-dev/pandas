from __future__ import annotations

import warnings

import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import dask
import dask.array as da
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")


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
    x = cupy.random.default_rng().random((10, 10, 10))
    a = da.from_array(x, chunks=(3, 4, 5))

    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))
    assert_eq(dfunc(a, 1), func(x, 1))
    assert_eq(dfunc(a, 2), func(x, 2))
    with dask.config.set(split_every=2):
        assert_eq(dfunc(a), func(x))
        assert_eq(dfunc(a, 0), func(x, 0))
        assert_eq(dfunc(a, 1), func(x, 1))
        assert_eq(dfunc(a, 2), func(x, 2))

    pytest.raises(ValueError, lambda: dfunc(a, 3))
    pytest.raises(TypeError, lambda: dfunc(a, (0, 1)))

    x2 = cupy.arange(10)
    a2 = da.from_array(x2, chunks=3)
    assert_eq(dfunc(a2), func(x2))
    assert_eq(dfunc(a2, 0), func(x2, 0))
    assert_eq(dfunc(a2, 0, split_every=2), func(x2, 0))


@pytest.mark.parametrize(
    ["dfunc", "func"], [(da.nanargmin, np.nanargmin), (da.nanargmax, np.nanargmax)]
)
def test_nanarg_reductions(dfunc, func):
    x = cupy.random.default_rng().random((10, 10, 10))
    x[5] = cupy.nan
    a = da.from_array(x, chunks=(3, 4, 5))
    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # All-NaN slice encountered
        with pytest.raises(ValueError):
            dfunc(a, 1).compute()

        with pytest.raises(ValueError):
            dfunc(a, 2).compute()

        x[:] = cupy.nan
        a = da.from_array(x, chunks=(3, 4, 5))
        with pytest.raises(ValueError):
            dfunc(a).compute()


@pytest.mark.parametrize("func", [np.cumsum, np.cumprod])
def test_cumreduction_with_cupy(func):
    a = cupy.ones((10, 10))
    b = da.from_array(a, chunks=(4, 4))
    result = func(b, axis=0)
    assert_eq(result, func(a, axis=0))
