from __future__ import annotations

import warnings

import pytest
from packaging.version import Version

scipy = pytest.importorskip("scipy")
import numpy as np

import dask.array as da
import dask.array.stats
from dask.array.utils import allclose, assert_eq
from dask.delayed import Delayed


@pytest.mark.parametrize(
    "kind, kwargs", [("skew", {}), ("kurtosis", {}), ("kurtosis", {"fisher": False})]
)
@pytest.mark.parametrize("single_dim", [True, False])
def test_measures(kind, kwargs, single_dim):
    np.random.seed(seed=1337)
    if single_dim:
        x = np.random.random(size=(30,))
    else:
        x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)
    dfunc = getattr(dask.array.stats, kind)
    sfunc = getattr(scipy.stats, kind)

    expected = sfunc(x, **kwargs)
    result = dfunc(y, **kwargs)
    if np.isscalar(expected):
        # make it an array to account for possible numeric errors
        expected = np.array(expected)
    assert_eq(result, expected)
    assert isinstance(result, da.Array)


def test_bias_raises():
    x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)

    with pytest.raises(NotImplementedError):
        dask.array.stats.skew(y, bias=False)

    with pytest.raises(NotImplementedError):
        dask.array.stats.kurtosis(y, bias=False)


@pytest.mark.parametrize(
    "kind", ["chisquare", "power_divergence", "normaltest", "skewtest", "kurtosistest"]
)
def test_one(kind):
    a = np.random.random(size=30)
    a_ = da.from_array(a, 3)

    dask_test = getattr(dask.array.stats, kind)
    scipy_test = getattr(scipy.stats, kind)

    result = dask_test(a_)
    expected = scipy_test(a)

    assert isinstance(result, Delayed)
    assert allclose(result.compute(), expected)


@pytest.mark.parametrize(
    "kind, kwargs",
    [
        ("ttest_ind", {}),
        ("ttest_ind", {"equal_var": False}),
        pytest.param(
            "ttest_1samp",
            {},
            marks=pytest.mark.xfail(
                # NOTE: using nested `Version` calls here to handle night scipy releases
                Version(Version(scipy.__version__).base_version) >= Version("1.10.0"),
                reason="https://github.com/dask/dask/issues/9499",
            ),
        ),
        ("ttest_rel", {}),
        ("chisquare", {}),
        ("power_divergence", {}),
        ("power_divergence", {"lambda_": 0}),
        ("power_divergence", {"lambda_": -1}),
        ("power_divergence", {"lambda_": "neyman"}),
    ],
)
def test_two(kind, kwargs):
    # The sums of observed and expected frequencies must match
    a = np.random.random(size=30)
    b = a[::-1]

    a_ = da.from_array(a, 3)
    b_ = da.from_array(b, 3)

    dask_test = getattr(dask.array.stats, kind)
    scipy_test = getattr(scipy.stats, kind)

    with warnings.catch_warnings():  # maybe overflow warning (power_divergence)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        result = dask_test(a_, b_, **kwargs)
        expected = scipy_test(a, b, **kwargs)

    assert isinstance(result, Delayed)
    assert allclose(result.compute(), expected)
    # fails occasionally. shouldn't this be exact?
    # assert dask.compute(*result) == expected


@pytest.mark.parametrize("k", range(5))
def test_moments(k):
    x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)

    expected = scipy.stats.moment(x, k)
    result = dask.array.stats.moment(y, k)
    assert_eq(result, expected)


def test_anova():
    np_args = [i * np.random.random(size=(30,)) for i in range(4)]
    da_args = [da.from_array(x, chunks=10) for x in np_args]

    result = dask.array.stats.f_oneway(*da_args)
    expected = scipy.stats.f_oneway(*np_args)

    assert allclose(result.compute(), expected)


@pytest.mark.parametrize(
    "func, nargs",
    [
        (dask.array.stats.ttest_1samp, 2),
        (dask.array.stats.ttest_rel, 2),
        (dask.array.stats.skewtest, 1),
        (dask.array.stats.kurtosis, 1),
        (dask.array.stats.kurtosistest, 1),
        (dask.array.stats.normaltest, 1),
        (dask.array.stats.moment, 1),
    ],
)
@pytest.mark.parametrize("nan_policy", ["omit", "raise"])
def test_nan_raises(func, nargs, nan_policy):
    with pytest.raises(NotImplementedError):
        func(*(None,) * nargs, nan_policy=nan_policy)


def test_power_divergence_invalid():
    a = np.random.random(size=30)
    a_ = da.from_array(a, 3)

    with pytest.raises(ValueError):
        dask.array.stats.power_divergence(a_, lambda_="wrong")


def test_skew_raises():
    a = da.ones((7,), chunks=(7,))
    with pytest.raises(ValueError, match="7 samples"):
        dask.array.stats.skewtest(a)


def test_skew_single_return_type():
    """This function tests the return type for the skew method for a 1d array."""
    numpy_array = np.random.random(size=(30,))
    dask_array = da.from_array(numpy_array, 3)
    result = dask.array.stats.skew(dask_array).compute()
    assert isinstance(result, np.float64)


def test_kurtosis_single_return_type():
    """This function tests the return type for the kurtosis method for a 1d array."""
    numpy_array = np.random.random(size=(30,))
    dask_array = da.from_array(numpy_array, 3)
    result = dask.array.stats.kurtosis(dask_array).compute()
    result_non_fisher = dask.array.stats.kurtosis(dask_array, fisher=False).compute()
    assert isinstance(result, np.float64)
    assert isinstance(result_non_fisher, np.float64)
