from datetime import (
    datetime,
    timedelta,
)

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
    bdate_range,
)


@pytest.fixture(name="raw", params=[True, False])
def fixture_raw(request):
    """raw keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    name="arithmetic_win_operators",
    params=[
        "sum",
        "mean",
        "median",
        "max",
        "min",
        "var",
        "std",
        "kurt",
        "skew",
        "count",
        "sem",
    ],
)
def fixture_arithmetic_win_operators(request):
    return request.param


@pytest.fixture(name="center", params=[True, False])
def fixture_center(request):
    return request.param


@pytest.fixture(name="min_periods", params=[None, 1])
def fixture_min_periods(request):
    return request.param


@pytest.fixture(name="parallel", params=[True, False])
def fixture_parallel(request):
    """parallel keyword argument for numba.jit"""
    return request.param


# Can parameterize nogil & nopython over True | False, but limiting per
# https://github.com/pandas-dev/pandas/pull/41971#issuecomment-860607472


@pytest.fixture(name="nogil", params=[False])
def fixture_nogil(request):
    """nogil keyword argument for numba.jit"""
    return request.param


@pytest.fixture(name="nopython", params=[True])
def fixture_nopython(request):
    """nopython keyword argument for numba.jit"""
    return request.param


@pytest.fixture(name="adjust", params=[True, False])
def fixture_adjust(request):
    """adjust keyword argument for ewm"""
    return request.param


@pytest.fixture(name="ignore_na", params=[True, False])
def fixture_ignore_na(request):
    """ignore_na keyword argument for ewm"""
    return request.param


@pytest.fixture(name="numeric_only", params=[True, False])
def fixture_numeric_only(request):
    """numeric_only keyword argument"""
    return request.param


@pytest.fixture(
    name="engine",
    params=[pytest.param("numba", marks=td.skip_if_no("numba")), "cython"],
)
def fixture_engine(request):
    """engine keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    name="engine_and_raw",
    params=[
        pytest.param(("numba", True), marks=td.skip_if_no("numba")),
        ("cython", True),
        ("cython", False),
    ],
)
def fixture_engine_and_raw(request):
    """engine and raw keyword arguments for rolling.apply"""
    return request.param


@pytest.fixture(
    name="halflife_with_times",
    params=["1 day", timedelta(days=1), np.timedelta64(1, "D")],
)
def fixture_halflife_with_times(request):
    """Halflife argument for EWM when times is specified."""
    return request.param


@pytest.fixture(name="series")
def fixture_series():
    """Make mocked series as fixture."""
    arr = np.random.randn(100)
    locs = np.arange(20, 40)
    arr[locs] = np.NaN
    series = Series(arr, index=bdate_range(datetime(2009, 1, 1), periods=100))
    return series


@pytest.fixture(name="frame")
def fixture_frame():
    """Make mocked frame as fixture."""
    return DataFrame(
        np.random.randn(100, 10),
        index=bdate_range(datetime(2009, 1, 1), periods=100),
        columns=np.arange(10),
    )


@pytest.fixture(name="step", params=[None, 1, 2, 5, 10])
def fixture_step(request):
    """step keyword argument for rolling window operations."""
    return request.param
