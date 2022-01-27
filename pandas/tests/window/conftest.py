from datetime import (
    datetime,
    timedelta,
)
from functools import partial

import numpy as np
import pytest

import pandas._libs.window.aggregations as window_aggregations
import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
    bdate_range,
    to_datetime,
)


@pytest.fixture(params=[True, False])
def raw(request):
    """raw keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    params=[
        "triang",
        "blackman",
        "hamming",
        "bartlett",
        "bohman",
        "blackmanharris",
        "nuttall",
        "barthann",
    ]
)
def win_types(request):
    return request.param


@pytest.fixture(params=["kaiser", "gaussian", "general_gaussian", "exponential"])
def win_types_special(request):
    return request.param


@pytest.fixture(
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
    ]
)
def arithmetic_win_operators(request):
    return request.param


@pytest.fixture(
    params=[
        ["sum", {}],
        ["mean", {}],
        ["median", {}],
        ["max", {}],
        ["min", {}],
        ["var", {}],
        ["var", {"ddof": 0}],
        ["std", {}],
        ["std", {"ddof": 0}],
    ]
)
def arithmetic_numba_supported_operators(request):
    return request.param


@pytest.fixture(params=["right", "left", "both", "neither"])
def closed(request):
    return request.param


@pytest.fixture(params=[True, False])
def center(request):
    return request.param


@pytest.fixture(params=[None, 1])
def min_periods(request):
    return request.param


@pytest.fixture(params=["single", "table"])
def method(request):
    """method keyword in rolling/expanding/ewm constructor"""
    return request.param


@pytest.fixture(params=[True, False])
def parallel(request):
    """parallel keyword argument for numba.jit"""
    return request.param


# Can parameterize nogil & nopython over True | False, but limiting per
# https://github.com/pandas-dev/pandas/pull/41971#issuecomment-860607472


@pytest.fixture(params=[False])
def nogil(request):
    """nogil keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True])
def nopython(request):
    """nopython keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True, False])
def adjust(request):
    """adjust keyword argument for ewm"""
    return request.param


@pytest.fixture(params=[True, False])
def ignore_na(request):
    """ignore_na keyword argument for ewm"""
    return request.param


@pytest.fixture(params=[pytest.param("numba", marks=td.skip_if_no("numba")), "cython"])
def engine(request):
    """engine keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    params=[
        pytest.param(("numba", True), marks=td.skip_if_no("numba")),
        ("cython", True),
        ("cython", False),
    ]
)
def engine_and_raw(request):
    """engine and raw keyword arguments for rolling.apply"""
    return request.param


@pytest.fixture
def times_frame():
    """Frame for testing times argument in EWM groupby."""
    return DataFrame(
        {
            "A": ["a", "b", "c", "a", "b", "c", "a", "b", "c", "a"],
            "B": [0, 0, 0, 1, 1, 1, 2, 2, 2, 3],
            "C": to_datetime(
                [
                    "2020-01-01",
                    "2020-01-01",
                    "2020-01-01",
                    "2020-01-02",
                    "2020-01-10",
                    "2020-01-22",
                    "2020-01-03",
                    "2020-01-23",
                    "2020-01-23",
                    "2020-01-04",
                ]
            ),
        }
    )


@pytest.fixture(params=["1 day", timedelta(days=1)])
def halflife_with_times(request):
    """Halflife argument for EWM when times is specified."""
    return request.param


@pytest.fixture(
    params=[
        "object",
        "category",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float16",
        "float32",
        "float64",
        "m8[ns]",
        "M8[ns]",
        "datetime64[ns, UTC]",
    ]
)
def dtypes(request):
    """Dtypes for window tests"""
    return request.param


@pytest.fixture(
    params=[
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1, 0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1, 1]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=["C", "C"]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1.0, 0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[0.0, 1]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=["C", 1]),
        DataFrame([[2.0, 4.0], [1.0, 2.0], [5.0, 2.0], [8.0, 1.0]], columns=[1, 0.0]),
        DataFrame([[2, 4.0], [1, 2.0], [5, 2.0], [8, 1.0]], columns=[0, 1.0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1.0]], columns=[1.0, "X"]),
    ]
)
def pairwise_frames(request):
    """Pairwise frames test_pairwise"""
    return request.param


@pytest.fixture
def pairwise_target_frame():
    """Pairwise target frame for test_pairwise"""
    return DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[0, 1])


@pytest.fixture
def pairwise_other_frame():
    """Pairwise other frame for test_pairwise"""
    return DataFrame(
        [[None, 1, 1], [None, 1, 2], [None, 3, 2], [None, 8, 1]],
        columns=["Y", "Z", "X"],
    )


@pytest.fixture
def series():
    """Make mocked series as fixture."""
    arr = np.random.randn(100)
    locs = np.arange(20, 40)
    arr[locs] = np.NaN
    series = Series(arr, index=bdate_range(datetime(2009, 1, 1), periods=100))
    return series


@pytest.fixture
def frame():
    """Make mocked frame as fixture."""
    return DataFrame(
        np.random.randn(100, 10),
        index=bdate_range(datetime(2009, 1, 1), periods=100),
        columns=np.arange(10),
    )


def _named_func(name_and_func):
    name, func = name_and_func
    if not hasattr(func, "func"):
        func = partial(func)
    func.__name__ = name
    return func


@pytest.fixture(
    params=[
        _named_func(x)
        for x in [
            ("roll_sum", window_aggregations.roll_sum),
            ("roll_mean", window_aggregations.roll_mean),
        ]
        + [
            (f"roll_var({ddof})", partial(window_aggregations.roll_var, ddof=ddof))
            for ddof in [0, 1]
        ]
        + [
            ("roll_skew", window_aggregations.roll_skew),
            ("roll_kurt", window_aggregations.roll_kurt),
            ("roll_mediac_c", window_aggregations.roll_median_c),
            ("roll_max", window_aggregations.roll_max),
            ("roll_min", window_aggregations.roll_min),
        ]
        + [
            (
                f"roll_quntile({quantile},{interpolation})",
                partial(
                    window_aggregations.roll_quantile,
                    quantile=quantile,
                    interpolation=interpolation,
                ),
            )
            for quantile in [0.25, 0.5, 0.75]
            for interpolation in window_aggregations.interpolation_types
        ]
        + [
            (
                f"roll_rank({percentile},{method},{ascending})",
                partial(
                    window_aggregations.roll_rank,
                    percentile=percentile,
                    method=method,
                    ascending=ascending,
                ),
            )
            for percentile in [True, False]
            for method in window_aggregations.rolling_rank_tiebreakers.keys()
            for ascending in [True, False]
        ]
    ]
)
def rolling_aggregation(request):
    """Make a named rolling aggregation function as fixture."""
    return request.param
