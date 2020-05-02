import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame, Series, notna


@pytest.fixture(params=[True, False])
def raw(request):
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
    params=["sum", "mean", "median", "max", "min", "var", "std", "kurt", "skew"]
)
def arithmetic_win_operators(request):
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


@pytest.fixture(params=[True, False])
def parallel(request):
    """parallel keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True, False])
def nogil(request):
    """nogil keyword argument for numba.jit"""
    return request.param


@pytest.fixture(params=[True, False])
def nopython(request):
    """nopython keyword argument for numba.jit"""
    return request.param


@pytest.fixture(
    params=[pytest.param("numba", marks=td.skip_if_no("numba", "0.46.0")), "cython"]
)
def engine(request):
    """engine keyword argument for rolling.apply"""
    return request.param


@pytest.fixture(
    params=[
        pytest.param(("numba", True), marks=td.skip_if_no("numba", "0.46.0")),
        ("cython", True),
        ("cython", False),
    ]
)
def engine_and_raw(request):
    """engine and raw keyword arguments for rolling.apply"""
    return request.param


# create the data only once as we are not setting it
def _create_consistency_data():
    def create_series():
        return [
            Series(dtype=object),
            Series([np.nan]),
            Series([np.nan, np.nan]),
            Series([3.0]),
            Series([np.nan, 3.0]),
            Series([3.0, np.nan]),
            Series([1.0, 3.0]),
            Series([2.0, 2.0]),
            Series([3.0, 1.0]),
            Series(
                [5.0, 5.0, 5.0, 5.0, np.nan, np.nan, np.nan, 5.0, 5.0, np.nan, np.nan]
            ),
            Series(
                [
                    np.nan,
                    5.0,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                ]
            ),
            Series(
                [
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    5.0,
                    5.0,
                    np.nan,
                    np.nan,
                ]
            ),
            Series(
                [
                    np.nan,
                    3.0,
                    np.nan,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    np.nan,
                    np.nan,
                    7.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    np.nan,
                    5.0,
                    np.nan,
                    2.0,
                    4.0,
                    0.0,
                    9.0,
                    np.nan,
                    np.nan,
                    3.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    2.0,
                    3.0,
                    np.nan,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    np.nan,
                    np.nan,
                    7.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(
                [
                    2.0,
                    5.0,
                    np.nan,
                    2.0,
                    4.0,
                    0.0,
                    9.0,
                    np.nan,
                    np.nan,
                    3.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                ]
            ),
            Series(range(10)),
            Series(range(20, 0, -2)),
        ]

    def create_dataframes():
        return [
            DataFrame(),
            DataFrame(columns=["a"]),
            DataFrame(columns=["a", "a"]),
            DataFrame(columns=["a", "b"]),
            DataFrame(np.arange(10).reshape((5, 2))),
            DataFrame(np.arange(25).reshape((5, 5))),
            DataFrame(np.arange(25).reshape((5, 5)), columns=["a", "b", 99, "d", "d"]),
        ] + [DataFrame(s) for s in create_series()]

    def is_constant(x):
        values = x.values.ravel()
        return len(set(values[notna(values)])) == 1

    def no_nans(x):
        return x.notna().all().all()

    # data is a tuple(object, is_constant, no_nans)
    data = create_series() + create_dataframes()

    return [(x, is_constant(x), no_nans(x)) for x in data]


@pytest.fixture(params=_create_consistency_data())
def consistency_data(request):
    """Create consistency data"""
    return request.param


base_functions_list = [
    (lambda v: Series(v).count(), None, "count"),
    (lambda v: Series(v).max(), None, "max"),
    (lambda v: Series(v).min(), None, "min"),
    (lambda v: Series(v).sum(), None, "sum"),
    (lambda v: Series(v).mean(), None, "mean"),
    (lambda v: Series(v).std(), 1, "std"),
    (lambda v: Series(v).cov(Series(v)), None, "cov"),
    (lambda v: Series(v).corr(Series(v)), None, "corr"),
    (lambda v: Series(v).var(), 1, "var"),
    # restore once GH 8086 is fixed
    # lambda v: Series(v).skew(), 3, 'skew'),
    # (lambda v: Series(v).kurt(), 4, 'kurt'),
    # restore once GH 8084 is fixed
    # lambda v: Series(v).quantile(0.3), None, 'quantile'),
    (lambda v: Series(v).median(), None, "median"),
    (np.nanmax, 1, "max"),
    (np.nanmin, 1, "min"),
    (np.nansum, 1, "sum"),
    (np.nanmean, 1, "mean"),
    (lambda v: np.nanstd(v, ddof=1), 1, "std"),
    (lambda v: np.nanvar(v, ddof=1), 1, "var"),
    (np.nanmedian, 1, "median"),
]

no_nan_functions_list = [
    (np.max, None, "max"),
    (np.min, None, "min"),
    (np.sum, None, "sum"),
    (np.mean, None, "mean"),
    (lambda v: np.std(v, ddof=1), 1, "std"),
    (lambda v: np.var(v, ddof=1), 1, "var"),
    (np.median, None, "median"),
]


@pytest.fixture()
def base_functions():
    """Fixture for base functions."""
    return base_functions_list


@pytest.fixture()
def no_nan_functions():
    """Fixture for no nan functions."""
    return no_nan_functions_list
