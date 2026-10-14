import itertools

import numpy as np
import pytest

import pandas as pd


def create_series():
    return [
        pd.Series(dtype=np.float64, name="a"),
        pd.Series([np.nan] * 5),
        pd.Series([1.0] * 5),
        pd.Series(range(5, 0, -1)),
        pd.Series(range(5)),
        pd.Series([np.nan, 1.0, np.nan, 1.0, 1.0]),
        pd.Series([np.nan, 1.0, np.nan, 2.0, 3.0]),
        pd.Series([np.nan, 1.0, np.nan, 3.0, 2.0]),
    ]


def create_dataframes():
    return [
        pd.DataFrame(columns=["a", "a"]),
        pd.DataFrame(np.arange(15).reshape((5, 3)), columns=["a", "a", 99]),
    ] + [pd.DataFrame(s) for s in create_series()]


def is_constant(x):
    values = x.values.ravel("K")
    return len(set(values[pd.notna(values)])) == 1


@pytest.fixture(
    params=(
        obj
        for obj in itertools.chain(create_series(), create_dataframes())
        if is_constant(obj)
    ),
)
def consistent_data(request):
    return request.param


@pytest.fixture(params=create_series())
def series_data(request):
    return request.param


@pytest.fixture(params=itertools.chain(create_series(), create_dataframes()))
def all_data(request):
    """
    Test:
        - Empty Series / DataFrame
        - All NaN
        - All consistent value
        - Monotonically decreasing
        - Monotonically increasing
        - Monotonically consistent with NaNs
        - Monotonically increasing with NaNs
        - Monotonically decreasing with NaNs
    """
    return request.param


@pytest.fixture(params=[0, 2])
def min_periods(request):
    return request.param
