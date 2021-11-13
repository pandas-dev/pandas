import itertools

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    notna,
)


# create the data only once as we are not setting it
def _create_consistency_data():
    def create_series():
        return [
            Series(dtype=np.float64, name="a"),
            Series([np.nan] * 5),
            Series([1.0] * 5),
            Series(range(5, 0, -1)),
            Series(range(5)),
            Series([np.nan, 1.0, np.nan, 1.0, 1.0]),
            Series([np.nan, 1.0, np.nan, 2.0, 3.0]),
            Series([np.nan, 1.0, np.nan, 3.0, 2.0]),
        ]

    def create_dataframes():
        return [
            DataFrame(columns=["a", "a"]),
            DataFrame(np.arange(15).reshape((5, 3)), columns=["a", "a", 99]),
        ] + [DataFrame(s) for s in create_series()]

    def is_constant(x):
        values = x.values.ravel("K")
        return len(set(values[notna(values)])) == 1

    def no_nans(x):
        return x.notna().all().all()

    return [
        (x, is_constant(x), no_nans(x))
        for x in itertools.chain(create_dataframes(), create_dataframes())
    ]


@pytest.fixture(params=_create_consistency_data())
def consistency_data(request):
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


@pytest.fixture(params=[(1, 0), (5, 1)])
def rolling_consistency_cases(request):
    """window, min_periods"""
    return request.param


@pytest.fixture(params=[0, 2])
def min_periods(request):
    return request.param
