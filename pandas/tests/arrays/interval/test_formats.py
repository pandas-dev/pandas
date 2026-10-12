import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import IntervalArray


def test_repr():
    # GH#25022
    arr = IntervalArray.from_tuples([(0, 1), (1, 2)])
    result = repr(arr)
    expected = (
        "<IntervalArray>\n[(0, 1], (1, 2]]\nLength: 2, dtype: interval[int64, right]"
    )
    assert result == expected


def test_series_repr_float_precision():
    # GH#25920
    arr = IntervalArray.from_arrays([0.55555555, np.nan], [1.55555555, np.nan])
    ser = pd.Series(arr)
    with pd.option_context("display.precision", 3):
        result = repr(ser)
    expected = "0    (0.556, 1.556]\n1               NaN\ndtype: interval"
    assert result == expected


@pytest.mark.parametrize(
    "breaks, expected",
    [
        ([1e-9, 2e-9], "0    (1.000000e-09, 2.000000e-09]\ndtype: interval"),
        ([1e20, 2e20], "0    (1.000000e+20, 2.000000e+20]\ndtype: interval"),
        ([0.5, 1.25, 2.0], "0    (0.50, 1.25]\n1    (1.25, 2.00]\ndtype: interval"),
    ],
)
def test_series_repr_float_matches_float64(breaks, expected):
    # GH#25920 endpoints are formatted like a float64 Series
    ser = pd.Series(IntervalArray.from_breaks(breaks))
    assert repr(ser) == expected


def test_series_repr_float_format():
    # GH#25920
    ser = pd.Series(IntervalArray.from_breaks([0.5, 1.5]))
    with pd.option_context("display.float_format", "{:.2e}".format):
        result = repr(ser)
    assert result == "0    (5.00e-01, 1.50e+00]\ndtype: interval"


@pytest.mark.parametrize(
    "ser, expected",
    [
        (
            pd.Series([1, 2], index=pd.IntervalIndex.from_breaks([0.5, 1.25, 2.0])),
            "(0.5, 1.25]    1\n(1.25, 2.0]    2\ndtype: int64",
        ),
        (
            pd.Series(pd.Categorical(pd.IntervalIndex.from_breaks([0.5, 1.25, 2.0]))),
            "0    (0.5, 1.25]\n1    (1.25, 2.0]\ndtype: category\n"
            "Categories (2, interval[float64, right]): [(0.5, 1.25], (1.25, 2.0]]",
        ),
    ],
)
def test_series_repr_float_interval_labels_and_categories_unchanged(ser, expected):
    # GH#25920 only interval-dtype values are formatted like float64
    assert repr(ser) == expected
