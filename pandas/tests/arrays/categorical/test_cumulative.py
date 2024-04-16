"""
Tests for Ordered Categorical Array cumulative operations.
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestAccumulator:
    @pytest.mark.parametrize(
        "method, order",
        [
            ["cummax", "abc"],
            ["cummin", "cab"],
        ],
    )
    def test_cummax_cummin_on_ordered_categorical(self, method, order):
        # GH#52335
        arr = pd.array(
            list("ababcab"), dtype=pd.CategoricalDtype(list(order), ordered=True)
        )
        result = getattr(arr, method)()
        tm.assert_series_equal(result, pd.Series(list("abbbccc")))

    @pytest.mark.parametrize(
        "method, order",
        [
            ["cummax", "abc"],
            ["cummin", "cab"],
        ],
    )
    def test_cummax_cummin_ordered_categorical_nan(self, method, order):
        # GH#52335
        arr = pd.array(
            ["a", np.nan, "b", "a", "c"],
            dtype=pd.CategoricalDtype(list(order), ordered=True),
        )
        result = getattr(arr, method)(skipna=True)
        tm.assert_series_equal(result, pd.array(["a", "NaN", "b", "b", "c"]))

        result = getattr(arr, method)(skipna=False)
        tm.assert_series_equal(result, pd.array(["a", "NaN", "NaN", "NaN", "NaN"]))
