"""
Tests for Ordered Categorical Array cumulative operations.
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestAccumulator:
    @pytest.mark.parametrize(
        "method, input, output",
        [
            ["cummax", [1, 2, 1, 2, 3, 3, 2, 1], [1, 2, 2, 2, 3, 3, 3, 3]],
            ["cummin", [3, 2, 3, 2, 1, 1, 2, 3], [3, 2, 2, 2, 1, 1, 1, 1]],
        ],
    )
    def test_cummax_cummin_on_ordered_categorical(self, method, input, output):
        # GH#52335
        result = pd.Categorical(input, ordered=True)._accumulate(method)
        tm.assert_extension_array_equal(result, pd.Categorical(output, ordered=True))

    @pytest.mark.parametrize(
        "method, skip, input, output",
        [
            ["cummax", True, [1, np.nan, 2, 1, 3], [1, np.nan, 2, 2, 3]],
            [
                "cummax",
                False,
                [1, np.nan, 2, 1, 3],
                [1, np.nan, np.nan, np.nan, np.nan],
            ],
            ["cummin", True, [3, np.nan, 2, 3, 1], [3, np.nan, 2, 2, 1]],
            [
                "cummin",
                False,
                [3, np.nan, 2, 3, 1],
                [3, np.nan, np.nan, np.nan, np.nan],
            ],
        ],
    )
    def test_cummax_cummin_ordered_categorical_nan(self, skip, method, input, output):
        # GH#52335
        result = pd.Categorical(input, ordered=True)._accumulate(method, skipna=skip)
        tm.assert_extension_array_equal(
            result, pd.Categorical(output, categories=[1, 2, 3], ordered=True)
        )
