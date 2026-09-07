import numpy as np

import pandas as pd
import pandas._testing as tm


class TestCategoricalComparisons:
    def test_categorical_nan_equality(self):
        cat = pd.Series(pd.Categorical(["a", "b", "c", np.nan]))
        expected = pd.Series([True, True, True, False])
        result = cat == cat
        tm.assert_series_equal(result, expected)

    def test_categorical_tuple_equality(self):
        # GH 18050
        ser = pd.Series([(0, 0), (0, 1), (0, 0), (1, 0), (1, 1)])
        expected = pd.Series([True, False, True, False, False])
        result = ser == (0, 0)
        tm.assert_series_equal(result, expected)

        result = ser.astype("category") == (0, 0)
        tm.assert_series_equal(result, expected)
