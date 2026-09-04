"""
We also test Series.notna in this file.
"""

import numpy as np

import pandas as pd
import pandas._testing as tm


class TestIsna:
    def test_isna_period_dtype(self):
        # GH#13737
        ser = pd.Series([pd.Period("2011-01", freq="M"), pd.Period("NaT", freq="M")])

        expected = pd.Series([False, True])

        result = ser.isna()
        tm.assert_series_equal(result, expected)

        result = ser.notna()
        tm.assert_series_equal(result, ~expected)

    def test_isna(self):
        ser = pd.Series([0, 5.4, 3, np.nan, -0.001])
        expected = pd.Series([False, False, False, True, False])
        tm.assert_series_equal(ser.isna(), expected)
        tm.assert_series_equal(ser.notna(), ~expected)

        ser = pd.Series(["hi", "", np.nan])
        expected = pd.Series([False, False, True])
        tm.assert_series_equal(ser.isna(), expected)
        tm.assert_series_equal(ser.notna(), ~expected)
