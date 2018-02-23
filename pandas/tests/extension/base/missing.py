import numpy as np

import pandas as pd
import pandas.util.testing as tm


class BaseMissingTests(object):
    def test_isna(self, data_missing):
        if data_missing._can_hold_na:
            expected = np.array([True, False])
        else:
            expected = np.array([False, False])

        result = pd.isna(data_missing)
        tm.assert_numpy_array_equal(result, expected)

        result = pd.Series(data_missing).isna()
        expected = pd.Series(expected)
        tm.assert_series_equal(result, expected)

    def test_dropna_series(self, data_missing):
        ser = pd.Series(data_missing)
        result = ser.dropna()
        expected = ser.iloc[[1]]
        tm.assert_series_equal(result, expected)

    def test_dropna_frame(self, data_missing):
        df = pd.DataFrame({"A": data_missing})

        # defaults
        result = df.dropna()
        expected = df.iloc[[1]]
        tm.assert_frame_equal(result, expected)

        # axis = 1
        result = df.dropna(axis='columns')
        expected = pd.DataFrame(index=[0, 1])
        tm.assert_frame_equal(result, expected)

        # multiple
        df = pd.DataFrame({"A": data_missing,
                           "B": [1, np.nan]})
        result = df.dropna()
        expected = df.iloc[:0]
        tm.assert_frame_equal(result, expected)
