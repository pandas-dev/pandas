import numpy as np
import pytest

from pandas import Series
import pandas.util.testing as tm


class TestSeriesRound:
    def test_round(self, datetime_series):
        datetime_series.index.name = "index_name"
        result = datetime_series.round(2)
        expected = Series(
            np.round(datetime_series.values, 2), index=datetime_series.index, name="ts"
        )
        tm.assert_series_equal(result, expected)
        assert result.name == datetime_series.name

    def test_round_numpy(self):
        # See GH#12600
        ser = Series([1.53, 1.36, 0.06])
        out = np.round(ser, decimals=0)
        expected = Series([2.0, 1.0, 0.0])
        tm.assert_series_equal(out, expected)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.round(ser, decimals=0, out=ser)

    def test_round_numpy_with_nan(self):
        # See GH#14197
        ser = Series([1.53, np.nan, 0.06])
        with tm.assert_produces_warning(None):
            result = ser.round()
        expected = Series([2.0, np.nan, 0.0])
        tm.assert_series_equal(result, expected)

    def test_round_builtin(self):
        ser = Series([1.123, 2.123, 3.123], index=range(3))
        result = round(ser)
        expected_rounded0 = Series([1.0, 2.0, 3.0], index=range(3))
        tm.assert_series_equal(result, expected_rounded0)

        decimals = 2
        expected_rounded = Series([1.12, 2.12, 3.12], index=range(3))
        result = round(ser, decimals)
        tm.assert_series_equal(result, expected_rounded)
