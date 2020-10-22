from datetime import date

import numpy as np
import pytest

from pandas import MultiIndex, NaT, Series, Timestamp, date_range, period_range
import pandas.testing as tm


class TestSetitemDT64Values:
    def test_setitem_none_nan(self):
        series = Series(date_range("1/1/2000", periods=10))
        series[3] = None
        assert series[3] is NaT

        series[3:5] = None
        assert series[4] is NaT

        series[5] = np.nan
        assert series[5] is NaT

        series[5:7] = np.nan
        assert series[6] is NaT

    def test_setitem_multiindex_empty_slice(self):
        # https://github.com/pandas-dev/pandas/issues/35878
        idx = MultiIndex.from_tuples([("a", 1), ("b", 2)])
        result = Series([1, 2], index=idx)
        expected = result.copy()
        result.loc[[]] = 0
        tm.assert_series_equal(result, expected)

    def test_setitem_with_string_index(self):
        # GH#23451
        ser = Series([1, 2, 3], index=["Date", "b", "other"])
        ser["Date"] = date.today()
        assert ser.Date == date.today()
        assert ser["Date"] == date.today()

    def test_setitem_with_different_tz_casts_to_object(self):
        # GH#24024
        ser = Series(date_range("2000", periods=2, tz="US/Central"))
        ser[0] = Timestamp("2000", tz="US/Eastern")
        expected = Series(
            [
                Timestamp("2000-01-01 00:00:00-05:00", tz="US/Eastern"),
                Timestamp("2000-01-02 00:00:00-06:00", tz="US/Central"),
            ],
            dtype=object,
        )
        tm.assert_series_equal(ser, expected)


class TestSetitemPeriodDtype:
    @pytest.mark.parametrize("na_val", [None, np.nan])
    def test_setitem_na_period_dtype_casts_to_nat(self, na_val):
        ser = Series(period_range("2000-01-01", periods=10, freq="D"))

        ser[3] = na_val
        assert ser[3] is NaT

        ser[3:5] = na_val
        assert ser[4] is NaT
