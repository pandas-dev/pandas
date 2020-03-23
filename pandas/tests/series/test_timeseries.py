import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, DatetimeIndex, Series, date_range, timedelta_range
import pandas._testing as tm


class TestTimeSeries:
    def test_mpl_compat_hack(self, datetime_series):

        # This is currently failing because the test was relying on
        # the DeprecationWarning coming through Index.__getitem__.
        # We want to implement a warning specifically for Series.__getitem__
        # at which point this will become a Deprecation/FutureWarning
        with tm.assert_produces_warning(None):
            # GH#30588 multi-dimensional indexing deprecated
            result = datetime_series[:, np.newaxis]
        expected = datetime_series.values[:, np.newaxis]
        tm.assert_almost_equal(result, expected)

    def test_timeseries_coercion(self):
        idx = tm.makeDateIndex(10000)
        ser = Series(np.random.randn(len(idx)), idx.astype(object))
        assert ser.index.is_all_dates
        assert isinstance(ser.index, DatetimeIndex)

    def test_contiguous_boolean_preserve_freq(self):
        rng = date_range("1/1/2000", "3/1/2000", freq="B")

        mask = np.zeros(len(rng), dtype=bool)
        mask[10:20] = True

        masked = rng[mask]
        expected = rng[10:20]
        assert expected.freq == rng.freq
        tm.assert_index_equal(masked, expected)

        mask[22] = True
        masked = rng[mask]
        assert masked.freq is None

    def test_promote_datetime_date(self):
        rng = date_range("1/1/2000", periods=20)
        ts = Series(np.random.randn(20), index=rng)

        ts_slice = ts[5:]
        ts2 = ts_slice.copy()
        ts2.index = [x.date() for x in ts2.index]

        result = ts + ts2
        result2 = ts2 + ts
        expected = ts + ts[5:]
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        # test asfreq
        result = ts2.asfreq("4H", method="ffill")
        expected = ts[5:].asfreq("4H", method="ffill")
        tm.assert_series_equal(result, expected)

        result = rng.get_indexer(ts2.index)
        expected = rng.get_indexer(ts_slice.index)
        tm.assert_numpy_array_equal(result, expected)

    def test_groupby_count_dateparseerror(self):
        dr = date_range(start="1/1/2012", freq="5min", periods=10)

        # BAD Example, datetimes first
        s = Series(np.arange(10), index=[dr, np.arange(10)])
        grouped = s.groupby(lambda x: x[1] % 2 == 0)
        result = grouped.count()

        s = Series(np.arange(10), index=[np.arange(10), dr])
        grouped = s.groupby(lambda x: x[0] % 2 == 0)
        expected = grouped.count()

        tm.assert_series_equal(result, expected)

    def test_series_map_box_timedelta(self):
        # GH 11349
        s = Series(timedelta_range("1 day 1 s", periods=5, freq="h"))

        def f(x):
            return x.total_seconds()

        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_asfreq_resample_set_correct_freq(self):
        # GH5613
        # we test if .asfreq() and .resample() set the correct value for .freq
        df = pd.DataFrame(
            {"date": ["2012-01-01", "2012-01-02", "2012-01-03"], "col": [1, 2, 3]}
        )
        df = df.set_index(pd.to_datetime(df.date))

        # testing the settings before calling .asfreq() and .resample()
        assert df.index.freq is None
        assert df.index.inferred_freq == "D"

        # does .asfreq() set .freq correctly?
        assert df.asfreq("D").index.freq == "D"

        # does .resample() set .freq correctly?
        assert df.resample("D").asfreq().index.freq == "D"

    def test_view_tz(self):
        # GH#24024
        ser = pd.Series(pd.date_range("2000", periods=4, tz="US/Central"))
        result = ser.view("i8")
        expected = pd.Series(
            [
                946706400000000000,
                946792800000000000,
                946879200000000000,
                946965600000000000,
            ]
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "US/Central"])
    def test_asarray_object_dt64(self, tz):
        ser = pd.Series(pd.date_range("2000", periods=2, tz=tz))

        with tm.assert_produces_warning(None):
            # Future behavior (for tzaware case) with no warning
            result = np.asarray(ser, dtype=object)

        expected = np.array(
            [pd.Timestamp("2000-01-01", tz=tz), pd.Timestamp("2000-01-02", tz=tz)]
        )
        tm.assert_numpy_array_equal(result, expected)

    def test_asarray_tz_naive(self):
        # This shouldn't produce a warning.
        ser = pd.Series(pd.date_range("2000", periods=2))
        expected = np.array(["2000-01-01", "2000-01-02"], dtype="M8[ns]")
        result = np.asarray(ser)

        tm.assert_numpy_array_equal(result, expected)

    def test_asarray_tz_aware(self):
        tz = "US/Central"
        ser = pd.Series(pd.date_range("2000", periods=2, tz=tz))
        expected = np.array(["2000-01-01T06", "2000-01-02T06"], dtype="M8[ns]")
        result = np.asarray(ser, dtype="datetime64[ns]")

        tm.assert_numpy_array_equal(result, expected)

        # Old behavior with no warning
        result = np.asarray(ser, dtype="M8[ns]")

        tm.assert_numpy_array_equal(result, expected)
