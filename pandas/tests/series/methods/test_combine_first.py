from datetime import datetime

import numpy as np

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestCombineFirst:
    def test_combine_first_period_datetime(self):
        # GH#3367
        didx = pd.date_range(start="1950-01-31", end="1950-07-31", freq="ME")
        pidx = pd.period_range(
            start=pd.Period("1950-1"), end=pd.Period("1950-7"), freq="M"
        )
        # check to be consistent with DatetimeIndex
        for idx in [didx, pidx]:
            a = pd.Series([1, np.nan, np.nan, 4, 5, np.nan, 7], index=idx)
            b = pd.Series([9, 9, 9, 9, 9, 9, 9], index=idx)

            result = a.combine_first(b)
            expected = pd.Series([1, 9, 9, 4, 5, 9, 7], index=idx, dtype=np.float64)
            tm.assert_series_equal(result, expected)

    def test_combine_first_name(self, datetime_series):
        result = datetime_series.combine_first(datetime_series[:5])
        assert result.name == datetime_series.name

    def test_combine_first(self, using_infer_string):
        values = np.arange(20, dtype=np.float64)
        series = pd.Series(values, index=np.arange(20, dtype=np.int64))

        series_copy = series * 2
        series_copy[::2] = np.nan

        # nothing used from the input
        combined = series.combine_first(series_copy)

        tm.assert_series_equal(combined, series)

        # Holes filled from input
        combined = series_copy.combine_first(series)
        assert np.isfinite(combined).all()

        tm.assert_series_equal(combined[::2], series[::2])
        tm.assert_series_equal(combined[1::2], series_copy[1::2])

        # mixed types
        index = pd.Index([str(i) for i in range(20)])
        floats = pd.Series(np.random.default_rng(2).standard_normal(20), index=index)
        strings = pd.Series([str(i) for i in range(10)], index=index[::2], dtype=object)

        combined = strings.combine_first(floats)

        tm.assert_series_equal(strings, combined.loc[index[::2]])
        tm.assert_series_equal(floats[1::2].astype(object), combined.loc[index[1::2]])

        # corner case
        ser = pd.Series([1.0, 2, 3], index=[0, 1, 2])
        empty = pd.Series([], index=[], dtype=object)
        result = ser.combine_first(empty)
        if not using_infer_string:
            ser.index = ser.index.astype("O")
        tm.assert_series_equal(result, ser.astype(object))

    def test_combine_first_dt64(self, unit):
        s0 = pd.to_datetime(pd.Series(["2010", np.nan])).dt.as_unit(unit)
        s1 = pd.to_datetime(pd.Series([np.nan, "2011"])).dt.as_unit(unit)
        rs = s0.combine_first(s1)
        xp = pd.to_datetime(pd.Series(["2010", "2011"])).dt.as_unit(unit)
        tm.assert_series_equal(rs, xp)

    def test_combine_first_dt64_casting_deprecation(self, unit):
        # GH#62931
        s0 = pd.to_datetime(pd.Series(["2010", np.nan])).dt.as_unit(unit)
        s1 = pd.Series([np.nan, "2011"])

        msg = "Silently casting non-datetime 'other' to datetime"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            rs = s0.combine_first(s1)

        xp = pd.Series([datetime(2010, 1, 1), "2011"], dtype=f"datetime64[{unit}]")
        if unit in ["s", "ms"]:
            # TODO: should _cast_pointwise_result attempt to preserve unit?
            xp = xp.dt.as_unit("us")
        tm.assert_series_equal(rs, xp)

    def test_combine_first_dt_tz_values(self, tz_naive_fixture):
        ser1 = pd.Series(
            pd.DatetimeIndex(["20150101", "20150102", "20150103"], tz=tz_naive_fixture),
            name="ser1",
        )
        ser2 = pd.Series(
            pd.DatetimeIndex(["20160514", "20160515", "20160516"], tz=tz_naive_fixture),
            index=[2, 3, 4],
            name="ser2",
        )
        result = ser1.combine_first(ser2)
        exp_vals = pd.DatetimeIndex(
            ["20150101", "20150102", "20150103", "20160515", "20160516"],
            tz=tz_naive_fixture,
        )
        exp = pd.Series(exp_vals, name="ser1")
        tm.assert_series_equal(exp, result)

    def test_combine_first_timezone_series_with_empty_series(self):
        # GH 41800
        time_index = pd.date_range(
            datetime(2021, 1, 1, 1),
            datetime(2021, 1, 1, 10),
            freq="h",
            tz="Europe/Rome",
        )
        s1 = pd.Series(range(10), index=time_index)
        s2 = pd.Series(index=time_index)
        result = s1.combine_first(s2)
        tm.assert_series_equal(result, s1.astype(np.float64))

    def test_combine_first_preserves_dtype(self):
        # GH51764
        s1 = pd.Series([1666880195890293744, 1666880195890293837])
        s2 = pd.Series([1, 2, 3])
        result = s1.combine_first(s2)
        expected = pd.Series([1666880195890293744, 1666880195890293837, 3])
        tm.assert_series_equal(result, expected)

    def test_combine_mixed_timezone(self):
        # GH 26283
        uniform_tz = pd.Series({pd.Timestamp("2019-05-01", tz="UTC"): 1.0})
        multi_tz = pd.Series(
            {
                pd.Timestamp("2019-05-01 01:00:00+0100", tz="Europe/London"): 2.0,
                pd.Timestamp("2019-05-02", tz="UTC"): 3.0,
            }
        )

        result = uniform_tz.combine_first(multi_tz)
        expected = pd.Series(
            [1.0, 3.0],
            index=pd.Index(
                [
                    pd.Timestamp("2019-05-01 00:00:00+00:00", tz="UTC"),
                    pd.Timestamp("2019-05-02 00:00:00+00:00", tz="UTC"),
                ],
                dtype="object",
            ),
        )
        tm.assert_series_equal(result, expected)

    def test_combine_first_none_not_nan(self):
        # GH#58977
        s1 = pd.Series([None, None, None], index=["a", "b", "c"])
        s2 = pd.Series([None, None, None], index=["b", "c", "d"])

        result = s1.combine_first(s2)
        expected = pd.Series([None] * 4, index=["a", "b", "c", "d"])
        tm.assert_series_equal(result, expected)


def test_combine_first_timestamp_names_anterior():
    # GH#65333
    s1 = pd.Series([0], name=pd.to_datetime("2026"))
    s3 = pd.Series([1, 3], name=pd.to_datetime("2025"))
    result = s1.combine_first(s3)
    expected = pd.Series([0, 3], index=[0, 1], name=pd.to_datetime("2026"))
    tm.assert_series_equal(result, expected)
