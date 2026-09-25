import datetime as dt
from datetime import datetime

import dateutil
import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestDatetimeConcat:
    def test_concat_datetime64_block(self):
        rng = pd.date_range("1/1/2000", periods=10)

        df = pd.DataFrame({"time": rng})

        result = pd.concat([df, df])
        assert (result.iloc[:10]["time"] == rng).all()
        assert (result.iloc[10:]["time"] == rng).all()

    def test_concat_datetime_datetime64_frame(self):
        # GH#2624
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), "hi"])

        df2_obj = pd.DataFrame.from_records(rows, columns=["date", "test"])

        ind = pd.date_range(start="2000/1/1", freq="D", periods=10)
        df1 = pd.DataFrame({"date": ind, "test": range(10)})

        # it works!
        pd.concat([df1, df2_obj])

    def test_concat_datetime_timezone(self):
        # GH 18523
        idx1 = pd.date_range(
            "2011-01-01", periods=3, freq="h", tz="Europe/Paris", unit="ns"
        )
        idx2 = pd.date_range(start=idx1[0], end=idx1[-1], freq="h", unit="ns")
        df1 = pd.DataFrame({"a": [1, 2, 3]}, index=idx1)
        df2 = pd.DataFrame({"b": [1, 2, 3]}, index=idx2)
        result = pd.concat([df1, df2], axis=1)

        exp_idx = pd.DatetimeIndex(
            [
                "2011-01-01 00:00:00+01:00",
                "2011-01-01 01:00:00+01:00",
                "2011-01-01 02:00:00+01:00",
            ],
            dtype="M8[ns, Europe/Paris]",
            freq="h",
        )
        expected = pd.DataFrame(
            [[1, 1], [2, 2], [3, 3]], index=exp_idx, columns=["a", "b"]
        )

        tm.assert_frame_equal(result, expected)

        idx3 = pd.date_range(
            "2011-01-01", periods=3, freq="h", tz="Asia/Tokyo", unit="ns"
        )
        df3 = pd.DataFrame({"b": [1, 2, 3]}, index=idx3)
        msg = "Sorting by default when concatenating all DatetimeIndex"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.concat([df1, df3], axis=1)

        exp_idx = pd.DatetimeIndex(
            [
                "2010-12-31 15:00:00+00:00",
                "2010-12-31 16:00:00+00:00",
                "2010-12-31 17:00:00+00:00",
                "2010-12-31 23:00:00+00:00",
                "2011-01-01 00:00:00+00:00",
                "2011-01-01 01:00:00+00:00",
            ]
        ).as_unit("ns")

        expected = pd.DataFrame(
            [
                [np.nan, 1],
                [np.nan, 2],
                [np.nan, 3],
                [1, np.nan],
                [2, np.nan],
                [3, np.nan],
            ],
            index=exp_idx,
            columns=["a", "b"],
        )

        tm.assert_frame_equal(result, expected)

        # GH 13783: Concat after resample
        result = pd.concat(
            [df1.resample("h").mean(), df2.resample("h").mean()], sort=True
        )
        expected = pd.DataFrame(
            {"a": [1, 2, 3] + [np.nan] * 3, "b": [np.nan] * 3 + [1, 2, 3]},
            index=idx1.append(idx1),
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_datetimeindex_freq(self):
        # GH 3232
        # Monotonic index result
        dr = pd.date_range("01-Jan-2013", periods=100, freq="50ms", tz="UTC")
        data = list(range(100))
        expected = pd.DataFrame(data, index=dr)
        result = pd.concat([expected[:50], expected[50:]])
        tm.assert_frame_equal(result, expected)

        # Non-monotonic index result
        result = pd.concat([expected[50:], expected[:50]])
        expected = pd.DataFrame(data[50:] + data[:50], index=dr[50:].append(dr[:50]))
        tm.assert_frame_equal(result, expected)

    def test_concat_datetimeindex_freq_mixed_unit(self):
        # GH#65920 - retain freq when concatting contiguous, evenly spaced
        # DatetimeIndexes that share a freq but differ in unit
        idx1 = pd.date_range("2024-01-01", periods=3, freq="s", unit="s")
        idx2 = pd.date_range("2024-01-01 00:00:03", periods=3, freq="s", unit="ns")
        result = idx1.append(idx2)
        expected = pd.date_range("2024-01-01", periods=6, freq="s", unit="ns")
        tm.assert_index_equal(result, expected)
        # Not checked by assert_index_equal
        assert result.freq == "s"

    def test_concat_datetimeindex_tz_convert_freq(self):
        # GH#41585 - concat after tz_convert should not raise when
        # the converted timestamps no longer conform to the original freq
        dti1 = pd.date_range(
            start="2020-01-01", end="2021-01-01", freq="MS", tz="CET", inclusive="left"
        ).tz_convert("UTC")
        df1 = pd.DataFrame({"full": [1] * len(dti1)}, index=dti1)

        dti2 = pd.date_range(
            start="2020-01-01", end="2021-02-01", freq="MS", tz="CET", inclusive="left"
        ).tz_convert("UTC")
        df2 = pd.DataFrame({"one_month_more": [1] * len(dti2)}, index=dti2)

        result = pd.concat([df1, df2], axis=1)
        expected_index = dti2.copy()
        expected_index.freq = result.index.freq
        expected = pd.DataFrame(
            {"full": [1.0] * 12 + [np.nan], "one_month_more": [1] * 13},
            index=expected_index,
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_multiindex_datetime_object_index(self):
        # https://github.com/pandas-dev/pandas/issues/11058
        idx = pd.Index(
            [dt.date(2013, 1, 1), dt.date(2014, 1, 1), dt.date(2015, 1, 1)],
            dtype="object",
        )

        s = pd.Series(
            ["a", "b"],
            index=pd.MultiIndex.from_arrays(
                [
                    [1, 2],
                    idx[:-1],
                ],
                names=["first", "second"],
            ),
        )
        s2 = pd.Series(
            ["a", "b"],
            index=pd.MultiIndex.from_arrays(
                [[1, 2], idx[::2]],
                names=["first", "second"],
            ),
        )
        mi = pd.MultiIndex.from_arrays(
            [[1, 2, 2], idx],
            names=["first", "second"],
        )
        assert mi.levels[1].dtype == object

        expected = pd.DataFrame(
            [["a", "a"], ["b", np.nan], [np.nan, "b"]],
            index=mi,
        )
        result = pd.concat([s, s2], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_concat_NaT_series(self):
        # GH 11693
        # test for merging NaT series with datetime series.
        x = pd.Series(
            pd.date_range(
                "20151124 08:00",
                "20151124 09:00",
                freq="1h",
                tz="US/Eastern",
                unit="ns",
            )
        )
        y = pd.Series(pd.NaT, index=[0, 1], dtype="datetime64[ns, US/Eastern]")
        expected = pd.Series([x[0], x[1], pd.NaT, pd.NaT])

        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT with tz
        expected = pd.Series(pd.NaT, index=range(4), dtype="datetime64[ns, US/Eastern]")
        result = pd.concat([y, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_NaT_series2(self):
        # without tz
        x = pd.Series(
            pd.date_range("20151124 08:00", "20151124 09:00", freq="1h", unit="ns")
        )
        y = pd.Series(
            pd.date_range("20151124 10:00", "20151124 11:00", freq="1h", unit="ns")
        )
        y[:] = pd.NaT
        expected = pd.Series([x[0], x[1], pd.NaT, pd.NaT])
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT without tz
        x[:] = pd.NaT
        expected = pd.Series(pd.NaT, index=range(4), dtype="datetime64[ns]")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "UTC"])
    def test_concat_NaT_dataframes(self, tz):
        # GH 12396

        dti = pd.DatetimeIndex([pd.NaT, pd.NaT], tz=tz)
        first = pd.DataFrame({0: dti})
        second = pd.DataFrame(
            [[pd.Timestamp("2015/01/01", tz=tz)], [pd.Timestamp("2016/01/01", tz=tz)]],
            index=[2, 3],
        )
        expected = pd.DataFrame(
            [
                pd.NaT,
                pd.NaT,
                pd.Timestamp("2015/01/01", tz=tz),
                pd.Timestamp("2016/01/01", tz=tz),
            ]
        )

        result = pd.concat([first, second], axis=0)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz1", [None, "UTC"])
    @pytest.mark.parametrize("tz2", [None, "UTC"])
    @pytest.mark.parametrize("item", [pd.NaT, pd.Timestamp("20150101").as_unit("ns")])
    def test_concat_NaT_dataframes_all_NaT_axis_0(self, tz1, tz2, item):
        # GH 12396

        # tz-naive
        first = pd.DataFrame([[pd.NaT], [pd.NaT]]).apply(
            lambda x: x.dt.tz_localize(tz1)
        )
        second = pd.DataFrame([item]).apply(lambda x: x.dt.tz_localize(tz2))

        result = pd.concat([first, second], axis=0)
        expected = pd.DataFrame(pd.Series([pd.NaT, pd.NaT, item], index=[0, 1, 0]))
        expected = expected.apply(lambda x: x.dt.tz_localize(tz2))
        if tz1 != tz2:
            expected = expected.astype(object)

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz1", [None, "UTC"])
    @pytest.mark.parametrize("tz2", [None, "UTC"])
    def test_concat_NaT_dataframes_all_NaT_axis_1(self, tz1, tz2):
        # GH 12396

        first = pd.DataFrame(pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1))
        second = pd.DataFrame(pd.Series([pd.NaT]).dt.tz_localize(tz2), columns=[1])
        expected = pd.DataFrame(
            {
                0: pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1),
                1: pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz2),
            }
        )
        result = pd.concat([first, second], axis=1)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz1", [None, "UTC"])
    @pytest.mark.parametrize("tz2", [None, "UTC"])
    def test_concat_NaT_series_dataframe_all_NaT(self, tz1, tz2):
        # GH 12396

        # tz-naive
        first = pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1)
        second = pd.DataFrame(
            [
                [pd.Timestamp("2015/01/01", tz=tz2)],
                [pd.Timestamp("2016/01/01", tz=tz2)],
            ],
            index=[2, 3],
        )

        expected = pd.DataFrame(
            [
                pd.NaT,
                pd.NaT,
                pd.Timestamp("2015/01/01", tz=tz2),
                pd.Timestamp("2016/01/01", tz=tz2),
            ]
        )
        if tz1 != tz2:
            expected = expected.astype(object)

        result = pd.concat([first, second])
        tm.assert_frame_equal(result, expected)

    def test_concat_compat_on_non_ns_datetime_EA(self):
        # GH#33331
        first = pd.Series(np.array([datetime(2010, 1, 1)], dtype="datetime64[D]"))
        second = pd.Series(pd.array(["a", "b"], dtype="category"))

        expected = pd.Series([pd.Timestamp("2010-01-01 00:00:00"), "a", "b"])

        result = pd.concat([first, second], ignore_index=True)
        tm.assert_series_equal(result, expected)


class TestTimezoneConcat:
    def test_concat_tz_series(self):
        # gh-11755: tz and no tz
        x = pd.Series(
            pd.date_range("20151124 08:00", "20151124 09:00", freq="1h", tz="UTC")
        )
        y = pd.Series(pd.date_range("2012-01-01", "2012-01-02"))
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_tz_series2(self):
        # gh-11887: concat tz and object
        x = pd.Series(
            pd.date_range("20151124 08:00", "20151124 09:00", freq="1h", tz="UTC")
        )
        y = pd.Series(["a", "b"])
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_tz_series3(self, unit, unit2):
        # see gh-12217 and gh-12306
        # Concatenating two UTC times
        first = pd.DataFrame([[datetime(2016, 1, 1)]], dtype=f"M8[{unit}]")
        first[0] = first[0].dt.tz_localize("UTC")

        second = pd.DataFrame([[datetime(2016, 1, 2)]], dtype=f"M8[{unit2}]")
        second[0] = second[0].dt.tz_localize("UTC")

        result = pd.concat([first, second])
        exp_unit = tm.get_finest_unit(unit, unit2)
        assert result[0].dtype == f"datetime64[{exp_unit}, UTC]"

    def test_concat_tz_series4(self, unit, unit2):
        # Concatenating two London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]], dtype=f"M8[{unit}]")
        first[0] = first[0].dt.tz_localize("Europe/London")

        second = pd.DataFrame([[datetime(2016, 1, 2)]], dtype=f"M8[{unit2}]")
        second[0] = second[0].dt.tz_localize("Europe/London")

        result = pd.concat([first, second])
        exp_unit = tm.get_finest_unit(unit, unit2)
        assert result[0].dtype == f"datetime64[{exp_unit}, Europe/London]"

    def test_concat_tz_series5(self, unit, unit2):
        # Concatenating 2+1 London times
        first = pd.DataFrame(
            [[datetime(2016, 1, 1)], [datetime(2016, 1, 2)]], dtype=f"M8[{unit}]"
        )
        first[0] = first[0].dt.tz_localize("Europe/London")

        second = pd.DataFrame([[datetime(2016, 1, 3)]], dtype=f"M8[{unit2}]")
        second[0] = second[0].dt.tz_localize("Europe/London")

        result = pd.concat([first, second])
        exp_unit = tm.get_finest_unit(unit, unit2)
        assert result[0].dtype == f"datetime64[{exp_unit}, Europe/London]"

    def test_concat_tz_series6(self, unit, unit2):
        # Concatenating 1+2 London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]], dtype=f"M8[{unit}]")
        first[0] = first[0].dt.tz_localize("Europe/London")

        second = pd.DataFrame(
            [[datetime(2016, 1, 2)], [datetime(2016, 1, 3)]], dtype=f"M8[{unit2}]"
        )
        second[0] = second[0].dt.tz_localize("Europe/London")

        result = pd.concat([first, second])
        exp_unit = tm.get_finest_unit(unit, unit2)
        assert result[0].dtype == f"datetime64[{exp_unit}, Europe/London]"

    def test_concat_tz_series_tzlocal(self):
        # see gh-13583
        x = [
            pd.Timestamp("2011-01-01", tz=dateutil.tz.tzlocal()),
            pd.Timestamp("2011-02-01", tz=dateutil.tz.tzlocal()),
        ]
        y = [
            pd.Timestamp("2012-01-01", tz=dateutil.tz.tzlocal()),
            pd.Timestamp("2012-02-01", tz=dateutil.tz.tzlocal()),
        ]

        result = pd.concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y))
        assert result.dtype == "datetime64[us, tzlocal()]"

    def test_concat_tz_series_with_datetimelike(self):
        # see gh-12620: tz and timedelta
        x = [
            pd.Timestamp("2011-01-01", tz="US/Eastern"),
            pd.Timestamp("2011-02-01", tz="US/Eastern"),
        ]
        y = [pd.Timedelta("1 day"), pd.Timedelta("2 day")]
        result = pd.concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype="object"))

        # tz and period
        y = [pd.Period("2011-03", freq="M"), pd.Period("2011-04", freq="M")]
        result = pd.concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype="object"))

    def test_concat_tz_frame(self):
        df2 = pd.DataFrame(
            {
                "A": pd.Timestamp("20130102", tz="US/Eastern"),
                "B": pd.Timestamp("20130603", tz="CET"),
            },
            index=range(5),
        )

        # concat
        df3 = pd.concat([df2.A.to_frame(), df2.B.to_frame()], axis=1)
        tm.assert_frame_equal(df2, df3)

    def test_concat_multiple_tzs(self):
        # GH#12467
        # combining datetime tz-aware and naive DataFrames
        ts1 = pd.Timestamp("2015-01-01", tz=None)
        ts2 = pd.Timestamp("2015-01-01", tz="UTC")
        ts3 = pd.Timestamp("2015-01-01", tz="EST")

        df1 = pd.DataFrame({"time": [ts1]})
        df2 = pd.DataFrame({"time": [ts2]})
        df3 = pd.DataFrame({"time": [ts3]})

        results = pd.concat([df1, df2]).reset_index(drop=True)
        expected = pd.DataFrame({"time": [ts1, ts2]}, dtype=object)
        tm.assert_frame_equal(results, expected)

        results = pd.concat([df1, df3]).reset_index(drop=True)
        expected = pd.DataFrame({"time": [ts1, ts3]}, dtype=object)
        tm.assert_frame_equal(results, expected)

        results = pd.concat([df2, df3]).reset_index(drop=True)
        expected = pd.DataFrame({"time": [ts2, ts3]})
        tm.assert_frame_equal(results, expected)

    def test_concat_multiindex_with_tz(self):
        # GH 6606
        df = pd.DataFrame(
            {
                "dt": pd.DatetimeIndex(
                    [
                        datetime(2014, 1, 1),
                        datetime(2014, 1, 2),
                        datetime(2014, 1, 3),
                    ],
                    dtype="M8[ns, US/Pacific]",
                ),
                "b": ["A", "B", "C"],
                "c": [1, 2, 3],
                "d": [4, 5, 6],
            }
        )
        df = df.set_index(["dt", "b"])

        exp_idx1 = pd.DatetimeIndex(
            ["2014-01-01", "2014-01-02", "2014-01-03"] * 2,
            dtype="M8[ns, US/Pacific]",
            name="dt",
        )
        exp_idx2 = pd.Index(["A", "B", "C"] * 2, name="b")
        exp_idx = pd.MultiIndex.from_arrays([exp_idx1, exp_idx2])
        expected = pd.DataFrame(
            {"c": [1, 2, 3] * 2, "d": [4, 5, 6] * 2}, index=exp_idx, columns=["c", "d"]
        )

        result = pd.concat([df, df])
        tm.assert_frame_equal(result, expected)

    def test_concat_tz_not_aligned(self):
        # GH#22796
        ts = pd.to_datetime([1, 2]).tz_localize("UTC")
        a = pd.DataFrame({"A": ts})
        b = pd.DataFrame({"A": ts, "B": ts})
        result = pd.concat([a, b], sort=True, ignore_index=True)
        expected = pd.DataFrame(
            {"A": list(ts) + list(ts), "B": [pd.NaT, pd.NaT, *list(ts)]}
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "t1",
        [
            "2015-01-01",
            pytest.param(
                pd.NaT,
                marks=pytest.mark.xfail(
                    reason="GH23037 incorrect dtype when concatenating"
                ),
            ),
        ],
    )
    def test_concat_tz_NaT(self, t1):
        # GH#22796
        # Concatenating tz-aware multicolumn DataFrames
        ts1 = pd.Timestamp(t1, tz="UTC")
        ts2 = pd.Timestamp("2015-01-01", tz="UTC")
        ts3 = pd.Timestamp("2015-01-01", tz="UTC")

        df1 = pd.DataFrame([[ts1, ts2]])
        df2 = pd.DataFrame([[ts3]])

        result = pd.concat([df1, df2])
        expected = pd.DataFrame([[ts1, ts2], [ts3, pd.NaT]], index=[0, 0])

        tm.assert_frame_equal(result, expected)

    def test_concat_tz_with_empty(self):
        # GH 9188
        result = pd.concat(
            [pd.DataFrame(pd.date_range("2000", periods=1, tz="UTC")), pd.DataFrame()]
        )
        expected = pd.DataFrame(pd.date_range("2000", periods=1, tz="UTC"))
        tm.assert_frame_equal(result, expected)


class TestPeriodConcat:
    def test_concat_period_series(self):
        x = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="D"))
        y = pd.Series(pd.PeriodIndex(["2015-10-01", "2016-01-01"], freq="D"))
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="Period[D]")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_period_multiple_freq_series(self):
        x = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="D"))
        y = pd.Series(pd.PeriodIndex(["2015-10-01", "2016-01-01"], freq="M"))
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == "object"

    def test_concat_period_other_series(self):
        x = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="D"))
        y = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="M"))
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == "object"

    def test_concat_period_other_series2(self):
        # non-period
        x = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="D"))
        y = pd.Series(pd.DatetimeIndex(["2015-11-01", "2015-12-01"]))
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == "object"

    def test_concat_period_other_series3(self):
        x = pd.Series(pd.PeriodIndex(["2015-11-01", "2015-12-01"], freq="D"))
        y = pd.Series(["A", "B"])
        expected = pd.Series([x[0], x[1], y[0], y[1]], dtype="object")
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == "object"

    def test_concat_keys_mixed_freq_period_columns(self):
        # GH#51489 concat with keys must not raise IncompatibleFrequency
        # when the inputs have PeriodIndex columns with different freq
        q_period = pd.period_range("2022-1-1", "2022-12-31", freq="Q")
        y_period = pd.period_range("2019-1-1", "2023-1-1", freq="Y")

        q_df = pd.DataFrame([range(len(q_period))], columns=q_period)
        y_df = pd.DataFrame([range(len(y_period))], columns=y_period)

        result = pd.concat([q_df, y_df], keys=["Quarterly", "Yearly"], axis=1)

        expected_inner = pd.Index(list(q_period) + list(y_period), dtype=object)
        expected_columns = pd.MultiIndex.from_arrays(
            [
                ["Quarterly"] * len(q_period) + ["Yearly"] * len(y_period),
                expected_inner,
            ]
        )
        expected = pd.DataFrame(
            [list(range(len(q_period))) + list(range(len(y_period)))],
            columns=expected_columns,
        )
        tm.assert_frame_equal(result, expected)


def test_concat_timedelta64_block():
    rng = pd.to_timedelta(np.arange(10), unit="s")

    df = pd.DataFrame({"time": rng})

    result = pd.concat([df, df])
    tm.assert_frame_equal(result.iloc[:10], df, check_index_type=False)
    tm.assert_frame_equal(result.iloc[10:], df, check_index_type=False)


def test_concat_multiindex_datetime_nat():
    # GH#44900
    left = pd.DataFrame({"a": 1}, index=pd.MultiIndex.from_tuples([(1, pd.NaT)]))
    right = pd.DataFrame(
        {"b": 2}, index=pd.MultiIndex.from_tuples([(1, pd.NaT), (2, pd.NaT)])
    )
    result = pd.concat([left, right], axis="columns")
    expected = pd.DataFrame(
        {"a": [1.0, np.nan], "b": 2},
        pd.MultiIndex.from_tuples([(1, pd.NaT), (2, pd.NaT)]),
    )
    tm.assert_frame_equal(result, expected)


def test_concat_float_datetime64():
    # GH#32934
    df_time = pd.DataFrame({"A": pd.array(["2000"], dtype="datetime64[ns]")})
    df_float = pd.DataFrame({"A": pd.array([1.0], dtype="float64")})

    expected = pd.DataFrame(
        {
            "A": [
                pd.array(["2000"], dtype="datetime64[ns]")[0],
                pd.array([1.0], dtype="float64")[0],
            ]
        },
        index=[0, 0],
    )
    result = pd.concat([df_time, df_float])
    tm.assert_frame_equal(result, expected)

    expected = pd.DataFrame({"A": pd.array([], dtype="object")})
    result = pd.concat([df_time.iloc[:0], df_float.iloc[:0]])
    tm.assert_frame_equal(result, expected)

    expected = pd.DataFrame({"A": pd.array([1.0], dtype="object")})
    result = pd.concat([df_time.iloc[:0], df_float])
    tm.assert_frame_equal(result, expected)

    expected = pd.DataFrame({"A": pd.array(["2000"], dtype="datetime64[ns]")}).astype(
        object
    )

    result = pd.concat([df_time, df_float.iloc[:0]])
    tm.assert_frame_equal(result, expected)


def test_concat_datetime64_different_resolutions():
    # GH#53307
    # Concatenating datetime64 columns with different resolutions
    # should preserve datetime64 dtype (not convert to object)
    df = pd.DataFrame(
        {
            "ints": range(2),
            "dates": pd.date_range("2000", periods=2, freq="min"),
        },
    )
    df2 = df.copy()
    df2["dates"] = df.dates.astype("M8[s]")

    combined = pd.concat([df, df2])

    # The result should be a datetime64 dtype, not object
    assert combined.dates.dtype.kind == "M"


def test_concat_non_ns_datetime_axis1(unit):
    # GH#58471 - concat with non-ns datetime unit on axis=1 should
    # preserve all data, matching the behavior of ns-resolution
    dti1 = pd.date_range("2024-01-01 00:00", periods=3, freq="5min", unit=unit)
    dti2 = pd.date_range("2024-01-01 00:15", periods=3, freq="5min", unit=unit)
    ser1 = pd.Series([1.0, 2.0, 3.0], index=dti1, name="a")
    ser2 = pd.Series([4.0, 5.0, 6.0], index=dti2, name="b")

    result = pd.concat([ser1, ser2], axis=1)

    expected_index = pd.date_range(
        "2024-01-01 00:00", periods=6, freq="5min", unit=unit
    )
    expected = pd.DataFrame(
        {
            "a": [1.0, 2.0, 3.0, np.nan, np.nan, np.nan],
            "b": [np.nan, np.nan, np.nan, 4.0, 5.0, 6.0],
        },
        index=expected_index,
    )
    tm.assert_frame_equal(result, expected)
