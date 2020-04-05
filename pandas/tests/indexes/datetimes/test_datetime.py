from datetime import date, timedelta

import dateutil
import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, DatetimeIndex, Index, NaT, Timestamp, date_range, offsets
import pandas._testing as tm

randn = np.random.randn


class TestDatetimeIndex:
    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = date_range("20130101", periods=3, tz="US/Eastern", name="foo")
        unpickled = tm.round_trip_pickle(index)
        tm.assert_index_equal(index, unpickled)

    def test_pickle(self):

        # GH#4606
        p = tm.round_trip_pickle(NaT)
        assert p is NaT

        idx = pd.to_datetime(["2013-01-01", NaT, "2014-01-06"])
        idx_p = tm.round_trip_pickle(idx)
        assert idx_p[0] == idx[0]
        assert idx_p[1] is NaT
        assert idx_p[2] == idx[2]

        # GH#11002
        # don't infer freq
        idx = date_range("1750-1-1", "2050-1-1", freq="7D")
        idx_p = tm.round_trip_pickle(idx)
        tm.assert_index_equal(idx, idx_p)

    def test_reindex_preserves_tz_if_target_is_empty_list_or_array(self):
        # GH7774
        index = date_range("20130101", periods=3, tz="US/Eastern")
        assert str(index.reindex([])[0].tz) == "US/Eastern"
        assert str(index.reindex(np.array([]))[0].tz) == "US/Eastern"

    def test_reindex_with_same_tz(self):
        # GH 32740
        rng_a = date_range("2010-01-01", "2010-01-02", periods=24, tz="utc")
        rng_b = date_range("2010-01-01", "2010-01-02", periods=23, tz="utc")
        result1, result2 = rng_a.reindex(
            rng_b, method="nearest", tolerance=timedelta(seconds=20)
        )
        expected_list1 = [
            "2010-01-01 00:00:00",
            "2010-01-01 01:05:27.272727272",
            "2010-01-01 02:10:54.545454545",
            "2010-01-01 03:16:21.818181818",
            "2010-01-01 04:21:49.090909090",
            "2010-01-01 05:27:16.363636363",
            "2010-01-01 06:32:43.636363636",
            "2010-01-01 07:38:10.909090909",
            "2010-01-01 08:43:38.181818181",
            "2010-01-01 09:49:05.454545454",
            "2010-01-01 10:54:32.727272727",
            "2010-01-01 12:00:00",
            "2010-01-01 13:05:27.272727272",
            "2010-01-01 14:10:54.545454545",
            "2010-01-01 15:16:21.818181818",
            "2010-01-01 16:21:49.090909090",
            "2010-01-01 17:27:16.363636363",
            "2010-01-01 18:32:43.636363636",
            "2010-01-01 19:38:10.909090909",
            "2010-01-01 20:43:38.181818181",
            "2010-01-01 21:49:05.454545454",
            "2010-01-01 22:54:32.727272727",
            "2010-01-02 00:00:00",
        ]
        expected1 = DatetimeIndex(
            expected_list1, dtype="datetime64[ns, UTC]", freq=None,
        )
        expected2 = np.array([0] + [-1] * 21 + [23], dtype=np.int64,)
        tm.assert_index_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    def test_time_loc(self):  # GH8667
        from datetime import time
        from pandas._libs.index import _SIZE_CUTOFF

        ns = _SIZE_CUTOFF + np.array([-100, 100], dtype=np.int64)
        key = time(15, 11, 30)
        start = key.hour * 3600 + key.minute * 60 + key.second
        step = 24 * 3600

        for n in ns:
            idx = pd.date_range("2014-11-26", periods=n, freq="S")
            ts = pd.Series(np.random.randn(n), index=idx)
            i = np.arange(start, n, step)

            tm.assert_numpy_array_equal(ts.index.get_loc(key), i, check_dtype=False)
            tm.assert_series_equal(ts[key], ts.iloc[i])

            left, right = ts.copy(), ts.copy()
            left[key] *= -10
            right.iloc[i] *= -10
            tm.assert_series_equal(left, right)

    def test_time_overflow_for_32bit_machines(self):
        # GH8943.  On some machines NumPy defaults to np.int32 (for example,
        # 32-bit Linux machines).  In the function _generate_regular_range
        # found in tseries/index.py, `periods` gets multiplied by `strides`
        # (which has value 1e9) and since the max value for np.int32 is ~2e9,
        # and since those machines won't promote np.int32 to np.int64, we get
        # overflow.
        periods = np.int_(1000)

        idx1 = pd.date_range(start="2000", periods=periods, freq="S")
        assert len(idx1) == periods

        idx2 = pd.date_range(end="2000", periods=periods, freq="S")
        assert len(idx2) == periods

    def test_nat(self):
        assert DatetimeIndex([np.nan])[0] is pd.NaT

    def test_week_of_month_frequency(self):
        # GH 5348: "ValueError: Could not evaluate WOM-1SUN" shouldn't raise
        d1 = date(2002, 9, 1)
        d2 = date(2013, 10, 27)
        d3 = date(2012, 9, 30)
        idx1 = DatetimeIndex([d1, d2])
        idx2 = DatetimeIndex([d3])
        result_append = idx1.append(idx2)
        expected = DatetimeIndex([d1, d2, d3])
        tm.assert_index_equal(result_append, expected)
        result_union = idx1.union(idx2)
        expected = DatetimeIndex([d1, d3, d2])
        tm.assert_index_equal(result_union, expected)

        # GH 5115
        result = date_range("2013-1-1", periods=4, freq="WOM-1SAT")
        dates = ["2013-01-05", "2013-02-02", "2013-03-02", "2013-04-06"]
        expected = DatetimeIndex(dates, freq="WOM-1SAT")
        tm.assert_index_equal(result, expected)

    def test_stringified_slice_with_tz(self):
        # GH#2658
        start = "2013-01-07"
        idx = date_range(start=start, freq="1d", periods=10, tz="US/Eastern")
        df = DataFrame(np.arange(10), index=idx)
        df["2013-01-14 23:44:34.437768-05:00":]  # no exception here

    def test_append_nondatetimeindex(self):
        rng = date_range("1/1/2000", periods=10)
        idx = Index(["a", "b", "c", "d"])

        result = rng.append(idx)
        assert isinstance(result[0], Timestamp)

    def test_map(self):
        rng = date_range("1/1/2000", periods=10)

        f = lambda x: x.strftime("%Y%m%d")
        result = rng.map(f)
        exp = Index([f(x) for x in rng], dtype="<U8")
        tm.assert_index_equal(result, exp)

    def test_map_fallthrough(self, capsys):
        # GH#22067, check we don't get warnings about silently ignored errors
        dti = date_range("2017-01-01", "2018-01-01", freq="B")

        dti.map(lambda x: pd.Period(year=x.year, month=x.month, freq="M"))

        captured = capsys.readouterr()
        assert captured.err == ""

    def test_iteration_preserves_tz(self):
        # see gh-8890
        index = date_range("2012-01-01", periods=3, freq="H", tz="US/Eastern")

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result == expected

        index = date_range(
            "2012-01-01", periods=3, freq="H", tz=dateutil.tz.tzoffset(None, -28800)
        )

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result._repr_base == expected._repr_base
            assert result == expected

        # 9100
        index = pd.DatetimeIndex(
            ["2014-12-01 03:32:39.987000-08:00", "2014-12-01 04:12:34.987000-08:00"]
        )
        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result._repr_base == expected._repr_base
            assert result == expected

    @pytest.mark.parametrize("periods", [0, 9999, 10000, 10001])
    def test_iteration_over_chunksize(self, periods):
        # GH21012

        index = date_range("2000-01-01 00:00:00", periods=periods, freq="min")
        num = 0
        for stamp in index:
            assert index[num] == stamp
            num += 1
        assert num == len(index)

    def test_misc_coverage(self):
        rng = date_range("1/1/2000", periods=5)
        result = rng.groupby(rng.day)
        assert isinstance(list(result.values())[0][0], Timestamp)

        idx = DatetimeIndex(["2000-01-03", "2000-01-01", "2000-01-02"])
        assert not idx.equals(list(idx))

        non_datetime = Index(list("abc"))
        assert not idx.equals(list(non_datetime))

    def test_string_index_series_name_converted(self):
        # #1644
        df = DataFrame(np.random.randn(10, 4), index=date_range("1/1/2000", periods=10))

        result = df.loc["1/3/2000"]
        assert result.name == df.index[2]

        result = df.T["1/3/2000"]
        assert result.name == df.index[2]

    def test_argmin_argmax(self):
        idx = DatetimeIndex(["2000-01-04", "2000-01-01", "2000-01-02"])
        assert idx.argmin() == 1
        assert idx.argmax() == 0

    def test_sort_values(self):
        idx = DatetimeIndex(["2000-01-04", "2000-01-01", "2000-01-02"])

        ordered = idx.sort_values()
        assert ordered.is_monotonic

        ordered = idx.sort_values(ascending=False)
        assert ordered[::-1].is_monotonic

        ordered, dexer = idx.sort_values(return_indexer=True)
        assert ordered.is_monotonic
        tm.assert_numpy_array_equal(dexer, np.array([1, 2, 0], dtype=np.intp))

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        assert ordered[::-1].is_monotonic
        tm.assert_numpy_array_equal(dexer, np.array([0, 2, 1], dtype=np.intp))

    def test_map_bug_1677(self):
        index = DatetimeIndex(["2012-04-25 09:30:00.393000"])
        f = index.asof

        result = index.map(f)
        expected = Index([f(index[0])])
        tm.assert_index_equal(result, expected)

    def test_groupby_function_tuple_1677(self):
        df = DataFrame(np.random.rand(100), index=date_range("1/1/2000", periods=100))
        monthly_group = df.groupby(lambda x: (x.year, x.month))

        result = monthly_group.mean()
        assert isinstance(result.index[0], tuple)

    def test_append_numpy_bug_1681(self):
        # another datetime64 bug
        dr = date_range("2011/1/1", "2012/1/1", freq="W-FRI")
        a = DataFrame()
        c = DataFrame({"A": "foo", "B": dr}, index=dr)

        result = a.append(c)
        assert (result["B"] == dr).all()

    def test_isin(self):
        index = tm.makeDateIndex(4)
        result = index.isin(index)
        assert result.all()

        result = index.isin(list(index))
        assert result.all()

        tm.assert_almost_equal(
            index.isin([index[2], 5]), np.array([False, False, True, False])
        )

    def assert_index_parameters(self, index):
        assert index.freq == "40960N"
        assert index.inferred_freq == "40960N"

    def test_ns_index(self):
        nsamples = 400
        ns = int(1e9 / 24414)
        dtstart = np.datetime64("2012-09-20T00:00:00")

        dt = dtstart + np.arange(nsamples) * np.timedelta64(ns, "ns")
        freq = ns * offsets.Nano()
        index = pd.DatetimeIndex(dt, freq=freq, name="time")
        self.assert_index_parameters(index)

        new_index = pd.date_range(start=index[0], end=index[-1], freq=index.freq)
        self.assert_index_parameters(new_index)

    def test_factorize(self):
        idx1 = DatetimeIndex(
            ["2014-01", "2014-01", "2014-02", "2014-02", "2014-03", "2014-03"]
        )

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(["2014-01", "2014-02", "2014-03"])

        arr, idx = idx1.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # tz must be preserved
        idx1 = idx1.tz_localize("Asia/Tokyo")
        exp_idx = exp_idx.tz_localize("Asia/Tokyo")

        arr, idx = idx1.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        idx2 = pd.DatetimeIndex(
            ["2014-03", "2014-03", "2014-02", "2014-01", "2014-03", "2014-01"]
        )

        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.intp)
        exp_idx = DatetimeIndex(["2014-01", "2014-02", "2014-03"])
        arr, idx = idx2.factorize(sort=True)
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        exp_arr = np.array([0, 0, 1, 2, 0, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(["2014-03", "2014-02", "2014-01"])
        arr, idx = idx2.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # freq must be preserved
        idx3 = date_range("2000-01", periods=4, freq="M", tz="Asia/Tokyo")
        exp_arr = np.array([0, 1, 2, 3], dtype=np.intp)
        arr, idx = idx3.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, idx3)

    def test_factorize_tz(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH#13750
        base = pd.date_range("2016-11-05", freq="H", periods=100, tz=tz)
        idx = base.repeat(5)

        exp_arr = np.arange(100, dtype=np.intp).repeat(5)

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, exp_arr)
            tm.assert_index_equal(res, base)

    def test_factorize_dst(self):
        # GH 13750
        idx = pd.date_range("2016-11-06", freq="H", periods=12, tz="US/Eastern")

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

        idx = pd.date_range("2016-06-13", freq="H", periods=12, tz="US/Eastern")

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

    @pytest.mark.parametrize(
        "arr, expected",
        [
            (pd.DatetimeIndex(["2017", "2017"]), pd.DatetimeIndex(["2017"])),
            (
                pd.DatetimeIndex(["2017", "2017"], tz="US/Eastern"),
                pd.DatetimeIndex(["2017"], tz="US/Eastern"),
            ),
        ],
    )
    def test_unique(self, arr, expected):
        result = arr.unique()
        tm.assert_index_equal(result, expected)
        # GH 21737
        # Ensure the underlying data is consistent
        assert result[0] == expected[0]

    def test_asarray_tz_naive(self):
        # This shouldn't produce a warning.
        idx = pd.date_range("2000", periods=2)
        # M8[ns] by default
        result = np.asarray(idx)

        expected = np.array(["2000-01-01", "2000-01-02"], dtype="M8[ns]")
        tm.assert_numpy_array_equal(result, expected)

        # optionally, object
        result = np.asarray(idx, dtype=object)

        expected = np.array([pd.Timestamp("2000-01-01"), pd.Timestamp("2000-01-02")])
        tm.assert_numpy_array_equal(result, expected)

    def test_asarray_tz_aware(self):
        tz = "US/Central"
        idx = pd.date_range("2000", periods=2, tz=tz)
        expected = np.array(["2000-01-01T06", "2000-01-02T06"], dtype="M8[ns]")
        result = np.asarray(idx, dtype="datetime64[ns]")

        tm.assert_numpy_array_equal(result, expected)

        # Old behavior with no warning
        result = np.asarray(idx, dtype="M8[ns]")

        tm.assert_numpy_array_equal(result, expected)

        # Future behavior with no warning
        expected = np.array(
            [pd.Timestamp("2000-01-01", tz=tz), pd.Timestamp("2000-01-02", tz=tz)]
        )
        result = np.asarray(idx, dtype=object)

        tm.assert_numpy_array_equal(result, expected)

    def test_to_frame_datetime_tz(self):
        # GH 25809
        idx = date_range(start="2019-01-01", end="2019-01-30", freq="D", tz="UTC")
        result = idx.to_frame()
        expected = DataFrame(idx, index=idx)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("name", [None, "name"])
    def test_index_map(self, name):
        # see GH20990
        count = 6
        index = pd.date_range("2018-01-01", periods=count, freq="M", name=name).map(
            lambda x: (x.year, x.month)
        )
        exp_index = pd.MultiIndex.from_product(
            ((2018,), range(1, 7)), names=[name, name]
        )
        tm.assert_index_equal(index, exp_index)

    def test_split_non_utc(self):
        # GH 14042
        indices = pd.date_range("2016-01-01 00:00:00+0200", freq="S", periods=10)
        result = np.split(indices, indices_or_sections=[])[0]
        expected = indices.copy()
        expected._set_freq(None)
        tm.assert_index_equal(result, expected)
