from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    Int64Index,
    Series,
    Timedelta,
    TimedeltaIndex,
    date_range,
    timedelta_range,
)
import pandas._testing as tm

from ..datetimelike import DatetimeLike

randn = np.random.randn


class TestTimedeltaIndex(DatetimeLike):
    _holder = TimedeltaIndex

    @pytest.fixture
    def indices(self):
        return tm.makeTimedeltaIndex(10)

    def create_index(self):
        return pd.to_timedelta(range(5), unit="d") + pd.offsets.Hour(1)

    def test_numeric_compat(self):
        # Dummy method to override super's version; this test is now done
        # in test_arithmetic.py
        pass

    def test_shift(self):
        pass  # this is handled in test_arithmetic.py

    def test_pickle_compat_construction(self):
        pass

    def test_fillna_timedelta(self):
        # GH 11343
        idx = pd.TimedeltaIndex(["1 day", pd.NaT, "3 day"])

        exp = pd.TimedeltaIndex(["1 day", "2 day", "3 day"])
        tm.assert_index_equal(idx.fillna(pd.Timedelta("2 day")), exp)

        exp = pd.TimedeltaIndex(["1 day", "3 hour", "3 day"])
        idx.fillna(pd.Timedelta("3 hour"))

        exp = pd.Index(
            [pd.Timedelta("1 day"), "x", pd.Timedelta("3 day")], dtype=object
        )
        tm.assert_index_equal(idx.fillna("x"), exp)

    def test_isin(self):

        index = tm.makeTimedeltaIndex(4)
        result = index.isin(index)
        assert result.all()

        result = index.isin(list(index))
        assert result.all()

        tm.assert_almost_equal(
            index.isin([index[2], 5]), np.array([False, False, True, False])
        )

    def test_factorize(self):
        idx1 = TimedeltaIndex(["1 day", "1 day", "2 day", "2 day", "3 day", "3 day"])

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = TimedeltaIndex(["1 day", "2 day", "3 day"])

        arr, idx = idx1.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # freq must be preserved
        idx3 = timedelta_range("1 day", periods=4, freq="s")
        exp_arr = np.array([0, 1, 2, 3], dtype=np.intp)
        arr, idx = idx3.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, idx3)

    def test_join_self(self, join_type):
        index = timedelta_range("1 day", periods=10)
        joined = index.join(index, how=join_type)
        tm.assert_index_equal(index, joined)

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(
            10,
            10,
            data_gen_f=lambda *args, **kwargs: randn(),
            r_idx_type="i",
            c_idx_type="td",
        )
        str(df)

        cols = df.columns.join(df.index, how="outer")
        joined = cols.join(df.columns)
        assert cols.dtype == np.dtype("O")
        assert cols.dtype == joined.dtype
        tm.assert_index_equal(cols, joined)

    def test_sort_values(self):

        idx = TimedeltaIndex(["4d", "1d", "2d"])

        ordered = idx.sort_values()
        assert ordered.is_monotonic

        ordered = idx.sort_values(ascending=False)
        assert ordered[::-1].is_monotonic

        ordered, dexer = idx.sort_values(return_indexer=True)
        assert ordered.is_monotonic

        tm.assert_numpy_array_equal(dexer, np.array([1, 2, 0]), check_dtype=False)

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        assert ordered[::-1].is_monotonic

        tm.assert_numpy_array_equal(dexer, np.array([0, 2, 1]), check_dtype=False)

    def test_argmin_argmax(self):
        idx = TimedeltaIndex(["1 day 00:00:05", "1 day 00:00:01", "1 day 00:00:02"])
        assert idx.argmin() == 1
        assert idx.argmax() == 0

    def test_misc_coverage(self):

        rng = timedelta_range("1 day", periods=5)
        result = rng.groupby(rng.days)
        assert isinstance(list(result.values())[0][0], Timedelta)

        idx = TimedeltaIndex(["3d", "1d", "2d"])
        assert not idx.equals(list(idx))

        non_td = Index(list("abc"))
        assert not idx.equals(list(non_td))

    def test_map(self):
        # test_map_dictlike generally tests

        rng = timedelta_range("1 day", periods=10)

        f = lambda x: x.days
        result = rng.map(f)
        exp = Int64Index([f(x) for x in rng])
        tm.assert_index_equal(result, exp)

    def test_pass_TimedeltaIndex_to_index(self):

        rng = timedelta_range("1 days", "10 days")
        idx = Index(rng, dtype=object)

        expected = Index(rng.to_pytimedelta(), dtype=object)

        tm.assert_numpy_array_equal(idx.values, expected.values)

    def test_pickle(self):

        rng = timedelta_range("1 days", periods=10)
        rng_p = tm.round_trip_pickle(rng)
        tm.assert_index_equal(rng, rng_p)

    def test_hash_error(self):
        index = timedelta_range("1 days", periods=10)
        with pytest.raises(
            TypeError, match=(f"unhashable type: {repr(type(index).__name__)}")
        ):
            hash(index)

    def test_append_join_nondatetimeindex(self):
        rng = timedelta_range("1 days", periods=10)
        idx = Index(["a", "b", "c", "d"])

        result = rng.append(idx)
        assert isinstance(result[0], Timedelta)

        # it works
        rng.join(idx, how="outer")

    def test_append_numpy_bug_1681(self):

        td = timedelta_range("1 days", "10 days", freq="2D")
        a = DataFrame()
        c = DataFrame({"A": "foo", "B": td}, index=td)
        str(c)

        result = a.append(c)
        assert (result["B"] == td).all()

    def test_fields(self):
        rng = timedelta_range("1 days, 10:11:12.100123456", periods=2, freq="s")
        tm.assert_index_equal(rng.days, Index([1, 1], dtype="int64"))
        tm.assert_index_equal(
            rng.seconds,
            Index([10 * 3600 + 11 * 60 + 12, 10 * 3600 + 11 * 60 + 13], dtype="int64"),
        )
        tm.assert_index_equal(
            rng.microseconds, Index([100 * 1000 + 123, 100 * 1000 + 123], dtype="int64")
        )
        tm.assert_index_equal(rng.nanoseconds, Index([456, 456], dtype="int64"))

        msg = "'TimedeltaIndex' object has no attribute '{}'"
        with pytest.raises(AttributeError, match=msg.format("hours")):
            rng.hours
        with pytest.raises(AttributeError, match=msg.format("minutes")):
            rng.minutes
        with pytest.raises(AttributeError, match=msg.format("milliseconds")):
            rng.milliseconds

        # with nat
        s = Series(rng)
        s[1] = np.nan

        tm.assert_series_equal(s.dt.days, Series([1, np.nan], index=[0, 1]))
        tm.assert_series_equal(
            s.dt.seconds, Series([10 * 3600 + 11 * 60 + 12, np.nan], index=[0, 1])
        )

        # preserve name (GH15589)
        rng.name = "name"
        assert rng.days.name == "name"

    def test_freq_conversion(self):

        # doc example

        # series
        td = Series(date_range("20130101", periods=4)) - Series(
            date_range("20121201", periods=4)
        )
        td[2] += timedelta(minutes=5, seconds=3)
        td[3] = np.nan

        result = td / np.timedelta64(1, "D")
        expected = Series([31, 31, (31 * 86400 + 5 * 60 + 3) / 86400.0, np.nan])
        tm.assert_series_equal(result, expected)

        result = td.astype("timedelta64[D]")
        expected = Series([31, 31, 31, np.nan])
        tm.assert_series_equal(result, expected)

        result = td / np.timedelta64(1, "s")
        expected = Series([31 * 86400, 31 * 86400, 31 * 86400 + 5 * 60 + 3, np.nan])
        tm.assert_series_equal(result, expected)

        result = td.astype("timedelta64[s]")
        tm.assert_series_equal(result, expected)

        # tdi
        td = TimedeltaIndex(td)

        result = td / np.timedelta64(1, "D")
        expected = Index([31, 31, (31 * 86400 + 5 * 60 + 3) / 86400.0, np.nan])
        tm.assert_index_equal(result, expected)

        result = td.astype("timedelta64[D]")
        expected = Index([31, 31, 31, np.nan])
        tm.assert_index_equal(result, expected)

        result = td / np.timedelta64(1, "s")
        expected = Index([31 * 86400, 31 * 86400, 31 * 86400 + 5 * 60 + 3, np.nan])
        tm.assert_index_equal(result, expected)

        result = td.astype("timedelta64[s]")
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("unit", ["Y", "y", "M"])
    def test_unit_m_y_raises(self, unit):
        msg = "Units 'M' and 'Y' are no longer supported"
        with pytest.raises(ValueError, match=msg):
            TimedeltaIndex([1, 3, 7], unit)


class TestTimeSeries:
    def test_series_box_timedelta(self):
        rng = timedelta_range("1 day 1 s", periods=5, freq="h")
        s = Series(rng)
        assert isinstance(s[1], Timedelta)
        assert isinstance(s.iat[2], Timedelta)
