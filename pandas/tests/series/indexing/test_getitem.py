"""
Series.__getitem__ test classes are organized by the type of key passed.
"""
from datetime import date, datetime, time

import numpy as np
import pytest

from pandas._libs.tslibs import conversion, timezones

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    DatetimeIndex,
    Index,
    Series,
    Timestamp,
    date_range,
    period_range,
)
import pandas._testing as tm
from pandas.core.indexing import IndexingError

from pandas.tseries.offsets import BDay


class TestSeriesGetitemScalars:
    def test_getitem_out_of_bounds_indexerror(self, datetime_series):
        # don't segfault, GH#495
        msg = r"index \d+ is out of bounds for axis 0 with size \d+"
        with pytest.raises(IndexError, match=msg):
            datetime_series[len(datetime_series)]

    def test_getitem_out_of_bounds_empty_rangeindex_keyerror(self):
        # GH#917
        # With a RangeIndex, an int key gives a KeyError
        ser = Series([], dtype=object)
        with pytest.raises(KeyError, match="-1"):
            ser[-1]

    def test_getitem_keyerror_with_int64index(self):
        ser = Series(np.random.randn(6), index=[0, 0, 1, 1, 2, 2])

        with pytest.raises(KeyError, match=r"^5$"):
            ser[5]

        with pytest.raises(KeyError, match=r"^'c'$"):
            ser["c"]

        # not monotonic
        ser = Series(np.random.randn(6), index=[2, 2, 0, 0, 1, 1])

        with pytest.raises(KeyError, match=r"^5$"):
            ser[5]

        with pytest.raises(KeyError, match=r"^'c'$"):
            ser["c"]

    def test_getitem_int64(self, datetime_series):
        idx = np.int64(5)
        assert datetime_series[idx] == datetime_series[5]

    # TODO: better name/GH ref?
    def test_getitem_regression(self):
        ser = Series(range(5), index=list(range(5)))
        result = ser[list(range(5))]
        tm.assert_series_equal(result, ser)

    # ------------------------------------------------------------------
    # Series with DatetimeIndex

    @pytest.mark.parametrize("tzstr", ["Europe/Berlin", "dateutil/Europe/Berlin"])
    def test_getitem_pydatetime_tz(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        index = date_range(
            start="2012-12-24 16:00", end="2012-12-24 18:00", freq="H", tz=tzstr
        )
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp("2012-12-24 17:00", tz=tzstr)

        dt = datetime(2012, 12, 24, 17, 0)
        time_datetime = conversion.localize_pydatetime(dt, tz)
        assert ts[time_pandas] == ts[time_datetime]

    @pytest.mark.parametrize("tz", ["US/Eastern", "dateutil/US/Eastern"])
    def test_string_index_alias_tz_aware(self, tz):
        rng = date_range("1/1/2000", periods=10, tz=tz)
        ser = Series(np.random.randn(len(rng)), index=rng)

        result = ser["1/3/2000"]
        tm.assert_almost_equal(result, ser[2])

    def test_getitem_time_object(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = Series(np.random.randn(len(rng)), index=rng)

        mask = (rng.hour == 9) & (rng.minute == 30)
        result = ts[time(9, 30)]
        expected = ts[mask]
        result.index = result.index._with_freq(None)
        tm.assert_series_equal(result, expected)

    # ------------------------------------------------------------------
    # Series with CategoricalIndex

    def test_getitem_scalar_categorical_index(self):
        cats = Categorical([Timestamp("12-31-1999"), Timestamp("12-31-2000")])

        ser = Series([1, 2], index=cats)

        expected = ser.iloc[0]
        result = ser[cats[0]]
        assert result == expected


class TestSeriesGetitemSlices:
    def test_getitem_partial_str_slice_with_datetimeindex(self):
        # GH#34860
        arr = date_range("1/1/2008", "1/1/2009")
        ser = arr.to_series()
        result = ser["2008"]

        rng = date_range(start="2008-01-01", end="2008-12-31")
        expected = Series(rng, index=rng)

        tm.assert_series_equal(result, expected)

    def test_getitem_slice_strings_with_datetimeindex(self):
        idx = DatetimeIndex(
            ["1/1/2000", "1/2/2000", "1/2/2000", "1/3/2000", "1/4/2000"]
        )

        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts["1/2/2000":]
        expected = ts[1:]
        tm.assert_series_equal(result, expected)

        result = ts["1/2/2000":"1/3/2000"]
        expected = ts[1:4]
        tm.assert_series_equal(result, expected)

    def test_getitem_slice_2d(self, datetime_series):
        # GH#30588 multi-dimensional indexing deprecated

        with tm.assert_produces_warning(FutureWarning):
            # GH#30867 Don't want to support this long-term, but
            # for now ensure that the warning from Index
            # doesn't comes through via Series.__getitem__.
            result = datetime_series[:, np.newaxis]
        expected = datetime_series.values[:, np.newaxis]
        tm.assert_almost_equal(result, expected)

    # FutureWarning from NumPy.
    @pytest.mark.filterwarnings("ignore:Using a non-tuple:FutureWarning")
    def test_getitem_median_slice_bug(self):
        index = date_range("20090415", "20090519", freq="2B")
        s = Series(np.random.randn(13), index=index)

        indexer = [slice(6, 7, None)]
        with tm.assert_produces_warning(FutureWarning):
            # GH#31299
            result = s[indexer]
        expected = s[indexer[0]]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "slc, positions",
        [
            [slice(date(2018, 1, 1), None), [0, 1, 2]],
            [slice(date(2019, 1, 2), None), [2]],
            [slice(date(2020, 1, 1), None), []],
            [slice(None, date(2020, 1, 1)), [0, 1, 2]],
            [slice(None, date(2019, 1, 1)), [0]],
        ],
    )
    def test_getitem_slice_date(self, slc, positions):
        # https://github.com/pandas-dev/pandas/issues/31501
        ser = Series(
            [0, 1, 2],
            DatetimeIndex(["2019-01-01", "2019-01-01T06:00:00", "2019-01-02"]),
        )
        result = ser[slc]
        expected = ser.take(positions)
        tm.assert_series_equal(result, expected)


class TestSeriesGetitemListLike:
    @pytest.mark.parametrize("box", [list, np.array, Index, pd.Series])
    def test_getitem_no_matches(self, box):
        # GH#33462 we expect the same behavior for list/ndarray/Index/Series
        ser = Series(["A", "B"])

        key = Series(["C"], dtype=object)
        key = box(key)

        msg = r"None of \[Index\(\['C'\], dtype='object'\)\] are in the \[index\]"
        with pytest.raises(KeyError, match=msg):
            ser[key]

    def test_getitem_intlist_intindex_periodvalues(self):
        ser = Series(period_range("2000-01-01", periods=10, freq="D"))

        result = ser[[2, 4]]
        exp = Series(
            [pd.Period("2000-01-03", freq="D"), pd.Period("2000-01-05", freq="D")],
            index=[2, 4],
            dtype="Period[D]",
        )
        tm.assert_series_equal(result, exp)
        assert result.dtype == "Period[D]"

    @pytest.mark.parametrize("box", [list, np.array, Index])
    def test_getitem_intlist_intervalindex_non_int(self, box):
        # GH#33404 fall back to positional since ints are unambiguous
        dti = date_range("2000-01-03", periods=3)._with_freq(None)
        ii = pd.IntervalIndex.from_breaks(dti)
        ser = Series(range(len(ii)), index=ii)

        expected = ser.iloc[:1]
        key = box([0])
        result = ser[key]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("box", [list, np.array, Index])
    @pytest.mark.parametrize("dtype", [np.int64, np.float64, np.uint64])
    def test_getitem_intlist_multiindex_numeric_level(self, dtype, box):
        # GH#33404 do _not_ fall back to positional since ints are ambiguous
        idx = Index(range(4)).astype(dtype)
        dti = date_range("2000-01-03", periods=3)
        mi = pd.MultiIndex.from_product([idx, dti])
        ser = Series(range(len(mi))[::-1], index=mi)

        key = box([5])
        with pytest.raises(KeyError, match="5"):
            ser[key]

    def test_getitem_uint_array_key(self, uint_dtype):
        # GH #37218
        ser = Series([1, 2, 3])
        key = np.array([4], dtype=uint_dtype)

        with pytest.raises(KeyError, match="4"):
            ser[key]
        with pytest.raises(KeyError, match="4"):
            ser.loc[key]


class TestGetitemBooleanMask:
    def test_getitem_boolean(self, string_series):
        ser = string_series
        mask = ser > ser.median()

        # passing list is OK
        result = ser[list(mask)]
        expected = ser[mask]
        tm.assert_series_equal(result, expected)
        tm.assert_index_equal(result.index, ser.index[mask])

    def test_getitem_boolean_empty(self):
        ser = Series([], dtype=np.int64)
        ser.index.name = "index_name"
        ser = ser[ser.isna()]
        assert ser.index.name == "index_name"
        assert ser.dtype == np.int64

        # GH#5877
        # indexing with empty series
        ser = Series(["A", "B"])
        expected = Series(dtype=object, index=Index([], dtype="int64"))
        result = ser[Series([], dtype=object)]
        tm.assert_series_equal(result, expected)

        # invalid because of the boolean indexer
        # that's empty or not-aligned
        msg = (
            r"Unalignable boolean Series provided as indexer \(index of "
            r"the boolean Series and of the indexed object do not match"
        )
        with pytest.raises(IndexingError, match=msg):
            ser[Series([], dtype=bool)]

        with pytest.raises(IndexingError, match=msg):
            ser[Series([True], dtype=bool)]

    def test_getitem_boolean_object(self, string_series):
        # using column from DataFrame

        ser = string_series
        mask = ser > ser.median()
        omask = mask.astype(object)

        # getitem
        result = ser[omask]
        expected = ser[mask]
        tm.assert_series_equal(result, expected)

        # setitem
        s2 = ser.copy()
        cop = ser.copy()
        cop[omask] = 5
        s2[mask] = 5
        tm.assert_series_equal(cop, s2)

        # nans raise exception
        omask[5:10] = np.nan
        msg = "Cannot mask with non-boolean array containing NA / NaN values"
        with pytest.raises(ValueError, match=msg):
            ser[omask]
        with pytest.raises(ValueError, match=msg):
            ser[omask] = 5

    def test_getitem_boolean_dt64_copies(self):
        # GH#36210
        dti = date_range("2016-01-01", periods=4, tz="US/Pacific")
        key = np.array([True, True, False, False])

        ser = Series(dti._data)

        res = ser[key]
        assert res._values._data.base is None

        # compare with numeric case for reference
        ser2 = Series(range(4))
        res2 = ser2[key]
        assert res2._values.base is None

    def test_getitem_boolean_corner(self, datetime_series):
        ts = datetime_series
        mask_shifted = ts.shift(1, freq=BDay()) > ts.median()

        msg = (
            r"Unalignable boolean Series provided as indexer \(index of "
            r"the boolean Series and of the indexed object do not match"
        )
        with pytest.raises(IndexingError, match=msg):
            ts[mask_shifted]

        with pytest.raises(IndexingError, match=msg):
            ts.loc[mask_shifted]

    def test_getitem_boolean_different_order(self, string_series):
        ordered = string_series.sort_values()

        sel = string_series[ordered > 0]
        exp = string_series[string_series > 0]
        tm.assert_series_equal(sel, exp)

    def test_getitem_boolean_contiguous_preserve_freq(self):
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


class TestGetitemCallable:
    def test_getitem_callable(self):
        # GH#12533
        ser = Series(4, index=list("ABCD"))
        result = ser[lambda x: "A"]
        assert result == ser.loc["A"]

        result = ser[lambda x: ["A", "B"]]
        expected = ser.loc[["A", "B"]]
        tm.assert_series_equal(result, expected)

        result = ser[lambda x: [True, False, True, True]]
        expected = ser.iloc[[0, 2, 3]]
        tm.assert_series_equal(result, expected)


def test_getitem_generator(string_series):
    gen = (x > 0 for x in string_series)
    result = string_series[gen]
    result2 = string_series[iter(string_series > 0)]
    expected = string_series[string_series > 0]
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result2, expected)


@pytest.mark.parametrize(
    "series",
    [
        Series([0, 1]),
        Series(date_range("2012-01-01", periods=2)),
        Series(date_range("2012-01-01", periods=2, tz="CET")),
    ],
)
def test_getitem_ndim_deprecated(series):
    with tm.assert_produces_warning(
        FutureWarning, match="Support for multi-dimensional indexing"
    ):
        result = series[:, None]

    expected = np.asarray(series)[:, None]
    tm.assert_numpy_array_equal(result, expected)


def test_getitem_multilevel_scalar_slice_not_implemented(
    multiindex_year_month_day_dataframe_random_data,
):
    # not implementing this for now
    df = multiindex_year_month_day_dataframe_random_data
    ser = df["A"]

    msg = r"\(2000, slice\(3, 4, None\)\)"
    with pytest.raises(TypeError, match=msg):
        ser[2000, 3:4]


def test_getitem_dataframe_raises():
    rng = list(range(10))
    ser = Series(10, index=rng)
    df = DataFrame(rng, index=rng)
    msg = (
        "Indexing a Series with DataFrame is not supported, "
        "use the appropriate DataFrame column"
    )
    with pytest.raises(TypeError, match=msg):
        ser[df > 5]


def test_getitem_assignment_series_aligment():
    # https://github.com/pandas-dev/pandas/issues/37427
    # with getitem, when assigning with a Series, it is not first aligned
    ser = Series(range(10))
    idx = np.array([2, 4, 9])
    ser[idx] = Series([10, 11, 12])
    expected = Series([0, 1, 10, 3, 11, 5, 6, 7, 8, 12])
    tm.assert_series_equal(ser, expected)


def test_getitem_duplicate_index_mistyped_key_raises_keyerror():
    # GH#29189 float_index.get_loc(None) should raise KeyError, not TypeError
    ser = Series([2, 5, 6, 8], index=[2.0, 4.0, 4.0, 5.0])
    with pytest.raises(KeyError, match="None"):
        ser[None]

    with pytest.raises(KeyError, match="None"):
        ser.index.get_loc(None)

    with pytest.raises(KeyError, match="None"):
        ser.index._engine.get_loc(None)


def test_getitem_1tuple_slice_without_multiindex():
    ser = Series(range(5))
    key = (slice(3),)

    result = ser[key]
    expected = ser[key[0]]
    tm.assert_series_equal(result, expected)
