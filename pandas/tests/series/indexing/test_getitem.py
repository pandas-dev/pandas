"""
Series.__getitem__ test classes are organized by the type of key passed.
"""
from datetime import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import conversion, timezones

import pandas as pd
from pandas import Series, Timestamp, date_range, period_range
import pandas._testing as tm


class TestSeriesGetitemScalars:

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


class TestSeriesGetitemSlices:
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


class TestSeriesGetitemListLike:
    @pytest.mark.parametrize("box", [list, np.array, pd.Index, pd.Series])
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
        exp = pd.Series(
            [pd.Period("2000-01-03", freq="D"), pd.Period("2000-01-05", freq="D")],
            index=[2, 4],
            dtype="Period[D]",
        )
        tm.assert_series_equal(result, exp)
        assert result.dtype == "Period[D]"

    @pytest.mark.parametrize("box", [list, np.array, pd.Index])
    def test_getitem_intlist_intervalindex_non_int(self, box):
        # GH#33404 fall back to positional since ints are unambiguous
        dti = date_range("2000-01-03", periods=3)
        ii = pd.IntervalIndex.from_breaks(dti)
        ser = Series(range(len(ii)), index=ii)

        expected = ser.iloc[:1]
        key = box([0])
        result = ser[key]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("box", [list, np.array, pd.Index])
    @pytest.mark.parametrize("dtype", [np.int64, np.float64, np.uint64])
    def test_getitem_intlist_multiindex_numeric_level(self, dtype, box):
        # GH#33404 do _not_ fall back to positional since ints are ambiguous
        idx = pd.Index(range(4)).astype(dtype)
        dti = date_range("2000-01-03", periods=3)
        mi = pd.MultiIndex.from_product([idx, dti])
        ser = Series(range(len(mi))[::-1], index=mi)

        key = box([5])
        with pytest.raises(KeyError, match="5"):
            ser[key]


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
    with tm.assert_produces_warning(FutureWarning):
        result = series[:, None]

    expected = np.asarray(series)[:, None]
    tm.assert_numpy_array_equal(result, expected)


def test_getitem_assignment_series_aligment():
    # https://github.com/pandas-dev/pandas/issues/37427
    # with getitem, when assigning with a Series, it is not first aligned
    s = Series(range(10))
    idx = np.array([2, 4, 9])
    s[idx] = Series([10, 11, 12])
    expected = Series([0, 1, 10, 3, 11, 5, 6, 7, 8, 12])
    tm.assert_series_equal(s, expected)
