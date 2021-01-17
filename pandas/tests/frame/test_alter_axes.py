from datetime import datetime

import numpy as np
import pytest
import pytz

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_interval_dtype,
    is_object_dtype,
)

from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    IntervalIndex,
    Series,
    Timestamp,
    cut,
    date_range,
)
import pandas._testing as tm


class TestDataFrameAlterAxes:
    @pytest.fixture
    def idx_expected(self):
        idx = DatetimeIndex(["2013-1-1 13:00", "2013-1-2 14:00"], name="B").tz_localize(
            "US/Pacific"
        )

        expected = Series(
            np.array(
                [
                    Timestamp("2013-01-01 13:00:00-0800", tz="US/Pacific"),
                    Timestamp("2013-01-02 14:00:00-0800", tz="US/Pacific"),
                ],
                dtype="object",
            ),
            name="B",
        )
        assert expected.dtype == idx.dtype
        return idx, expected

    def test_to_series_keep_tz_deprecated_true(self, idx_expected):
        # convert to series while keeping the timezone
        idx, expected = idx_expected

        msg = "stop passing 'keep_tz'"
        with tm.assert_produces_warning(FutureWarning) as m:
            result = idx.to_series(keep_tz=True, index=[0, 1])
        assert msg in str(m[0].message)

        tm.assert_series_equal(result, expected)

    def test_to_series_keep_tz_deprecated_false(self, idx_expected):
        idx, expected = idx_expected

        with tm.assert_produces_warning(FutureWarning) as m:
            result = idx.to_series(keep_tz=False, index=[0, 1])
        tm.assert_series_equal(result, expected.dt.tz_convert(None))
        msg = "do 'idx.tz_convert(None)' before calling"
        assert msg in str(m[0].message)

    def test_setitem_dt64series(self, idx_expected):
        # convert to utc
        idx, expected = idx_expected
        df = DataFrame(np.random.randn(2, 1), columns=["A"])
        df["B"] = idx

        with tm.assert_produces_warning(FutureWarning) as m:
            df["B"] = idx.to_series(keep_tz=False, index=[0, 1])
        msg = "do 'idx.tz_convert(None)' before calling"
        assert msg in str(m[0].message)

        result = df["B"]
        comp = Series(idx.tz_convert("UTC").tz_localize(None), name="B")
        tm.assert_series_equal(result, comp)

    def test_setitem_datetimeindex(self, idx_expected):
        # setting a DataFrame column with a tzaware DTI retains the dtype
        idx, expected = idx_expected
        df = DataFrame(np.random.randn(2, 1), columns=["A"])

        # assign to frame
        df["B"] = idx
        result = df["B"]
        tm.assert_series_equal(result, expected)

    def test_setitem_object_array_of_tzaware_datetimes(self, idx_expected):
        # setting a DataFrame column with a tzaware DTI retains the dtype
        idx, expected = idx_expected
        df = DataFrame(np.random.randn(2, 1), columns=["A"])

        # object array of datetimes with a tz
        df["B"] = idx.to_pydatetime()
        result = df["B"]
        tm.assert_series_equal(result, expected)

    def test_constructor_from_tzaware_datetimeindex(self, idx_expected):
        # don't cast a DatetimeIndex WITH a tz, leave as object
        # GH 6032
        idx, expected = idx_expected

        # convert index to series
        result = Series(idx)
        tm.assert_series_equal(result, expected)

    def test_set_axis_setattr_index(self):
        # GH 6785
        # set the index manually

        df = DataFrame([{"ts": datetime(2014, 4, 1, tzinfo=pytz.utc), "foo": 1}])
        expected = df.set_index("ts")
        df.index = df["ts"]
        df.pop("ts")
        tm.assert_frame_equal(df, expected)

    def test_dti_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = date_range("2011/01/01", periods=6, freq="M", tz="US/Eastern")
        idx2 = date_range("2013", periods=6, freq="A", tz="Asia/Tokyo")

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.reindex(idx2)
        tm.assert_index_equal(df.index, idx2)

    def test_dti_set_index_reindex_with_tz(self):
        # GH 11314
        # with tz
        index = date_range(
            datetime(2015, 10, 1), datetime(2015, 10, 1, 23), freq="H", tz="US/Eastern"
        )
        df = DataFrame(np.random.randn(24, 1), columns=["a"], index=index)
        new_index = date_range(
            datetime(2015, 10, 2), datetime(2015, 10, 2, 23), freq="H", tz="US/Eastern"
        )

        result = df.set_index(new_index)
        assert result.index.freq == index.freq

    # Renaming

    def test_assign_columns(self, float_frame):
        float_frame["hi"] = "there"

        df = float_frame.copy()
        df.columns = ["foo", "bar", "baz", "quux", "foo2"]
        tm.assert_series_equal(float_frame["C"], df["baz"], check_names=False)
        tm.assert_series_equal(float_frame["hi"], df["foo2"], check_names=False)


class TestIntervalIndex:
    def test_setitem(self):

        df = DataFrame({"A": range(10)})
        ser = cut(df["A"], 5)
        assert isinstance(ser.cat.categories, IntervalIndex)

        # B & D end up as Categoricals
        # the remainer are converted to in-line objects
        # contining an IntervalIndex.values
        df["B"] = ser
        df["C"] = np.array(ser)
        df["D"] = ser.values
        df["E"] = np.array(ser.values)

        assert is_categorical_dtype(df["B"].dtype)
        assert is_interval_dtype(df["B"].cat.categories)
        assert is_categorical_dtype(df["D"].dtype)
        assert is_interval_dtype(df["D"].cat.categories)

        assert is_object_dtype(df["C"])
        assert is_object_dtype(df["E"])

        # they compare equal as Index
        # when converted to numpy objects
        c = lambda x: Index(np.array(x))
        tm.assert_index_equal(c(df.B), c(df.B))
        tm.assert_index_equal(c(df.B), c(df.C), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.D), check_names=False)
        tm.assert_index_equal(c(df.C), c(df.D), check_names=False)

        # B & D are the same Series
        tm.assert_series_equal(df["B"], df["B"])
        tm.assert_series_equal(df["B"], df["D"], check_names=False)

        # C & E are the same Series
        tm.assert_series_equal(df["C"], df["C"])
        tm.assert_series_equal(df["C"], df["E"], check_names=False)

    def test_set_reset_index(self):

        df = DataFrame({"A": range(10)})
        s = cut(df.A, 5)
        df["B"] = s
        df = df.set_index("B")

        df = df.reset_index()
