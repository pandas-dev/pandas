from datetime import datetime
import inspect

import numpy as np
import pytest

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_interval_dtype,
    is_object_dtype,
)

from pandas import (
    Categorical,
    DataFrame,
    DatetimeIndex,
    Index,
    IntervalIndex,
    Series,
    Timestamp,
    cut,
    date_range,
    to_datetime,
)
import pandas._testing as tm


class TestDataFrameAlterAxes:
    def test_set_index_directly(self, float_string_frame):
        df = float_string_frame
        idx = Index(np.arange(len(df))[::-1])

        df.index = idx
        tm.assert_index_equal(df.index, idx)
        with pytest.raises(ValueError, match="Length mismatch"):
            df.index = idx[::2]

    def test_convert_dti_to_series(self):
        # don't cast a DatetimeIndex WITH a tz, leave as object
        # GH 6032
        idx = DatetimeIndex(
            to_datetime(["2013-1-1 13:00", "2013-1-2 14:00"]), name="B"
        ).tz_localize("US/Pacific")
        df = DataFrame(np.random.randn(2, 1), columns=["A"])

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

        # convert index to series
        result = Series(idx)
        tm.assert_series_equal(result, expected)

        # assign to frame
        df["B"] = idx
        result = df["B"]
        tm.assert_series_equal(result, expected)

        # convert to series while keeping the timezone
        msg = "stop passing 'keep_tz'"
        with tm.assert_produces_warning(FutureWarning) as m:
            result = idx.to_series(keep_tz=True, index=[0, 1])
        tm.assert_series_equal(result, expected)
        assert msg in str(m[0].message)

        # convert to utc
        with tm.assert_produces_warning(FutureWarning) as m:
            df["B"] = idx.to_series(keep_tz=False, index=[0, 1])
        result = df["B"]
        comp = Series(DatetimeIndex(expected.values).tz_localize(None), name="B")
        tm.assert_series_equal(result, comp)
        msg = "do 'idx.tz_convert(None)' before calling"
        assert msg in str(m[0].message)

        result = idx.to_series(index=[0, 1])
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = idx.to_series(keep_tz=False, index=[0, 1])
        tm.assert_series_equal(result, expected.dt.tz_convert(None))
        msg = "do 'idx.tz_convert(None)' before calling"
        assert msg in str(m[0].message)

        # list of datetimes with a tz
        df["B"] = idx.to_pydatetime()
        result = df["B"]
        tm.assert_series_equal(result, expected)

        # GH 6785
        # set the index manually
        import pytz

        df = DataFrame([{"ts": datetime(2014, 4, 1, tzinfo=pytz.utc), "foo": 1}])
        expected = df.set_index("ts")
        df.index = df["ts"]
        df.pop("ts")
        tm.assert_frame_equal(df, expected)

    def test_set_columns(self, float_string_frame):
        cols = Index(np.arange(len(float_string_frame.columns)))
        float_string_frame.columns = cols
        with pytest.raises(ValueError, match="Length mismatch"):
            float_string_frame.columns = cols[::2]

    def test_dti_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = date_range("2011/01/01", periods=6, freq="M", tz="US/Eastern")
        idx2 = date_range("2013", periods=6, freq="A", tz="Asia/Tokyo")

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.reindex(idx2)
        tm.assert_index_equal(df.index, idx2)

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

    def test_reindex_api_equivalence(self):
        # equivalence of the labels/axis and index/columns API's
        df = DataFrame(
            [[1, 2, 3], [3, 4, 5], [5, 6, 7]],
            index=["a", "b", "c"],
            columns=["d", "e", "f"],
        )

        res1 = df.reindex(["b", "a"])
        res2 = df.reindex(index=["b", "a"])
        res3 = df.reindex(labels=["b", "a"])
        res4 = df.reindex(labels=["b", "a"], axis=0)
        res5 = df.reindex(["b", "a"], axis=0)
        for res in [res2, res3, res4, res5]:
            tm.assert_frame_equal(res1, res)

        res1 = df.reindex(columns=["e", "d"])
        res2 = df.reindex(["e", "d"], axis=1)
        res3 = df.reindex(labels=["e", "d"], axis=1)
        for res in [res2, res3]:
            tm.assert_frame_equal(res1, res)

        res1 = df.reindex(index=["b", "a"], columns=["e", "d"])
        res2 = df.reindex(columns=["e", "d"], index=["b", "a"])
        res3 = df.reindex(labels=["b", "a"], axis=0).reindex(labels=["e", "d"], axis=1)
        for res in [res2, res3]:
            tm.assert_frame_equal(res1, res)

    def test_assign_columns(self, float_frame):
        float_frame["hi"] = "there"

        df = float_frame.copy()
        df.columns = ["foo", "bar", "baz", "quux", "foo2"]
        tm.assert_series_equal(float_frame["C"], df["baz"], check_names=False)
        tm.assert_series_equal(float_frame["hi"], df["foo2"], check_names=False)

    def test_set_index_preserve_categorical_dtype(self):
        # GH13743, GH13854
        df = DataFrame(
            {
                "A": [1, 2, 1, 1, 2],
                "B": [10, 16, 22, 28, 34],
                "C1": Categorical(list("abaab"), categories=list("bac"), ordered=False),
                "C2": Categorical(list("abaab"), categories=list("bac"), ordered=True),
            }
        )
        for cols in ["C1", "C2", ["A", "C1"], ["A", "C2"], ["C1", "C2"]]:
            result = df.set_index(cols).reset_index()
            result = result.reindex(columns=df.columns)
            tm.assert_frame_equal(result, df)

    def test_rename_signature(self):
        sig = inspect.signature(DataFrame.rename)
        parameters = set(sig.parameters)
        assert parameters == {
            "self",
            "mapper",
            "index",
            "columns",
            "axis",
            "inplace",
            "copy",
            "level",
            "errors",
        }

    def test_reindex_signature(self):
        sig = inspect.signature(DataFrame.reindex)
        parameters = set(sig.parameters)
        assert parameters == {
            "self",
            "labels",
            "index",
            "columns",
            "axis",
            "limit",
            "copy",
            "level",
            "method",
            "fill_value",
            "tolerance",
        }


class TestIntervalIndex:
    def test_setitem(self):

        df = DataFrame({"A": range(10)})
        s = cut(df.A, 5)
        assert isinstance(s.cat.categories, IntervalIndex)

        # B & D end up as Categoricals
        # the remainer are converted to in-line objects
        # contining an IntervalIndex.values
        df["B"] = s
        df["C"] = np.array(s)
        df["D"] = s.values
        df["E"] = np.array(s.values)

        assert is_categorical_dtype(df["B"].dtype)
        assert is_interval_dtype(df["B"].cat.categories)
        assert is_categorical_dtype(df["D"].dtype)
        assert is_interval_dtype(df["D"].cat.categories)

        assert is_object_dtype(df["C"])
        assert is_object_dtype(df["E"])

        # they compare equal as Index
        # when converted to numpy objects
        c = lambda x: Index(np.array(x))
        tm.assert_index_equal(c(df.B), c(df.B), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.C), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.D), check_names=False)
        tm.assert_index_equal(c(df.B), c(df.D), check_names=False)

        # B & D are the same Series
        tm.assert_series_equal(df["B"], df["B"], check_names=False)
        tm.assert_series_equal(df["B"], df["D"], check_names=False)

        # C & E are the same Series
        tm.assert_series_equal(df["C"], df["C"], check_names=False)
        tm.assert_series_equal(df["C"], df["E"], check_names=False)

    def test_set_reset_index(self):

        df = DataFrame({"A": range(10)})
        s = cut(df.A, 5)
        df["B"] = s
        df = df.set_index("B")

        df = df.reset_index()
