from collections import OrderedDict
from datetime import date, datetime, timedelta
import random
import re

import numpy as np
import pytest

from pandas.core.dtypes.common import is_categorical_dtype, is_object_dtype
from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    DatetimeIndex,
    Float64Index,
    Int64Index,
    IntervalIndex,
    MultiIndex,
    PeriodIndex,
    RangeIndex,
    Series,
    TimedeltaIndex,
    UInt64Index,
)
import pandas._testing as tm
from pandas.api.types import CategoricalDtype as CDT
from pandas.core.reshape.concat import concat
from pandas.core.reshape.merge import MergeError, merge

N = 50
NGROUPS = 8


def get_test_data(ngroups=NGROUPS, n=N):
    unique_groups = list(range(ngroups))
    arr = np.asarray(np.tile(unique_groups, n // ngroups))

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[: n - len(arr)])

    random.shuffle(arr)
    return arr


def get_series():
    return [
        pd.Series([1], dtype="int64"),
        pd.Series([1], dtype="Int64"),
        pd.Series([1.23]),
        pd.Series(["foo"]),
        pd.Series([True]),
        pd.Series([pd.Timestamp("2018-01-01")]),
        pd.Series([pd.Timestamp("2018-01-01", tz="US/Eastern")]),
    ]


def get_series_na():
    return [
        pd.Series([np.nan], dtype="Int64"),
        pd.Series([np.nan], dtype="float"),
        pd.Series([np.nan], dtype="object"),
        pd.Series([pd.NaT]),
    ]


@pytest.fixture(params=get_series(), ids=lambda x: x.dtype.name)
def series_of_dtype(request):
    """
    A parametrized fixture returning a variety of Series of different
    dtypes
    """
    return request.param


@pytest.fixture(params=get_series(), ids=lambda x: x.dtype.name)
def series_of_dtype2(request):
    """
    A duplicate of the series_of_dtype fixture, so that it can be used
    twice by a single function
    """
    return request.param


@pytest.fixture(params=get_series_na(), ids=lambda x: x.dtype.name)
def series_of_dtype_all_na(request):
    """
    A parametrized fixture returning a variety of Series with all NA
    values
    """
    return request.param


class TestMerge:
    def setup_method(self, method):
        # aggregate multiple columns
        self.df = DataFrame(
            {
                "key1": get_test_data(),
                "key2": get_test_data(),
                "data1": np.random.randn(N),
                "data2": np.random.randn(N),
            }
        )

        # exclude a couple keys for fun
        self.df = self.df[self.df["key2"] > 1]

        self.df2 = DataFrame(
            {
                "key1": get_test_data(n=N // 5),
                "key2": get_test_data(ngroups=NGROUPS // 2, n=N // 5),
                "value": np.random.randn(N // 5),
            }
        )

        self.left = DataFrame(
            {"key": ["a", "b", "c", "d", "e", "e", "a"], "v1": np.random.randn(7)}
        )
        self.right = DataFrame({"v2": np.random.randn(4)}, index=["d", "b", "c", "a"])

    def test_merge_inner_join_empty(self):
        # GH 15328
        df_empty = pd.DataFrame()
        df_a = pd.DataFrame({"a": [1, 2]}, index=[0, 1], dtype="int64")
        result = pd.merge(df_empty, df_a, left_index=True, right_index=True)
        expected = pd.DataFrame({"a": []}, index=[], dtype="int64")
        tm.assert_frame_equal(result, expected)

    def test_merge_common(self):
        joined = merge(self.df, self.df2)
        exp = merge(self.df, self.df2, on=["key1", "key2"])
        tm.assert_frame_equal(joined, exp)

    def test_merge_non_string_columns(self):
        # https://github.com/pandas-dev/pandas/issues/17962
        # Checks that method runs for non string column names
        left = pd.DataFrame(
            {0: [1, 0, 1, 0], 1: [0, 1, 0, 0], 2: [0, 0, 2, 0], 3: [1, 0, 0, 3]}
        )

        right = left.astype(float)
        expected = left
        result = pd.merge(left, right)
        tm.assert_frame_equal(expected, result)

    def test_merge_index_as_on_arg(self):
        # GH14355

        left = self.df.set_index("key1")
        right = self.df2.set_index("key1")
        result = merge(left, right, on="key1")
        expected = merge(self.df, self.df2, on="key1").set_index("key1")
        tm.assert_frame_equal(result, expected)

    def test_merge_index_singlekey_right_vs_left(self):
        left = DataFrame(
            {"key": ["a", "b", "c", "d", "e", "e", "a"], "v1": np.random.randn(7)}
        )
        right = DataFrame({"v2": np.random.randn(4)}, index=["d", "b", "c", "a"])

        merged1 = merge(
            left, right, left_on="key", right_index=True, how="left", sort=False
        )
        merged2 = merge(
            right, left, right_on="key", left_index=True, how="right", sort=False
        )
        tm.assert_frame_equal(merged1, merged2.loc[:, merged1.columns])

        merged1 = merge(
            left, right, left_on="key", right_index=True, how="left", sort=True
        )
        merged2 = merge(
            right, left, right_on="key", left_index=True, how="right", sort=True
        )
        tm.assert_frame_equal(merged1, merged2.loc[:, merged1.columns])

    def test_merge_index_singlekey_inner(self):
        left = DataFrame(
            {"key": ["a", "b", "c", "d", "e", "e", "a"], "v1": np.random.randn(7)}
        )
        right = DataFrame({"v2": np.random.randn(4)}, index=["d", "b", "c", "a"])

        # inner join
        result = merge(left, right, left_on="key", right_index=True, how="inner")
        expected = left.join(right, on="key").loc[result.index]
        tm.assert_frame_equal(result, expected)

        result = merge(right, left, right_on="key", left_index=True, how="inner")
        expected = left.join(right, on="key").loc[result.index]
        tm.assert_frame_equal(result, expected.loc[:, result.columns])

    def test_merge_misspecified(self):
        msg = "Must pass right_on or right_index=True"
        with pytest.raises(pd.errors.MergeError, match=msg):
            merge(self.left, self.right, left_index=True)
        msg = "Must pass left_on or left_index=True"
        with pytest.raises(pd.errors.MergeError, match=msg):
            merge(self.left, self.right, right_index=True)

        msg = (
            'Can only pass argument "on" OR "left_on" and "right_on", not '
            "a combination of both"
        )
        with pytest.raises(pd.errors.MergeError, match=msg):
            merge(self.left, self.left, left_on="key", on="key")

        msg = r"len\(right_on\) must equal len\(left_on\)"
        with pytest.raises(ValueError, match=msg):
            merge(self.df, self.df2, left_on=["key1"], right_on=["key1", "key2"])

    def test_index_and_on_parameters_confusion(self):
        msg = "right_index parameter must be of type bool, not <class 'list'>"
        with pytest.raises(ValueError, match=msg):
            merge(
                self.df,
                self.df2,
                how="left",
                left_index=False,
                right_index=["key1", "key2"],
            )
        msg = "left_index parameter must be of type bool, not <class 'list'>"
        with pytest.raises(ValueError, match=msg):
            merge(
                self.df,
                self.df2,
                how="left",
                left_index=["key1", "key2"],
                right_index=False,
            )
        with pytest.raises(ValueError, match=msg):
            merge(
                self.df,
                self.df2,
                how="left",
                left_index=["key1", "key2"],
                right_index=["key1", "key2"],
            )

    def test_merge_overlap(self):
        merged = merge(self.left, self.left, on="key")
        exp_len = (self.left["key"].value_counts() ** 2).sum()
        assert len(merged) == exp_len
        assert "v1_x" in merged
        assert "v1_y" in merged

    def test_merge_different_column_key_names(self):
        left = DataFrame({"lkey": ["foo", "bar", "baz", "foo"], "value": [1, 2, 3, 4]})
        right = DataFrame({"rkey": ["foo", "bar", "qux", "foo"], "value": [5, 6, 7, 8]})

        merged = left.merge(
            right, left_on="lkey", right_on="rkey", how="outer", sort=True
        )

        exp = pd.Series(["bar", "baz", "foo", "foo", "foo", "foo", np.nan], name="lkey")
        tm.assert_series_equal(merged["lkey"], exp)

        exp = pd.Series(["bar", np.nan, "foo", "foo", "foo", "foo", "qux"], name="rkey")
        tm.assert_series_equal(merged["rkey"], exp)

        exp = pd.Series([2, 3, 1, 1, 4, 4, np.nan], name="value_x")
        tm.assert_series_equal(merged["value_x"], exp)

        exp = pd.Series([6, np.nan, 5, 8, 5, 8, 7], name="value_y")
        tm.assert_series_equal(merged["value_y"], exp)

    def test_merge_copy(self):
        left = DataFrame({"a": 0, "b": 1}, index=range(10))
        right = DataFrame({"c": "foo", "d": "bar"}, index=range(10))

        merged = merge(left, right, left_index=True, right_index=True, copy=True)

        merged["a"] = 6
        assert (left["a"] == 0).all()

        merged["d"] = "peekaboo"
        assert (right["d"] == "bar").all()

    def test_merge_nocopy(self):
        left = DataFrame({"a": 0, "b": 1}, index=range(10))
        right = DataFrame({"c": "foo", "d": "bar"}, index=range(10))

        merged = merge(left, right, left_index=True, right_index=True, copy=False)

        merged["a"] = 6
        assert (left["a"] == 6).all()

        merged["d"] = "peekaboo"
        assert (right["d"] == "peekaboo").all()

    def test_intelligently_handle_join_key(self):
        # #733, be a bit more 1337 about not returning unconsolidated DataFrame

        left = DataFrame(
            {"key": [1, 1, 2, 2, 3], "value": list(range(5))}, columns=["value", "key"]
        )
        right = DataFrame({"key": [1, 1, 2, 3, 4, 5], "rvalue": list(range(6))})

        joined = merge(left, right, on="key", how="outer")
        expected = DataFrame(
            {
                "key": [1, 1, 1, 1, 2, 2, 3, 4, 5],
                "value": np.array([0, 0, 1, 1, 2, 3, 4, np.nan, np.nan]),
                "rvalue": [0, 1, 0, 1, 2, 2, 3, 4, 5],
            },
            columns=["value", "key", "rvalue"],
        )
        tm.assert_frame_equal(joined, expected)

    def test_merge_join_key_dtype_cast(self):
        # #8596

        df1 = DataFrame({"key": [1], "v1": [10]})
        df2 = DataFrame({"key": [2], "v1": [20]})
        df = merge(df1, df2, how="outer")
        assert df["key"].dtype == "int64"

        df1 = DataFrame({"key": [True], "v1": [1]})
        df2 = DataFrame({"key": [False], "v1": [0]})
        df = merge(df1, df2, how="outer")

        # GH13169
        # this really should be bool
        assert df["key"].dtype == "object"

        df1 = DataFrame({"val": [1]})
        df2 = DataFrame({"val": [2]})
        lkey = np.array([1])
        rkey = np.array([2])
        df = merge(df1, df2, left_on=lkey, right_on=rkey, how="outer")
        assert df["key_0"].dtype == "int64"

    def test_handle_join_key_pass_array(self):
        left = DataFrame(
            {"key": [1, 1, 2, 2, 3], "value": np.arange(5)}, columns=["value", "key"]
        )
        right = DataFrame({"rvalue": np.arange(6)})
        key = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on="key", right_on=key, how="outer")
        merged2 = merge(right, left, left_on=key, right_on="key", how="outer")

        tm.assert_series_equal(merged["key"], merged2["key"])
        assert merged["key"].notna().all()
        assert merged2["key"].notna().all()

        left = DataFrame({"value": np.arange(5)}, columns=["value"])
        right = DataFrame({"rvalue": np.arange(6)})
        lkey = np.array([1, 1, 2, 2, 3])
        rkey = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on=lkey, right_on=rkey, how="outer")
        tm.assert_series_equal(
            merged["key_0"], Series([1, 1, 1, 1, 2, 2, 3, 4, 5], name="key_0")
        )

        left = DataFrame({"value": np.arange(3)})
        right = DataFrame({"rvalue": np.arange(6)})

        key = np.array([0, 1, 1, 2, 2, 3], dtype=np.int64)
        merged = merge(left, right, left_index=True, right_on=key, how="outer")
        tm.assert_series_equal(merged["key_0"], Series(key, name="key_0"))

    def test_no_overlap_more_informative_error(self):
        dt = datetime.now()
        df1 = DataFrame({"x": ["a"]}, index=[dt])

        df2 = DataFrame({"y": ["b", "c"]}, index=[dt, dt])

        msg = (
            "No common columns to perform merge on. "
            f"Merge options: left_on={None}, right_on={None}, "
            f"left_index={False}, right_index={False}"
        )

        with pytest.raises(MergeError, match=msg):
            merge(df1, df2)

    def test_merge_non_unique_indexes(self):

        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        dt4 = datetime(2012, 5, 4)

        df1 = DataFrame({"x": ["a"]}, index=[dt])
        df2 = DataFrame({"y": ["b", "c"]}, index=[dt, dt])
        _check_merge(df1, df2)

        # Not monotonic
        df1 = DataFrame({"x": ["a", "b", "q"]}, index=[dt2, dt, dt4])
        df2 = DataFrame(
            {"y": ["c", "d", "e", "f", "g", "h"]}, index=[dt3, dt3, dt2, dt2, dt, dt]
        )
        _check_merge(df1, df2)

        df1 = DataFrame({"x": ["a", "b"]}, index=[dt, dt])
        df2 = DataFrame({"y": ["c", "d"]}, index=[dt, dt])
        _check_merge(df1, df2)

    def test_merge_non_unique_index_many_to_many(self):
        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        df1 = DataFrame({"x": ["a", "b", "c", "d"]}, index=[dt2, dt2, dt, dt])
        df2 = DataFrame(
            {"y": ["e", "f", "g", " h", "i"]}, index=[dt2, dt2, dt3, dt, dt]
        )
        _check_merge(df1, df2)

    def test_left_merge_empty_dataframe(self):
        left = DataFrame({"key": [1], "value": [2]})
        right = DataFrame({"key": []})

        result = merge(left, right, on="key", how="left")
        tm.assert_frame_equal(result, left)

        result = merge(right, left, on="key", how="right")
        tm.assert_frame_equal(result, left)

    @pytest.mark.parametrize(
        "kwarg",
        [
            dict(left_index=True, right_index=True),
            dict(left_index=True, right_on="x"),
            dict(left_on="a", right_index=True),
            dict(left_on="a", right_on="x"),
        ],
    )
    def test_merge_left_empty_right_empty(self, join_type, kwarg):
        # GH 10824
        left = pd.DataFrame(columns=["a", "b", "c"])
        right = pd.DataFrame(columns=["x", "y", "z"])

        exp_in = pd.DataFrame(
            columns=["a", "b", "c", "x", "y", "z"],
            index=pd.Index([], dtype=object),
            dtype=object,
        )

        result = pd.merge(left, right, how=join_type, **kwarg)
        tm.assert_frame_equal(result, exp_in)

    def test_merge_left_empty_right_notempty(self):
        # GH 10824
        left = pd.DataFrame(columns=["a", "b", "c"])
        right = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["x", "y", "z"])

        exp_out = pd.DataFrame(
            {
                "a": np.array([np.nan] * 3, dtype=object),
                "b": np.array([np.nan] * 3, dtype=object),
                "c": np.array([np.nan] * 3, dtype=object),
                "x": [1, 4, 7],
                "y": [2, 5, 8],
                "z": [3, 6, 9],
            },
            columns=["a", "b", "c", "x", "y", "z"],
        )
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how="inner", **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how="left", **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how="right", **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how="outer", **kwarg)
            tm.assert_frame_equal(result, exp)

        for kwarg in [
            dict(left_index=True, right_index=True),
            dict(left_index=True, right_on="x"),
        ]:
            check1(exp_in, kwarg)
            check2(exp_out, kwarg)

        kwarg = dict(left_on="a", right_index=True)
        check1(exp_in, kwarg)
        exp_out["a"] = [0, 1, 2]
        check2(exp_out, kwarg)

        kwarg = dict(left_on="a", right_on="x")
        check1(exp_in, kwarg)
        exp_out["a"] = np.array([np.nan] * 3, dtype=object)
        check2(exp_out, kwarg)

    def test_merge_left_notempty_right_empty(self):
        # GH 10824
        left = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["a", "b", "c"])
        right = pd.DataFrame(columns=["x", "y", "z"])

        exp_out = pd.DataFrame(
            {
                "a": [1, 4, 7],
                "b": [2, 5, 8],
                "c": [3, 6, 9],
                "x": np.array([np.nan] * 3, dtype=object),
                "y": np.array([np.nan] * 3, dtype=object),
                "z": np.array([np.nan] * 3, dtype=object),
            },
            columns=["a", "b", "c", "x", "y", "z"],
        )
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how="inner", **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how="right", **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how="left", **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how="outer", **kwarg)
            tm.assert_frame_equal(result, exp)

            for kwarg in [
                dict(left_index=True, right_index=True),
                dict(left_index=True, right_on="x"),
                dict(left_on="a", right_index=True),
                dict(left_on="a", right_on="x"),
            ]:
                check1(exp_in, kwarg)
                check2(exp_out, kwarg)

    def test_merge_empty_frame(self, series_of_dtype, series_of_dtype2):
        # GH 25183
        df = pd.DataFrame(
            {"key": series_of_dtype, "value": series_of_dtype2},
            columns=["key", "value"],
        )
        df_empty = df[:0]
        expected = pd.DataFrame(
            {
                "value_x": pd.Series(dtype=df.dtypes["value"]),
                "key": pd.Series(dtype=df.dtypes["key"]),
                "value_y": pd.Series(dtype=df.dtypes["value"]),
            },
            columns=["value_x", "key", "value_y"],
        )
        actual = df_empty.merge(df, on="key")
        tm.assert_frame_equal(actual, expected)

    def test_merge_all_na_column(self, series_of_dtype, series_of_dtype_all_na):
        # GH 25183
        df_left = pd.DataFrame(
            {"key": series_of_dtype, "value": series_of_dtype_all_na},
            columns=["key", "value"],
        )
        df_right = pd.DataFrame(
            {"key": series_of_dtype, "value": series_of_dtype_all_na},
            columns=["key", "value"],
        )
        expected = pd.DataFrame(
            {
                "key": series_of_dtype,
                "value_x": series_of_dtype_all_na,
                "value_y": series_of_dtype_all_na,
            },
            columns=["key", "value_x", "value_y"],
        )
        actual = df_left.merge(df_right, on="key")
        tm.assert_frame_equal(actual, expected)

    def test_merge_nosort(self):
        # GH#2098, TODO: anything to do?

        d = {
            "var1": np.random.randint(0, 10, size=10),
            "var2": np.random.randint(0, 10, size=10),
            "var3": [
                datetime(2012, 1, 12),
                datetime(2011, 2, 4),
                datetime(2010, 2, 3),
                datetime(2012, 1, 12),
                datetime(2011, 2, 4),
                datetime(2012, 4, 3),
                datetime(2012, 3, 4),
                datetime(2008, 5, 1),
                datetime(2010, 2, 3),
                datetime(2012, 2, 3),
            ],
        }
        df = DataFrame.from_dict(d)
        var3 = df.var3.unique()
        var3.sort()
        new = DataFrame.from_dict({"var3": var3, "var8": np.random.random(7)})

        result = df.merge(new, on="var3", sort=False)
        exp = merge(df, new, on="var3", sort=False)
        tm.assert_frame_equal(result, exp)

        assert (df.var3.unique() == result.var3.unique()).all()

    def test_merge_nan_right(self):
        df1 = DataFrame({"i1": [0, 1], "i2": [0, 1]})
        df2 = DataFrame({"i1": [0], "i3": [0]})
        result = df1.join(df2, on="i1", rsuffix="_")
        expected = (
            DataFrame(
                {
                    "i1": {0: 0.0, 1: 1},
                    "i2": {0: 0, 1: 1},
                    "i1_": {0: 0, 1: np.nan},
                    "i3": {0: 0.0, 1: np.nan},
                    None: {0: 0, 1: 0},
                }
            )
            .set_index(None)
            .reset_index()[["i1", "i2", "i1_", "i3"]]
        )
        tm.assert_frame_equal(result, expected, check_dtype=False)

        df1 = DataFrame({"i1": [0, 1], "i2": [0.5, 1.5]})
        df2 = DataFrame({"i1": [0], "i3": [0.7]})
        result = df1.join(df2, rsuffix="_", on="i1")
        expected = DataFrame(
            {
                "i1": {0: 0, 1: 1},
                "i1_": {0: 0.0, 1: np.nan},
                "i2": {0: 0.5, 1: 1.5},
                "i3": {0: 0.69999999999999996, 1: np.nan},
            }
        )[["i1", "i2", "i1_", "i3"]]
        tm.assert_frame_equal(result, expected)

    def test_merge_type(self):
        class NotADataFrame(DataFrame):
            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.df)
        result = nad.merge(self.df2, on="key1")

        assert isinstance(result, NotADataFrame)

    def test_join_append_timedeltas(self):
        # timedelta64 issues with join/merge
        # GH 5695

        d = {"d": datetime(2013, 11, 5, 5, 56), "t": timedelta(0, 22500)}
        df = DataFrame(columns=list("dt"))
        df = df.append(d, ignore_index=True)
        result = df.append(d, ignore_index=True)
        expected = DataFrame(
            {
                "d": [datetime(2013, 11, 5, 5, 56), datetime(2013, 11, 5, 5, 56)],
                "t": [timedelta(0, 22500), timedelta(0, 22500)],
            }
        )
        tm.assert_frame_equal(result, expected)

        td = np.timedelta64(300000000)
        lhs = DataFrame(Series([td, td], index=["A", "B"]))
        rhs = DataFrame(Series([td], index=["A"]))

        result = lhs.join(rhs, rsuffix="r", how="left")
        expected = DataFrame(
            {
                "0": Series([td, td], index=list("AB")),
                "0r": Series([td, pd.NaT], index=list("AB")),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_other_datetime_unit(self):
        # GH 13389
        df1 = pd.DataFrame({"entity_id": [101, 102]})
        s = pd.Series([None, None], index=[101, 102], name="days")

        for dtype in [
            "datetime64[D]",
            "datetime64[h]",
            "datetime64[m]",
            "datetime64[s]",
            "datetime64[ms]",
            "datetime64[us]",
            "datetime64[ns]",
        ]:

            df2 = s.astype(dtype).to_frame("days")
            # coerces to datetime64[ns], thus should not be affected
            assert df2["days"].dtype == "datetime64[ns]"

            result = df1.merge(df2, left_on="entity_id", right_index=True)

            exp = pd.DataFrame(
                {
                    "entity_id": [101, 102],
                    "days": np.array(["nat", "nat"], dtype="datetime64[ns]"),
                },
                columns=["entity_id", "days"],
            )
            tm.assert_frame_equal(result, exp)

    @pytest.mark.parametrize("unit", ["D", "h", "m", "s", "ms", "us", "ns"])
    def test_other_timedelta_unit(self, unit):
        # GH 13389
        df1 = pd.DataFrame({"entity_id": [101, 102]})
        s = pd.Series([None, None], index=[101, 102], name="days")

        dtype = f"m8[{unit}]"
        df2 = s.astype(dtype).to_frame("days")
        assert df2["days"].dtype == "m8[ns]"

        result = df1.merge(df2, left_on="entity_id", right_index=True)

        exp = pd.DataFrame(
            {"entity_id": [101, 102], "days": np.array(["nat", "nat"], dtype=dtype)},
            columns=["entity_id", "days"],
        )
        tm.assert_frame_equal(result, exp)

    def test_overlapping_columns_error_message(self):
        df = DataFrame({"key": [1, 2, 3], "v1": [4, 5, 6], "v2": [7, 8, 9]})
        df2 = DataFrame({"key": [1, 2, 3], "v1": [4, 5, 6], "v2": [7, 8, 9]})

        df.columns = ["key", "foo", "foo"]
        df2.columns = ["key", "bar", "bar"]
        expected = DataFrame(
            {
                "key": [1, 2, 3],
                "v1": [4, 5, 6],
                "v2": [7, 8, 9],
                "v3": [4, 5, 6],
                "v4": [7, 8, 9],
            }
        )
        expected.columns = ["key", "foo", "foo", "bar", "bar"]
        tm.assert_frame_equal(merge(df, df2), expected)

        # #2649, #10639
        df2.columns = ["key1", "foo", "foo"]
        msg = r"Data columns not unique: Index\(\['foo', 'foo'\], dtype='object'\)"
        with pytest.raises(MergeError, match=msg):
            merge(df, df2)

    def test_merge_on_datetime64tz(self):

        # GH11405
        left = pd.DataFrame(
            {
                "key": pd.date_range("20151010", periods=2, tz="US/Eastern"),
                "value": [1, 2],
            }
        )
        right = pd.DataFrame(
            {
                "key": pd.date_range("20151011", periods=3, tz="US/Eastern"),
                "value": [1, 2, 3],
            }
        )

        expected = DataFrame(
            {
                "key": pd.date_range("20151010", periods=4, tz="US/Eastern"),
                "value_x": [1, 2, np.nan, np.nan],
                "value_y": [np.nan, 1, 2, 3],
            }
        )
        result = pd.merge(left, right, on="key", how="outer")
        tm.assert_frame_equal(result, expected)

        left = pd.DataFrame(
            {
                "key": [1, 2],
                "value": pd.date_range("20151010", periods=2, tz="US/Eastern"),
            }
        )
        right = pd.DataFrame(
            {
                "key": [2, 3],
                "value": pd.date_range("20151011", periods=2, tz="US/Eastern"),
            }
        )
        expected = DataFrame(
            {
                "key": [1, 2, 3],
                "value_x": list(pd.date_range("20151010", periods=2, tz="US/Eastern"))
                + [pd.NaT],
                "value_y": [pd.NaT]
                + list(pd.date_range("20151011", periods=2, tz="US/Eastern")),
            }
        )
        result = pd.merge(left, right, on="key", how="outer")
        tm.assert_frame_equal(result, expected)
        assert result["value_x"].dtype == "datetime64[ns, US/Eastern]"
        assert result["value_y"].dtype == "datetime64[ns, US/Eastern]"

    def test_merge_on_datetime64tz_empty(self):
        # https://github.com/pandas-dev/pandas/issues/25014
        dtz = pd.DatetimeTZDtype(tz="UTC")
        right = pd.DataFrame(
            {
                "date": [pd.Timestamp("2018", tz=dtz.tz)],
                "value": [4.0],
                "date2": [pd.Timestamp("2019", tz=dtz.tz)],
            },
            columns=["date", "value", "date2"],
        )
        left = right[:0]
        result = left.merge(right, on="date")
        expected = pd.DataFrame(
            {
                "value_x": pd.Series(dtype=float),
                "date2_x": pd.Series(dtype=dtz),
                "date": pd.Series(dtype=dtz),
                "value_y": pd.Series(dtype=float),
                "date2_y": pd.Series(dtype=dtz),
            },
            columns=["value_x", "date2_x", "date", "value_y", "date2_y"],
        )
        tm.assert_frame_equal(result, expected)

    def test_merge_datetime64tz_with_dst_transition(self):
        # GH 18885
        df1 = pd.DataFrame(
            pd.date_range("2017-10-29 01:00", periods=4, freq="H", tz="Europe/Madrid"),
            columns=["date"],
        )
        df1["value"] = 1
        df2 = pd.DataFrame(
            {
                "date": pd.to_datetime(
                    [
                        "2017-10-29 03:00:00",
                        "2017-10-29 04:00:00",
                        "2017-10-29 05:00:00",
                    ]
                ),
                "value": 2,
            }
        )
        df2["date"] = df2["date"].dt.tz_localize("UTC").dt.tz_convert("Europe/Madrid")
        result = pd.merge(df1, df2, how="outer", on="date")
        expected = pd.DataFrame(
            {
                "date": pd.date_range(
                    "2017-10-29 01:00", periods=7, freq="H", tz="Europe/Madrid"
                ),
                "value_x": [1] * 4 + [np.nan] * 3,
                "value_y": [np.nan] * 4 + [2] * 3,
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_merge_non_unique_period_index(self):
        # GH #16871
        index = pd.period_range("2016-01-01", periods=16, freq="M")
        df = DataFrame(list(range(len(index))), index=index, columns=["pnum"])
        df2 = concat([df, df])
        result = df.merge(df2, left_index=True, right_index=True, how="inner")
        expected = DataFrame(
            np.tile(np.arange(16, dtype=np.int64).repeat(2).reshape(-1, 1), 2),
            columns=["pnum_x", "pnum_y"],
            index=df2.sort_index().index,
        )
        tm.assert_frame_equal(result, expected)

    def test_merge_on_periods(self):
        left = pd.DataFrame(
            {"key": pd.period_range("20151010", periods=2, freq="D"), "value": [1, 2]}
        )
        right = pd.DataFrame(
            {
                "key": pd.period_range("20151011", periods=3, freq="D"),
                "value": [1, 2, 3],
            }
        )

        expected = DataFrame(
            {
                "key": pd.period_range("20151010", periods=4, freq="D"),
                "value_x": [1, 2, np.nan, np.nan],
                "value_y": [np.nan, 1, 2, 3],
            }
        )
        result = pd.merge(left, right, on="key", how="outer")
        tm.assert_frame_equal(result, expected)

        left = pd.DataFrame(
            {"key": [1, 2], "value": pd.period_range("20151010", periods=2, freq="D")}
        )
        right = pd.DataFrame(
            {"key": [2, 3], "value": pd.period_range("20151011", periods=2, freq="D")}
        )

        exp_x = pd.period_range("20151010", periods=2, freq="D")
        exp_y = pd.period_range("20151011", periods=2, freq="D")
        expected = DataFrame(
            {
                "key": [1, 2, 3],
                "value_x": list(exp_x) + [pd.NaT],
                "value_y": [pd.NaT] + list(exp_y),
            }
        )
        result = pd.merge(left, right, on="key", how="outer")
        tm.assert_frame_equal(result, expected)
        assert result["value_x"].dtype == "Period[D]"
        assert result["value_y"].dtype == "Period[D]"

    def test_indicator(self):
        # PR #10054. xref #7412 and closes #8790.
        df1 = DataFrame(
            {"col1": [0, 1], "col_conflict": [1, 2], "col_left": ["a", "b"]}
        )
        df1_copy = df1.copy()

        df2 = DataFrame(
            {
                "col1": [1, 2, 3, 4, 5],
                "col_conflict": [1, 2, 3, 4, 5],
                "col_right": [2, 2, 2, 2, 2],
            }
        )
        df2_copy = df2.copy()

        df_result = DataFrame(
            {
                "col1": [0, 1, 2, 3, 4, 5],
                "col_conflict_x": [1, 2, np.nan, np.nan, np.nan, np.nan],
                "col_left": ["a", "b", np.nan, np.nan, np.nan, np.nan],
                "col_conflict_y": [np.nan, 1, 2, 3, 4, 5],
                "col_right": [np.nan, 2, 2, 2, 2, 2],
            }
        )
        df_result["_merge"] = Categorical(
            [
                "left_only",
                "both",
                "right_only",
                "right_only",
                "right_only",
                "right_only",
            ],
            categories=["left_only", "right_only", "both"],
        )

        df_result = df_result[
            [
                "col1",
                "col_conflict_x",
                "col_left",
                "col_conflict_y",
                "col_right",
                "_merge",
            ]
        ]

        test = merge(df1, df2, on="col1", how="outer", indicator=True)
        tm.assert_frame_equal(test, df_result)
        test = df1.merge(df2, on="col1", how="outer", indicator=True)
        tm.assert_frame_equal(test, df_result)

        # No side effects
        tm.assert_frame_equal(df1, df1_copy)
        tm.assert_frame_equal(df2, df2_copy)

        # Check with custom name
        df_result_custom_name = df_result
        df_result_custom_name = df_result_custom_name.rename(
            columns={"_merge": "custom_name"}
        )

        test_custom_name = merge(
            df1, df2, on="col1", how="outer", indicator="custom_name"
        )
        tm.assert_frame_equal(test_custom_name, df_result_custom_name)
        test_custom_name = df1.merge(
            df2, on="col1", how="outer", indicator="custom_name"
        )
        tm.assert_frame_equal(test_custom_name, df_result_custom_name)

        # Check only accepts strings and booleans
        msg = "indicator option can only accept boolean or string arguments"
        with pytest.raises(ValueError, match=msg):
            merge(df1, df2, on="col1", how="outer", indicator=5)
        with pytest.raises(ValueError, match=msg):
            df1.merge(df2, on="col1", how="outer", indicator=5)

        # Check result integrity

        test2 = merge(df1, df2, on="col1", how="left", indicator=True)
        assert (test2._merge != "right_only").all()
        test2 = df1.merge(df2, on="col1", how="left", indicator=True)
        assert (test2._merge != "right_only").all()

        test3 = merge(df1, df2, on="col1", how="right", indicator=True)
        assert (test3._merge != "left_only").all()
        test3 = df1.merge(df2, on="col1", how="right", indicator=True)
        assert (test3._merge != "left_only").all()

        test4 = merge(df1, df2, on="col1", how="inner", indicator=True)
        assert (test4._merge == "both").all()
        test4 = df1.merge(df2, on="col1", how="inner", indicator=True)
        assert (test4._merge == "both").all()

        # Check if working name in df
        for i in ["_right_indicator", "_left_indicator", "_merge"]:
            df_badcolumn = DataFrame({"col1": [1, 2], i: [2, 2]})

            msg = (
                "Cannot use `indicator=True` option when data contains a "
                f"column named {i}|"
                "Cannot use name of an existing column for indicator column"
            )
            with pytest.raises(ValueError, match=msg):
                merge(df1, df_badcolumn, on="col1", how="outer", indicator=True)
            with pytest.raises(ValueError, match=msg):
                df1.merge(df_badcolumn, on="col1", how="outer", indicator=True)

        # Check for name conflict with custom name
        df_badcolumn = DataFrame({"col1": [1, 2], "custom_column_name": [2, 2]})

        msg = "Cannot use name of an existing column for indicator column"
        with pytest.raises(ValueError, match=msg):
            merge(
                df1,
                df_badcolumn,
                on="col1",
                how="outer",
                indicator="custom_column_name",
            )
        with pytest.raises(ValueError, match=msg):
            df1.merge(
                df_badcolumn, on="col1", how="outer", indicator="custom_column_name"
            )

        # Merge on multiple columns
        df3 = DataFrame({"col1": [0, 1], "col2": ["a", "b"]})

        df4 = DataFrame({"col1": [1, 1, 3], "col2": ["b", "x", "y"]})

        hand_coded_result = DataFrame(
            {"col1": [0, 1, 1, 3], "col2": ["a", "b", "x", "y"]}
        )
        hand_coded_result["_merge"] = Categorical(
            ["left_only", "both", "right_only", "right_only"],
            categories=["left_only", "right_only", "both"],
        )

        test5 = merge(df3, df4, on=["col1", "col2"], how="outer", indicator=True)
        tm.assert_frame_equal(test5, hand_coded_result)
        test5 = df3.merge(df4, on=["col1", "col2"], how="outer", indicator=True)
        tm.assert_frame_equal(test5, hand_coded_result)

    def test_validation(self):
        left = DataFrame(
            {"a": ["a", "b", "c", "d"], "b": ["cat", "dog", "weasel", "horse"]},
            index=range(4),
        )

        right = DataFrame(
            {
                "a": ["a", "b", "c", "d", "e"],
                "c": ["meow", "bark", "um... weasel noise?", "nay", "chirp"],
            },
            index=range(5),
        )

        # Make sure no side effects.
        left_copy = left.copy()
        right_copy = right.copy()

        result = merge(left, right, left_index=True, right_index=True, validate="1:1")
        tm.assert_frame_equal(left, left_copy)
        tm.assert_frame_equal(right, right_copy)

        # make sure merge still correct
        expected = DataFrame(
            {
                "a_x": ["a", "b", "c", "d"],
                "b": ["cat", "dog", "weasel", "horse"],
                "a_y": ["a", "b", "c", "d"],
                "c": ["meow", "bark", "um... weasel noise?", "nay"],
            },
            index=range(4),
            columns=["a_x", "b", "a_y", "c"],
        )

        result = merge(
            left, right, left_index=True, right_index=True, validate="one_to_one"
        )
        tm.assert_frame_equal(result, expected)

        expected_2 = DataFrame(
            {
                "a": ["a", "b", "c", "d"],
                "b": ["cat", "dog", "weasel", "horse"],
                "c": ["meow", "bark", "um... weasel noise?", "nay"],
            },
            index=range(4),
        )

        result = merge(left, right, on="a", validate="1:1")
        tm.assert_frame_equal(left, left_copy)
        tm.assert_frame_equal(right, right_copy)
        tm.assert_frame_equal(result, expected_2)

        result = merge(left, right, on="a", validate="one_to_one")
        tm.assert_frame_equal(result, expected_2)

        # One index, one column
        expected_3 = DataFrame(
            {
                "b": ["cat", "dog", "weasel", "horse"],
                "a": ["a", "b", "c", "d"],
                "c": ["meow", "bark", "um... weasel noise?", "nay"],
            },
            columns=["b", "a", "c"],
            index=range(4),
        )

        left_index_reset = left.set_index("a")
        result = merge(
            left_index_reset,
            right,
            left_index=True,
            right_on="a",
            validate="one_to_one",
        )
        tm.assert_frame_equal(result, expected_3)

        # Dups on right
        right_w_dups = right.append(pd.DataFrame({"a": ["e"], "c": ["moo"]}, index=[4]))
        merge(
            left,
            right_w_dups,
            left_index=True,
            right_index=True,
            validate="one_to_many",
        )

        msg = "Merge keys are not unique in right dataset; not a one-to-one merge"
        with pytest.raises(MergeError, match=msg):
            merge(
                left,
                right_w_dups,
                left_index=True,
                right_index=True,
                validate="one_to_one",
            )

        with pytest.raises(MergeError, match=msg):
            merge(left, right_w_dups, on="a", validate="one_to_one")

        # Dups on left
        left_w_dups = left.append(
            pd.DataFrame({"a": ["a"], "c": ["cow"]}, index=[3]), sort=True
        )
        merge(
            left_w_dups,
            right,
            left_index=True,
            right_index=True,
            validate="many_to_one",
        )

        msg = "Merge keys are not unique in left dataset; not a one-to-one merge"
        with pytest.raises(MergeError, match=msg):
            merge(
                left_w_dups,
                right,
                left_index=True,
                right_index=True,
                validate="one_to_one",
            )

        with pytest.raises(MergeError, match=msg):
            merge(left_w_dups, right, on="a", validate="one_to_one")

        # Dups on both
        merge(left_w_dups, right_w_dups, on="a", validate="many_to_many")

        msg = "Merge keys are not unique in right dataset; not a many-to-one merge"
        with pytest.raises(MergeError, match=msg):
            merge(
                left_w_dups,
                right_w_dups,
                left_index=True,
                right_index=True,
                validate="many_to_one",
            )

        msg = "Merge keys are not unique in left dataset; not a one-to-many merge"
        with pytest.raises(MergeError, match=msg):
            merge(left_w_dups, right_w_dups, on="a", validate="one_to_many")

        # Check invalid arguments
        msg = "Not a valid argument for validate"
        with pytest.raises(ValueError, match=msg):
            merge(left, right, on="a", validate="jibberish")

        # Two column merge, dups in both, but jointly no dups.
        left = DataFrame(
            {
                "a": ["a", "a", "b", "b"],
                "b": [0, 1, 0, 1],
                "c": ["cat", "dog", "weasel", "horse"],
            },
            index=range(4),
        )

        right = DataFrame(
            {
                "a": ["a", "a", "b"],
                "b": [0, 1, 0],
                "d": ["meow", "bark", "um... weasel noise?"],
            },
            index=range(3),
        )

        expected_multi = DataFrame(
            {
                "a": ["a", "a", "b"],
                "b": [0, 1, 0],
                "c": ["cat", "dog", "weasel"],
                "d": ["meow", "bark", "um... weasel noise?"],
            },
            index=range(3),
        )

        msg = (
            "Merge keys are not unique in either left or right dataset; "
            "not a one-to-one merge"
        )
        with pytest.raises(MergeError, match=msg):
            merge(left, right, on="a", validate="1:1")

        result = merge(left, right, on=["a", "b"], validate="1:1")
        tm.assert_frame_equal(result, expected_multi)

    def test_merge_two_empty_df_no_division_error(self):
        # GH17776, PR #17846
        a = pd.DataFrame({"a": [], "b": [], "c": []})
        with np.errstate(divide="raise"):
            merge(a, a, on=("a", "b"))

    @pytest.mark.parametrize("how", ["right", "outer"])
    @pytest.mark.parametrize(
        "index,expected_index",
        [
            (
                CategoricalIndex([1, 2, 4]),
                CategoricalIndex([1, 2, 4, None, None, None]),
            ),
            (
                DatetimeIndex(["2001-01-01", "2002-02-02", "2003-03-03"]),
                DatetimeIndex(
                    ["2001-01-01", "2002-02-02", "2003-03-03", pd.NaT, pd.NaT, pd.NaT]
                ),
            ),
            (Float64Index([1, 2, 3]), Float64Index([1, 2, 3, None, None, None])),
            (Int64Index([1, 2, 3]), Float64Index([1, 2, 3, None, None, None])),
            (
                IntervalIndex.from_tuples([(1, 2), (2, 3), (3, 4)]),
                IntervalIndex.from_tuples(
                    [(1, 2), (2, 3), (3, 4), np.nan, np.nan, np.nan]
                ),
            ),
            (
                PeriodIndex(["2001-01-01", "2001-01-02", "2001-01-03"], freq="D"),
                PeriodIndex(
                    ["2001-01-01", "2001-01-02", "2001-01-03", pd.NaT, pd.NaT, pd.NaT],
                    freq="D",
                ),
            ),
            (
                TimedeltaIndex(["1d", "2d", "3d"]),
                TimedeltaIndex(["1d", "2d", "3d", pd.NaT, pd.NaT, pd.NaT]),
            ),
        ],
    )
    def test_merge_on_index_with_more_values(self, how, index, expected_index):
        # GH 24212
        # pd.merge gets [0, 1, 2, -1, -1, -1] as left_indexer, ensure that
        # -1 is interpreted as a missing value instead of the last element
        df1 = pd.DataFrame({"a": [1, 2, 3], "key": [0, 2, 2]}, index=index)
        df2 = pd.DataFrame({"b": [1, 2, 3, 4, 5]})
        result = df1.merge(df2, left_on="key", right_index=True, how=how)
        expected = pd.DataFrame(
            [
                [1.0, 0, 1],
                [2.0, 2, 3],
                [3.0, 2, 3],
                [np.nan, 1, 2],
                [np.nan, 3, 4],
                [np.nan, 4, 5],
            ],
            columns=["a", "key", "b"],
        )
        expected.set_index(expected_index, inplace=True)
        tm.assert_frame_equal(result, expected)

    def test_merge_right_index_right(self):
        # Note: the expected output here is probably incorrect.
        # See https://github.com/pandas-dev/pandas/issues/17257 for more.
        # We include this as a regression test for GH-24897.
        left = pd.DataFrame({"a": [1, 2, 3], "key": [0, 1, 1]})
        right = pd.DataFrame({"b": [1, 2, 3]})

        expected = pd.DataFrame(
            {"a": [1, 2, 3, None], "key": [0, 1, 1, 2], "b": [1, 2, 2, 3]},
            columns=["a", "key", "b"],
            index=[0, 1, 2, np.nan],
        )
        result = left.merge(right, left_on="key", right_index=True, how="right")
        tm.assert_frame_equal(result, expected)

    def test_merge_take_missing_values_from_index_of_other_dtype(self):
        # GH 24212
        left = pd.DataFrame(
            {
                "a": [1, 2, 3],
                "key": pd.Categorical(["a", "a", "b"], categories=list("abc")),
            }
        )
        right = pd.DataFrame(
            {"b": [1, 2, 3]}, index=pd.CategoricalIndex(["a", "b", "c"])
        )
        result = left.merge(right, left_on="key", right_index=True, how="right")
        expected = pd.DataFrame(
            {
                "a": [1, 2, 3, None],
                "key": pd.Categorical(["a", "a", "b", "c"]),
                "b": [1, 1, 2, 3],
            },
            index=[0, 1, 2, np.nan],
        )
        expected = expected.reindex(columns=["a", "key", "b"])
        tm.assert_frame_equal(result, expected)

    def test_merge_readonly(self):
        # https://github.com/pandas-dev/pandas/issues/27943
        data1 = pd.DataFrame(
            np.arange(20).reshape((4, 5)) + 1, columns=["a", "b", "c", "d", "e"]
        )
        data2 = pd.DataFrame(
            np.arange(20).reshape((5, 4)) + 1, columns=["a", "b", "x", "y"]
        )

        data1._data.blocks[0].values.flags.writeable = False
        data1.merge(data2)  # no error


def _check_merge(x, y):
    for how in ["inner", "left", "outer"]:
        result = x.join(y, how=how)

        expected = merge(x.reset_index(), y.reset_index(), how=how, sort=True)
        expected = expected.set_index("index")

        # TODO check_names on merge?
        tm.assert_frame_equal(result, expected, check_names=False)


class TestMergeDtypes:
    @pytest.mark.parametrize(
        "right_vals", [["foo", "bar"], Series(["foo", "bar"]).astype("category")]
    )
    def test_different(self, right_vals):

        left = DataFrame(
            {
                "A": ["foo", "bar"],
                "B": Series(["foo", "bar"]).astype("category"),
                "C": [1, 2],
                "D": [1.0, 2.0],
                "E": Series([1, 2], dtype="uint64"),
                "F": Series([1, 2], dtype="int32"),
            }
        )
        right = DataFrame({"A": right_vals})

        # GH 9780
        # We allow merging on object and categorical cols and cast
        # categorical cols to object
        result = pd.merge(left, right, on="A")
        assert is_object_dtype(result.A.dtype)

    @pytest.mark.parametrize("d1", [np.int64, np.int32, np.int16, np.int8, np.uint8])
    @pytest.mark.parametrize("d2", [np.int64, np.float64, np.float32, np.float16])
    def test_join_multi_dtypes(self, d1, d2):

        dtype1 = np.dtype(d1)
        dtype2 = np.dtype(d2)

        left = DataFrame(
            {
                "k1": np.array([0, 1, 2] * 8, dtype=dtype1),
                "k2": ["foo", "bar"] * 12,
                "v": np.array(np.arange(24), dtype=np.int64),
            }
        )

        index = MultiIndex.from_tuples([(2, "bar"), (1, "foo")])
        right = DataFrame({"v2": np.array([5, 7], dtype=dtype2)}, index=index)

        result = left.join(right, on=["k1", "k2"])

        expected = left.copy()

        if dtype2.kind == "i":
            dtype2 = np.dtype("float64")
        expected["v2"] = np.array(np.nan, dtype=dtype2)
        expected.loc[(expected.k1 == 2) & (expected.k2 == "bar"), "v2"] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == "foo"), "v2"] = 7

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=["k1", "k2"], sort=True)
        expected.sort_values(["k1", "k2"], kind="mergesort", inplace=True)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "int_vals, float_vals, exp_vals",
        [
            ([1, 2, 3], [1.0, 2.0, 3.0], {"X": [1, 2, 3], "Y": [1.0, 2.0, 3.0]}),
            ([1, 2, 3], [1.0, 3.0], {"X": [1, 3], "Y": [1.0, 3.0]}),
            ([1, 2], [1.0, 2.0, 3.0], {"X": [1, 2], "Y": [1.0, 2.0]}),
        ],
    )
    def test_merge_on_ints_floats(self, int_vals, float_vals, exp_vals):
        # GH 16572
        # Check that float column is not cast to object if
        # merging on float and int columns
        A = DataFrame({"X": int_vals})
        B = DataFrame({"Y": float_vals})
        expected = DataFrame(exp_vals)

        result = A.merge(B, left_on="X", right_on="Y")
        tm.assert_frame_equal(result, expected)

        result = B.merge(A, left_on="Y", right_on="X")
        tm.assert_frame_equal(result, expected[["Y", "X"]])

    def test_merge_key_dtype_cast(self):
        # GH 17044
        df1 = DataFrame({"key": [1.0, 2.0], "v1": [10, 20]}, columns=["key", "v1"])
        df2 = DataFrame({"key": [2], "v2": [200]}, columns=["key", "v2"])
        result = df1.merge(df2, on="key", how="left")
        expected = DataFrame(
            {"key": [1.0, 2.0], "v1": [10, 20], "v2": [np.nan, 200.0]},
            columns=["key", "v1", "v2"],
        )
        tm.assert_frame_equal(result, expected)

    def test_merge_on_ints_floats_warning(self):
        # GH 16572
        # merge will produce a warning when merging on int and
        # float columns where the float values are not exactly
        # equal to their int representation
        A = DataFrame({"X": [1, 2, 3]})
        B = DataFrame({"Y": [1.1, 2.5, 3.0]})
        expected = DataFrame({"X": [3], "Y": [3.0]})

        with tm.assert_produces_warning(UserWarning):
            result = A.merge(B, left_on="X", right_on="Y")
            tm.assert_frame_equal(result, expected)

        with tm.assert_produces_warning(UserWarning):
            result = B.merge(A, left_on="Y", right_on="X")
            tm.assert_frame_equal(result, expected[["Y", "X"]])

        # test no warning if float has NaNs
        B = DataFrame({"Y": [np.nan, np.nan, 3.0]})

        with tm.assert_produces_warning(None):
            result = B.merge(A, left_on="Y", right_on="X")
            tm.assert_frame_equal(result, expected[["Y", "X"]])

    def test_merge_incompat_infer_boolean_object(self):
        # GH21119: bool + object bool merge OK
        df1 = DataFrame({"key": Series([True, False], dtype=object)})
        df2 = DataFrame({"key": [True, False]})

        expected = DataFrame({"key": [True, False]}, dtype=object)
        result = pd.merge(df1, df2, on="key")
        tm.assert_frame_equal(result, expected)
        result = pd.merge(df2, df1, on="key")
        tm.assert_frame_equal(result, expected)

        # with missing value
        df1 = DataFrame({"key": Series([True, False, np.nan], dtype=object)})
        df2 = DataFrame({"key": [True, False]})

        expected = DataFrame({"key": [True, False]}, dtype=object)
        result = pd.merge(df1, df2, on="key")
        tm.assert_frame_equal(result, expected)
        result = pd.merge(df2, df1, on="key")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "df1_vals, df2_vals",
        [
            # merge on category coerces to object
            ([0, 1, 2], Series(["a", "b", "a"]).astype("category")),
            ([0.0, 1.0, 2.0], Series(["a", "b", "a"]).astype("category")),
            # no not infer
            ([0, 1], pd.Series([False, True], dtype=object)),
            ([0, 1], pd.Series([False, True], dtype=bool)),
        ],
    )
    def test_merge_incompat_dtypes_are_ok(self, df1_vals, df2_vals):
        # these are explicitly allowed incompat merges, that pass thru
        # the result type is dependent on if the values on the rhs are
        # inferred, otherwise these will be coerced to object

        df1 = DataFrame({"A": df1_vals})
        df2 = DataFrame({"A": df2_vals})

        result = pd.merge(df1, df2, on=["A"])
        assert is_object_dtype(result.A.dtype)
        result = pd.merge(df2, df1, on=["A"])
        assert is_object_dtype(result.A.dtype)

    @pytest.mark.parametrize(
        "df1_vals, df2_vals",
        [
            # do not infer to numeric
            (Series([1, 2], dtype="uint64"), ["a", "b", "c"]),
            (Series([1, 2], dtype="int32"), ["a", "b", "c"]),
            ([0, 1, 2], ["0", "1", "2"]),
            ([0.0, 1.0, 2.0], ["0", "1", "2"]),
            ([0, 1, 2], ["0", "1", "2"]),
            (
                pd.date_range("1/1/2011", periods=2, freq="D"),
                ["2011-01-01", "2011-01-02"],
            ),
            (pd.date_range("1/1/2011", periods=2, freq="D"), [0, 1]),
            (pd.date_range("1/1/2011", periods=2, freq="D"), [0.0, 1.0]),
            (
                pd.date_range("20130101", periods=3),
                pd.date_range("20130101", periods=3, tz="US/Eastern"),
            ),
        ],
    )
    def test_merge_incompat_dtypes_error(self, df1_vals, df2_vals):
        # GH 9780, GH 15800
        # Raise a ValueError when a user tries to merge on
        # dtypes that are incompatible (e.g., obj and int/float)

        df1 = DataFrame({"A": df1_vals})
        df2 = DataFrame({"A": df2_vals})

        msg = (
            f"You are trying to merge on {df1['A'].dtype} and "
            f"{df2['A'].dtype} columns. If you wish to proceed "
            "you should use pd.concat"
        )
        msg = re.escape(msg)
        with pytest.raises(ValueError, match=msg):
            pd.merge(df1, df2, on=["A"])

        # Check that error still raised when swapping order of dataframes
        msg = (
            f"You are trying to merge on {df2['A'].dtype} and "
            f"{df1['A'].dtype} columns. If you wish to proceed "
            "you should use pd.concat"
        )
        msg = re.escape(msg)
        with pytest.raises(ValueError, match=msg):
            pd.merge(df2, df1, on=["A"])


@pytest.fixture
def left():
    np.random.seed(1234)
    return DataFrame(
        {
            "X": Series(np.random.choice(["foo", "bar"], size=(10,))).astype(
                CDT(["foo", "bar"])
            ),
            "Y": np.random.choice(["one", "two", "three"], size=(10,)),
        }
    )


@pytest.fixture
def right():
    np.random.seed(1234)
    return DataFrame(
        {"X": Series(["foo", "bar"]).astype(CDT(["foo", "bar"])), "Z": [1, 2]}
    )


class TestMergeCategorical:
    def test_identical(self, left):
        # merging on the same, should preserve dtypes
        merged = pd.merge(left, left, on="X")
        result = merged.dtypes.sort_index()
        expected = Series(
            [CategoricalDtype(), np.dtype("O"), np.dtype("O")],
            index=["X", "Y_x", "Y_y"],
        )
        tm.assert_series_equal(result, expected)

    def test_basic(self, left, right):
        # we have matching Categorical dtypes in X
        # so should preserve the merged column
        merged = pd.merge(left, right, on="X")
        result = merged.dtypes.sort_index()
        expected = Series(
            [CategoricalDtype(), np.dtype("O"), np.dtype("int64")],
            index=["X", "Y", "Z"],
        )
        tm.assert_series_equal(result, expected)

    def test_merge_categorical(self):
        # GH 9426

        right = DataFrame(
            {
                "c": {0: "a", 1: "b", 2: "c", 3: "d", 4: "e"},
                "d": {0: "null", 1: "null", 2: "null", 3: "null", 4: "null"},
            }
        )
        left = DataFrame(
            {
                "a": {0: "f", 1: "f", 2: "f", 3: "f", 4: "f"},
                "b": {0: "g", 1: "g", 2: "g", 3: "g", 4: "g"},
            }
        )
        df = pd.merge(left, right, how="left", left_on="b", right_on="c")

        # object-object
        expected = df.copy()

        # object-cat
        # note that we propagate the category
        # because we don't have any matching rows
        cright = right.copy()
        cright["d"] = cright["d"].astype("category")
        result = pd.merge(left, cright, how="left", left_on="b", right_on="c")
        expected["d"] = expected["d"].astype(CategoricalDtype(["null"]))
        tm.assert_frame_equal(result, expected)

        # cat-object
        cleft = left.copy()
        cleft["b"] = cleft["b"].astype("category")
        result = pd.merge(cleft, cright, how="left", left_on="b", right_on="c")
        tm.assert_frame_equal(result, expected)

        # cat-cat
        cright = right.copy()
        cright["d"] = cright["d"].astype("category")
        cleft = left.copy()
        cleft["b"] = cleft["b"].astype("category")
        result = pd.merge(cleft, cright, how="left", left_on="b", right_on="c")
        tm.assert_frame_equal(result, expected)

    def tests_merge_categorical_unordered_equal(self):
        # GH-19551
        df1 = DataFrame(
            {
                "Foo": Categorical(["A", "B", "C"], categories=["A", "B", "C"]),
                "Left": ["A0", "B0", "C0"],
            }
        )

        df2 = DataFrame(
            {
                "Foo": Categorical(["C", "B", "A"], categories=["C", "B", "A"]),
                "Right": ["C1", "B1", "A1"],
            }
        )
        result = pd.merge(df1, df2, on=["Foo"])
        expected = DataFrame(
            {
                "Foo": pd.Categorical(["A", "B", "C"]),
                "Left": ["A0", "B0", "C0"],
                "Right": ["A1", "B1", "C1"],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_other_columns(self, left, right):
        # non-merge columns should preserve if possible
        right = right.assign(Z=right.Z.astype("category"))

        merged = pd.merge(left, right, on="X")
        result = merged.dtypes.sort_index()
        expected = Series(
            [CategoricalDtype(), np.dtype("O"), CategoricalDtype()],
            index=["X", "Y", "Z"],
        )
        tm.assert_series_equal(result, expected)

        # categories are preserved
        assert left.X.values.is_dtype_equal(merged.X.values)
        assert right.Z.values.is_dtype_equal(merged.Z.values)

    @pytest.mark.parametrize(
        "change",
        [
            lambda x: x,
            lambda x: x.astype(CDT(["foo", "bar", "bah"])),
            lambda x: x.astype(CDT(ordered=True)),
        ],
    )
    def test_dtype_on_merged_different(self, change, join_type, left, right):
        # our merging columns, X now has 2 different dtypes
        # so we must be object as a result

        X = change(right.X.astype("object"))
        right = right.assign(X=X)
        assert is_categorical_dtype(left.X.values)
        # assert not left.X.values.is_dtype_equal(right.X.values)

        merged = pd.merge(left, right, on="X", how=join_type)

        result = merged.dtypes.sort_index()
        expected = Series(
            [np.dtype("O"), np.dtype("O"), np.dtype("int64")], index=["X", "Y", "Z"]
        )
        tm.assert_series_equal(result, expected)

    def test_self_join_multiple_categories(self):
        # GH 16767
        # non-duplicates should work with multiple categories
        m = 5
        df = pd.DataFrame(
            {
                "a": ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"] * m,
                "b": ["t", "w", "x", "y", "z"] * 2 * m,
                "c": [
                    letter
                    for each in ["m", "n", "u", "p", "o"]
                    for letter in [each] * 2 * m
                ],
                "d": [
                    letter
                    for each in [
                        "aa",
                        "bb",
                        "cc",
                        "dd",
                        "ee",
                        "ff",
                        "gg",
                        "hh",
                        "ii",
                        "jj",
                    ]
                    for letter in [each] * m
                ],
            }
        )

        # change them all to categorical variables
        df = df.apply(lambda x: x.astype("category"))

        # self-join should equal ourselves
        result = pd.merge(df, df, on=list(df.columns))

        tm.assert_frame_equal(result, df)

    def test_dtype_on_categorical_dates(self):
        # GH 16900
        # dates should not be coerced to ints

        df = pd.DataFrame(
            [[date(2001, 1, 1), 1.1], [date(2001, 1, 2), 1.3]], columns=["date", "num2"]
        )
        df["date"] = df["date"].astype("category")

        df2 = pd.DataFrame(
            [[date(2001, 1, 1), 1.3], [date(2001, 1, 3), 1.4]], columns=["date", "num4"]
        )
        df2["date"] = df2["date"].astype("category")

        expected_outer = pd.DataFrame(
            [
                [pd.Timestamp("2001-01-01"), 1.1, 1.3],
                [pd.Timestamp("2001-01-02"), 1.3, np.nan],
                [pd.Timestamp("2001-01-03"), np.nan, 1.4],
            ],
            columns=["date", "num2", "num4"],
        )
        result_outer = pd.merge(df, df2, how="outer", on=["date"])
        tm.assert_frame_equal(result_outer, expected_outer)

        expected_inner = pd.DataFrame(
            [[pd.Timestamp("2001-01-01"), 1.1, 1.3]], columns=["date", "num2", "num4"]
        )
        result_inner = pd.merge(df, df2, how="inner", on=["date"])
        tm.assert_frame_equal(result_inner, expected_inner)

    @pytest.mark.parametrize("ordered", [True, False])
    @pytest.mark.parametrize(
        "category_column,categories,expected_categories",
        [
            ([False, True, True, False], [True, False], [True, False]),
            ([2, 1, 1, 2], [1, 2], [1, 2]),
            (["False", "True", "True", "False"], ["True", "False"], ["True", "False"]),
        ],
    )
    def test_merging_with_bool_or_int_cateorical_column(
        self, category_column, categories, expected_categories, ordered
    ):
        # GH 17187
        # merging with a boolean/int categorical column
        df1 = pd.DataFrame({"id": [1, 2, 3, 4], "cat": category_column})
        df1["cat"] = df1["cat"].astype(CDT(categories, ordered=ordered))
        df2 = pd.DataFrame({"id": [2, 4], "num": [1, 9]})
        result = df1.merge(df2)
        expected = pd.DataFrame(
            {"id": [2, 4], "cat": expected_categories, "num": [1, 9]}
        )
        expected["cat"] = expected["cat"].astype(CDT(categories, ordered=ordered))
        tm.assert_frame_equal(expected, result)

    def test_merge_on_int_array(self):
        # GH 23020
        df = pd.DataFrame({"A": pd.Series([1, 2, np.nan], dtype="Int64"), "B": 1})
        result = pd.merge(df, df, on="A")
        expected = pd.DataFrame(
            {"A": pd.Series([1, 2, np.nan], dtype="Int64"), "B_x": 1, "B_y": 1}
        )
        tm.assert_frame_equal(result, expected)


@pytest.fixture
def left_df():
    return DataFrame({"a": [20, 10, 0]}, index=[2, 1, 0])


@pytest.fixture
def right_df():
    return DataFrame({"b": [300, 100, 200]}, index=[3, 1, 2])


class TestMergeOnIndexes:
    @pytest.mark.parametrize(
        "how, sort, expected",
        [
            ("inner", False, DataFrame({"a": [20, 10], "b": [200, 100]}, index=[2, 1])),
            ("inner", True, DataFrame({"a": [10, 20], "b": [100, 200]}, index=[1, 2])),
            (
                "left",
                False,
                DataFrame({"a": [20, 10, 0], "b": [200, 100, np.nan]}, index=[2, 1, 0]),
            ),
            (
                "left",
                True,
                DataFrame({"a": [0, 10, 20], "b": [np.nan, 100, 200]}, index=[0, 1, 2]),
            ),
            (
                "right",
                False,
                DataFrame(
                    {"a": [np.nan, 10, 20], "b": [300, 100, 200]}, index=[3, 1, 2]
                ),
            ),
            (
                "right",
                True,
                DataFrame(
                    {"a": [10, 20, np.nan], "b": [100, 200, 300]}, index=[1, 2, 3]
                ),
            ),
            (
                "outer",
                False,
                DataFrame(
                    {"a": [0, 10, 20, np.nan], "b": [np.nan, 100, 200, 300]},
                    index=[0, 1, 2, 3],
                ),
            ),
            (
                "outer",
                True,
                DataFrame(
                    {"a": [0, 10, 20, np.nan], "b": [np.nan, 100, 200, 300]},
                    index=[0, 1, 2, 3],
                ),
            ),
        ],
    )
    def test_merge_on_indexes(self, left_df, right_df, how, sort, expected):
        result = pd.merge(
            left_df, right_df, left_index=True, right_index=True, how=how, sort=sort
        )
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        CategoricalIndex(["A", "B"], categories=["A", "B"], name="index_col"),
        Float64Index([1.0, 2.0], name="index_col"),
        Int64Index([1, 2], name="index_col"),
        UInt64Index([1, 2], name="index_col"),
        RangeIndex(start=0, stop=2, name="index_col"),
        DatetimeIndex(["2018-01-01", "2018-01-02"], name="index_col"),
    ],
    ids=lambda x: type(x).__name__,
)
def test_merge_index_types(index):
    # gh-20777
    # assert key access is consistent across index types
    left = DataFrame({"left_data": [1, 2]}, index=index)
    right = DataFrame({"right_data": [1.0, 2.0]}, index=index)

    result = left.merge(right, on=["index_col"])

    expected = DataFrame(
        OrderedDict([("left_data", [1, 2]), ("right_data", [1.0, 2.0])]), index=index
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "on,left_on,right_on,left_index,right_index,nm",
    [
        (["outer", "inner"], None, None, False, False, "B"),
        (None, None, None, True, True, "B"),
        (None, ["outer", "inner"], None, False, True, "B"),
        (None, None, ["outer", "inner"], True, False, "B"),
        (["outer", "inner"], None, None, False, False, None),
        (None, None, None, True, True, None),
        (None, ["outer", "inner"], None, False, True, None),
        (None, None, ["outer", "inner"], True, False, None),
    ],
)
def test_merge_series(on, left_on, right_on, left_index, right_index, nm):
    # GH 21220
    a = pd.DataFrame(
        {"A": [1, 2, 3, 4]},
        index=pd.MultiIndex.from_product(
            [["a", "b"], [0, 1]], names=["outer", "inner"]
        ),
    )
    b = pd.Series(
        [1, 2, 3, 4],
        index=pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=["outer", "inner"]
        ),
        name=nm,
    )
    expected = pd.DataFrame(
        {"A": [2, 4], "B": [1, 3]},
        index=pd.MultiIndex.from_product([["a", "b"], [1]], names=["outer", "inner"]),
    )
    if nm is not None:
        result = pd.merge(
            a,
            b,
            on=on,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
        )
        tm.assert_frame_equal(result, expected)
    else:
        msg = "Cannot merge a Series without a name"
        with pytest.raises(ValueError, match=msg):
            result = pd.merge(
                a,
                b,
                on=on,
                left_on=left_on,
                right_on=right_on,
                left_index=left_index,
                right_index=right_index,
            )


@pytest.mark.parametrize(
    "col1, col2, kwargs, expected_cols",
    [
        (0, 0, dict(suffixes=("", "_dup")), ["0", "0_dup"]),
        (0, 0, dict(suffixes=(None, "_dup")), [0, "0_dup"]),
        (0, 0, dict(suffixes=("_x", "_y")), ["0_x", "0_y"]),
        ("a", 0, dict(suffixes=(None, "_y")), ["a", 0]),
        (0.0, 0.0, dict(suffixes=("_x", None)), ["0.0_x", 0.0]),
        ("b", "b", dict(suffixes=(None, "_y")), ["b", "b_y"]),
        ("a", "a", dict(suffixes=("_x", None)), ["a_x", "a"]),
        ("a", "b", dict(suffixes=("_x", None)), ["a", "b"]),
        ("a", "a", dict(suffixes=[None, "_x"]), ["a", "a_x"]),
        (0, 0, dict(suffixes=["_a", None]), ["0_a", 0]),
        ("a", "a", dict(), ["a_x", "a_y"]),
        (0, 0, dict(), ["0_x", "0_y"]),
    ],
)
def test_merge_suffix(col1, col2, kwargs, expected_cols):
    # issue: 24782
    a = pd.DataFrame({col1: [1, 2, 3]})
    b = pd.DataFrame({col2: [4, 5, 6]})

    expected = pd.DataFrame([[1, 4], [2, 5], [3, 6]], columns=expected_cols)

    result = a.merge(b, left_index=True, right_index=True, **kwargs)
    tm.assert_frame_equal(result, expected)

    result = pd.merge(a, b, left_index=True, right_index=True, **kwargs)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "col1, col2, suffixes",
    [
        ("a", "a", [None, None]),
        ("a", "a", (None, None)),
        ("a", "a", ("", None)),
        (0, 0, [None, None]),
        (0, 0, (None, "")),
    ],
)
def test_merge_suffix_error(col1, col2, suffixes):
    # issue: 24782
    a = pd.DataFrame({col1: [1, 2, 3]})
    b = pd.DataFrame({col2: [3, 4, 5]})

    # TODO: might reconsider current raise behaviour, see issue 24782
    msg = "columns overlap but no suffix specified"
    with pytest.raises(ValueError, match=msg):
        pd.merge(a, b, left_index=True, right_index=True, suffixes=suffixes)


@pytest.mark.parametrize("col1, col2, suffixes", [("a", "a", None), (0, 0, None)])
def test_merge_suffix_none_error(col1, col2, suffixes):
    # issue: 24782
    a = pd.DataFrame({col1: [1, 2, 3]})
    b = pd.DataFrame({col2: [3, 4, 5]})

    # TODO: might reconsider current raise behaviour, see GH24782
    msg = "iterable"
    with pytest.raises(TypeError, match=msg):
        pd.merge(a, b, left_index=True, right_index=True, suffixes=suffixes)


@pytest.mark.parametrize("cat_dtype", ["one", "two"])
@pytest.mark.parametrize("reverse", [True, False])
def test_merge_equal_cat_dtypes(cat_dtype, reverse):
    # see gh-22501
    cat_dtypes = {
        "one": CategoricalDtype(categories=["a", "b", "c"], ordered=False),
        "two": CategoricalDtype(categories=["a", "b", "c"], ordered=False),
    }

    df1 = DataFrame(
        {"foo": Series(["a", "b", "c"]).astype(cat_dtypes["one"]), "left": [1, 2, 3]}
    ).set_index("foo")

    data_foo = ["a", "b", "c"]
    data_right = [1, 2, 3]

    if reverse:
        data_foo.reverse()
        data_right.reverse()

    df2 = DataFrame(
        {"foo": Series(data_foo).astype(cat_dtypes[cat_dtype]), "right": data_right}
    ).set_index("foo")

    result = df1.merge(df2, left_index=True, right_index=True)

    expected = DataFrame(
        {
            "left": [1, 2, 3],
            "right": [1, 2, 3],
            "foo": Series(["a", "b", "c"]).astype(cat_dtypes["one"]),
        }
    ).set_index("foo")

    # Categorical is unordered, so don't check ordering.
    tm.assert_frame_equal(result, expected, check_categorical=False)


def test_merge_equal_cat_dtypes2():
    # see gh-22501
    cat_dtype = CategoricalDtype(categories=["a", "b", "c"], ordered=False)

    # Test Data
    df1 = DataFrame(
        {"foo": Series(["a", "b"]).astype(cat_dtype), "left": [1, 2]}
    ).set_index("foo")

    df2 = DataFrame(
        {"foo": Series(["a", "b", "c"]).astype(cat_dtype), "right": [3, 2, 1]}
    ).set_index("foo")

    result = df1.merge(df2, left_index=True, right_index=True)

    expected = DataFrame(
        {"left": [1, 2], "right": [3, 2], "foo": Series(["a", "b"]).astype(cat_dtype)}
    ).set_index("foo")

    # Categorical is unordered, so don't check ordering.
    tm.assert_frame_equal(result, expected, check_categorical=False)


def test_merge_on_cat_and_ext_array():
    # GH 28668
    right = DataFrame(
        {"a": Series([pd.Interval(0, 1), pd.Interval(1, 2)], dtype="interval")}
    )
    left = right.copy()
    left["a"] = left["a"].astype("category")

    result = pd.merge(left, right, how="inner", on="a")
    expected = right.copy()

    tm.assert_frame_equal(result, expected)


def test_merge_multiindex_columns():
    # Issue #28518
    # Verify that merging two dataframes give the expected labels
    # The original cause of this issue come from a bug lexsort_depth and is tested in
    # test_lexsort_depth

    letters = ["a", "b", "c", "d"]
    numbers = ["1", "2", "3"]
    index = pd.MultiIndex.from_product((letters, numbers), names=["outer", "inner"])

    frame_x = pd.DataFrame(columns=index)
    frame_x["id"] = ""
    frame_y = pd.DataFrame(columns=index)
    frame_y["id"] = ""

    l_suf = "_x"
    r_suf = "_y"
    result = frame_x.merge(frame_y, on="id", suffixes=((l_suf, r_suf)))

    # Constructing the expected results
    expected_labels = [l + l_suf for l in letters] + [l + r_suf for l in letters]
    expected_index = pd.MultiIndex.from_product(
        [expected_labels, numbers], names=["outer", "inner"]
    )
    expected = pd.DataFrame(columns=expected_index)
    expected["id"] = ""

    tm.assert_frame_equal(result, expected)


def test_merge_datetime_upcast_dtype():
    # https://github.com/pandas-dev/pandas/issues/31208
    df1 = pd.DataFrame({"x": ["a", "b", "c"], "y": ["1", "2", "4"]})
    df2 = pd.DataFrame(
        {"y": ["1", "2", "3"], "z": pd.to_datetime(["2000", "2001", "2002"])}
    )
    result = pd.merge(df1, df2, how="left", on="y")
    expected = pd.DataFrame(
        {
            "x": ["a", "b", "c"],
            "y": ["1", "2", "4"],
            "z": pd.to_datetime(["2000", "2001", "NaT"]),
        }
    )
    tm.assert_frame_equal(result, expected)
