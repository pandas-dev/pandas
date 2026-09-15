from collections import (
    abc,
    deque,
)
from collections.abc import Iterator
from datetime import datetime
from decimal import Decimal
import itertools

import numpy as np
import pytest

from pandas.errors import (
    InvalidIndexError,
    Pandas4Warning,
)

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import SparseArray
from pandas.tests.extension.decimal import to_decimal


class TestConcatenate:
    def test_append_concat(self):
        # GH#1815
        d1 = pd.date_range("12/31/1990", "12/31/1999", freq="YE-DEC")
        d2 = pd.date_range("12/31/2000", "12/31/2009", freq="YE-DEC")

        s1 = pd.Series(np.random.default_rng(2).standard_normal(10), d1)
        s2 = pd.Series(np.random.default_rng(2).standard_normal(10), d2)

        s1 = s1.to_period()
        s2 = s2.to_period()

        # drops index
        result = pd.concat([s1, s2])
        assert isinstance(result.index, pd.PeriodIndex)
        assert result.index[0] == s1.index[0]

    def test_concat_copy(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3)))
        df2 = pd.DataFrame(
            np.random.default_rng(2).integers(0, 10, size=4).reshape(4, 1)
        )
        df3 = pd.DataFrame({5: "foo"}, index=range(4))

        result = pd.concat([df, df2, df3], axis=1)

        for block in result._mgr.blocks:
            arr = block.values
            if arr.dtype.kind == "f":
                assert arr.base is df._mgr.blocks[0].values.base
            elif arr.dtype.kind in ["i", "u"]:
                assert arr.base is df2._mgr.blocks[0].values.base
            elif arr.dtype == object:
                assert arr.base is not None
            elif arr.dtype == "string":
                tm.shares_memory(arr, df3._mgr.blocks[0].values)

        # Float block was consolidated.
        df4 = pd.DataFrame(np.random.default_rng(2).standard_normal((4, 1)))
        result = pd.concat([df, df2, df3, df4], axis=1)
        for blocks in result._mgr.blocks:
            arr = blocks.values
            if arr.dtype.kind == "f":
                # this is a view on some array in either df or df4
                assert any(
                    np.shares_memory(arr, block.values)
                    for block in itertools.chain(df._mgr.blocks, df4._mgr.blocks)
                )
            elif arr.dtype.kind in ["i", "u"]:
                assert arr.base is df2._mgr.blocks[0].values.base
            elif arr.dtype == object:
                # this is a view on df3
                assert any(
                    np.shares_memory(arr, block.values) for block in df3._mgr.blocks
                )

    def test_concat_with_group_keys(self):
        # axis=0
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((3, 4)))
        df2 = pd.DataFrame(np.random.default_rng(2).standard_normal((4, 4)))

        result = pd.concat([df, df2], keys=[0, 1])
        exp_index = pd.MultiIndex.from_arrays(
            [[0, 0, 0, 1, 1, 1, 1], [0, 1, 2, 0, 1, 2, 3]]
        )
        expected = pd.DataFrame(np.r_[df.values, df2.values], index=exp_index)
        tm.assert_frame_equal(result, expected)

        result = pd.concat([df, df], keys=[0, 1])
        exp_index2 = pd.MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]])
        expected = pd.DataFrame(np.r_[df.values, df.values], index=exp_index2)
        tm.assert_frame_equal(result, expected)

        # axis=1
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3)))
        df2 = pd.DataFrame(np.random.default_rng(2).standard_normal((4, 4)))

        result = pd.concat([df, df2], keys=[0, 1], axis=1)
        expected = pd.DataFrame(np.c_[df.values, df2.values], columns=exp_index)
        tm.assert_frame_equal(result, expected)

        result = pd.concat([df, df], keys=[0, 1], axis=1)
        expected = pd.DataFrame(np.c_[df.values, df.values], columns=exp_index2)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_specific_levels(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((10, 4)))
        pieces = [df.iloc[:, [0, 1]], df.iloc[:, [2]], df.iloc[:, [3]]]
        level = ["three", "two", "one", "zero"]
        result = pd.concat(
            pieces,
            axis=1,
            keys=["one", "two", "three"],
            levels=[level],
            names=["group_key"],
        )

        tm.assert_index_equal(
            result.columns.levels[0], pd.Index(level, name="group_key")
        )
        tm.assert_index_equal(
            result.columns.levels[1], pd.RangeIndex(start=0, stop=4, step=1), exact=True
        )

        assert result.columns.names == ["group_key", None]

    @pytest.mark.parametrize("mapping", ["mapping", "dict"])
    def test_concat_mapping(self, mapping, non_dict_mapping_subclass):
        constructor = dict if mapping == "dict" else non_dict_mapping_subclass
        frames = constructor(
            {
                "foo": pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3))),
                "bar": pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3))),
                "baz": pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3))),
                "qux": pd.DataFrame(np.random.default_rng(2).standard_normal((4, 3))),
            }
        )

        sorted_keys = list(frames.keys())

        result = pd.concat(frames)
        expected = pd.concat([frames[k] for k in sorted_keys], keys=sorted_keys)
        tm.assert_frame_equal(result, expected)

        result = pd.concat(frames, axis=1)
        expected = pd.concat([frames[k] for k in sorted_keys], keys=sorted_keys, axis=1)
        tm.assert_frame_equal(result, expected)

        keys = ["baz", "foo", "bar"]
        result = pd.concat(frames, keys=keys)
        expected = pd.concat([frames[k] for k in keys], keys=keys)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_and_levels(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 3)))
        df2 = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 4)))

        levels = [["foo", "baz"], ["one", "two"]]
        names = ["first", "second"]
        result = pd.concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            levels=levels,
            names=names,
        )
        expected = pd.concat([df, df2, df, df2])
        exp_index = pd.MultiIndex(
            levels=[*levels, [0]],
            codes=[[0, 0, 1, 1], [0, 1, 0, 1], [0, 0, 0, 0]],
            names=[*names, None],
        )
        expected.index = exp_index

        tm.assert_frame_equal(result, expected)

        # no names
        result = pd.concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            levels=levels,
        )
        assert result.index.names == (None,) * 3

        # no levels
        result = pd.concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            names=["first", "second"],
        )
        assert result.index.names == ("first", "second", None)
        tm.assert_index_equal(
            result.index.levels[0], pd.Index(["baz", "foo"], name="first")
        )

    def test_concat_keys_levels_no_overlap(self):
        # GH #1406
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 3)), index=["a"])
        df2 = pd.DataFrame(
            np.random.default_rng(2).standard_normal((1, 4)), index=["b"]
        )

        msg = "Values not found in passed level"
        with pytest.raises(ValueError, match=msg):
            pd.concat([df, df], keys=["one", "two"], levels=[["foo", "bar", "baz"]])

        msg = "Key one not in level"
        with pytest.raises(ValueError, match=msg):
            pd.concat([df, df2], keys=["one", "two"], levels=[["foo", "bar", "baz"]])

    def test_crossed_dtypes_weird_corner(self):
        columns = ["A", "B", "C", "D"]
        df1 = pd.DataFrame(
            {
                "A": np.array([1, 2, 3, 4], dtype="f8"),
                "B": np.array([1, 2, 3, 4], dtype="i8"),
                "C": np.array([1, 2, 3, 4], dtype="f8"),
                "D": np.array([1, 2, 3, 4], dtype="i8"),
            },
            columns=columns,
        )

        df2 = pd.DataFrame(
            {
                "A": np.array([1, 2, 3, 4], dtype="i8"),
                "B": np.array([1, 2, 3, 4], dtype="f8"),
                "C": np.array([1, 2, 3, 4], dtype="i8"),
                "D": np.array([1, 2, 3, 4], dtype="f8"),
            },
            columns=columns,
        )

        appended = pd.concat([df1, df2], ignore_index=True)
        expected = pd.DataFrame(
            np.concatenate([df1.values, df2.values], axis=0), columns=columns
        )
        tm.assert_frame_equal(appended, expected)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 3)), index=["a"])
        df2 = pd.DataFrame(
            np.random.default_rng(2).standard_normal((1, 4)), index=["b"]
        )
        result = pd.concat([df, df2], keys=["one", "two"], names=["first", "second"])
        assert result.index.names == ("first", "second")

    def test_with_mixed_tuples(self, sort):
        # 10697
        # columns have mixed tuples, so handle properly
        df1 = pd.DataFrame({"A": "foo", ("B", 1): "bar"}, index=range(2))
        df2 = pd.DataFrame({"B": "foo", ("B", 1): "bar"}, index=range(2))

        # it works
        pd.concat([df1, df2], sort=sort)

    def test_concat_mixed_objs_columns(self):
        # Test column-wise concat for mixed series/frames (axis=1)
        # G2385

        index = pd.date_range("01-Jan-2013", periods=10, freq="h")
        arr = np.arange(10, dtype="int64")
        s1 = pd.Series(arr, index=index)
        s2 = pd.Series(arr, index=index)
        df = pd.DataFrame(arr.reshape(-1, 1), index=index)

        expected = pd.DataFrame(
            np.repeat(arr, 2).reshape(-1, 2), index=index, columns=[0, 0]
        )
        result = pd.concat([df, df], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(
            np.repeat(arr, 2).reshape(-1, 2), index=index, columns=[0, 1]
        )
        result = pd.concat([s1, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=[0, 1, 2]
        )
        result = pd.concat([s1, s2, s1], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(
            np.repeat(arr, 5).reshape(-1, 5), index=index, columns=[0, 0, 1, 2, 3]
        )
        result = pd.concat([s1, df, s2, s2, s1], axis=1)
        tm.assert_frame_equal(result, expected)

        # with names
        s1.name = "foo"
        expected = pd.DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=["foo", 0, 0]
        )
        result = pd.concat([s1, df, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        s2.name = "bar"
        expected = pd.DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=["foo", 0, "bar"]
        )
        result = pd.concat([s1, df, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        # ignore index
        expected = pd.DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=[0, 1, 2]
        )
        result = pd.concat([s1, df, s2], axis=1, ignore_index=True)
        tm.assert_frame_equal(result, expected)

    def test_concat_mixed_objs_index(self):
        # Test row-wise concat for mixed series/frames with a common name
        # GH2385, GH15047

        index = pd.date_range("01-Jan-2013", periods=10, freq="h")
        arr = np.arange(10, dtype="int64")
        s1 = pd.Series(arr, index=index)
        s2 = pd.Series(arr, index=index)
        df = pd.DataFrame(arr.reshape(-1, 1), index=index)

        expected = pd.DataFrame(
            np.tile(arr, 3).reshape(-1, 1), index=index.tolist() * 3, columns=[0]
        )
        result = pd.concat([s1, df, s2])
        tm.assert_frame_equal(result, expected)

    def test_concat_mixed_objs_index_names(self):
        # Test row-wise concat for mixed series/frames with distinct names
        # GH2385, GH15047
        # GH #60723 & GH #56257 (Updated the test case,
        # as the above GH PR ones were incorrect)

        index = pd.date_range("01-Jan-2013", periods=10, freq="h")
        arr = np.arange(10, dtype="int64")
        s1 = pd.Series(arr, index=index, name="foo")
        s2 = pd.Series(arr, index=index, name="bar")
        df = pd.DataFrame(arr.reshape(-1, 1), index=index)

        expected = pd.DataFrame(
            np.kron(np.where(np.identity(3) == 1, 1, np.nan), arr).T,
            index=index.tolist() * 3,
            columns=["foo", 0, "bar"],
        )
        result = pd.concat([s1, df, s2])
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(
            np.kron(np.where(np.identity(3) == 1, 1, np.nan), arr).T,
            index=np.arange(30, dtype=np.int64),
            columns=["foo", 0, "bar"],
        )
        result = pd.concat([s1, df, s2], ignore_index=True)
        tm.assert_frame_equal(result, expected)

    def test_dtype_coercion(self):
        # 12411
        df = pd.DataFrame(
            {"date": [pd.Timestamp("20130101").tz_localize("UTC"), pd.NaT]}
        )

        result = pd.concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 12045
        df = pd.DataFrame({"date": [datetime(2012, 1, 1), datetime(1012, 1, 2)]})
        result = pd.concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 11594
        df = pd.DataFrame({"text": ["some words"] + [None] * 9})
        result = pd.concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

    def test_concat_single_with_key(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((10, 4)))

        result = pd.concat([df], keys=["foo"])
        expected = pd.concat([df, df], keys=["foo", "bar"])
        tm.assert_frame_equal(result, expected[:10])

    def test_concat_no_items_raises(self):
        with pytest.raises(ValueError, match="No objects to concatenate"):
            pd.concat([])

    def test_concat_exclude_none(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((10, 4)))

        pieces = [df[:5], None, None, df[5:]]
        result = pd.concat(pieces)
        tm.assert_frame_equal(result, df)
        with pytest.raises(ValueError, match="All objects passed were None"):
            pd.concat([None, None])

    def test_concat_keys_with_none(self):
        # #1649
        df0 = pd.DataFrame([[10, 20, 30], [10, 20, 30], [10, 20, 30]])

        result = pd.concat({"a": None, "b": df0, "c": df0[:2], "d": df0[:1], "e": df0})
        expected = pd.concat({"b": df0, "c": df0[:2], "d": df0[:1], "e": df0})
        tm.assert_frame_equal(result, expected)

        result = pd.concat(
            [None, df0, df0[:2], df0[:1], df0], keys=["a", "b", "c", "d", "e"]
        )
        expected = pd.concat([df0, df0[:2], df0[:1], df0], keys=["b", "c", "d", "e"])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("klass", [range, pd.RangeIndex])
    @pytest.mark.parametrize("include_none", [True, False])
    def test_concat_preserves_rangeindex(self, klass, include_none):
        df = pd.DataFrame([1, 2])
        df2 = pd.DataFrame([3, 4])
        data = [df, None, df2, None] if include_none else [df, df2]
        keys_length = 4 if include_none else 2
        result = pd.concat(data, keys=klass(keys_length))
        expected = pd.DataFrame(
            [1, 2, 3, 4],
            index=pd.MultiIndex(
                levels=(
                    pd.RangeIndex(start=0, stop=keys_length, step=keys_length / 2),
                    pd.RangeIndex(start=0, stop=2, step=1),
                ),
                codes=(
                    np.array([0, 0, 1, 1], dtype=np.int8),
                    np.array([0, 1, 0, 1], dtype=np.int8),
                ),
            ),
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_bug_1719(self):
        ts1 = pd.Series(
            np.arange(10, dtype=np.float64),
            index=pd.date_range("2020-01-01", periods=10),
        )
        ts2 = ts1[::2]

        # to join with union
        # these two are of different length!
        left = pd.concat([ts1, ts2], join="outer", axis=1)
        msg = "Sorting by default when concatenating all DatetimeIndex is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            right = pd.concat([ts2, ts1], join="outer", axis=1)

        assert len(left) == len(right)

    def test_concat_bug_2972(self):
        ts0 = pd.Series(np.zeros(5))
        ts1 = pd.Series(np.ones(5))
        ts0.name = ts1.name = "same name"
        result = pd.concat([ts0, ts1], axis=1)

        expected = pd.DataFrame({0: ts0, 1: ts1})
        expected.columns = ["same name", "same name"]
        tm.assert_frame_equal(result, expected)

    def test_concat_bug_3602(self):
        # GH 3602, duplicate columns
        df1 = pd.DataFrame(
            {
                "firmNo": [0, 0, 0, 0],
                "prc": [6, 6, 6, 6],
                "stringvar": ["rrr", "rrr", "rrr", "rrr"],
            }
        )
        df2 = pd.DataFrame(
            {"C": [9, 10, 11, 12], "misc": [1, 2, 3, 4], "prc": [6, 6, 6, 6]}
        )
        expected = pd.DataFrame(
            [
                [0, 6, "rrr", 9, 1, 6],
                [0, 6, "rrr", 10, 2, 6],
                [0, 6, "rrr", 11, 3, 6],
                [0, 6, "rrr", 12, 4, 6],
            ]
        )
        expected.columns = ["firmNo", "prc", "stringvar", "C", "misc", "prc"]

        result = pd.concat([df1, df2], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_concat_iterables(self):
        # GH8645 check concat works with tuples, list, generators, and weird
        # stuff like deque and custom iterables
        df1 = pd.DataFrame([1, 2, 3])
        df2 = pd.DataFrame([4, 5, 6])
        expected = pd.DataFrame([1, 2, 3, 4, 5, 6])
        tm.assert_frame_equal(pd.concat((df1, df2), ignore_index=True), expected)
        tm.assert_frame_equal(pd.concat([df1, df2], ignore_index=True), expected)
        tm.assert_frame_equal(
            pd.concat((df for df in (df1, df2)), ignore_index=True), expected
        )
        tm.assert_frame_equal(pd.concat(deque((df1, df2)), ignore_index=True), expected)

        class CustomIterator1:
            def __len__(self) -> int:
                return 2

            def __getitem__(self, index):
                try:
                    return {0: df1, 1: df2}[index]
                except KeyError as err:
                    raise IndexError from err

        tm.assert_frame_equal(pd.concat(CustomIterator1(), ignore_index=True), expected)

        class CustomIterator2(abc.Iterable):
            def __iter__(self) -> Iterator:
                yield df1
                yield df2

        tm.assert_frame_equal(pd.concat(CustomIterator2(), ignore_index=True), expected)

    def test_concat_order(self):
        # GH 17344, GH#47331
        dfs = [pd.DataFrame(index=range(3), columns=["a", 1, None])]
        dfs += [
            pd.DataFrame(index=range(3), columns=[None, 1, "a"]) for _ in range(100)
        ]

        result = pd.concat(dfs, sort=True).columns
        expected = pd.Index([1, "a", None])
        tm.assert_index_equal(result, expected)

    def test_concat_different_extension_dtypes_upcasts(self):
        a = pd.Series(pd.array([1, 2], dtype="Int64"))
        b = pd.Series(to_decimal([1, 2]))

        result = pd.concat([a, b], ignore_index=True)
        expected = pd.Series([1, 2, Decimal(1), Decimal(2)], dtype=object)
        tm.assert_series_equal(result, expected)

    def test_concat_ordered_dict(self):
        # GH 21510
        expected = pd.concat(
            [pd.Series(range(3)), pd.Series(range(4))], keys=["First", "Another"]
        )
        result = pd.concat(
            {"First": pd.Series(range(3)), "Another": pd.Series(range(4))}
        )
        tm.assert_series_equal(result, expected)

    def test_concat_duplicate_indices_raise(self):
        # GH 45888: test raise for concat DataFrames with duplicate indices
        # https://github.com/pandas-dev/pandas/issues/36263
        df1 = pd.DataFrame(
            np.random.default_rng(2).standard_normal(5),
            index=[0, 1, 2, 3, 3],
            columns=["a"],
        )
        df2 = pd.DataFrame(
            np.random.default_rng(2).standard_normal(5),
            index=[0, 1, 2, 2, 4],
            columns=["b"],
        )
        msg = "Reindexing only valid with uniquely valued Index objects"
        with pytest.raises(InvalidIndexError, match=msg):
            pd.concat([df1, df2], axis=1)


def test_concat_no_unnecessary_upcast(float_numpy_dtype, frame_or_series):
    # GH 13247
    dims = frame_or_series(dtype=object).ndim
    dt = float_numpy_dtype

    dfs = [
        frame_or_series(np.array([1], dtype=dt, ndmin=dims)),
        frame_or_series(np.array([np.nan], dtype=dt, ndmin=dims)),
        frame_or_series(np.array([5], dtype=dt, ndmin=dims)),
    ]
    x = pd.concat(dfs)
    assert x.values.dtype == dt


def test_concat_will_upcast(frame_or_series, any_signed_int_numpy_dtype):
    dt = any_signed_int_numpy_dtype
    dims = frame_or_series().ndim
    dfs = [
        frame_or_series(np.array([1], dtype=dt, ndmin=dims)),
        frame_or_series(np.array([np.nan], ndmin=dims)),
        frame_or_series(np.array([5], dtype=dt, ndmin=dims)),
    ]
    x = pd.concat(dfs)
    assert x.values.dtype == "float64"


def test_concat_empty_and_non_empty_frame_regression():
    # GH 18178 regression test
    df1 = pd.DataFrame({"foo": [1]})
    df2 = pd.DataFrame({"foo": []})
    expected = pd.DataFrame({"foo": [1.0]})
    result = pd.concat([df1, df2])
    tm.assert_frame_equal(result, expected)


def test_concat_sparse():
    # GH 23557
    a = pd.Series(SparseArray([0, 1, 2]))
    expected = pd.DataFrame(data=[[0, 0], [1, 1], [2, 2]]).astype(
        pd.SparseDtype(np.int64, 0)
    )
    result = pd.concat([a, a], axis=1)
    tm.assert_frame_equal(result, expected)


def test_concat_dense_sparse():
    # GH 30668
    dtype = pd.SparseDtype(np.float64, None)
    a = pd.Series(pd.arrays.SparseArray([1, None]), dtype=dtype)
    b = pd.Series([1], dtype=float)
    expected = pd.Series(data=[1, None, 1], index=[0, 1, 0]).astype(dtype)
    result = pd.concat([a, b], axis=0)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("keys", [["e", "f", "f"], ["f", "e", "f"]])
def test_duplicate_keys(keys):
    # GH 33654
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    s1 = pd.Series([7, 8, 9], name="c")
    s2 = pd.Series([10, 11, 12], name="d")
    result = pd.concat([df, s1, s2], axis=1, keys=keys)
    expected_values = [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
    expected_columns = pd.MultiIndex.from_tuples(
        [(keys[0], "a"), (keys[0], "b"), (keys[1], "c"), (keys[2], "d")]
    )
    expected = pd.DataFrame(expected_values, columns=expected_columns)
    tm.assert_frame_equal(result, expected)


def test_duplicate_keys_same_frame():
    # GH 43595
    keys = ["e", "e"]
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    result = pd.concat([df, df], axis=1, keys=keys)
    expected_values = [[1, 4, 1, 4], [2, 5, 2, 5], [3, 6, 3, 6]]
    expected_columns = pd.MultiIndex.from_tuples(
        [(keys[0], "a"), (keys[0], "b"), (keys[1], "a"), (keys[1], "b")]
    )
    expected = pd.DataFrame(expected_values, columns=expected_columns)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "obj",
    [
        tm.SubclassedDataFrame({"A": np.arange(0, 10)}),
        tm.SubclassedSeries(np.arange(0, 10), name="A"),
    ],
)
def test_concat_preserves_subclass(obj):
    # GH28330 -- preserve subclass

    result = pd.concat([obj, obj])
    assert isinstance(result, type(obj))


def test_concat_frame_axis0_extension_dtypes():
    # preserve extension dtype (through common_dtype mechanism)
    df1 = pd.DataFrame({"a": pd.array([1, 2, 3], dtype="Int64")})
    df2 = pd.DataFrame({"a": np.array([4, 5, 6])})

    result = pd.concat([df1, df2], ignore_index=True)
    expected = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)

    result = pd.concat([df2, df1], ignore_index=True)
    expected = pd.DataFrame({"a": [4, 5, 6, 1, 2, 3]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)


def test_concat_preserves_extension_int64_dtype():
    # GH 24768
    df_a = pd.DataFrame({"a": [-1]}, dtype="Int64")
    df_b = pd.DataFrame({"b": [1]}, dtype="Int64")
    result = pd.concat([df_a, df_b], ignore_index=True)
    expected = pd.DataFrame({"a": [-1, None], "b": [None, 1]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype1,dtype2,expected_dtype",
    [
        ("bool", "bool", "bool"),
        ("boolean", "bool", "boolean"),
        ("bool", "boolean", "boolean"),
        ("boolean", "boolean", "boolean"),
    ],
)
def test_concat_bool_types(dtype1, dtype2, expected_dtype):
    # GH 42800
    ser1 = pd.Series([True, False], dtype=dtype1)
    ser2 = pd.Series([False, True], dtype=dtype2)
    result = pd.concat([ser1, ser2], ignore_index=True)
    expected = pd.Series([True, False, False, True], dtype=expected_dtype)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    ("keys", "integrity"),
    [
        (["red"] * 3, True),
        (["red"] * 3, False),
        (["red", "blue", "red"], False),
        (["red", "blue", "red"], True),
    ],
)
def test_concat_repeated_keys(keys, integrity):
    # GH: 20816
    series_list = [pd.Series({"a": 1}), pd.Series({"b": 2}), pd.Series({"c": 3})]
    result = pd.concat(series_list, keys=keys, verify_integrity=integrity)
    tuples = list(zip(keys, ["a", "b", "c"], strict=True))
    expected = pd.Series([1, 2, 3], index=pd.MultiIndex.from_tuples(tuples))
    tm.assert_series_equal(result, expected)


def test_concat_null_object_with_dti():
    # GH#40841
    dti = pd.DatetimeIndex(
        ["2021-04-08 21:21:14+00:00"], dtype="datetime64[ns, UTC]", name="Time (UTC)"
    )
    right = pd.DataFrame(data={"C": [0.5274]}, index=dti)

    idx = pd.Index([None], dtype="object", name="Maybe Time (UTC)")
    left = pd.DataFrame(data={"A": [None], "B": [np.nan]}, index=idx)

    result = pd.concat([left, right], axis="columns")

    exp_index = pd.Index([None, dti[0]], dtype=object)
    expected = pd.DataFrame(
        {
            "A": np.array([None, np.nan], dtype=object),
            "B": [np.nan, np.nan],
            "C": [np.nan, 0.5274],
        },
        index=exp_index,
    )
    tm.assert_frame_equal(result, expected)


def test_concat_multiindex_with_empty_rangeindex():
    # GH#41234
    mi = pd.MultiIndex.from_tuples([("B", 1), ("C", 1)])
    df1 = pd.DataFrame([[1, 2]], columns=mi)
    df2 = pd.DataFrame(index=[1], columns=pd.RangeIndex(0))

    result = pd.concat([df1, df2])
    expected = pd.DataFrame([[1, 2], [np.nan, np.nan]], columns=mi)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data",
    [
        pd.Series(data=[1, 2]),
        pd.DataFrame(
            data={
                "col1": [1, 2],
            }
        ),
        pd.DataFrame(dtype=float),
        pd.Series(dtype=float),
    ],
)
def test_concat_drop_attrs(data):
    # GH#41828
    df1 = data.copy()
    df1.attrs = {1: 1}
    df2 = data.copy()
    df2.attrs = {1: 2}
    df = pd.concat([df1, df2])
    assert len(df.attrs) == 0


@pytest.mark.parametrize(
    "data",
    [
        pd.Series(data=[1, 2]),
        pd.DataFrame(
            data={
                "col1": [1, 2],
            }
        ),
        pd.DataFrame(dtype=float),
        pd.Series(dtype=float),
    ],
)
def test_concat_retain_attrs(data):
    # GH#41828
    df1 = data.copy()
    df1.attrs = {1: 1}
    df2 = data.copy()
    df2.attrs = {1: 1}
    df = pd.concat([df1, df2])
    assert df.attrs[1] == 1


@pytest.mark.parametrize("df_dtype", ["float64", "int64", "datetime64[ns]"])
@pytest.mark.parametrize("empty_dtype", [None, "float64", "object"])
def test_concat_ignore_empty_object_float(empty_dtype, df_dtype):
    # https://github.com/pandas-dev/pandas/issues/45637
    df = pd.DataFrame({"foo": [1, 2], "bar": [1, 2]}, dtype=df_dtype)
    empty = pd.DataFrame(columns=["foo", "bar"], dtype=empty_dtype)

    needs_update = False
    if df_dtype == "datetime64[ns]" or (
        df_dtype == "float64" and empty_dtype != "float64"
    ):
        needs_update = True

    result = pd.concat([empty, df])
    expected = df
    if df_dtype == "int64":
        # TODO what exact behaviour do we want for integer eventually?
        if empty_dtype == "float64":
            expected = df.astype("float64")
        else:
            expected = df.astype("object")

    if needs_update:
        # GH#40893 changed the expected here to retain dependence on empty
        expected = expected.astype(object)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("df_dtype", ["float64", "int64", "datetime64[ns]"])
@pytest.mark.parametrize("empty_dtype", [None, "float64", "object"])
def test_concat_ignore_all_na_object_float(empty_dtype, df_dtype):
    df = pd.DataFrame({"foo": [1, 2], "bar": [1, 2]}, dtype=df_dtype)
    empty = pd.DataFrame({"foo": [np.nan], "bar": [np.nan]}, dtype=empty_dtype)

    if df_dtype == "int64":
        # TODO what exact behaviour do we want for integer eventually?
        if empty_dtype == "object":
            df_dtype = "object"
        else:
            df_dtype = "float64"

    needs_update = False
    if empty_dtype != df_dtype and empty_dtype is not None:
        needs_update = True
    elif df_dtype == "datetime64[ns]":
        needs_update = True

    result = pd.concat([empty, df], ignore_index=True)

    expected = pd.DataFrame(
        {"foo": [np.nan, 1, 2], "bar": [np.nan, 1, 2]}, dtype=df_dtype
    )
    if needs_update:
        # GH#40893 changed the expected here to retain dependence on empty
        expected = expected.astype(object)
        expected.iloc[0] = np.nan
    tm.assert_frame_equal(result, expected)


def test_concat_ignore_empty_from_reindex():
    # https://github.com/pandas-dev/pandas/pull/43507#issuecomment-920375856
    df1 = pd.DataFrame({"a": [1], "b": [pd.Timestamp("2012-01-01")]})
    df2 = pd.DataFrame({"a": [2]})

    aligned = df2.reindex(columns=df1.columns)

    result = pd.concat([df1, aligned], ignore_index=True)

    expected = pd.DataFrame(
        {
            "a": [1, 2],
            "b": pd.array([pd.Timestamp("2012-01-01"), np.nan], dtype=object),
        },
        dtype=object,
    )
    expected["a"] = expected["a"].astype("int64")
    tm.assert_frame_equal(result, expected)


def test_concat_mismatched_keys_length():
    # GH#43485
    ser = pd.Series(range(5))
    sers = [ser + n for n in range(4)]
    keys = ["A", "B", "C"]

    msg = r"The length of the keys"
    with pytest.raises(ValueError, match=msg):
        pd.concat(sers, keys=keys, axis=1)
    with pytest.raises(ValueError, match=msg):
        pd.concat(sers, keys=keys, axis=0)
    with pytest.raises(ValueError, match=msg):
        pd.concat((x for x in sers), keys=(y for y in keys), axis=1)
    with pytest.raises(ValueError, match=msg):
        pd.concat((x for x in sers), keys=(y for y in keys), axis=0)


def test_concat_multiindex_with_category():
    df1 = pd.DataFrame(
        {
            "c1": pd.Series(list("abc"), dtype="category"),
            "c2": pd.Series(list("eee"), dtype="category"),
            "i2": pd.Series([1, 2, 3]),
        }
    )
    df1 = df1.set_index(["c1", "c2"])
    df2 = pd.DataFrame(
        {
            "c1": pd.Series(list("abc"), dtype="category"),
            "c2": pd.Series(list("eee"), dtype="category"),
            "i2": pd.Series([4, 5, 6]),
        }
    )
    df2 = df2.set_index(["c1", "c2"])
    result = pd.concat([df1, df2])
    expected = pd.DataFrame(
        {
            "c1": pd.Series(list("abcabc"), dtype="category"),
            "c2": pd.Series(list("eeeeee"), dtype="category"),
            "i2": pd.Series([1, 2, 3, 4, 5, 6]),
        }
    )
    expected = expected.set_index(["c1", "c2"])
    tm.assert_frame_equal(result, expected)


def test_concat_ea_upcast():
    # GH#54848
    df1 = pd.DataFrame(["a"], dtype="string")
    df2 = pd.DataFrame([1], dtype="Int64")
    result = pd.concat([df1, df2])
    expected = pd.DataFrame(["a", 1], index=[0, 0])
    tm.assert_frame_equal(result, expected)


def test_concat_none_with_timezone_timestamp():
    # GH#52093
    df1 = pd.DataFrame([{"A": None}])
    df2 = pd.DataFrame([{"A": pd.Timestamp("1990-12-20 00:00:00+00:00")}])
    result = pd.concat([df1, df2], ignore_index=True)
    expected = pd.DataFrame(
        {"A": [None, pd.Timestamp("1990-12-20 00:00:00+00:00")]}, dtype=object
    )
    tm.assert_frame_equal(result, expected)


def test_concat_with_series_and_frame_returns_rangeindex_columns():
    ser = pd.Series([0])
    df = pd.DataFrame([1, 2])
    result = pd.concat([ser, df])
    expected = pd.DataFrame([0, 1, 2], index=[0, 0, 1])
    tm.assert_frame_equal(result, expected, check_column_type=True)


def test_concat_with_moot_ignore_index_and_keys():
    df1 = pd.DataFrame([[0]])
    df2 = pd.DataFrame([[42]])

    ignore_index = True
    keys = ["df1", "df2"]
    msg = f"Cannot set {ignore_index=} and specify keys. Either should be used."
    with pytest.raises(ValueError, match=msg):
        pd.concat([df1, df2], keys=keys, ignore_index=ignore_index)


@pytest.mark.parametrize(
    "inputs, ignore_index, axis, expected",
    [
        # Concatenating DataFrame and named Series without ignore_index
        (
            [pd.DataFrame({"a": [0, 1], "b": [2, 3]}), pd.Series([4, 5], name="c")],
            False,
            0,
            pd.DataFrame(
                {
                    "a": [0, 1, None, None],
                    "b": [2, 3, None, None],
                    "c": [None, None, 4, 5],
                },
                index=[0, 1, 0, 1],
            ),
        ),
        # Concatenating DataFrame and named Series with ignore_index
        (
            [pd.DataFrame({"a": [0, 1], "b": [2, 3]}), pd.Series([4, 5], name="c")],
            True,
            0,
            pd.DataFrame(
                {
                    "a": [0, 1, None, None],
                    "b": [2, 3, None, None],
                    "c": [None, None, 4, 5],
                },
                index=[0, 1, 2, 3],
            ),
        ),
        # Concatenating DataFrame and unnamed Series along columns
        (
            [
                pd.DataFrame({"a": [0, 1], "b": [2, 3]}),
                pd.Series([4, 5]),
                pd.Series([4, 5]),
            ],
            False,
            1,
            pd.DataFrame(
                {"a": [0, 1], "b": [2, 3], 0: [4, 5], 1: [4, 5]}, index=[0, 1]
            ),
        ),
        # Concatenating DataFrame and unnamed Series along columns with ignore_index
        (
            [
                pd.DataFrame({"a": [0, 1], "b": [2, 3]}),
                pd.Series([4, 5]),
                pd.Series([4, 5]),
            ],
            True,
            1,
            pd.DataFrame({0: [0, 1], 1: [2, 3], 2: [4, 5], 3: [4, 5]}, index=[0, 1]),
        ),
    ],
)
def test_concat_of_series_and_frame(inputs, ignore_index, axis, expected):
    # GH #60723 and #56257
    result = pd.concat(inputs, ignore_index=ignore_index, axis=axis)
    tm.assert_frame_equal(result, expected)


def test_concat_keys_overlapping_intervalindex_level():
    # GH#64825 an overlapping IntervalIndex passed as levels is not
    #  unique-as-index, but its values are distinct, so concat should resolve it
    value_index = pd.IntervalIndex.from_tuples([(0.0, 1.0), (1.0, 2.0)], name="foo")
    level_index = pd.IntervalIndex.from_tuples([(0.0, 10.0), (0.0, 20.0)], name="bar")
    values = [
        pd.Series([1.0, 3.0], name=pd.Interval(0.0, 10.0), index=value_index),
        pd.Series([5.0, 7.0], name=pd.Interval(0.0, 20.0), index=value_index),
    ]
    result = pd.concat(values, keys=level_index, levels=[level_index], names=["bar"])
    expected = pd.Series(
        [1.0, 3.0, 5.0, 7.0],
        index=pd.MultiIndex.from_product([level_index, value_index]),
    )
    tm.assert_series_equal(result, expected)


def test_concat_mismatched_nlevels_with_keys_gh25413():
    # GH#25413
    df1 = pd.DataFrame(np.arange(12).reshape(4, 3), columns=["A", "B", "C"])
    df2 = pd.DataFrame(
        np.arange(12).reshape(4, 3),
        columns=["A", "B", "C"],
        index=pd.MultiIndex.from_product([["a", "b"], [1, 2]]),
    )
    msg = "Cannot concat indices that do not have the same number of levels"
    with pytest.raises(ValueError, match=msg):
        pd.concat([df1, df2], keys=[1, 2])


def test_concat_keys_overlapping_intervalindex_value_index():
    # GH#64825 the per-Series index is itself an overlapping IntervalIndex
    value_index = pd.IntervalIndex.from_tuples([(0.0, 10.0), (0.0, 20.0)], name="foo")
    level_index = pd.IntervalIndex.from_tuples([(0.0, 10.0), (0.0, 20.0)], name="bar")
    values = [
        pd.Series([1.0, 3.0], name=pd.Interval(0.0, 10.0), index=value_index),
        pd.Series([5.0, 7.0], name=pd.Interval(0.0, 20.0), index=value_index),
    ]
    result = pd.concat(values, keys=level_index, levels=[level_index], names=["bar"])
    expected = pd.Series(
        [1.0, 3.0, 5.0, 7.0],
        index=pd.MultiIndex.from_product([level_index, value_index]),
    )
    tm.assert_series_equal(result, expected)


def test_concat_keys_duplicate_index_equal_lengths():
    # GH#20565 concat with keys where each frame shares an identical duplicate
    #  index must dedupe the shared level rather than repeat it per position
    df1 = pd.DataFrame(np.arange(6).reshape(3, 2), columns=["A", "B"], index=["Z1"] * 3)
    df2 = pd.DataFrame(np.arange(6).reshape(3, 2), columns=["A", "B"], index=["Z1"] * 3)
    result = pd.concat([df1, df2], keys=["Key1", "Key2"], names=["KEY", "ID"])
    expected_index = pd.MultiIndex(
        levels=[["Key1", "Key2"], ["Z1"]],
        codes=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0]],
        names=["KEY", "ID"],
    )
    tm.assert_index_equal(result.index, expected_index)


def test_concat_tuple_index_series_axis1_not_multiindex():
    # GH#24783 concatenating Series whose index labels are tuples must keep a
    #  flat Index of those tuples, not explode them into a MultiIndex
    s1 = pd.Series([1.0, 2.0], index=[("a", "b"), ("x", "y", "z")], name="s1")
    s2 = pd.Series([3.0, 4.0], index=[("a", "b"), ("j", "k", "l")], name="s2")
    result = pd.concat([s1, s2], axis=1, sort=False)
    assert result.index.nlevels == 1
    index_arr = np.empty(3, dtype=object)
    index_arr[:] = [("a", "b"), ("x", "y", "z"), ("j", "k", "l")]
    expected = pd.DataFrame(
        {"s1": [1.0, 2.0, np.nan], "s2": [3.0, np.nan, 4.0]},
        index=pd.Index(index_arr),
    )
    tm.assert_frame_equal(result, expected)
