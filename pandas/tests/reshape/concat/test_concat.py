from collections import abc, deque
from decimal import Decimal
from warnings import catch_warnings

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, concat, date_range
import pandas._testing as tm
from pandas.core.arrays import SparseArray
from pandas.core.construction import create_series_with_explicit_dtype
from pandas.tests.extension.decimal import to_decimal


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param


class TestConcatenate:
    def test_concat_copy(self):
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randint(0, 10, size=4).reshape(4, 1))
        df3 = DataFrame({5: "foo"}, index=range(4))

        # These are actual copies.
        result = concat([df, df2, df3], axis=1, copy=True)

        for b in result._mgr.blocks:
            assert b.values.base is None

        # These are the same.
        result = concat([df, df2, df3], axis=1, copy=False)

        for b in result._mgr.blocks:
            if b.is_float:
                assert b.values.base is df._mgr.blocks[0].values.base
            elif b.is_integer:
                assert b.values.base is df2._mgr.blocks[0].values.base
            elif b.is_object:
                assert b.values.base is not None

        # Float block was consolidated.
        df4 = DataFrame(np.random.randn(4, 1))
        result = concat([df, df2, df3, df4], axis=1, copy=False)
        for b in result._mgr.blocks:
            if b.is_float:
                assert b.values.base is None
            elif b.is_integer:
                assert b.values.base is df2._mgr.blocks[0].values.base
            elif b.is_object:
                assert b.values.base is not None

    def test_concat_with_group_keys(self):
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        # axis=0
        df = DataFrame(np.random.randn(3, 4))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1])
        exp_index = MultiIndex.from_arrays(
            [[0, 0, 0, 1, 1, 1, 1], [0, 1, 2, 0, 1, 2, 3]]
        )
        expected = DataFrame(np.r_[df.values, df2.values], index=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1])
        exp_index2 = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]])
        expected = DataFrame(np.r_[df.values, df.values], index=exp_index2)
        tm.assert_frame_equal(result, expected)

        # axis=1
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df2.values], columns=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df.values], columns=exp_index2)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_specific_levels(self):
        df = DataFrame(np.random.randn(10, 4))
        pieces = [df.iloc[:, [0, 1]], df.iloc[:, [2]], df.iloc[:, [3]]]
        level = ["three", "two", "one", "zero"]
        result = concat(
            pieces,
            axis=1,
            keys=["one", "two", "three"],
            levels=[level],
            names=["group_key"],
        )

        tm.assert_index_equal(result.columns.levels[0], Index(level, name="group_key"))
        tm.assert_index_equal(result.columns.levels[1], Index([0, 1, 2, 3]))

        assert result.columns.names == ["group_key", None]

    def test_concat_dataframe_keys_bug(self, sort):
        t1 = DataFrame(
            {"value": Series([1, 2, 3], index=Index(["a", "b", "c"], name="id"))}
        )
        t2 = DataFrame({"value": Series([7, 8], index=Index(["a", "b"], name="id"))})

        # it works
        result = concat([t1, t2], axis=1, keys=["t1", "t2"], sort=sort)
        assert list(result.columns) == [("t1", "value"), ("t2", "value")]

    @pytest.mark.parametrize("mapping", ["mapping", "dict"])
    def test_concat_mapping(self, mapping, non_dict_mapping_subclass):
        constructor = dict if mapping == "dict" else non_dict_mapping_subclass
        frames = constructor(
            {
                "foo": DataFrame(np.random.randn(4, 3)),
                "bar": DataFrame(np.random.randn(4, 3)),
                "baz": DataFrame(np.random.randn(4, 3)),
                "qux": DataFrame(np.random.randn(4, 3)),
            }
        )

        sorted_keys = list(frames.keys())

        result = concat(frames)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys)
        tm.assert_frame_equal(result, expected)

        result = concat(frames, axis=1)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys, axis=1)
        tm.assert_frame_equal(result, expected)

        keys = ["baz", "foo", "bar"]
        result = concat(frames, keys=keys)
        expected = concat([frames[k] for k in keys], keys=keys)
        tm.assert_frame_equal(result, expected)

    def test_concat_ignore_index(self, sort):
        frame1 = DataFrame(
            {"test1": ["a", "b", "c"], "test2": [1, 2, 3], "test3": [4.5, 3.2, 1.2]}
        )
        frame2 = DataFrame({"test3": [5.2, 2.2, 4.3]})
        frame1.index = Index(["x", "y", "z"])
        frame2.index = Index(["x", "y", "q"])

        v1 = concat([frame1, frame2], axis=1, ignore_index=True, sort=sort)

        nan = np.nan
        expected = DataFrame(
            [
                [nan, nan, nan, 4.3],
                ["a", 1, 4.5, 5.2],
                ["b", 2, 3.2, 2.2],
                ["c", 3, 1.2, nan],
            ],
            index=Index(["q", "x", "y", "z"]),
        )
        if not sort:
            expected = expected.loc[["x", "y", "z", "q"]]

        tm.assert_frame_equal(v1, expected)

    @pytest.mark.parametrize(
        "name_in1,name_in2,name_in3,name_out",
        [
            ("idx", "idx", "idx", "idx"),
            ("idx", "idx", None, None),
            ("idx", None, None, None),
            ("idx1", "idx2", None, None),
            ("idx1", "idx1", "idx2", None),
            ("idx1", "idx2", "idx3", None),
            (None, None, None, None),
        ],
    )
    def test_concat_same_index_names(self, name_in1, name_in2, name_in3, name_out):
        # GH13475
        indices = [
            Index(["a", "b", "c"], name=name_in1),
            Index(["b", "c", "d"], name=name_in2),
            Index(["c", "d", "e"], name=name_in3),
        ]
        frames = [
            DataFrame({c: [0, 1, 2]}, index=i) for i, c in zip(indices, ["x", "y", "z"])
        ]
        result = pd.concat(frames, axis=1)

        exp_ind = Index(["a", "b", "c", "d", "e"], name=name_out)
        expected = DataFrame(
            {
                "x": [0, 1, 2, np.nan, np.nan],
                "y": [np.nan, 0, 1, 2, np.nan],
                "z": [np.nan, np.nan, 0, 1, 2],
            },
            index=exp_ind,
        )

        tm.assert_frame_equal(result, expected)

    def test_concat_multiindex_with_keys(self):
        index = MultiIndex(
            levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
            codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
            names=["first", "second"],
        )
        frame = DataFrame(
            np.random.randn(10, 3),
            index=index,
            columns=Index(["A", "B", "C"], name="exp"),
        )
        result = concat([frame, frame], keys=[0, 1], names=["iteration"])

        assert result.index.names == ("iteration",) + index.names
        tm.assert_frame_equal(result.loc[0], frame)
        tm.assert_frame_equal(result.loc[1], frame)
        assert result.index.nlevels == 3

    def test_concat_multiindex_with_none_in_index_names(self):
        # GH 15787
        index = pd.MultiIndex.from_product([[1], range(5)], names=["level1", None])
        df = DataFrame({"col": range(5)}, index=index, dtype=np.int32)

        result = concat([df, df], keys=[1, 2], names=["level2"])
        index = pd.MultiIndex.from_product(
            [[1, 2], [1], range(5)], names=["level2", "level1", None]
        )
        expected = DataFrame({"col": list(range(5)) * 2}, index=index, dtype=np.int32)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df[:2]], keys=[1, 2], names=["level2"])
        level2 = [1] * 5 + [2] * 2
        level1 = [1] * 7
        no_name = list(range(5)) + list(range(2))
        tuples = list(zip(level2, level1, no_name))
        index = pd.MultiIndex.from_tuples(tuples, names=["level2", "level1", None])
        expected = DataFrame({"col": no_name}, index=index, dtype=np.int32)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_and_levels(self):
        df = DataFrame(np.random.randn(1, 3))
        df2 = DataFrame(np.random.randn(1, 4))

        levels = [["foo", "baz"], ["one", "two"]]
        names = ["first", "second"]
        result = concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            levels=levels,
            names=names,
        )
        expected = concat([df, df2, df, df2])
        exp_index = MultiIndex(
            levels=levels + [[0]],
            codes=[[0, 0, 1, 1], [0, 1, 0, 1], [0, 0, 0, 0]],
            names=names + [None],
        )
        expected.index = exp_index

        tm.assert_frame_equal(result, expected)

        # no names
        result = concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            levels=levels,
        )
        assert result.index.names == (None,) * 3

        # no levels
        result = concat(
            [df, df2, df, df2],
            keys=[("foo", "one"), ("foo", "two"), ("baz", "one"), ("baz", "two")],
            names=["first", "second"],
        )
        assert result.index.names == ("first", "second", None)
        tm.assert_index_equal(
            result.index.levels[0], Index(["baz", "foo"], name="first")
        )

    def test_concat_keys_levels_no_overlap(self):
        # GH #1406
        df = DataFrame(np.random.randn(1, 3), index=["a"])
        df2 = DataFrame(np.random.randn(1, 4), index=["b"])

        msg = "Values not found in passed level"
        with pytest.raises(ValueError, match=msg):
            concat([df, df], keys=["one", "two"], levels=[["foo", "bar", "baz"]])

        msg = "Key one not in level"
        with pytest.raises(ValueError, match=msg):
            concat([df, df2], keys=["one", "two"], levels=[["foo", "bar", "baz"]])

    def test_concat_rename_index(self):
        a = DataFrame(
            np.random.rand(3, 3),
            columns=list("ABC"),
            index=Index(list("abc"), name="index_a"),
        )
        b = DataFrame(
            np.random.rand(3, 3),
            columns=list("ABC"),
            index=Index(list("abc"), name="index_b"),
        )

        result = concat([a, b], keys=["key0", "key1"], names=["lvl0", "lvl1"])

        exp = concat([a, b], keys=["key0", "key1"], names=["lvl0"])
        names = list(exp.index.names)
        names[1] = "lvl1"
        exp.index.set_names(names, inplace=True)

        tm.assert_frame_equal(result, exp)
        assert result.index.names == exp.index.names

    def test_crossed_dtypes_weird_corner(self):
        columns = ["A", "B", "C", "D"]
        df1 = DataFrame(
            {
                "A": np.array([1, 2, 3, 4], dtype="f8"),
                "B": np.array([1, 2, 3, 4], dtype="i8"),
                "C": np.array([1, 2, 3, 4], dtype="f8"),
                "D": np.array([1, 2, 3, 4], dtype="i8"),
            },
            columns=columns,
        )

        df2 = DataFrame(
            {
                "A": np.array([1, 2, 3, 4], dtype="i8"),
                "B": np.array([1, 2, 3, 4], dtype="f8"),
                "C": np.array([1, 2, 3, 4], dtype="i8"),
                "D": np.array([1, 2, 3, 4], dtype="f8"),
            },
            columns=columns,
        )

        appended = df1.append(df2, ignore_index=True)
        expected = DataFrame(
            np.concatenate([df1.values, df2.values], axis=0), columns=columns
        )
        tm.assert_frame_equal(appended, expected)

        df = DataFrame(np.random.randn(1, 3), index=["a"])
        df2 = DataFrame(np.random.randn(1, 4), index=["b"])
        result = concat([df, df2], keys=["one", "two"], names=["first", "second"])
        assert result.index.names == ("first", "second")

    def test_dups_index(self):
        # GH 4771

        # single dtypes
        df = DataFrame(
            np.random.randint(0, 10, size=40).reshape(10, 4),
            columns=["A", "A", "C", "C"],
        )

        result = concat([df, df], axis=1)
        tm.assert_frame_equal(result.iloc[:, :4], df)
        tm.assert_frame_equal(result.iloc[:, 4:], df)

        result = concat([df, df], axis=0)
        tm.assert_frame_equal(result.iloc[:10], df)
        tm.assert_frame_equal(result.iloc[10:], df)

        # multi dtypes
        df = concat(
            [
                DataFrame(np.random.randn(10, 4), columns=["A", "A", "B", "B"]),
                DataFrame(
                    np.random.randint(0, 10, size=20).reshape(10, 2), columns=["A", "C"]
                ),
            ],
            axis=1,
        )

        result = concat([df, df], axis=1)
        tm.assert_frame_equal(result.iloc[:, :6], df)
        tm.assert_frame_equal(result.iloc[:, 6:], df)

        result = concat([df, df], axis=0)
        tm.assert_frame_equal(result.iloc[:10], df)
        tm.assert_frame_equal(result.iloc[10:], df)

        # append
        result = df.iloc[0:8, :].append(df.iloc[8:])
        tm.assert_frame_equal(result, df)

        result = df.iloc[0:8, :].append(df.iloc[8:9]).append(df.iloc[9:10])
        tm.assert_frame_equal(result, df)

        expected = concat([df, df], axis=0)
        result = df.append(df)
        tm.assert_frame_equal(result, expected)

    def test_with_mixed_tuples(self, sort):
        # 10697
        # columns have mixed tuples, so handle properly
        df1 = DataFrame({"A": "foo", ("B", 1): "bar"}, index=range(2))
        df2 = DataFrame({"B": "foo", ("B", 1): "bar"}, index=range(2))

        # it works
        concat([df1, df2], sort=sort)

    def test_concat_mixed_objs(self):

        # concat mixed series/frames
        # G2385

        # axis 1
        index = date_range("01-Jan-2013", periods=10, freq="H")
        arr = np.arange(10, dtype="int64")
        s1 = Series(arr, index=index)
        s2 = Series(arr, index=index)
        df = DataFrame(arr.reshape(-1, 1), index=index)

        expected = DataFrame(
            np.repeat(arr, 2).reshape(-1, 2), index=index, columns=[0, 0]
        )
        result = concat([df, df], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            np.repeat(arr, 2).reshape(-1, 2), index=index, columns=[0, 1]
        )
        result = concat([s1, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=[0, 1, 2]
        )
        result = concat([s1, s2, s1], axis=1)
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            np.repeat(arr, 5).reshape(-1, 5), index=index, columns=[0, 0, 1, 2, 3]
        )
        result = concat([s1, df, s2, s2, s1], axis=1)
        tm.assert_frame_equal(result, expected)

        # with names
        s1.name = "foo"
        expected = DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=["foo", 0, 0]
        )
        result = concat([s1, df, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        s2.name = "bar"
        expected = DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=["foo", 0, "bar"]
        )
        result = concat([s1, df, s2], axis=1)
        tm.assert_frame_equal(result, expected)

        # ignore index
        expected = DataFrame(
            np.repeat(arr, 3).reshape(-1, 3), index=index, columns=[0, 1, 2]
        )
        result = concat([s1, df, s2], axis=1, ignore_index=True)
        tm.assert_frame_equal(result, expected)

        # axis 0
        expected = DataFrame(
            np.tile(arr, 3).reshape(-1, 1), index=index.tolist() * 3, columns=[0]
        )
        result = concat([s1, df, s2])
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(np.tile(arr, 3).reshape(-1, 1), columns=[0])
        result = concat([s1, df, s2], ignore_index=True)
        tm.assert_frame_equal(result, expected)

    def test_dtype_coerceion(self):

        # 12411
        df = DataFrame({"date": [pd.Timestamp("20130101").tz_localize("UTC"), pd.NaT]})

        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 12045
        import datetime

        df = DataFrame(
            {"date": [datetime.datetime(2012, 1, 1), datetime.datetime(1012, 1, 2)]}
        )
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 11594
        df = DataFrame({"text": ["some words"] + [None] * 9})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

    def test_concat_single_with_key(self):
        df = DataFrame(np.random.randn(10, 4))

        result = concat([df], keys=["foo"])
        expected = concat([df, df], keys=["foo", "bar"])
        tm.assert_frame_equal(result, expected[:10])

    def test_concat_exclude_none(self):
        df = DataFrame(np.random.randn(10, 4))

        pieces = [df[:5], None, None, df[5:]]
        result = concat(pieces)
        tm.assert_frame_equal(result, df)
        with pytest.raises(ValueError, match="All objects passed were None"):
            concat([None, None])

    def test_concat_timedelta64_block(self):
        from pandas import to_timedelta

        rng = to_timedelta(np.arange(10), unit="s")

        df = DataFrame({"time": rng})

        result = concat([df, df])
        assert (result.iloc[:10]["time"] == rng).all()
        assert (result.iloc[10:]["time"] == rng).all()

    def test_concat_keys_with_none(self):
        # #1649
        df0 = DataFrame([[10, 20, 30], [10, 20, 30], [10, 20, 30]])

        result = concat(dict(a=None, b=df0, c=df0[:2], d=df0[:1], e=df0))
        expected = concat(dict(b=df0, c=df0[:2], d=df0[:1], e=df0))
        tm.assert_frame_equal(result, expected)

        result = concat(
            [None, df0, df0[:2], df0[:1], df0], keys=["a", "b", "c", "d", "e"]
        )
        expected = concat([df0, df0[:2], df0[:1], df0], keys=["b", "c", "d", "e"])
        tm.assert_frame_equal(result, expected)

    def test_concat_bug_1719(self):
        ts1 = tm.makeTimeSeries()
        ts2 = tm.makeTimeSeries()[::2]

        # to join with union
        # these two are of different length!
        left = concat([ts1, ts2], join="outer", axis=1)
        right = concat([ts2, ts1], join="outer", axis=1)

        assert len(left) == len(right)

    def test_concat_bug_2972(self):
        ts0 = Series(np.zeros(5))
        ts1 = Series(np.ones(5))
        ts0.name = ts1.name = "same name"
        result = concat([ts0, ts1], axis=1)

        expected = DataFrame({0: ts0, 1: ts1})
        expected.columns = ["same name", "same name"]
        tm.assert_frame_equal(result, expected)

    def test_concat_bug_3602(self):

        # GH 3602, duplicate columns
        df1 = DataFrame(
            {
                "firmNo": [0, 0, 0, 0],
                "prc": [6, 6, 6, 6],
                "stringvar": ["rrr", "rrr", "rrr", "rrr"],
            }
        )
        df2 = DataFrame(
            {"C": [9, 10, 11, 12], "misc": [1, 2, 3, 4], "prc": [6, 6, 6, 6]}
        )
        expected = DataFrame(
            [
                [0, 6, "rrr", 9, 1, 6],
                [0, 6, "rrr", 10, 2, 6],
                [0, 6, "rrr", 11, 3, 6],
                [0, 6, "rrr", 12, 4, 6],
            ]
        )
        expected.columns = ["firmNo", "prc", "stringvar", "C", "misc", "prc"]

        result = concat([df1, df2], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_concat_iterables(self):
        # GH8645 check concat works with tuples, list, generators, and weird
        # stuff like deque and custom iterables
        df1 = DataFrame([1, 2, 3])
        df2 = DataFrame([4, 5, 6])
        expected = DataFrame([1, 2, 3, 4, 5, 6])
        tm.assert_frame_equal(concat((df1, df2), ignore_index=True), expected)
        tm.assert_frame_equal(concat([df1, df2], ignore_index=True), expected)
        tm.assert_frame_equal(
            concat((df for df in (df1, df2)), ignore_index=True), expected
        )
        tm.assert_frame_equal(concat(deque((df1, df2)), ignore_index=True), expected)

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
            def __iter__(self):
                yield df1
                yield df2

        tm.assert_frame_equal(pd.concat(CustomIterator2(), ignore_index=True), expected)

    def test_default_index(self):
        # is_series and ignore_index
        s1 = Series([1, 2, 3], name="x")
        s2 = Series([4, 5, 6], name="y")
        res = pd.concat([s1, s2], axis=1, ignore_index=True)
        assert isinstance(res.columns, pd.RangeIndex)
        exp = DataFrame([[1, 4], [2, 5], [3, 6]])
        # use check_index_type=True to check the result have
        # RangeIndex (default index)
        tm.assert_frame_equal(res, exp, check_index_type=True, check_column_type=True)

        # is_series and all inputs have no names
        s1 = Series([1, 2, 3])
        s2 = Series([4, 5, 6])
        res = pd.concat([s1, s2], axis=1, ignore_index=False)
        assert isinstance(res.columns, pd.RangeIndex)
        exp = DataFrame([[1, 4], [2, 5], [3, 6]])
        exp.columns = pd.RangeIndex(2)
        tm.assert_frame_equal(res, exp, check_index_type=True, check_column_type=True)

        # is_dataframe and ignore_index
        df1 = DataFrame({"A": [1, 2], "B": [5, 6]})
        df2 = DataFrame({"A": [3, 4], "B": [7, 8]})

        res = pd.concat([df1, df2], axis=0, ignore_index=True)
        exp = DataFrame([[1, 5], [2, 6], [3, 7], [4, 8]], columns=["A", "B"])
        tm.assert_frame_equal(res, exp, check_index_type=True, check_column_type=True)

        res = pd.concat([df1, df2], axis=1, ignore_index=True)
        exp = DataFrame([[1, 5, 3, 7], [2, 6, 4, 8]])
        tm.assert_frame_equal(res, exp, check_index_type=True, check_column_type=True)

    def test_concat_multiindex_rangeindex(self):
        # GH13542
        # when multi-index levels are RangeIndex objects
        # there is a bug in concat with objects of len 1

        df = DataFrame(np.random.randn(9, 2))
        df.index = MultiIndex(
            levels=[pd.RangeIndex(3), pd.RangeIndex(3)],
            codes=[np.repeat(np.arange(3), 3), np.tile(np.arange(3), 3)],
        )

        res = concat([df.iloc[[2, 3, 4], :], df.iloc[[5], :]])
        exp = df.iloc[[2, 3, 4, 5], :]
        tm.assert_frame_equal(res, exp)

    def test_concat_multiindex_dfs_with_deepcopy(self):
        # GH 9967
        from copy import deepcopy

        example_multiindex1 = pd.MultiIndex.from_product([["a"], ["b"]])
        example_dataframe1 = DataFrame([0], index=example_multiindex1)

        example_multiindex2 = pd.MultiIndex.from_product([["a"], ["c"]])
        example_dataframe2 = DataFrame([1], index=example_multiindex2)

        example_dict = {"s1": example_dataframe1, "s2": example_dataframe2}
        expected_index = pd.MultiIndex(
            levels=[["s1", "s2"], ["a"], ["b", "c"]],
            codes=[[0, 1], [0, 0], [0, 1]],
            names=["testname", None, None],
        )
        expected = DataFrame([[0], [1]], index=expected_index)
        result_copy = pd.concat(deepcopy(example_dict), names=["testname"])
        tm.assert_frame_equal(result_copy, expected)
        result_no_copy = pd.concat(example_dict, names=["testname"])
        tm.assert_frame_equal(result_no_copy, expected)

    def test_concat_order(self):
        # GH 17344
        dfs = [DataFrame(index=range(3), columns=["a", 1, None])]
        dfs += [DataFrame(index=range(3), columns=[None, 1, "a"]) for i in range(100)]

        result = pd.concat(dfs, sort=True).columns
        expected = dfs[0].columns
        tm.assert_index_equal(result, expected)

    def test_concat_different_extension_dtypes_upcasts(self):
        a = Series(pd.core.arrays.integer_array([1, 2]))
        b = Series(to_decimal([1, 2]))

        result = pd.concat([a, b], ignore_index=True)
        expected = Series([1, 2, Decimal(1), Decimal(2)], dtype=object)
        tm.assert_series_equal(result, expected)

    def test_concat_odered_dict(self):
        # GH 21510
        expected = pd.concat(
            [Series(range(3)), Series(range(4))], keys=["First", "Another"]
        )
        result = pd.concat(
            dict([("First", Series(range(3))), ("Another", Series(range(4)))])
        )
        tm.assert_series_equal(result, expected)

    def test_concat_empty_dataframe_dtypes(self):
        df = DataFrame(columns=list("abc"))
        df["a"] = df["a"].astype(np.bool_)
        df["b"] = df["b"].astype(np.int32)
        df["c"] = df["c"].astype(np.float64)

        result = pd.concat([df, df])
        assert result["a"].dtype == np.bool_
        assert result["b"].dtype == np.int32
        assert result["c"].dtype == np.float64

        result = pd.concat([df, df.astype(np.float64)])
        assert result["a"].dtype == np.object_
        assert result["b"].dtype == np.float64
        assert result["c"].dtype == np.float64


@pytest.mark.parametrize("pdt", [Series, pd.DataFrame])
@pytest.mark.parametrize("dt", np.sctypes["float"])
def test_concat_no_unnecessary_upcast(dt, pdt):
    # GH 13247
    dims = pdt(dtype=object).ndim

    dfs = [
        pdt(np.array([1], dtype=dt, ndmin=dims)),
        pdt(np.array([np.nan], dtype=dt, ndmin=dims)),
        pdt(np.array([5], dtype=dt, ndmin=dims)),
    ]
    x = pd.concat(dfs)
    assert x.values.dtype == dt


@pytest.mark.parametrize("pdt", [create_series_with_explicit_dtype, pd.DataFrame])
@pytest.mark.parametrize("dt", np.sctypes["int"])
def test_concat_will_upcast(dt, pdt):
    with catch_warnings(record=True):
        dims = pdt().ndim
        dfs = [
            pdt(np.array([1], dtype=dt, ndmin=dims)),
            pdt(np.array([np.nan], ndmin=dims)),
            pdt(np.array([5], dtype=dt, ndmin=dims)),
        ]
        x = pd.concat(dfs)
        assert x.values.dtype == "float64"


def test_concat_empty_and_non_empty_frame_regression():
    # GH 18178 regression test
    df1 = DataFrame({"foo": [1]})
    df2 = DataFrame({"foo": []})
    expected = DataFrame({"foo": [1.0]})
    result = pd.concat([df1, df2])
    tm.assert_frame_equal(result, expected)


def test_concat_sparse():
    # GH 23557
    a = Series(SparseArray([0, 1, 2]))
    expected = DataFrame(data=[[0, 0], [1, 1], [2, 2]]).astype(
        pd.SparseDtype(np.int64, 0)
    )
    result = pd.concat([a, a], axis=1)
    tm.assert_frame_equal(result, expected)


def test_concat_dense_sparse():
    # GH 30668
    a = Series(pd.arrays.SparseArray([1, None]), dtype=float)
    b = Series([1], dtype=float)
    expected = Series(data=[1, None, 1], index=[0, 1, 0]).astype(
        pd.SparseDtype(np.float64, None)
    )
    result = pd.concat([a, b], axis=0)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("test_series", [True, False])
def test_concat_copy_index(test_series, axis):
    # GH 29879
    if test_series:
        ser = Series([1, 2])
        comb = concat([ser, ser], axis=axis, copy=True)
        assert comb.index is not ser.index
    else:
        df = DataFrame([[1, 2], [3, 4]], columns=["a", "b"])
        comb = concat([df, df], axis=axis, copy=True)
        assert comb.index is not df.index
        assert comb.columns is not df.columns


@pytest.mark.parametrize("keys", [["e", "f", "f"], ["f", "e", "f"]])
def test_duplicate_keys(keys):
    # GH 33654
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    s1 = Series([7, 8, 9], name="c")
    s2 = Series([10, 11, 12], name="d")
    result = concat([df, s1, s2], axis=1, keys=keys)
    expected_values = [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
    expected_columns = pd.MultiIndex.from_tuples(
        [(keys[0], "a"), (keys[0], "b"), (keys[1], "c"), (keys[2], "d")]
    )
    expected = DataFrame(expected_values, columns=expected_columns)
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

    result = concat([obj, obj])
    assert isinstance(result, type(obj))


def test_concat_frame_axis0_extension_dtypes():
    # preserve extension dtype (through common_dtype mechanism)
    df1 = DataFrame({"a": pd.array([1, 2, 3], dtype="Int64")})
    df2 = DataFrame({"a": np.array([4, 5, 6])})

    result = pd.concat([df1, df2], ignore_index=True)
    expected = DataFrame({"a": [1, 2, 3, 4, 5, 6]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)

    result = pd.concat([df2, df1], ignore_index=True)
    expected = DataFrame({"a": [4, 5, 6, 1, 2, 3]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)


def test_concat_preserves_extension_int64_dtype():
    # GH 24768
    df_a = DataFrame({"a": [-1]}, dtype="Int64")
    df_b = DataFrame({"b": [1]}, dtype="Int64")
    result = pd.concat([df_a, df_b], ignore_index=True)
    expected = DataFrame({"a": [-1, None], "b": [None, 1]}, dtype="Int64")
    tm.assert_frame_equal(result, expected)
