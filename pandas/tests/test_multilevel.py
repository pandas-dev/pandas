import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import (
    ArrowDtype,
    DataFrame,
    MultiIndex,
    Series,
)
import pandas._testing as tm


class TestMultiLevel:
    def test_reindex_level(self, multiindex_year_month_day_dataframe_random_data):
        # axis=0
        ymd = multiindex_year_month_day_dataframe_random_data

        month_sums = ymd.groupby("month").sum()
        result = month_sums.reindex(ymd.index, level=1)
        expected = ymd.groupby(level="month").transform("sum")

        tm.assert_frame_equal(result, expected)

        # Series
        result = month_sums["A"].reindex(ymd.index, level=1)
        expected = ymd["A"].groupby(level="month").transform("sum")
        tm.assert_series_equal(result, expected, check_names=False)

    def test_reindex(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data

        expected = frame.iloc[[0, 3]]
        reindexed = frame.loc[[("foo", "one"), ("bar", "one")]]
        tm.assert_frame_equal(reindexed, expected)

    def test_reindex_preserve_levels(
        self, multiindex_year_month_day_dataframe_random_data
    ):
        ymd = multiindex_year_month_day_dataframe_random_data

        new_index = ymd.index[::10]
        chunk = ymd.reindex(new_index)
        assert chunk.index.is_(new_index)

        chunk = ymd.loc[new_index]
        assert chunk.index.equals(new_index)

        ymdT = ymd.T
        chunk = ymdT.reindex(columns=new_index)
        assert chunk.columns.is_(new_index)

        chunk = ymdT.loc[:, new_index]
        assert chunk.columns.equals(new_index)

    def test_groupby_transform(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data

        s = frame["A"]
        grouper = s.index.get_level_values(0)

        grouped = s.groupby(grouper, group_keys=False)

        applied = grouped.apply(lambda x: x * 2)
        expected = grouped.transform(lambda x: x * 2)
        result = applied.reindex(expected.index)
        tm.assert_series_equal(result, expected, check_names=False)

    def test_groupby_corner(self):
        midx = MultiIndex(
            levels=[["foo"], ["bar"], ["baz"]],
            codes=[[0], [0], [0]],
            names=["one", "two", "three"],
        )
        df = DataFrame(
            [np.random.default_rng(2).random(4)],
            columns=["a", "b", "c", "d"],
            index=midx,
        )
        # should work
        df.groupby(level="three")

    def test_setitem_with_expansion_multiindex_columns(
        self, multiindex_year_month_day_dataframe_random_data
    ):
        ymd = multiindex_year_month_day_dataframe_random_data

        df = ymd[:5].T
        df[2000, 1, 10] = df[2000, 1, 7]
        assert isinstance(df.columns, MultiIndex)
        assert (df[2000, 1, 10] == df[2000, 1, 7]).all()

    def test_alignment(self):
        x = Series(
            data=[1, 2, 3], index=MultiIndex.from_tuples([("A", 1), ("A", 2), ("B", 3)])
        )

        y = Series(
            data=[4, 5, 6], index=MultiIndex.from_tuples([("Z", 1), ("Z", 2), ("B", 3)])
        )

        res = x - y
        exp_index = x.index.union(y.index)
        exp = x.reindex(exp_index) - y.reindex(exp_index)
        tm.assert_series_equal(res, exp)

        # hit non-monotonic code path
        res = x[::-1] - y[::-1]
        exp_index = x.index.union(y.index)
        exp = x.reindex(exp_index) - y.reindex(exp_index)
        tm.assert_series_equal(res, exp)

    def test_groupby_multilevel(self, multiindex_year_month_day_dataframe_random_data):
        ymd = multiindex_year_month_day_dataframe_random_data

        result = ymd.groupby(level=[0, 1]).mean()

        k1 = ymd.index.get_level_values(0)
        k2 = ymd.index.get_level_values(1)

        expected = ymd.groupby([k1, k2]).mean()

        tm.assert_frame_equal(result, expected)
        assert result.index.names == ymd.index.names[:2]

        result2 = ymd.groupby(level=ymd.index.names[:2]).mean()
        tm.assert_frame_equal(result, result2)

    def test_multilevel_consolidate(self):
        index = MultiIndex.from_tuples(
            [("foo", "one"), ("foo", "two"), ("bar", "one"), ("bar", "two")]
        )
        df = DataFrame(
            np.random.default_rng(2).standard_normal((4, 4)), index=index, columns=index
        )
        df["Totals", ""] = df.sum(axis=1)
        df = df._consolidate()

    def test_level_with_tuples(self):
        index = MultiIndex(
            levels=[[("foo", "bar", 0), ("foo", "baz", 0), ("foo", "qux", 0)], [0, 1]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )

        series = Series(np.random.default_rng(2).standard_normal(6), index=index)
        frame = DataFrame(np.random.default_rng(2).standard_normal((6, 4)), index=index)

        result = series[("foo", "bar", 0)]
        result2 = series.loc[("foo", "bar", 0)]
        expected = series[:2]
        expected.index = expected.index.droplevel(0)
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        with pytest.raises(KeyError, match=r"^\(\('foo', 'bar', 0\), 2\)$"):
            series[("foo", "bar", 0), 2]

        result = frame.loc[("foo", "bar", 0)]
        result2 = frame.xs(("foo", "bar", 0))
        expected = frame[:2]
        expected.index = expected.index.droplevel(0)
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result2, expected)

        index = MultiIndex(
            levels=[[("foo", "bar"), ("foo", "baz"), ("foo", "qux")], [0, 1]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )

        series = Series(np.random.default_rng(2).standard_normal(6), index=index)
        frame = DataFrame(np.random.default_rng(2).standard_normal((6, 4)), index=index)

        result = series[("foo", "bar")]
        result2 = series.loc[("foo", "bar")]
        expected = series[:2]
        expected.index = expected.index.droplevel(0)
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        result = frame.loc[("foo", "bar")]
        result2 = frame.xs(("foo", "bar"))
        expected = frame[:2]
        expected.index = expected.index.droplevel(0)
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result2, expected)

    def test_reindex_level_partial_selection(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data

        result = frame.reindex(["foo", "qux"], level=0)
        expected = frame.iloc[[0, 1, 2, 7, 8, 9]]
        tm.assert_frame_equal(result, expected)

        result = frame.T.reindex(["foo", "qux"], axis=1, level=0)
        tm.assert_frame_equal(result, expected.T)

        result = frame.loc[["foo", "qux"]]
        tm.assert_frame_equal(result, expected)

        result = frame["A"].loc[["foo", "qux"]]
        tm.assert_series_equal(result, expected["A"])

        result = frame.T.loc[:, ["foo", "qux"]]
        tm.assert_frame_equal(result, expected.T)

    @pytest.mark.parametrize("d", [4, "d"])
    def test_empty_frame_groupby_dtypes_consistency(self, d):
        # GH 20888
        group_keys = ["a", "b", "c"]
        df = DataFrame({"a": [1], "b": [2], "c": [3], "d": [d]})

        g = df[df.a == 2].groupby(group_keys)
        result = g.first().index
        expected = MultiIndex(
            levels=[[1], [2], [3]], codes=[[], [], []], names=["a", "b", "c"]
        )

        tm.assert_index_equal(result, expected)

    def test_duplicate_groupby_issues(self):
        idx_tp = [
            ("600809", "20061231"),
            ("600809", "20070331"),
            ("600809", "20070630"),
            ("600809", "20070331"),
        ]
        dt = ["demo", "demo", "demo", "demo"]

        idx = MultiIndex.from_tuples(idx_tp, names=["STK_ID", "RPT_Date"])
        s = Series(dt, index=idx)

        result = s.groupby(s.index).first()
        assert len(result) == 3

    def test_subsets_multiindex_dtype(self):
        # GH 20757
        data = [["x", 1]]
        columns = [("a", "b", np.nan), ("a", "c", 0.0)]
        df = DataFrame(data, columns=MultiIndex.from_tuples(columns))
        expected = df.dtypes.a.b
        result = df.a.b.dtypes
        tm.assert_series_equal(result, expected)

    def test_datetime_object_multiindex(self):
        data_dic = {
            (0, datetime.date(2018, 3, 3)): {"A": 1, "B": 10},
            (0, datetime.date(2018, 3, 4)): {"A": 2, "B": 11},
            (1, datetime.date(2018, 3, 3)): {"A": 3, "B": 12},
            (1, datetime.date(2018, 3, 4)): {"A": 4, "B": 13},
        }
        result = DataFrame.from_dict(data_dic, orient="index")
        data = {"A": [1, 2, 3, 4], "B": [10, 11, 12, 13]}
        index = [
            [0, 0, 1, 1],
            [
                datetime.date(2018, 3, 3),
                datetime.date(2018, 3, 4),
                datetime.date(2018, 3, 3),
                datetime.date(2018, 3, 4),
            ],
        ]
        expected = DataFrame(data=data, index=index)

        tm.assert_frame_equal(result, expected)

    def test_multiindex_with_na(self):
        df = DataFrame(
            [
                ["A", np.nan, 1.23, 4.56],
                ["A", "G", 1.23, 4.56],
                ["A", "D", 9.87, 10.54],
            ],
            columns=["pivot_0", "pivot_1", "col_1", "col_2"],
        ).set_index(["pivot_0", "pivot_1"])

        df.at[("A", "F"), "col_2"] = 0.0

        expected = DataFrame(
            [
                ["A", np.nan, 1.23, 4.56],
                ["A", "G", 1.23, 4.56],
                ["A", "D", 9.87, 10.54],
                ["A", "F", np.nan, 0.0],
            ],
            columns=["pivot_0", "pivot_1", "col_1", "col_2"],
        ).set_index(["pivot_0", "pivot_1"])

        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("na", [None, np.nan])
    def test_multiindex_insert_level_with_na(self, na):
        # GH 59003
        df = DataFrame([0], columns=[["A"], ["B"]])
        df[na, "B"] = 1
        tm.assert_frame_equal(df[na], DataFrame([1], columns=["B"]))

    def test_multiindex_dt_with_nan(self):
        # GH#60388
        df = DataFrame(
            [
                [1, np.nan, 5, np.nan],
                [2, np.nan, 6, np.nan],
                [np.nan, 3, np.nan, 7],
                [np.nan, 4, np.nan, 8],
            ],
            index=Series(["a", "b", "c", "d"], dtype=object, name="sub"),
            columns=MultiIndex.from_product(
                [
                    ["value1", "value2"],
                    [datetime.datetime(2024, 11, 1), datetime.datetime(2024, 11, 2)],
                ],
                names=[None, "Date"],
            ),
        )
        df = df.reset_index()
        result = df[df.columns[0]]
        expected = Series(["a", "b", "c", "d"], name=("sub", np.nan))
        tm.assert_series_equal(result, expected)

    @pytest.mark.filterwarnings("ignore:Passing a BlockManager:DeprecationWarning")
    def test_multiindex_with_pyarrow_categorical(self):
        # GH#53051
        pa = pytest.importorskip("pyarrow")

        df = DataFrame(
            {"string_column": ["A", "B", "C"], "number_column": [1, 2, 3]}
        ).astype(
            {
                "string_column": ArrowDtype(pa.dictionary(pa.int32(), pa.string())),
                "number_column": "float[pyarrow]",
            }
        )

        df = df.set_index(["string_column", "number_column"])

        df_expected = DataFrame(
            index=MultiIndex.from_arrays(
                [["A", "B", "C"], [1, 2, 3]], names=["string_column", "number_column"]
            )
        )
        tm.assert_frame_equal(
            df,
            df_expected,
            check_index_type=False,
            check_column_type=False,
        )


class TestSorted:
    """everything you wanted to test about sorting"""

    def test_sort_non_lexsorted(self):
        # degenerate case where we sort but don't
        # have a satisfying result :<
        # GH 15797
        idx = MultiIndex(
            [["A", "B", "C"], ["c", "b", "a"]], [[0, 1, 2, 0, 1, 2], [0, 2, 1, 1, 0, 2]]
        )

        df = DataFrame({"col": range(len(idx))}, index=idx, dtype="int64")
        assert df.index.is_monotonic_increasing is False

        sorted = df.sort_index()
        assert sorted.index.is_monotonic_increasing is True

        expected = DataFrame(
            {"col": [1, 4, 5, 2]},
            index=MultiIndex.from_tuples(
                [("B", "a"), ("B", "c"), ("C", "a"), ("C", "b")]
            ),
            dtype="int64",
        )
        result = sorted.loc[pd.IndexSlice["B":"C", "a":"c"], :]
        tm.assert_frame_equal(result, expected)
