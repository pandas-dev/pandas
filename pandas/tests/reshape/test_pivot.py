from datetime import date, datetime, timedelta
from itertools import product

import numpy as np
import pytest

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Grouper,
    Index,
    MultiIndex,
    Series,
    concat,
    date_range,
)
import pandas._testing as tm
from pandas.api.types import CategoricalDtype as CDT
from pandas.core.reshape.pivot import crosstab, pivot_table


@pytest.fixture(params=[True, False])
def dropna(request):
    return request.param


@pytest.fixture(params=[([0] * 4, [1] * 4), (range(0, 3), range(1, 4))])
def interval_values(request, closed):
    left, right = request.param
    return Categorical(pd.IntervalIndex.from_arrays(left, right, closed))


class TestPivotTable:
    def setup_method(self, method):
        self.data = DataFrame(
            {
                "A": [
                    "foo",
                    "foo",
                    "foo",
                    "foo",
                    "bar",
                    "bar",
                    "bar",
                    "bar",
                    "foo",
                    "foo",
                    "foo",
                ],
                "B": [
                    "one",
                    "one",
                    "one",
                    "two",
                    "one",
                    "one",
                    "one",
                    "two",
                    "two",
                    "two",
                    "one",
                ],
                "C": [
                    "dull",
                    "dull",
                    "shiny",
                    "dull",
                    "dull",
                    "shiny",
                    "shiny",
                    "dull",
                    "shiny",
                    "shiny",
                    "shiny",
                ],
                "D": np.random.randn(11),
                "E": np.random.randn(11),
                "F": np.random.randn(11),
            }
        )

    def test_pivot_table(self, observed):
        index = ["A", "B"]
        columns = "C"
        table = pivot_table(
            self.data, values="D", index=index, columns=columns, observed=observed
        )

        table2 = self.data.pivot_table(
            values="D", index=index, columns=columns, observed=observed
        )
        tm.assert_frame_equal(table, table2)

        # this works
        pivot_table(self.data, values="D", index=index, observed=observed)

        if len(index) > 1:
            assert table.index.names == tuple(index)
        else:
            assert table.index.name == index[0]

        if len(columns) > 1:
            assert table.columns.names == columns
        else:
            assert table.columns.name == columns[0]

        expected = self.data.groupby(index + [columns])["D"].agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_table_categorical_observed_equal(self, observed):
        # issue #24923
        df = pd.DataFrame(
            {"col1": list("abcde"), "col2": list("fghij"), "col3": [1, 2, 3, 4, 5]}
        )

        expected = df.pivot_table(
            index="col1", values="col3", columns="col2", aggfunc=np.sum, fill_value=0
        )

        expected.index = expected.index.astype("category")
        expected.columns = expected.columns.astype("category")

        df.col1 = df.col1.astype("category")
        df.col2 = df.col2.astype("category")

        result = df.pivot_table(
            index="col1",
            values="col3",
            columns="col2",
            aggfunc=np.sum,
            fill_value=0,
            observed=observed,
        )

        tm.assert_frame_equal(result, expected)

    def test_pivot_table_nocols(self):
        df = DataFrame(
            {"rows": ["a", "b", "c"], "cols": ["x", "y", "z"], "values": [1, 2, 3]}
        )
        rs = df.pivot_table(columns="cols", aggfunc=np.sum)
        xp = df.pivot_table(index="cols", aggfunc=np.sum).T
        tm.assert_frame_equal(rs, xp)

        rs = df.pivot_table(columns="cols", aggfunc={"values": "mean"})
        xp = df.pivot_table(index="cols", aggfunc={"values": "mean"}).T
        tm.assert_frame_equal(rs, xp)

    def test_pivot_table_dropna(self):
        df = DataFrame(
            {
                "amount": {0: 60000, 1: 100000, 2: 50000, 3: 30000},
                "customer": {0: "A", 1: "A", 2: "B", 3: "C"},
                "month": {0: 201307, 1: 201309, 2: 201308, 3: 201310},
                "product": {0: "a", 1: "b", 2: "c", 3: "d"},
                "quantity": {0: 2000000, 1: 500000, 2: 1000000, 3: 1000000},
            }
        )
        pv_col = df.pivot_table(
            "quantity", "month", ["customer", "product"], dropna=False
        )
        pv_ind = df.pivot_table(
            "quantity", ["customer", "product"], "month", dropna=False
        )

        m = MultiIndex.from_tuples(
            [
                ("A", "a"),
                ("A", "b"),
                ("A", "c"),
                ("A", "d"),
                ("B", "a"),
                ("B", "b"),
                ("B", "c"),
                ("B", "d"),
                ("C", "a"),
                ("C", "b"),
                ("C", "c"),
                ("C", "d"),
            ],
            names=["customer", "product"],
        )
        tm.assert_index_equal(pv_col.columns, m)
        tm.assert_index_equal(pv_ind.index, m)

    def test_pivot_table_categorical(self):

        cat1 = Categorical(
            ["a", "a", "b", "b"], categories=["a", "b", "z"], ordered=True
        )
        cat2 = Categorical(
            ["c", "d", "c", "d"], categories=["c", "d", "y"], ordered=True
        )
        df = DataFrame({"A": cat1, "B": cat2, "values": [1, 2, 3, 4]})
        result = pd.pivot_table(df, values="values", index=["A", "B"], dropna=True)

        exp_index = pd.MultiIndex.from_arrays([cat1, cat2], names=["A", "B"])
        expected = DataFrame({"values": [1, 2, 3, 4]}, index=exp_index)
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_dropna_categoricals(self, dropna):
        # GH 15193
        categories = ["a", "b", "c", "d"]

        df = DataFrame(
            {
                "A": ["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                "B": [1, 2, 3, 1, 2, 3, 1, 2, 3],
                "C": range(0, 9),
            }
        )

        df["A"] = df["A"].astype(CDT(categories, ordered=False))
        result = df.pivot_table(index="B", columns="A", values="C", dropna=dropna)
        expected_columns = Series(["a", "b", "c"], name="A")
        expected_columns = expected_columns.astype(CDT(categories, ordered=False))
        expected_index = Series([1, 2, 3], name="B")
        expected = DataFrame(
            [[0, 3, 6], [1, 4, 7], [2, 5, 8]],
            index=expected_index,
            columns=expected_columns,
        )
        if not dropna:
            # add back the non observed to compare
            expected = expected.reindex(columns=Categorical(categories)).astype("float")

        tm.assert_frame_equal(result, expected)

    def test_pivot_with_non_observable_dropna(self, dropna):
        # gh-21133
        df = pd.DataFrame(
            {
                "A": pd.Categorical(
                    [np.nan, "low", "high", "low", "high"],
                    categories=["low", "high"],
                    ordered=True,
                ),
                "B": range(5),
            }
        )

        result = df.pivot_table(index="A", values="B", dropna=dropna)
        expected = pd.DataFrame(
            {"B": [2, 3]},
            index=pd.Index(
                pd.Categorical.from_codes(
                    [0, 1], categories=["low", "high"], ordered=True
                ),
                name="A",
            ),
        )

        tm.assert_frame_equal(result, expected)

        # gh-21378
        df = pd.DataFrame(
            {
                "A": pd.Categorical(
                    ["left", "low", "high", "low", "high"],
                    categories=["low", "high", "left"],
                    ordered=True,
                ),
                "B": range(5),
            }
        )

        result = df.pivot_table(index="A", values="B", dropna=dropna)
        expected = pd.DataFrame(
            {"B": [2, 3, 0]},
            index=pd.Index(
                pd.Categorical.from_codes(
                    [0, 1, 2], categories=["low", "high", "left"], ordered=True
                ),
                name="A",
            ),
        )

        tm.assert_frame_equal(result, expected)

    def test_pivot_with_interval_index(self, interval_values, dropna):
        # GH 25814
        df = DataFrame({"A": interval_values, "B": 1})
        result = df.pivot_table(index="A", values="B", dropna=dropna)
        expected = DataFrame({"B": 1}, index=Index(interval_values.unique(), name="A"))
        tm.assert_frame_equal(result, expected)

    def test_pivot_with_interval_index_margins(self):
        # GH 25815
        ordered_cat = pd.IntervalIndex.from_arrays([0, 0, 1, 1], [1, 1, 2, 2])
        df = DataFrame(
            {
                "A": np.arange(4, 0, -1, dtype=np.intp),
                "B": ["a", "b", "a", "b"],
                "C": pd.Categorical(ordered_cat, ordered=True).sort_values(
                    ascending=False
                ),
            }
        )

        pivot_tab = pd.pivot_table(
            df, index="C", columns="B", values="A", aggfunc="sum", margins=True
        )

        result = pivot_tab["All"]
        expected = Series(
            [3, 7, 10],
            index=Index([pd.Interval(0, 1), pd.Interval(1, 2), "All"], name="C"),
            name="All",
            dtype=np.intp,
        )
        tm.assert_series_equal(result, expected)

    def test_pass_array(self):
        result = self.data.pivot_table("D", index=self.data.A, columns=self.data.C)
        expected = self.data.pivot_table("D", index="A", columns="C")
        tm.assert_frame_equal(result, expected)

    def test_pass_function(self):
        result = self.data.pivot_table("D", index=lambda x: x // 5, columns=self.data.C)
        expected = self.data.pivot_table("D", index=self.data.index // 5, columns="C")
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_multiple(self):
        index = ["A", "B"]
        columns = "C"
        table = pivot_table(self.data, index=index, columns=columns)
        expected = self.data.groupby(index + [columns]).agg(np.mean).unstack()
        tm.assert_frame_equal(table, expected)

    def test_pivot_dtypes(self):

        # can convert dtypes
        f = DataFrame(
            {
                "a": ["cat", "bat", "cat", "bat"],
                "v": [1, 2, 3, 4],
                "i": ["a", "b", "a", "b"],
            }
        )
        assert f.dtypes["v"] == "int64"

        z = pivot_table(
            f, values="v", index=["a"], columns=["i"], fill_value=0, aggfunc=np.sum
        )
        result = z.dtypes
        expected = Series([np.dtype("int64")] * 2, index=Index(list("ab"), name="i"))
        tm.assert_series_equal(result, expected)

        # cannot convert dtypes
        f = DataFrame(
            {
                "a": ["cat", "bat", "cat", "bat"],
                "v": [1.5, 2.5, 3.5, 4.5],
                "i": ["a", "b", "a", "b"],
            }
        )
        assert f.dtypes["v"] == "float64"

        z = pivot_table(
            f, values="v", index=["a"], columns=["i"], fill_value=0, aggfunc=np.mean
        )
        result = z.dtypes
        expected = Series([np.dtype("float64")] * 2, index=Index(list("ab"), name="i"))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "columns,values",
        [
            ("bool1", ["float1", "float2"]),
            ("bool1", ["float1", "float2", "bool1"]),
            ("bool2", ["float1", "float2", "bool1"]),
        ],
    )
    def test_pivot_preserve_dtypes(self, columns, values):
        # GH 7142 regression test
        v = np.arange(5, dtype=np.float64)
        df = DataFrame(
            {"float1": v, "float2": v + 2.0, "bool1": v <= 2, "bool2": v <= 3}
        )

        df_res = df.reset_index().pivot_table(
            index="index", columns=columns, values=values
        )

        result = dict(df_res.dtypes)
        expected = {
            col: np.dtype("O") if col[0].startswith("b") else np.dtype("float64")
            for col in df_res
        }
        assert result == expected

    def test_pivot_no_values(self):
        # GH 14380
        idx = pd.DatetimeIndex(
            ["2011-01-01", "2011-02-01", "2011-01-02", "2011-01-01", "2011-01-02"]
        )
        df = pd.DataFrame({"A": [1, 2, 3, 4, 5]}, index=idx)
        res = df.pivot_table(index=df.index.month, columns=df.index.day)

        exp_columns = pd.MultiIndex.from_tuples([("A", 1), ("A", 2)])
        exp = pd.DataFrame(
            [[2.5, 4.0], [2.0, np.nan]], index=[1, 2], columns=exp_columns
        )
        tm.assert_frame_equal(res, exp)

        df = pd.DataFrame(
            {
                "A": [1, 2, 3, 4, 5],
                "dt": pd.date_range("2011-01-01", freq="D", periods=5),
            },
            index=idx,
        )
        res = df.pivot_table(
            index=df.index.month, columns=pd.Grouper(key="dt", freq="M")
        )
        exp_columns = pd.MultiIndex.from_tuples([("A", pd.Timestamp("2011-01-31"))])
        exp_columns.names = [None, "dt"]
        exp = pd.DataFrame([3.25, 2.0], index=[1, 2], columns=exp_columns)
        tm.assert_frame_equal(res, exp)

        res = df.pivot_table(
            index=pd.Grouper(freq="A"), columns=pd.Grouper(key="dt", freq="M")
        )
        exp = pd.DataFrame(
            [3], index=pd.DatetimeIndex(["2011-12-31"]), columns=exp_columns
        )
        tm.assert_frame_equal(res, exp)

    def test_pivot_multi_values(self):
        result = pivot_table(
            self.data, values=["D", "E"], index="A", columns=["B", "C"], fill_value=0
        )
        expected = pivot_table(
            self.data.drop(["F"], axis=1), index="A", columns=["B", "C"], fill_value=0
        )
        tm.assert_frame_equal(result, expected)

    def test_pivot_multi_functions(self):
        f = lambda func: pivot_table(
            self.data, values=["D", "E"], index=["A", "B"], columns="C", aggfunc=func
        )
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=["mean", "std"], axis=1)
        tm.assert_frame_equal(result, expected)

        # margins not supported??
        f = lambda func: pivot_table(
            self.data,
            values=["D", "E"],
            index=["A", "B"],
            columns="C",
            aggfunc=func,
            margins=True,
        )
        result = f([np.mean, np.std])
        means = f(np.mean)
        stds = f(np.std)
        expected = concat([means, stds], keys=["mean", "std"], axis=1)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_index_with_nan(self, method):
        # GH 3588
        nan = np.nan
        df = DataFrame(
            {
                "a": ["R1", "R2", nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, 17, 20],
            }
        )
        if method:
            result = df.pivot("a", "b", "c")
        else:
            result = pd.pivot(df, "a", "b", "c")
        expected = DataFrame(
            [
                [nan, nan, 17, nan],
                [10, nan, nan, nan],
                [nan, 15, nan, nan],
                [nan, nan, nan, 20],
            ],
            index=Index([nan, "R1", "R2", "R4"], name="a"),
            columns=Index(["C1", "C2", "C3", "C4"], name="b"),
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df.pivot("b", "a", "c"), expected.T)

        # GH9491
        df = DataFrame(
            {
                "a": pd.date_range("2014-02-01", periods=6, freq="D"),
                "c": 100 + np.arange(6),
            }
        )
        df["b"] = df["a"] - pd.Timestamp("2014-02-02")
        df.loc[1, "a"] = df.loc[3, "a"] = nan
        df.loc[1, "b"] = df.loc[4, "b"] = nan

        if method:
            pv = df.pivot("a", "b", "c")
        else:
            pv = pd.pivot(df, "a", "b", "c")
        assert pv.notna().values.sum() == len(df)

        for _, row in df.iterrows():
            assert pv.loc[row["a"], row["b"]] == row["c"]

        if method:
            result = df.pivot("b", "a", "c")
        else:
            result = pd.pivot(df, "b", "a", "c")
        tm.assert_frame_equal(result, pv.T)

    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_with_tz(self, method):
        # GH 5878
        df = DataFrame(
            {
                "dt1": [
                    datetime(2013, 1, 1, 9, 0),
                    datetime(2013, 1, 2, 9, 0),
                    datetime(2013, 1, 1, 9, 0),
                    datetime(2013, 1, 2, 9, 0),
                ],
                "dt2": [
                    datetime(2014, 1, 1, 9, 0),
                    datetime(2014, 1, 1, 9, 0),
                    datetime(2014, 1, 2, 9, 0),
                    datetime(2014, 1, 2, 9, 0),
                ],
                "data1": np.arange(4, dtype="int64"),
                "data2": np.arange(4, dtype="int64"),
            }
        )

        df["dt1"] = df["dt1"].apply(lambda d: pd.Timestamp(d, tz="US/Pacific"))
        df["dt2"] = df["dt2"].apply(lambda d: pd.Timestamp(d, tz="Asia/Tokyo"))

        exp_col1 = Index(["data1", "data1", "data2", "data2"])
        exp_col2 = pd.DatetimeIndex(
            ["2014/01/01 09:00", "2014/01/02 09:00"] * 2, name="dt2", tz="Asia/Tokyo"
        )
        exp_col = pd.MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame(
            [[0, 2, 0, 2], [1, 3, 1, 3]],
            index=pd.DatetimeIndex(
                ["2013/01/01 09:00", "2013/01/02 09:00"], name="dt1", tz="US/Pacific"
            ),
            columns=exp_col,
        )

        if method:
            pv = df.pivot(index="dt1", columns="dt2")
        else:
            pv = pd.pivot(df, index="dt1", columns="dt2")
        tm.assert_frame_equal(pv, expected)

        expected = DataFrame(
            [[0, 2], [1, 3]],
            index=pd.DatetimeIndex(
                ["2013/01/01 09:00", "2013/01/02 09:00"], name="dt1", tz="US/Pacific"
            ),
            columns=pd.DatetimeIndex(
                ["2014/01/01 09:00", "2014/01/02 09:00"], name="dt2", tz="Asia/Tokyo"
            ),
        )

        if method:
            pv = df.pivot(index="dt1", columns="dt2", values="data1")
        else:
            pv = pd.pivot(df, index="dt1", columns="dt2", values="data1")
        tm.assert_frame_equal(pv, expected)

    def test_pivot_tz_in_values(self):
        # GH 14948
        df = pd.DataFrame(
            [
                {
                    "uid": "aa",
                    "ts": pd.Timestamp("2016-08-12 13:00:00-0700", tz="US/Pacific"),
                },
                {
                    "uid": "aa",
                    "ts": pd.Timestamp("2016-08-12 08:00:00-0700", tz="US/Pacific"),
                },
                {
                    "uid": "aa",
                    "ts": pd.Timestamp("2016-08-12 14:00:00-0700", tz="US/Pacific"),
                },
                {
                    "uid": "aa",
                    "ts": pd.Timestamp("2016-08-25 11:00:00-0700", tz="US/Pacific"),
                },
                {
                    "uid": "aa",
                    "ts": pd.Timestamp("2016-08-25 13:00:00-0700", tz="US/Pacific"),
                },
            ]
        )

        df = df.set_index("ts").reset_index()
        mins = df.ts.map(lambda x: x.replace(hour=0, minute=0, second=0, microsecond=0))

        result = pd.pivot_table(
            df.set_index("ts").reset_index(),
            values="ts",
            index=["uid"],
            columns=[mins],
            aggfunc=np.min,
        )
        expected = pd.DataFrame(
            [
                [
                    pd.Timestamp("2016-08-12 08:00:00-0700", tz="US/Pacific"),
                    pd.Timestamp("2016-08-25 11:00:00-0700", tz="US/Pacific"),
                ]
            ],
            index=pd.Index(["aa"], name="uid"),
            columns=pd.DatetimeIndex(
                [
                    pd.Timestamp("2016-08-12 00:00:00", tz="US/Pacific"),
                    pd.Timestamp("2016-08-25 00:00:00", tz="US/Pacific"),
                ],
                name="ts",
            ),
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_periods(self, method):
        df = DataFrame(
            {
                "p1": [
                    pd.Period("2013-01-01", "D"),
                    pd.Period("2013-01-02", "D"),
                    pd.Period("2013-01-01", "D"),
                    pd.Period("2013-01-02", "D"),
                ],
                "p2": [
                    pd.Period("2013-01", "M"),
                    pd.Period("2013-01", "M"),
                    pd.Period("2013-02", "M"),
                    pd.Period("2013-02", "M"),
                ],
                "data1": np.arange(4, dtype="int64"),
                "data2": np.arange(4, dtype="int64"),
            }
        )

        exp_col1 = Index(["data1", "data1", "data2", "data2"])
        exp_col2 = pd.PeriodIndex(["2013-01", "2013-02"] * 2, name="p2", freq="M")
        exp_col = pd.MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame(
            [[0, 2, 0, 2], [1, 3, 1, 3]],
            index=pd.PeriodIndex(["2013-01-01", "2013-01-02"], name="p1", freq="D"),
            columns=exp_col,
        )
        if method:
            pv = df.pivot(index="p1", columns="p2")
        else:
            pv = pd.pivot(df, index="p1", columns="p2")
        tm.assert_frame_equal(pv, expected)

        expected = DataFrame(
            [[0, 2], [1, 3]],
            index=pd.PeriodIndex(["2013-01-01", "2013-01-02"], name="p1", freq="D"),
            columns=pd.PeriodIndex(["2013-01", "2013-02"], name="p2", freq="M"),
        )
        if method:
            pv = df.pivot(index="p1", columns="p2", values="data1")
        else:
            pv = pd.pivot(df, index="p1", columns="p2", values="data1")
        tm.assert_frame_equal(pv, expected)

    def test_pivot_periods_with_margins(self):
        # GH 28323
        df = DataFrame(
            {
                "a": [1, 1, 2, 2],
                "b": [
                    pd.Period("2019Q1"),
                    pd.Period("2019Q2"),
                    pd.Period("2019Q1"),
                    pd.Period("2019Q2"),
                ],
                "x": 1.0,
            }
        )

        expected = DataFrame(
            data=1.0,
            index=pd.Index([1, 2, "All"], name="a"),
            columns=pd.Index(
                [pd.Period("2019Q1"), pd.Period("2019Q2"), "All"], name="b"
            ),
        )

        result = df.pivot_table(index="a", columns="b", values="x", margins=True)
        tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize(
        "values",
        [
            ["baz", "zoo"],
            np.array(["baz", "zoo"]),
            pd.Series(["baz", "zoo"]),
            pd.Index(["baz", "zoo"]),
        ],
    )
    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_with_list_like_values(self, values, method):
        # issue #17160
        df = pd.DataFrame(
            {
                "foo": ["one", "one", "one", "two", "two", "two"],
                "bar": ["A", "B", "C", "A", "B", "C"],
                "baz": [1, 2, 3, 4, 5, 6],
                "zoo": ["x", "y", "z", "q", "w", "t"],
            }
        )

        if method:
            result = df.pivot(index="foo", columns="bar", values=values)
        else:
            result = pd.pivot(df, index="foo", columns="bar", values=values)

        data = [[1, 2, 3, "x", "y", "z"], [4, 5, 6, "q", "w", "t"]]
        index = Index(data=["one", "two"], name="foo")
        columns = MultiIndex(
            levels=[["baz", "zoo"], ["A", "B", "C"]],
            codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]],
            names=[None, "bar"],
        )
        expected = DataFrame(data=data, index=index, columns=columns, dtype="object")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "values",
        [
            ["bar", "baz"],
            np.array(["bar", "baz"]),
            pd.Series(["bar", "baz"]),
            pd.Index(["bar", "baz"]),
        ],
    )
    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_with_list_like_values_nans(self, values, method):
        # issue #17160
        df = pd.DataFrame(
            {
                "foo": ["one", "one", "one", "two", "two", "two"],
                "bar": ["A", "B", "C", "A", "B", "C"],
                "baz": [1, 2, 3, 4, 5, 6],
                "zoo": ["x", "y", "z", "q", "w", "t"],
            }
        )

        if method:
            result = df.pivot(index="zoo", columns="foo", values=values)
        else:
            result = pd.pivot(df, index="zoo", columns="foo", values=values)

        data = [
            [np.nan, "A", np.nan, 4],
            [np.nan, "C", np.nan, 6],
            [np.nan, "B", np.nan, 5],
            ["A", np.nan, 1, np.nan],
            ["B", np.nan, 2, np.nan],
            ["C", np.nan, 3, np.nan],
        ]
        index = Index(data=["q", "t", "w", "x", "y", "z"], name="zoo")
        columns = MultiIndex(
            levels=[["bar", "baz"], ["one", "two"]],
            codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
            names=[None, "foo"],
        )
        expected = DataFrame(data=data, index=index, columns=columns, dtype="object")
        tm.assert_frame_equal(result, expected)

    def test_pivot_columns_none_raise_error(self):
        # GH 30924
        df = pd.DataFrame(
            {"col1": ["a", "b", "c"], "col2": [1, 2, 3], "col3": [1, 2, 3]}
        )
        msg = r"pivot\(\) missing 1 required argument: 'columns'"
        with pytest.raises(TypeError, match=msg):
            df.pivot(index="col1", values="col3")

    @pytest.mark.xfail(
        reason="MultiIndexed unstack with tuple names fails with KeyError GH#19966"
    )
    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_with_multiindex(self, method):
        # issue #17160
        index = Index(data=[0, 1, 2, 3, 4, 5])
        data = [
            ["one", "A", 1, "x"],
            ["one", "B", 2, "y"],
            ["one", "C", 3, "z"],
            ["two", "A", 4, "q"],
            ["two", "B", 5, "w"],
            ["two", "C", 6, "t"],
        ]
        columns = MultiIndex(
            levels=[["bar", "baz"], ["first", "second"]],
            codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
        )
        df = DataFrame(data=data, index=index, columns=columns, dtype="object")
        if method:
            result = df.pivot(
                index=("bar", "first"),
                columns=("bar", "second"),
                values=("baz", "first"),
            )
        else:
            result = pd.pivot(
                df,
                index=("bar", "first"),
                columns=("bar", "second"),
                values=("baz", "first"),
            )

        data = {
            "A": Series([1, 4], index=["one", "two"]),
            "B": Series([2, 5], index=["one", "two"]),
            "C": Series([3, 6], index=["one", "two"]),
        }
        expected = DataFrame(data)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("method", [True, False])
    def test_pivot_with_tuple_of_values(self, method):
        # issue #17160
        df = pd.DataFrame(
            {
                "foo": ["one", "one", "one", "two", "two", "two"],
                "bar": ["A", "B", "C", "A", "B", "C"],
                "baz": [1, 2, 3, 4, 5, 6],
                "zoo": ["x", "y", "z", "q", "w", "t"],
            }
        )
        with pytest.raises(KeyError, match=r"^\('bar', 'baz'\)$"):
            # tuple is seen as a single column name
            if method:
                df.pivot(index="zoo", columns="foo", values=("bar", "baz"))
            else:
                pd.pivot(df, index="zoo", columns="foo", values=("bar", "baz"))

    def test_margins(self):
        def _check_output(
            result, values_col, index=["A", "B"], columns=["C"], margins_col="All"
        ):
            col_margins = result.loc[result.index[:-1], margins_col]
            expected_col_margins = self.data.groupby(index)[values_col].mean()
            tm.assert_series_equal(col_margins, expected_col_margins, check_names=False)
            assert col_margins.name == margins_col

            result = result.sort_index()
            index_margins = result.loc[(margins_col, "")].iloc[:-1]

            expected_ix_margins = self.data.groupby(columns)[values_col].mean()
            tm.assert_series_equal(
                index_margins, expected_ix_margins, check_names=False
            )
            assert index_margins.name == (margins_col, "")

            grand_total_margins = result.loc[(margins_col, ""), margins_col]
            expected_total_margins = self.data[values_col].mean()
            assert grand_total_margins == expected_total_margins

        # column specified
        result = self.data.pivot_table(
            values="D", index=["A", "B"], columns="C", margins=True, aggfunc=np.mean
        )
        _check_output(result, "D")

        # Set a different margins_name (not 'All')
        result = self.data.pivot_table(
            values="D",
            index=["A", "B"],
            columns="C",
            margins=True,
            aggfunc=np.mean,
            margins_name="Totals",
        )
        _check_output(result, "D", margins_col="Totals")

        # no column specified
        table = self.data.pivot_table(
            index=["A", "B"], columns="C", margins=True, aggfunc=np.mean
        )
        for value_col in table.columns.levels[0]:
            _check_output(table[value_col], value_col)

        # no col

        # to help with a buglet
        self.data.columns = [k * 2 for k in self.data.columns]
        table = self.data.pivot_table(index=["AA", "BB"], margins=True, aggfunc=np.mean)
        for value_col in table.columns:
            totals = table.loc[("All", ""), value_col]
            assert totals == self.data[value_col].mean()

        table = self.data.pivot_table(index=["AA", "BB"], margins=True, aggfunc="mean")
        for item in ["DD", "EE", "FF"]:
            totals = table.loc[("All", ""), item]
            assert totals == self.data[item].mean()

    def test_margins_dtype(self):
        # GH 17013

        df = self.data.copy()
        df[["D", "E", "F"]] = np.arange(len(df) * 3).reshape(len(df), 3)

        mi_val = list(product(["bar", "foo"], ["one", "two"])) + [("All", "")]
        mi = MultiIndex.from_tuples(mi_val, names=("A", "B"))
        expected = DataFrame(
            {"dull": [12, 21, 3, 9, 45], "shiny": [33, 0, 36, 51, 120]}, index=mi
        ).rename_axis("C", axis=1)
        expected["All"] = expected["dull"] + expected["shiny"]

        result = df.pivot_table(
            values="D",
            index=["A", "B"],
            columns="C",
            margins=True,
            aggfunc=np.sum,
            fill_value=0,
        )

        tm.assert_frame_equal(expected, result)

    @pytest.mark.xfail(reason="GH#17035 (len of floats is casted back to floats)")
    def test_margins_dtype_len(self):
        mi_val = list(product(["bar", "foo"], ["one", "two"])) + [("All", "")]
        mi = MultiIndex.from_tuples(mi_val, names=("A", "B"))
        expected = DataFrame(
            {"dull": [1, 1, 2, 1, 5], "shiny": [2, 0, 2, 2, 6]}, index=mi
        ).rename_axis("C", axis=1)
        expected["All"] = expected["dull"] + expected["shiny"]

        result = self.data.pivot_table(
            values="D",
            index=["A", "B"],
            columns="C",
            margins=True,
            aggfunc=len,
            fill_value=0,
        )

        tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize("cols", [(1, 2), ("a", "b"), (1, "b"), ("a", 1)])
    def test_pivot_table_multiindex_only(self, cols):
        # GH 17038
        df2 = DataFrame({cols[0]: [1, 2, 3], cols[1]: [1, 2, 3], "v": [4, 5, 6]})

        result = df2.pivot_table(values="v", columns=cols)
        expected = DataFrame(
            [[4, 5, 6]],
            columns=MultiIndex.from_tuples([(1, 1), (2, 2), (3, 3)], names=cols),
            index=Index(["v"]),
        )

        tm.assert_frame_equal(result, expected)

    def test_pivot_integer_columns(self):
        # caused by upstream bug in unstack

        d = date.min
        data = list(
            product(
                ["foo", "bar"],
                ["A", "B", "C"],
                ["x1", "x2"],
                [d + timedelta(i) for i in range(20)],
                [1.0],
            )
        )
        df = DataFrame(data)
        table = df.pivot_table(values=4, index=[0, 1, 3], columns=[2])

        df2 = df.rename(columns=str)
        table2 = df2.pivot_table(values="4", index=["0", "1", "3"], columns=["2"])

        tm.assert_frame_equal(table, table2, check_names=False)

    def test_pivot_no_level_overlap(self):
        # GH #1181

        data = DataFrame(
            {
                "a": ["a", "a", "a", "a", "b", "b", "b", "b"] * 2,
                "b": [0, 0, 0, 0, 1, 1, 1, 1] * 2,
                "c": (["foo"] * 4 + ["bar"] * 4) * 2,
                "value": np.random.randn(16),
            }
        )

        table = data.pivot_table("value", index="a", columns=["b", "c"])

        grouped = data.groupby(["a", "b", "c"])["value"].mean()
        expected = grouped.unstack("b").unstack("c").dropna(axis=1, how="all")
        tm.assert_frame_equal(table, expected)

    def test_pivot_columns_lexsorted(self):

        n = 10000

        dtype = np.dtype(
            [
                ("Index", object),
                ("Symbol", object),
                ("Year", int),
                ("Month", int),
                ("Day", int),
                ("Quantity", int),
                ("Price", float),
            ]
        )

        products = np.array(
            [
                ("SP500", "ADBE"),
                ("SP500", "NVDA"),
                ("SP500", "ORCL"),
                ("NDQ100", "AAPL"),
                ("NDQ100", "MSFT"),
                ("NDQ100", "GOOG"),
                ("FTSE", "DGE.L"),
                ("FTSE", "TSCO.L"),
                ("FTSE", "GSK.L"),
            ],
            dtype=[("Index", object), ("Symbol", object)],
        )
        items = np.empty(n, dtype=dtype)
        iproduct = np.random.randint(0, len(products), n)
        items["Index"] = products["Index"][iproduct]
        items["Symbol"] = products["Symbol"][iproduct]
        dr = pd.date_range(date(2000, 1, 1), date(2010, 12, 31))
        dates = dr[np.random.randint(0, len(dr), n)]
        items["Year"] = dates.year
        items["Month"] = dates.month
        items["Day"] = dates.day
        items["Price"] = np.random.lognormal(4.0, 2.0, n)

        df = DataFrame(items)

        pivoted = df.pivot_table(
            "Price",
            index=["Month", "Day"],
            columns=["Index", "Symbol", "Year"],
            aggfunc="mean",
        )

        assert pivoted.columns.is_monotonic

    def test_pivot_complex_aggfunc(self):
        f = {"D": ["std"], "E": ["sum"]}
        expected = self.data.groupby(["A", "B"]).agg(f).unstack("B")
        result = self.data.pivot_table(index="A", columns="B", aggfunc=f)

        tm.assert_frame_equal(result, expected)

    def test_margins_no_values_no_cols(self):
        # Regression test on pivot table: no values or cols passed.
        result = self.data[["A", "B"]].pivot_table(
            index=["A", "B"], aggfunc=len, margins=True
        )
        result_list = result.tolist()
        assert sum(result_list[:-1]) == result_list[-1]

    def test_margins_no_values_two_rows(self):
        # Regression test on pivot table: no values passed but rows are a
        # multi-index
        result = self.data[["A", "B", "C"]].pivot_table(
            index=["A", "B"], columns="C", aggfunc=len, margins=True
        )
        assert result.All.tolist() == [3.0, 1.0, 4.0, 3.0, 11.0]

    def test_margins_no_values_one_row_one_col(self):
        # Regression test on pivot table: no values passed but row and col
        # defined
        result = self.data[["A", "B"]].pivot_table(
            index="A", columns="B", aggfunc=len, margins=True
        )
        assert result.All.tolist() == [4.0, 7.0, 11.0]

    def test_margins_no_values_two_row_two_cols(self):
        # Regression test on pivot table: no values passed but rows and cols
        # are multi-indexed
        self.data["D"] = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"]
        result = self.data[["A", "B", "C", "D"]].pivot_table(
            index=["A", "B"], columns=["C", "D"], aggfunc=len, margins=True
        )
        assert result.All.tolist() == [3.0, 1.0, 4.0, 3.0, 11.0]

    @pytest.mark.parametrize("margin_name", ["foo", "one", 666, None, ["a", "b"]])
    def test_pivot_table_with_margins_set_margin_name(self, margin_name):
        # see gh-3335
        msg = (
            r'Conflicting name "{}" in margins|'
            "margins_name argument must be a string"
        ).format(margin_name)
        with pytest.raises(ValueError, match=msg):
            # multi-index index
            pivot_table(
                self.data,
                values="D",
                index=["A", "B"],
                columns=["C"],
                margins=True,
                margins_name=margin_name,
            )
        with pytest.raises(ValueError, match=msg):
            # multi-index column
            pivot_table(
                self.data,
                values="D",
                index=["C"],
                columns=["A", "B"],
                margins=True,
                margins_name=margin_name,
            )
        with pytest.raises(ValueError, match=msg):
            # non-multi-index index/column
            pivot_table(
                self.data,
                values="D",
                index=["A"],
                columns=["B"],
                margins=True,
                margins_name=margin_name,
            )

    def test_pivot_timegrouper(self):
        df = DataFrame(
            {
                "Branch": "A A A A A A A B".split(),
                "Buyer": "Carl Mark Carl Carl Joe Joe Joe Carl".split(),
                "Quantity": [1, 3, 5, 1, 8, 1, 9, 3],
                "Date": [
                    datetime(2013, 1, 1),
                    datetime(2013, 1, 1),
                    datetime(2013, 10, 1),
                    datetime(2013, 10, 2),
                    datetime(2013, 10, 1),
                    datetime(2013, 10, 2),
                    datetime(2013, 12, 2),
                    datetime(2013, 12, 2),
                ],
            }
        ).set_index("Date")

        expected = DataFrame(
            np.array([10, 18, 3], dtype="int64").reshape(1, 3),
            index=[datetime(2013, 12, 31)],
            columns="Carl Joe Mark".split(),
        )
        expected.index.name = "Date"
        expected.columns.name = "Buyer"

        result = pivot_table(
            df,
            index=Grouper(freq="A"),
            columns="Buyer",
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index="Buyer",
            columns=Grouper(freq="A"),
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

        expected = DataFrame(
            np.array([1, np.nan, 3, 9, 18, np.nan]).reshape(2, 3),
            index=[datetime(2013, 1, 1), datetime(2013, 7, 1)],
            columns="Carl Joe Mark".split(),
        )
        expected.index.name = "Date"
        expected.columns.name = "Buyer"

        result = pivot_table(
            df,
            index=Grouper(freq="6MS"),
            columns="Buyer",
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index="Buyer",
            columns=Grouper(freq="6MS"),
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

        # passing the name
        df = df.reset_index()
        result = pivot_table(
            df,
            index=Grouper(freq="6MS", key="Date"),
            columns="Buyer",
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index="Buyer",
            columns=Grouper(freq="6MS", key="Date"),
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

        msg = "'The grouper name foo is not found'"
        with pytest.raises(KeyError, match=msg):
            pivot_table(
                df,
                index=Grouper(freq="6MS", key="foo"),
                columns="Buyer",
                values="Quantity",
                aggfunc=np.sum,
            )
        with pytest.raises(KeyError, match=msg):
            pivot_table(
                df,
                index="Buyer",
                columns=Grouper(freq="6MS", key="foo"),
                values="Quantity",
                aggfunc=np.sum,
            )

        # passing the level
        df = df.set_index("Date")
        result = pivot_table(
            df,
            index=Grouper(freq="6MS", level="Date"),
            columns="Buyer",
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index="Buyer",
            columns=Grouper(freq="6MS", level="Date"),
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

        msg = "The level foo is not valid"
        with pytest.raises(ValueError, match=msg):
            pivot_table(
                df,
                index=Grouper(freq="6MS", level="foo"),
                columns="Buyer",
                values="Quantity",
                aggfunc=np.sum,
            )
        with pytest.raises(ValueError, match=msg):
            pivot_table(
                df,
                index="Buyer",
                columns=Grouper(freq="6MS", level="foo"),
                values="Quantity",
                aggfunc=np.sum,
            )

        # double grouper
        df = DataFrame(
            {
                "Branch": "A A A A A A A B".split(),
                "Buyer": "Carl Mark Carl Carl Joe Joe Joe Carl".split(),
                "Quantity": [1, 3, 5, 1, 8, 1, 9, 3],
                "Date": [
                    datetime(2013, 11, 1, 13, 0),
                    datetime(2013, 9, 1, 13, 5),
                    datetime(2013, 10, 1, 20, 0),
                    datetime(2013, 10, 2, 10, 0),
                    datetime(2013, 11, 1, 20, 0),
                    datetime(2013, 10, 2, 10, 0),
                    datetime(2013, 10, 2, 12, 0),
                    datetime(2013, 12, 5, 14, 0),
                ],
                "PayDay": [
                    datetime(2013, 10, 4, 0, 0),
                    datetime(2013, 10, 15, 13, 5),
                    datetime(2013, 9, 5, 20, 0),
                    datetime(2013, 11, 2, 10, 0),
                    datetime(2013, 10, 7, 20, 0),
                    datetime(2013, 9, 5, 10, 0),
                    datetime(2013, 12, 30, 12, 0),
                    datetime(2013, 11, 20, 14, 0),
                ],
            }
        )

        result = pivot_table(
            df,
            index=Grouper(freq="M", key="Date"),
            columns=Grouper(freq="M", key="PayDay"),
            values="Quantity",
            aggfunc=np.sum,
        )
        expected = DataFrame(
            np.array(
                [
                    np.nan,
                    3,
                    np.nan,
                    np.nan,
                    6,
                    np.nan,
                    1,
                    9,
                    np.nan,
                    9,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    3,
                    np.nan,
                ]
            ).reshape(4, 4),
            index=[
                datetime(2013, 9, 30),
                datetime(2013, 10, 31),
                datetime(2013, 11, 30),
                datetime(2013, 12, 31),
            ],
            columns=[
                datetime(2013, 9, 30),
                datetime(2013, 10, 31),
                datetime(2013, 11, 30),
                datetime(2013, 12, 31),
            ],
        )
        expected.index.name = "Date"
        expected.columns.name = "PayDay"

        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index=Grouper(freq="M", key="PayDay"),
            columns=Grouper(freq="M", key="Date"),
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

        tuples = [
            (datetime(2013, 9, 30), datetime(2013, 10, 31)),
            (datetime(2013, 10, 31), datetime(2013, 9, 30)),
            (datetime(2013, 10, 31), datetime(2013, 11, 30)),
            (datetime(2013, 10, 31), datetime(2013, 12, 31)),
            (datetime(2013, 11, 30), datetime(2013, 10, 31)),
            (datetime(2013, 12, 31), datetime(2013, 11, 30)),
        ]
        idx = MultiIndex.from_tuples(tuples, names=["Date", "PayDay"])
        expected = DataFrame(
            np.array(
                [3, np.nan, 6, np.nan, 1, np.nan, 9, np.nan, 9, np.nan, np.nan, 3]
            ).reshape(6, 2),
            index=idx,
            columns=["A", "B"],
        )
        expected.columns.name = "Branch"

        result = pivot_table(
            df,
            index=[Grouper(freq="M", key="Date"), Grouper(freq="M", key="PayDay")],
            columns=["Branch"],
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index=["Branch"],
            columns=[Grouper(freq="M", key="Date"), Grouper(freq="M", key="PayDay")],
            values="Quantity",
            aggfunc=np.sum,
        )
        tm.assert_frame_equal(result, expected.T)

    def test_pivot_datetime_tz(self):
        dates1 = [
            "2011-07-19 07:00:00",
            "2011-07-19 08:00:00",
            "2011-07-19 09:00:00",
            "2011-07-19 07:00:00",
            "2011-07-19 08:00:00",
            "2011-07-19 09:00:00",
        ]
        dates2 = [
            "2013-01-01 15:00:00",
            "2013-01-01 15:00:00",
            "2013-01-01 15:00:00",
            "2013-02-01 15:00:00",
            "2013-02-01 15:00:00",
            "2013-02-01 15:00:00",
        ]
        df = DataFrame(
            {
                "label": ["a", "a", "a", "b", "b", "b"],
                "dt1": dates1,
                "dt2": dates2,
                "value1": np.arange(6, dtype="int64"),
                "value2": [1, 2] * 3,
            }
        )
        df["dt1"] = df["dt1"].apply(lambda d: pd.Timestamp(d, tz="US/Pacific"))
        df["dt2"] = df["dt2"].apply(lambda d: pd.Timestamp(d, tz="Asia/Tokyo"))

        exp_idx = pd.DatetimeIndex(
            ["2011-07-19 07:00:00", "2011-07-19 08:00:00", "2011-07-19 09:00:00"],
            tz="US/Pacific",
            name="dt1",
        )
        exp_col1 = Index(["value1", "value1"])
        exp_col2 = Index(["a", "b"], name="label")
        exp_col = MultiIndex.from_arrays([exp_col1, exp_col2])
        expected = DataFrame([[0, 3], [1, 4], [2, 5]], index=exp_idx, columns=exp_col)
        result = pivot_table(df, index=["dt1"], columns=["label"], values=["value1"])
        tm.assert_frame_equal(result, expected)

        exp_col1 = Index(["sum", "sum", "sum", "sum", "mean", "mean", "mean", "mean"])
        exp_col2 = Index(["value1", "value1", "value2", "value2"] * 2)
        exp_col3 = pd.DatetimeIndex(
            ["2013-01-01 15:00:00", "2013-02-01 15:00:00"] * 4,
            tz="Asia/Tokyo",
            name="dt2",
        )
        exp_col = MultiIndex.from_arrays([exp_col1, exp_col2, exp_col3])
        expected = DataFrame(
            np.array(
                [
                    [0, 3, 1, 2, 0, 3, 1, 2],
                    [1, 4, 2, 1, 1, 4, 2, 1],
                    [2, 5, 1, 2, 2, 5, 1, 2],
                ],
                dtype="int64",
            ),
            index=exp_idx,
            columns=exp_col,
        )

        result = pivot_table(
            df,
            index=["dt1"],
            columns=["dt2"],
            values=["value1", "value2"],
            aggfunc=[np.sum, np.mean],
        )
        tm.assert_frame_equal(result, expected)

    def test_pivot_dtaccessor(self):
        # GH 8103
        dates1 = [
            "2011-07-19 07:00:00",
            "2011-07-19 08:00:00",
            "2011-07-19 09:00:00",
            "2011-07-19 07:00:00",
            "2011-07-19 08:00:00",
            "2011-07-19 09:00:00",
        ]
        dates2 = [
            "2013-01-01 15:00:00",
            "2013-01-01 15:00:00",
            "2013-01-01 15:00:00",
            "2013-02-01 15:00:00",
            "2013-02-01 15:00:00",
            "2013-02-01 15:00:00",
        ]
        df = DataFrame(
            {
                "label": ["a", "a", "a", "b", "b", "b"],
                "dt1": dates1,
                "dt2": dates2,
                "value1": np.arange(6, dtype="int64"),
                "value2": [1, 2] * 3,
            }
        )
        df["dt1"] = df["dt1"].apply(lambda d: pd.Timestamp(d))
        df["dt2"] = df["dt2"].apply(lambda d: pd.Timestamp(d))

        result = pivot_table(
            df, index="label", columns=df["dt1"].dt.hour, values="value1"
        )

        exp_idx = Index(["a", "b"], name="label")
        expected = DataFrame(
            {7: [0, 3], 8: [1, 4], 9: [2, 5]},
            index=exp_idx,
            columns=Index([7, 8, 9], name="dt1"),
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df, index=df["dt2"].dt.month, columns=df["dt1"].dt.hour, values="value1"
        )

        expected = DataFrame(
            {7: [0, 3], 8: [1, 4], 9: [2, 5]},
            index=Index([1, 2], name="dt2"),
            columns=Index([7, 8, 9], name="dt1"),
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index=df["dt2"].dt.year.values,
            columns=[df["dt1"].dt.hour, df["dt2"].dt.month],
            values="value1",
        )

        exp_col = MultiIndex.from_arrays(
            [[7, 7, 8, 8, 9, 9], [1, 2] * 3], names=["dt1", "dt2"]
        )
        expected = DataFrame(
            np.array([[0, 3, 1, 4, 2, 5]], dtype="int64"), index=[2013], columns=exp_col
        )
        tm.assert_frame_equal(result, expected)

        result = pivot_table(
            df,
            index=np.array(["X", "X", "X", "X", "Y", "Y"]),
            columns=[df["dt1"].dt.hour, df["dt2"].dt.month],
            values="value1",
        )
        expected = DataFrame(
            np.array(
                [[0, 3, 1, np.nan, 2, np.nan], [np.nan, np.nan, np.nan, 4, np.nan, 5]]
            ),
            index=["X", "Y"],
            columns=exp_col,
        )
        tm.assert_frame_equal(result, expected)

    def test_daily(self):
        rng = date_range("1/1/2000", "12/31/2004", freq="D")
        ts = Series(np.random.randn(len(rng)), index=rng)

        annual = pivot_table(
            DataFrame(ts), index=ts.index.year, columns=ts.index.dayofyear
        )
        annual.columns = annual.columns.droplevel(0)

        doy = np.asarray(ts.index.dayofyear)

        for i in range(1, 367):
            subset = ts[doy == i]
            subset.index = subset.index.year

            result = annual[i].dropna()
            tm.assert_series_equal(result, subset, check_names=False)
            assert result.name == i

    def test_monthly(self):
        rng = date_range("1/1/2000", "12/31/2004", freq="M")
        ts = Series(np.random.randn(len(rng)), index=rng)

        annual = pivot_table(
            pd.DataFrame(ts), index=ts.index.year, columns=ts.index.month
        )
        annual.columns = annual.columns.droplevel(0)

        month = ts.index.month
        for i in range(1, 13):
            subset = ts[month == i]
            subset.index = subset.index.year
            result = annual[i].dropna()
            tm.assert_series_equal(result, subset, check_names=False)
            assert result.name == i

    def test_pivot_table_with_iterator_values(self):
        # GH 12017
        aggs = {"D": "sum", "E": "mean"}

        pivot_values_list = pd.pivot_table(
            self.data, index=["A"], values=list(aggs.keys()), aggfunc=aggs
        )

        pivot_values_keys = pd.pivot_table(
            self.data, index=["A"], values=aggs.keys(), aggfunc=aggs
        )
        tm.assert_frame_equal(pivot_values_keys, pivot_values_list)

        agg_values_gen = (value for value in aggs.keys())
        pivot_values_gen = pd.pivot_table(
            self.data, index=["A"], values=agg_values_gen, aggfunc=aggs
        )
        tm.assert_frame_equal(pivot_values_gen, pivot_values_list)

    def test_pivot_table_margins_name_with_aggfunc_list(self):
        # GH 13354
        margins_name = "Weekly"
        costs = pd.DataFrame(
            {
                "item": ["bacon", "cheese", "bacon", "cheese"],
                "cost": [2.5, 4.5, 3.2, 3.3],
                "day": ["M", "M", "T", "T"],
            }
        )
        table = costs.pivot_table(
            index="item",
            columns="day",
            margins=True,
            margins_name=margins_name,
            aggfunc=[np.mean, max],
        )
        ix = pd.Index(["bacon", "cheese", margins_name], dtype="object", name="item")
        tups = [
            ("mean", "cost", "M"),
            ("mean", "cost", "T"),
            ("mean", "cost", margins_name),
            ("max", "cost", "M"),
            ("max", "cost", "T"),
            ("max", "cost", margins_name),
        ]
        cols = pd.MultiIndex.from_tuples(tups, names=[None, None, "day"])
        expected = pd.DataFrame(table.values, index=ix, columns=cols)
        tm.assert_frame_equal(table, expected)

    @pytest.mark.xfail(reason="GH#17035 (np.mean of ints is casted back to ints)")
    def test_categorical_margins(self, observed):
        # GH 10989
        df = pd.DataFrame(
            {"x": np.arange(8), "y": np.arange(8) // 4, "z": np.arange(8) % 2}
        )

        expected = pd.DataFrame([[1.0, 2.0, 1.5], [5, 6, 5.5], [3, 4, 3.5]])
        expected.index = Index([0, 1, "All"], name="y")
        expected.columns = Index([0, 1, "All"], name="z")

        table = df.pivot_table("x", "y", "z", dropna=observed, margins=True)
        tm.assert_frame_equal(table, expected)

    @pytest.mark.xfail(reason="GH#17035 (np.mean of ints is casted back to ints)")
    def test_categorical_margins_category(self, observed):
        df = pd.DataFrame(
            {"x": np.arange(8), "y": np.arange(8) // 4, "z": np.arange(8) % 2}
        )

        expected = pd.DataFrame([[1.0, 2.0, 1.5], [5, 6, 5.5], [3, 4, 3.5]])
        expected.index = Index([0, 1, "All"], name="y")
        expected.columns = Index([0, 1, "All"], name="z")

        df.y = df.y.astype("category")
        df.z = df.z.astype("category")
        table = df.pivot_table("x", "y", "z", dropna=observed, margins=True)
        tm.assert_frame_equal(table, expected)

    def test_margins_casted_to_float(self, observed):
        # GH 24893
        df = pd.DataFrame(
            {
                "A": [2, 4, 6, 8],
                "B": [1, 4, 5, 8],
                "C": [1, 3, 4, 6],
                "D": ["X", "X", "Y", "Y"],
            }
        )

        result = pd.pivot_table(df, index="D", margins=True)
        expected = pd.DataFrame(
            {"A": [3, 7, 5], "B": [2.5, 6.5, 4.5], "C": [2, 5, 3.5]},
            index=pd.Index(["X", "Y", "All"], name="D"),
        )
        tm.assert_frame_equal(result, expected)

    def test_pivot_with_categorical(self, observed, ordered_fixture):
        # gh-21370
        idx = [np.nan, "low", "high", "low", np.nan]
        col = [np.nan, "A", "B", np.nan, "A"]
        df = pd.DataFrame(
            {
                "In": pd.Categorical(
                    idx, categories=["low", "high"], ordered=ordered_fixture
                ),
                "Col": pd.Categorical(
                    col, categories=["A", "B"], ordered=ordered_fixture
                ),
                "Val": range(1, 6),
            }
        )
        # case with index/columns/value
        result = df.pivot_table(
            index="In", columns="Col", values="Val", observed=observed
        )

        expected_cols = pd.CategoricalIndex(
            ["A", "B"], ordered=ordered_fixture, name="Col"
        )

        expected = pd.DataFrame(
            data=[[2.0, np.nan], [np.nan, 3.0]], columns=expected_cols
        )
        expected.index = Index(
            pd.Categorical(
                ["low", "high"], categories=["low", "high"], ordered=ordered_fixture
            ),
            name="In",
        )

        tm.assert_frame_equal(result, expected)

        # case with columns/value
        result = df.pivot_table(columns="Col", values="Val", observed=observed)

        expected = pd.DataFrame(
            data=[[3.5, 3.0]], columns=expected_cols, index=Index(["Val"])
        )

        tm.assert_frame_equal(result, expected)

    def test_categorical_aggfunc(self, observed):
        # GH 9534
        df = pd.DataFrame(
            {"C1": ["A", "B", "C", "C"], "C2": ["a", "a", "b", "b"], "V": [1, 2, 3, 4]}
        )
        df["C1"] = df["C1"].astype("category")
        result = df.pivot_table(
            "V", index="C1", columns="C2", dropna=observed, aggfunc="count"
        )

        expected_index = pd.CategoricalIndex(
            ["A", "B", "C"], categories=["A", "B", "C"], ordered=False, name="C1"
        )
        expected_columns = pd.Index(["a", "b"], name="C2")
        expected_data = np.array([[1.0, np.nan], [1.0, np.nan], [np.nan, 2.0]])
        expected = pd.DataFrame(
            expected_data, index=expected_index, columns=expected_columns
        )
        tm.assert_frame_equal(result, expected)

    def test_categorical_pivot_index_ordering(self, observed):
        # GH 8731
        df = pd.DataFrame(
            {
                "Sales": [100, 120, 220],
                "Month": ["January", "January", "January"],
                "Year": [2013, 2014, 2013],
            }
        )
        months = [
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]
        df["Month"] = df["Month"].astype("category").cat.set_categories(months)
        result = df.pivot_table(
            values="Sales",
            index="Month",
            columns="Year",
            dropna=observed,
            aggfunc="sum",
        )
        expected_columns = pd.Int64Index([2013, 2014], name="Year")
        expected_index = pd.CategoricalIndex(
            ["January"], categories=months, ordered=False, name="Month"
        )
        expected = pd.DataFrame(
            [[320, 120]], index=expected_index, columns=expected_columns
        )
        if not observed:
            result = result.dropna().astype(np.int64)

        tm.assert_frame_equal(result, expected)

    def test_pivot_table_not_series(self):
        # GH 4386
        # pivot_table always returns a DataFrame
        # when values is not list like and columns is None
        # and aggfunc is not instance of list
        df = DataFrame({"col1": [3, 4, 5], "col2": ["C", "D", "E"], "col3": [1, 3, 9]})

        result = df.pivot_table("col1", index=["col3", "col2"], aggfunc=np.sum)
        m = MultiIndex.from_arrays([[1, 3, 9], ["C", "D", "E"]], names=["col3", "col2"])
        expected = DataFrame([3, 4, 5], index=m, columns=["col1"])

        tm.assert_frame_equal(result, expected)

        result = df.pivot_table("col1", index="col3", columns="col2", aggfunc=np.sum)
        expected = DataFrame(
            [[3, np.NaN, np.NaN], [np.NaN, 4, np.NaN], [np.NaN, np.NaN, 5]],
            index=Index([1, 3, 9], name="col3"),
            columns=Index(["C", "D", "E"], name="col2"),
        )

        tm.assert_frame_equal(result, expected)

        result = df.pivot_table("col1", index="col3", aggfunc=[np.sum])
        m = MultiIndex.from_arrays([["sum"], ["col1"]])
        expected = DataFrame([3, 4, 5], index=Index([1, 3, 9], name="col3"), columns=m)

        tm.assert_frame_equal(result, expected)

    def test_pivot_margins_name_unicode(self):
        # issue #13292
        greek = "\u0394\u03bf\u03ba\u03b9\u03bc\u03ae"
        frame = pd.DataFrame({"foo": [1, 2, 3]})
        table = pd.pivot_table(
            frame, index=["foo"], aggfunc=len, margins=True, margins_name=greek
        )
        index = pd.Index([1, 2, 3, greek], dtype="object", name="foo")
        expected = pd.DataFrame(index=index)
        tm.assert_frame_equal(table, expected)

    def test_pivot_string_as_func(self):
        # GH #18713
        # for correctness purposes
        data = DataFrame(
            {
                "A": [
                    "foo",
                    "foo",
                    "foo",
                    "foo",
                    "bar",
                    "bar",
                    "bar",
                    "bar",
                    "foo",
                    "foo",
                    "foo",
                ],
                "B": [
                    "one",
                    "one",
                    "one",
                    "two",
                    "one",
                    "one",
                    "one",
                    "two",
                    "two",
                    "two",
                    "one",
                ],
                "C": range(11),
            }
        )

        result = pivot_table(data, index="A", columns="B", aggfunc="sum")
        mi = MultiIndex(
            levels=[["C"], ["one", "two"]], codes=[[0, 0], [0, 1]], names=[None, "B"]
        )
        expected = DataFrame(
            {("C", "one"): {"bar": 15, "foo": 13}, ("C", "two"): {"bar": 7, "foo": 20}},
            columns=mi,
        ).rename_axis("A")
        tm.assert_frame_equal(result, expected)

        result = pivot_table(data, index="A", columns="B", aggfunc=["sum", "mean"])
        mi = MultiIndex(
            levels=[["sum", "mean"], ["C"], ["one", "two"]],
            codes=[[0, 0, 1, 1], [0, 0, 0, 0], [0, 1, 0, 1]],
            names=[None, None, "B"],
        )
        expected = DataFrame(
            {
                ("mean", "C", "one"): {"bar": 5.0, "foo": 3.25},
                ("mean", "C", "two"): {"bar": 7.0, "foo": 6.666666666666667},
                ("sum", "C", "one"): {"bar": 15, "foo": 13},
                ("sum", "C", "two"): {"bar": 7, "foo": 20},
            },
            columns=mi,
        ).rename_axis("A")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "f, f_numpy",
        [
            ("sum", np.sum),
            ("mean", np.mean),
            ("std", np.std),
            (["sum", "mean"], [np.sum, np.mean]),
            (["sum", "std"], [np.sum, np.std]),
            (["std", "mean"], [np.std, np.mean]),
        ],
    )
    def test_pivot_string_func_vs_func(self, f, f_numpy):
        # GH #18713
        # for consistency purposes
        result = pivot_table(self.data, index="A", columns="B", aggfunc=f)
        expected = pivot_table(self.data, index="A", columns="B", aggfunc=f_numpy)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.slow
    def test_pivot_number_of_levels_larger_than_int32(self):
        # GH 20601
        df = DataFrame(
            {"ind1": np.arange(2 ** 16), "ind2": np.arange(2 ** 16), "count": 0}
        )

        msg = "Unstacked DataFrame is too big, causing int32 overflow"
        with pytest.raises(ValueError, match=msg):
            df.pivot_table(
                index="ind1", columns="ind2", values="count", aggfunc="count"
            )

    def test_pivot_table_aggfunc_dropna(self, dropna):
        # GH 22159
        df = pd.DataFrame(
            {
                "fruit": ["apple", "peach", "apple"],
                "size": [1, 1, 2],
                "taste": [7, 6, 6],
            }
        )

        def ret_one(x):
            return 1

        def ret_sum(x):
            return sum(x)

        def ret_none(x):
            return np.nan

        result = pd.pivot_table(
            df, columns="fruit", aggfunc=[ret_sum, ret_none, ret_one], dropna=dropna
        )

        data = [[3, 1, np.nan, np.nan, 1, 1], [13, 6, np.nan, np.nan, 1, 1]]
        col = pd.MultiIndex.from_product(
            [["ret_sum", "ret_none", "ret_one"], ["apple", "peach"]],
            names=[None, "fruit"],
        )
        expected = pd.DataFrame(data, index=["size", "taste"], columns=col)

        if dropna:
            expected = expected.dropna(axis="columns")

        tm.assert_frame_equal(result, expected)

    def test_pivot_table_aggfunc_scalar_dropna(self, dropna):
        # GH 22159
        df = pd.DataFrame(
            {"A": ["one", "two", "one"], "x": [3, np.nan, 2], "y": [1, np.nan, np.nan]}
        )

        result = pd.pivot_table(df, columns="A", aggfunc=np.mean, dropna=dropna)

        data = [[2.5, np.nan], [1, np.nan]]
        col = pd.Index(["one", "two"], name="A")
        expected = pd.DataFrame(data, index=["x", "y"], columns=col)

        if dropna:
            expected = expected.dropna(axis="columns")

        tm.assert_frame_equal(result, expected)

    def test_pivot_table_empty_aggfunc(self):
        # GH 9186
        df = pd.DataFrame(
            {
                "A": [2, 2, 3, 3, 2],
                "id": [5, 6, 7, 8, 9],
                "C": ["p", "q", "q", "p", "q"],
                "D": [None, None, None, None, None],
            }
        )
        result = df.pivot_table(index="A", columns="D", values="id", aggfunc=np.size)
        expected = pd.DataFrame()
        tm.assert_frame_equal(result, expected)

    def test_pivot_table_no_column_raises(self):
        # GH 10326
        def agg(l):
            return np.mean(l)

        foo = pd.DataFrame(
            {"X": [0, 0, 1, 1], "Y": [0, 1, 0, 1], "Z": [10, 20, 30, 40]}
        )
        with pytest.raises(KeyError, match="notpresent"):
            foo.pivot_table("notpresent", "X", "Y", aggfunc=agg)


class TestCrosstab:
    def setup_method(self, method):
        df = DataFrame(
            {
                "A": [
                    "foo",
                    "foo",
                    "foo",
                    "foo",
                    "bar",
                    "bar",
                    "bar",
                    "bar",
                    "foo",
                    "foo",
                    "foo",
                ],
                "B": [
                    "one",
                    "one",
                    "one",
                    "two",
                    "one",
                    "one",
                    "one",
                    "two",
                    "two",
                    "two",
                    "one",
                ],
                "C": [
                    "dull",
                    "dull",
                    "shiny",
                    "dull",
                    "dull",
                    "shiny",
                    "shiny",
                    "dull",
                    "shiny",
                    "shiny",
                    "shiny",
                ],
                "D": np.random.randn(11),
                "E": np.random.randn(11),
                "F": np.random.randn(11),
            }
        )

        self.df = df.append(df, ignore_index=True)

    def test_crosstab_single(self):
        df = self.df
        result = crosstab(df["A"], df["C"])
        expected = df.groupby(["A", "C"]).size().unstack()
        tm.assert_frame_equal(result, expected.fillna(0).astype(np.int64))

    def test_crosstab_multiple(self):
        df = self.df

        result = crosstab(df["A"], [df["B"], df["C"]])
        expected = df.groupby(["A", "B", "C"]).size()
        expected = expected.unstack("B").unstack("C").fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

        result = crosstab([df["B"], df["C"]], df["A"])
        expected = df.groupby(["B", "C", "A"]).size()
        expected = expected.unstack("A").fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_ndarray(self):
        a = np.random.randint(0, 5, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 10, size=100)

        df = DataFrame({"a": a, "b": b, "c": c})

        result = crosstab(a, [b, c], rownames=["a"], colnames=("b", "c"))
        expected = crosstab(df["a"], [df["b"], df["c"]])
        tm.assert_frame_equal(result, expected)

        result = crosstab([b, c], a, colnames=["a"], rownames=("b", "c"))
        expected = crosstab([df["b"], df["c"]], df["a"])
        tm.assert_frame_equal(result, expected)

        # assign arbitrary names
        result = crosstab(self.df["A"].values, self.df["C"].values)
        assert result.index.name == "row_0"
        assert result.columns.name == "col_0"

    def test_crosstab_non_aligned(self):
        # GH 17005
        a = pd.Series([0, 1, 1], index=["a", "b", "c"])
        b = pd.Series([3, 4, 3, 4, 3], index=["a", "b", "c", "d", "f"])
        c = np.array([3, 4, 3])

        expected = pd.DataFrame(
            [[1, 0], [1, 1]],
            index=Index([0, 1], name="row_0"),
            columns=Index([3, 4], name="col_0"),
        )

        result = crosstab(a, b)
        tm.assert_frame_equal(result, expected)

        result = crosstab(a, c)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_margins(self):
        a = np.random.randint(0, 7, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 5, size=100)

        df = DataFrame({"a": a, "b": b, "c": c})

        result = crosstab(a, [b, c], rownames=["a"], colnames=("b", "c"), margins=True)

        assert result.index.names == ("a",)
        assert result.columns.names == ["b", "c"]

        all_cols = result["All", ""]
        exp_cols = df.groupby(["a"]).size().astype("i8")
        # to keep index.name
        exp_margin = Series([len(df)], index=Index(["All"], name="a"))
        exp_cols = exp_cols.append(exp_margin)
        exp_cols.name = ("All", "")

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.loc["All"]
        exp_rows = df.groupby(["b", "c"]).size().astype("i8")
        exp_rows = exp_rows.append(Series([len(df)], index=[("All", "")]))
        exp_rows.name = "All"

        exp_rows = exp_rows.reindex(all_rows.index)
        exp_rows = exp_rows.fillna(0).astype(np.int64)
        tm.assert_series_equal(all_rows, exp_rows)

    def test_crosstab_margins_set_margin_name(self):
        # GH 15972
        a = np.random.randint(0, 7, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 5, size=100)

        df = DataFrame({"a": a, "b": b, "c": c})

        result = crosstab(
            a,
            [b, c],
            rownames=["a"],
            colnames=("b", "c"),
            margins=True,
            margins_name="TOTAL",
        )

        assert result.index.names == ("a",)
        assert result.columns.names == ["b", "c"]

        all_cols = result["TOTAL", ""]
        exp_cols = df.groupby(["a"]).size().astype("i8")
        # to keep index.name
        exp_margin = Series([len(df)], index=Index(["TOTAL"], name="a"))
        exp_cols = exp_cols.append(exp_margin)
        exp_cols.name = ("TOTAL", "")

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.loc["TOTAL"]
        exp_rows = df.groupby(["b", "c"]).size().astype("i8")
        exp_rows = exp_rows.append(Series([len(df)], index=[("TOTAL", "")]))
        exp_rows.name = "TOTAL"

        exp_rows = exp_rows.reindex(all_rows.index)
        exp_rows = exp_rows.fillna(0).astype(np.int64)
        tm.assert_series_equal(all_rows, exp_rows)

        msg = "margins_name argument must be a string"
        for margins_name in [666, None, ["a", "b"]]:
            with pytest.raises(ValueError, match=msg):
                crosstab(
                    a,
                    [b, c],
                    rownames=["a"],
                    colnames=("b", "c"),
                    margins=True,
                    margins_name=margins_name,
                )

    def test_crosstab_pass_values(self):
        a = np.random.randint(0, 7, size=100)
        b = np.random.randint(0, 3, size=100)
        c = np.random.randint(0, 5, size=100)
        values = np.random.randn(100)

        table = crosstab(
            [a, b], c, values, aggfunc=np.sum, rownames=["foo", "bar"], colnames=["baz"]
        )

        df = DataFrame({"foo": a, "bar": b, "baz": c, "values": values})

        expected = df.pivot_table(
            "values", index=["foo", "bar"], columns="baz", aggfunc=np.sum
        )
        tm.assert_frame_equal(table, expected)

    def test_crosstab_dropna(self):
        # GH 3820
        a = np.array(["foo", "foo", "foo", "bar", "bar", "foo", "foo"], dtype=object)
        b = np.array(["one", "one", "two", "one", "two", "two", "two"], dtype=object)
        c = np.array(
            ["dull", "dull", "dull", "dull", "dull", "shiny", "shiny"], dtype=object
        )
        res = pd.crosstab(a, [b, c], rownames=["a"], colnames=["b", "c"], dropna=False)
        m = MultiIndex.from_tuples(
            [("one", "dull"), ("one", "shiny"), ("two", "dull"), ("two", "shiny")],
            names=["b", "c"],
        )
        tm.assert_index_equal(res.columns, m)

    def test_crosstab_no_overlap(self):
        # GS 10291

        s1 = pd.Series([1, 2, 3], index=[1, 2, 3])
        s2 = pd.Series([4, 5, 6], index=[4, 5, 6])

        actual = crosstab(s1, s2)
        expected = pd.DataFrame()

        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna(self):
        # GH 12577
        # pivot_table counts null into margin ('All')
        # when margins=true and dropna=true

        df = pd.DataFrame({"a": [1, 2, 2, 2, 2, np.nan], "b": [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [2, 3, 5]])
        expected.index = Index([1.0, 2.0, "All"], name="a")
        expected.columns = Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

        df = DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, 2, np.nan], "b": [3, np.nan, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = Index([1.0, 2.0, "All"], name="a")
        expected.columns = Index([3.0, 4.0, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

        df = DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, np.nan, 2], "b": [3, 3, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = Index([1.0, 2.0, "All"], name="a")
        expected.columns = Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

        # GH 12642
        # _add_margins raises KeyError: Level None not found
        # when margins=True and dropna=False
        df = pd.DataFrame({"a": [1, 2, 2, 2, 2, np.nan], "b": [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [2, 4, 6]])
        expected.index = Index([1.0, 2.0, "All"], name="a")
        expected.columns = Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

        df = DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, 2, np.nan], "b": [3, np.nan, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 4, 6]])
        expected.index = Index([1.0, 2.0, "All"], name="a")
        expected.columns = Index([3.0, 4.0, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

        a = np.array(["foo", "foo", "foo", "bar", "bar", "foo", "foo"], dtype=object)
        b = np.array(["one", "one", "two", "one", "two", np.nan, "two"], dtype=object)
        c = np.array(
            ["dull", "dull", "dull", "dull", "dull", "shiny", "shiny"], dtype=object
        )

        actual = pd.crosstab(
            a, [b, c], rownames=["a"], colnames=["b", "c"], margins=True, dropna=False
        )
        m = MultiIndex.from_arrays(
            [
                ["one", "one", "two", "two", "All"],
                ["dull", "shiny", "dull", "shiny", ""],
            ],
            names=["b", "c"],
        )
        expected = DataFrame(
            [[1, 0, 1, 0, 2], [2, 0, 1, 1, 5], [3, 0, 2, 1, 7]], columns=m
        )
        expected.index = Index(["bar", "foo", "All"], name="a")
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab(
            [a, b], c, rownames=["a", "b"], colnames=["c"], margins=True, dropna=False
        )
        m = MultiIndex.from_arrays(
            [["bar", "bar", "foo", "foo", "All"], ["one", "two", "one", "two", ""]],
            names=["a", "b"],
        )
        expected = DataFrame(
            [[1, 0, 1], [1, 0, 1], [2, 0, 2], [1, 1, 2], [5, 2, 7]], index=m
        )
        expected.columns = Index(["dull", "shiny", "All"], name="c")
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab(
            [a, b], c, rownames=["a", "b"], colnames=["c"], margins=True, dropna=True
        )
        m = MultiIndex.from_arrays(
            [["bar", "bar", "foo", "foo", "All"], ["one", "two", "one", "two", ""]],
            names=["a", "b"],
        )
        expected = DataFrame(
            [[1, 0, 1], [1, 0, 1], [2, 0, 2], [1, 1, 2], [5, 1, 6]], index=m
        )
        expected.columns = Index(["dull", "shiny", "All"], name="c")
        tm.assert_frame_equal(actual, expected)

    def test_crosstab_normalize(self):
        # Issue 12578
        df = pd.DataFrame(
            {"a": [1, 2, 2, 2, 2], "b": [3, 3, 4, 4, 4], "c": [1, 1, np.nan, 1, 1]}
        )

        rindex = pd.Index([1, 2], name="a")
        cindex = pd.Index([3, 4], name="b")
        full_normal = pd.DataFrame([[0.2, 0], [0.2, 0.6]], index=rindex, columns=cindex)
        row_normal = pd.DataFrame(
            [[1.0, 0], [0.25, 0.75]], index=rindex, columns=cindex
        )
        col_normal = pd.DataFrame([[0.5, 0], [0.5, 1.0]], index=rindex, columns=cindex)

        # Check all normalize args
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize="all"), full_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize=True), full_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize="index"), row_normal)
        tm.assert_frame_equal(pd.crosstab(df.a, df.b, normalize="columns"), col_normal)
        tm.assert_frame_equal(
            pd.crosstab(df.a, df.b, normalize=1),
            pd.crosstab(df.a, df.b, normalize="columns"),
        )
        tm.assert_frame_equal(
            pd.crosstab(df.a, df.b, normalize=0),
            pd.crosstab(df.a, df.b, normalize="index"),
        )

        row_normal_margins = pd.DataFrame(
            [[1.0, 0], [0.25, 0.75], [0.4, 0.6]],
            index=pd.Index([1, 2, "All"], name="a", dtype="object"),
            columns=pd.Index([3, 4], name="b", dtype="object"),
        )
        col_normal_margins = pd.DataFrame(
            [[0.5, 0, 0.2], [0.5, 1.0, 0.8]],
            index=pd.Index([1, 2], name="a", dtype="object"),
            columns=pd.Index([3, 4, "All"], name="b", dtype="object"),
        )

        all_normal_margins = pd.DataFrame(
            [[0.2, 0, 0.2], [0.2, 0.6, 0.8], [0.4, 0.6, 1]],
            index=pd.Index([1, 2, "All"], name="a", dtype="object"),
            columns=pd.Index([3, 4, "All"], name="b", dtype="object"),
        )
        tm.assert_frame_equal(
            pd.crosstab(df.a, df.b, normalize="index", margins=True), row_normal_margins
        )
        tm.assert_frame_equal(
            pd.crosstab(df.a, df.b, normalize="columns", margins=True),
            col_normal_margins,
        )
        tm.assert_frame_equal(
            pd.crosstab(df.a, df.b, normalize=True, margins=True), all_normal_margins
        )

        # Test arrays
        pd.crosstab(
            [np.array([1, 1, 2, 2]), np.array([1, 2, 1, 2])], np.array([1, 2, 1, 2])
        )

        # Test with aggfunc
        norm_counts = pd.DataFrame(
            [[0.25, 0, 0.25], [0.25, 0.5, 0.75], [0.5, 0.5, 1]],
            index=pd.Index([1, 2, "All"], name="a", dtype="object"),
            columns=pd.Index([3, 4, "All"], name="b"),
        )
        test_case = pd.crosstab(
            df.a, df.b, df.c, aggfunc="count", normalize="all", margins=True
        )
        tm.assert_frame_equal(test_case, norm_counts)

        df = pd.DataFrame(
            {"a": [1, 2, 2, 2, 2], "b": [3, 3, 4, 4, 4], "c": [0, 4, np.nan, 3, 3]}
        )

        norm_sum = pd.DataFrame(
            [[0, 0, 0.0], [0.4, 0.6, 1], [0.4, 0.6, 1]],
            index=pd.Index([1, 2, "All"], name="a", dtype="object"),
            columns=pd.Index([3, 4, "All"], name="b", dtype="object"),
        )
        test_case = pd.crosstab(
            df.a, df.b, df.c, aggfunc=np.sum, normalize="all", margins=True
        )
        tm.assert_frame_equal(test_case, norm_sum)

    def test_crosstab_with_empties(self):
        # Check handling of empties
        df = pd.DataFrame(
            {
                "a": [1, 2, 2, 2, 2],
                "b": [3, 3, 4, 4, 4],
                "c": [np.nan, np.nan, np.nan, np.nan, np.nan],
            }
        )

        empty = pd.DataFrame(
            [[0.0, 0.0], [0.0, 0.0]],
            index=pd.Index([1, 2], name="a", dtype="int64"),
            columns=pd.Index([3, 4], name="b"),
        )

        for i in [True, "index", "columns"]:
            calculated = pd.crosstab(
                df.a, df.b, values=df.c, aggfunc="count", normalize=i
            )
            tm.assert_frame_equal(empty, calculated)

        nans = pd.DataFrame(
            [[0.0, np.nan], [0.0, 0.0]],
            index=pd.Index([1, 2], name="a", dtype="int64"),
            columns=pd.Index([3, 4], name="b"),
        )

        calculated = pd.crosstab(
            df.a, df.b, values=df.c, aggfunc="count", normalize=False
        )
        tm.assert_frame_equal(nans, calculated)

    def test_crosstab_errors(self):
        # Issue 12578

        df = pd.DataFrame(
            {"a": [1, 2, 2, 2, 2], "b": [3, 3, 4, 4, 4], "c": [1, 1, np.nan, 1, 1]}
        )

        error = "values cannot be used without an aggfunc."
        with pytest.raises(ValueError, match=error):
            pd.crosstab(df.a, df.b, values=df.c)

        error = "aggfunc cannot be used without values"
        with pytest.raises(ValueError, match=error):
            pd.crosstab(df.a, df.b, aggfunc=np.mean)

        error = "Not a valid normalize argument"
        with pytest.raises(ValueError, match=error):
            pd.crosstab(df.a, df.b, normalize="42")

        with pytest.raises(ValueError, match=error):
            pd.crosstab(df.a, df.b, normalize=42)

        error = "Not a valid margins argument"
        with pytest.raises(ValueError, match=error):
            pd.crosstab(df.a, df.b, normalize="all", margins=42)

    def test_crosstab_with_categorial_columns(self):
        # GH 8860
        df = pd.DataFrame(
            {
                "MAKE": ["Honda", "Acura", "Tesla", "Honda", "Honda", "Acura"],
                "MODEL": ["Sedan", "Sedan", "Electric", "Pickup", "Sedan", "Sedan"],
            }
        )
        categories = ["Sedan", "Electric", "Pickup"]
        df["MODEL"] = df["MODEL"].astype("category").cat.set_categories(categories)
        result = pd.crosstab(df["MAKE"], df["MODEL"])

        expected_index = pd.Index(["Acura", "Honda", "Tesla"], name="MAKE")
        expected_columns = pd.CategoricalIndex(
            categories, categories=categories, ordered=False, name="MODEL"
        )
        expected_data = [[2, 0, 0], [2, 0, 1], [0, 1, 0]]
        expected = pd.DataFrame(
            expected_data, index=expected_index, columns=expected_columns
        )
        tm.assert_frame_equal(result, expected)

    def test_crosstab_with_numpy_size(self):
        # GH 4003
        df = pd.DataFrame(
            {
                "A": ["one", "one", "two", "three"] * 6,
                "B": ["A", "B", "C"] * 8,
                "C": ["foo", "foo", "foo", "bar", "bar", "bar"] * 4,
                "D": np.random.randn(24),
                "E": np.random.randn(24),
            }
        )
        result = pd.crosstab(
            index=[df["A"], df["B"]],
            columns=[df["C"]],
            margins=True,
            aggfunc=np.size,
            values=df["D"],
        )
        expected_index = pd.MultiIndex(
            levels=[["All", "one", "three", "two"], ["", "A", "B", "C"]],
            codes=[[1, 1, 1, 2, 2, 2, 3, 3, 3, 0], [1, 2, 3, 1, 2, 3, 1, 2, 3, 0]],
            names=["A", "B"],
        )
        expected_column = pd.Index(["bar", "foo", "All"], dtype="object", name="C")
        expected_data = np.array(
            [
                [2.0, 2.0, 4.0],
                [2.0, 2.0, 4.0],
                [2.0, 2.0, 4.0],
                [2.0, np.nan, 2.0],
                [np.nan, 2.0, 2.0],
                [2.0, np.nan, 2.0],
                [np.nan, 2.0, 2.0],
                [2.0, np.nan, 2.0],
                [np.nan, 2.0, 2.0],
                [12.0, 12.0, 24.0],
            ]
        )
        expected = pd.DataFrame(
            expected_data, index=expected_index, columns=expected_column
        )
        tm.assert_frame_equal(result, expected)

    def test_crosstab_dup_index_names(self):
        # GH 13279
        s = pd.Series(range(3), name="foo")

        result = pd.crosstab(s, s)
        expected_index = pd.Index(range(3), name="foo")
        expected = pd.DataFrame(
            np.eye(3, dtype=np.int64), index=expected_index, columns=expected_index
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("names", [["a", ("b", "c")], [("a", "b"), "c"]])
    def test_crosstab_tuple_name(self, names):
        s1 = pd.Series(range(3), name=names[0])
        s2 = pd.Series(range(1, 4), name=names[1])

        mi = pd.MultiIndex.from_arrays([range(3), range(1, 4)], names=names)
        expected = pd.Series(1, index=mi).unstack(1, fill_value=0)

        result = pd.crosstab(s1, s2)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_both_tuple_names(self):
        # GH 18321
        s1 = pd.Series(range(3), name=("a", "b"))
        s2 = pd.Series(range(3), name=("c", "d"))

        expected = pd.DataFrame(
            np.eye(3, dtype="int64"),
            index=pd.Index(range(3), name=("a", "b")),
            columns=pd.Index(range(3), name=("c", "d")),
        )
        result = crosstab(s1, s2)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_unsorted_order(self):
        df = pd.DataFrame({"b": [3, 1, 2], "a": [5, 4, 6]}, index=["C", "A", "B"])
        result = pd.crosstab(df.index, [df.b, df.a])
        e_idx = pd.Index(["A", "B", "C"], name="row_0")
        e_columns = pd.MultiIndex.from_tuples(
            [(1, 4), (2, 6), (3, 5)], names=["b", "a"]
        )
        expected = pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=e_idx, columns=e_columns
        )
        tm.assert_frame_equal(result, expected)

    def test_margin_normalize(self):
        # GH 27500
        df = pd.DataFrame(
            {
                "A": ["foo", "foo", "foo", "foo", "foo", "bar", "bar", "bar", "bar"],
                "B": ["one", "one", "one", "two", "two", "one", "one", "two", "two"],
                "C": [
                    "small",
                    "large",
                    "large",
                    "small",
                    "small",
                    "large",
                    "small",
                    "small",
                    "large",
                ],
                "D": [1, 2, 2, 3, 3, 4, 5, 6, 7],
                "E": [2, 4, 5, 5, 6, 6, 8, 9, 9],
            }
        )
        # normalize on index
        result = pd.crosstab(
            [df.A, df.B], df.C, margins=True, margins_name="Sub-Total", normalize=0
        )
        expected = pd.DataFrame(
            [[0.5, 0.5], [0.5, 0.5], [0.666667, 0.333333], [0, 1], [0.444444, 0.555556]]
        )
        expected.index = MultiIndex(
            levels=[["Sub-Total", "bar", "foo"], ["", "one", "two"]],
            codes=[[1, 1, 2, 2, 0], [1, 2, 1, 2, 0]],
            names=["A", "B"],
        )
        expected.columns = Index(["large", "small"], dtype="object", name="C")
        tm.assert_frame_equal(result, expected)

        # normalize on columns
        result = pd.crosstab(
            [df.A, df.B], df.C, margins=True, margins_name="Sub-Total", normalize=1
        )
        expected = pd.DataFrame(
            [
                [0.25, 0.2, 0.222222],
                [0.25, 0.2, 0.222222],
                [0.5, 0.2, 0.333333],
                [0, 0.4, 0.222222],
            ]
        )
        expected.columns = Index(
            ["large", "small", "Sub-Total"], dtype="object", name="C"
        )
        expected.index = MultiIndex(
            levels=[["bar", "foo"], ["one", "two"]],
            codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
            names=["A", "B"],
        )
        tm.assert_frame_equal(result, expected)

        # normalize on both index and column
        result = pd.crosstab(
            [df.A, df.B], df.C, margins=True, margins_name="Sub-Total", normalize=True
        )
        expected = pd.DataFrame(
            [
                [0.111111, 0.111111, 0.222222],
                [0.111111, 0.111111, 0.222222],
                [0.222222, 0.111111, 0.333333],
                [0.000000, 0.222222, 0.222222],
                [0.444444, 0.555555, 1],
            ]
        )
        expected.columns = Index(
            ["large", "small", "Sub-Total"], dtype="object", name="C"
        )
        expected.index = MultiIndex(
            levels=[["Sub-Total", "bar", "foo"], ["", "one", "two"]],
            codes=[[1, 1, 2, 2, 0], [1, 2, 1, 2, 0]],
            names=["A", "B"],
        )
        tm.assert_frame_equal(result, expected)
