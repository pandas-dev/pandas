import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def df():
    df = pd.DataFrame(
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
            "D": np.random.default_rng(2).standard_normal(11),
            "E": np.random.default_rng(2).standard_normal(11),
            "F": np.random.default_rng(2).standard_normal(11),
        }
    )

    return pd.concat([df, df], ignore_index=True)


class TestCrosstab:
    def test_crosstab_single(self, df):
        result = pd.crosstab(df["A"], df["C"])
        expected = df.groupby(["A", "C"]).size().unstack()
        tm.assert_frame_equal(result, expected.fillna(0).astype(np.int64))

    def test_crosstab_multiple(self, df):
        result = pd.crosstab(df["A"], [df["B"], df["C"]])
        expected = df.groupby(["A", "B", "C"]).size()
        expected = expected.unstack("B").unstack("C").fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

        result = pd.crosstab([df["B"], df["C"]], df["A"])
        expected = df.groupby(["B", "C", "A"]).size()
        expected = expected.unstack("A").fillna(0).astype(np.int64)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("box", [np.array, list, tuple])
    def test_crosstab_ndarray(self, box):
        # GH 44076
        a = box(np.random.default_rng(2).integers(0, 5, size=100))
        b = box(np.random.default_rng(2).integers(0, 3, size=100))
        c = box(np.random.default_rng(2).integers(0, 10, size=100))

        df = pd.DataFrame({"a": a, "b": b, "c": c})

        result = pd.crosstab(a, [b, c], rownames=["a"], colnames=("b", "c"))
        expected = pd.crosstab(df["a"], [df["b"], df["c"]])
        tm.assert_frame_equal(result, expected)

        result = pd.crosstab([b, c], a, colnames=["a"], rownames=("b", "c"))
        expected = pd.crosstab([df["b"], df["c"]], df["a"])
        tm.assert_frame_equal(result, expected)

        # assign arbitrary names
        result = pd.crosstab(a, c)
        expected = pd.crosstab(df["a"], df["c"])
        expected.index.names = ["row_0"]
        expected.columns.names = ["col_0"]
        tm.assert_frame_equal(result, expected)

    def test_crosstab_non_aligned(self):
        # GH 17005
        a = pd.Series([0, 1, 1], index=["a", "b", "c"])
        b = pd.Series([3, 4, 3, 4, 3], index=["a", "b", "c", "d", "f"])
        c = np.array([3, 4, 3], dtype=np.int64)

        expected = pd.DataFrame(
            [[1, 0], [1, 1]],
            index=pd.Index([0, 1], name="row_0"),
            columns=pd.Index([3, 4], name="col_0"),
        )

        result = pd.crosstab(a, b)
        tm.assert_frame_equal(result, expected)

        result = pd.crosstab(a, c)
        tm.assert_frame_equal(result, expected)

    def test_crosstab_margins(self):
        a = np.random.default_rng(2).integers(0, 7, size=100)
        b = np.random.default_rng(2).integers(0, 3, size=100)
        c = np.random.default_rng(2).integers(0, 5, size=100)

        df = pd.DataFrame({"a": a, "b": b, "c": c})

        result = pd.crosstab(
            a, [b, c], rownames=["a"], colnames=("b", "c"), margins=True
        )

        assert result.index.names == ("a",)
        assert result.columns.names == ["b", "c"]

        all_cols = result["All", ""]
        exp_cols = df.groupby(["a"]).size().astype("i8")
        # to keep index.name
        exp_margin = pd.Series([len(df)], index=pd.Index(["All"], name="a"))
        exp_cols = pd.concat([exp_cols, exp_margin])
        exp_cols.name = ("All", "")

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.loc["All"]
        exp_rows = df.groupby(["b", "c"]).size().astype("i8")
        exp_rows = pd.concat([exp_rows, pd.Series([len(df)], index=[("All", "")])])
        exp_rows.name = "All"

        exp_rows = exp_rows.reindex(all_rows.index)
        exp_rows = exp_rows.fillna(0).astype(np.int64)
        tm.assert_series_equal(all_rows, exp_rows)

    def test_crosstab_margins_set_margin_name(self):
        # GH 15972
        a = np.random.default_rng(2).integers(0, 7, size=100)
        b = np.random.default_rng(2).integers(0, 3, size=100)
        c = np.random.default_rng(2).integers(0, 5, size=100)

        df = pd.DataFrame({"a": a, "b": b, "c": c})

        result = pd.crosstab(
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
        exp_margin = pd.Series([len(df)], index=pd.Index(["TOTAL"], name="a"))
        exp_cols = pd.concat([exp_cols, exp_margin])
        exp_cols.name = ("TOTAL", "")

        tm.assert_series_equal(all_cols, exp_cols)

        all_rows = result.loc["TOTAL"]
        exp_rows = df.groupby(["b", "c"]).size().astype("i8")
        exp_rows = pd.concat([exp_rows, pd.Series([len(df)], index=[("TOTAL", "")])])
        exp_rows.name = "TOTAL"

        exp_rows = exp_rows.reindex(all_rows.index)
        exp_rows = exp_rows.fillna(0).astype(np.int64)
        tm.assert_series_equal(all_rows, exp_rows)

        msg = "margins_name argument must be a string"
        for margins_name in [666, None, ["a", "b"]]:
            with pytest.raises(ValueError, match=msg):
                pd.crosstab(
                    a,
                    [b, c],
                    rownames=["a"],
                    colnames=("b", "c"),
                    margins=True,
                    margins_name=margins_name,
                )

    def test_crosstab_pass_values(self):
        a = np.random.default_rng(2).integers(0, 7, size=100)
        b = np.random.default_rng(2).integers(0, 3, size=100)
        c = np.random.default_rng(2).integers(0, 5, size=100)
        values = np.random.default_rng(2).standard_normal(100)

        table = pd.crosstab(
            [a, b], c, values, aggfunc="sum", rownames=["foo", "bar"], colnames=["baz"]
        )

        df = pd.DataFrame({"foo": a, "bar": b, "baz": c, "values": values})

        expected = df.pivot_table(
            "values", index=["foo", "bar"], columns="baz", aggfunc="sum"
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
        m = pd.MultiIndex.from_tuples(
            [("one", "dull"), ("one", "shiny"), ("two", "dull"), ("two", "shiny")],
            names=["b", "c"],
        )
        tm.assert_index_equal(res.columns, m)

    def test_crosstab_no_overlap(self):
        # GS 10291

        s1 = pd.Series([1, 2, 3], index=[1, 2, 3])
        s2 = pd.Series([4, 5, 6], index=[4, 5, 6])

        actual = pd.crosstab(s1, s2)
        expected = pd.DataFrame(
            index=pd.Index([], dtype="int64", name="row_0"),
            columns=pd.Index([], dtype="int64", name="col_0"),
        )

        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna(self):
        # GH 12577
        # pivot_table counts null into margin ('All')
        # when margins=true and dropna=true

        df = pd.DataFrame({"a": [1, 2, 2, 2, 2, np.nan], "b": [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [2, 3, 5]])
        expected.index = pd.Index([1.0, 2.0, "All"], name="a")
        expected.columns = pd.Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna2(self):
        df = pd.DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, 2, np.nan], "b": [3, np.nan, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = pd.Index([1.0, 2.0, "All"], name="a")
        expected.columns = pd.Index([3.0, 4.0, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna3(self):
        df = pd.DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, np.nan, 2], "b": [3, 3, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=True)
        expected = pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]])
        expected.index = pd.Index([1.0, 2.0, "All"], name="a")
        expected.columns = pd.Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna4(self):
        # GH 12642
        # _add_margins raises KeyError: Level None not found
        # when margins=True and dropna=False
        # GH: 10772: Keep np.nan in result with dropna=False
        df = pd.DataFrame({"a": [1, 2, 2, 2, 2, np.nan], "b": [3, 3, 4, 4, 4, 4]})
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame([[1, 0, 1], [1, 3, 4], [0, 1, 1], [2, 4, 6]])
        expected.index = pd.Index([1.0, 2.0, np.nan, "All"], name="a")
        expected.columns = pd.Index([3, 4, "All"], name="b")
        tm.assert_frame_equal(actual, expected)

    def test_margin_dropna5(self):
        # GH: 10772: Keep np.nan in result with dropna=False
        df = pd.DataFrame(
            {"a": [1, np.nan, np.nan, np.nan, 2, np.nan], "b": [3, np.nan, 4, 4, 4, 4]}
        )
        actual = pd.crosstab(df.a, df.b, margins=True, dropna=False)
        expected = pd.DataFrame(
            [[1, 0, 0, 1.0], [0, 1, 0, 1.0], [0, 3, 1, 4.0], [1, 4, 1, 6.0]]
        )
        expected.index = pd.Index([1.0, 2.0, np.nan, "All"], name="a")
        expected.columns = pd.Index([3.0, 4.0, np.nan, "All"], name="b")
        tm.assert_frame_equal(actual, expected, check_dtype=False)

    def test_margin_dropna6(self):
        # GH: 10772: Keep np.nan in result with dropna=False
        a = np.array(["foo", "foo", "foo", "bar", "bar", "foo", "foo"], dtype=object)
        b = np.array(["one", "one", "two", "one", "two", np.nan, "two"], dtype=object)
        c = np.array(
            ["dull", "dull", "dull", "dull", "dull", "shiny", "shiny"], dtype=object
        )

        actual = pd.crosstab(
            a, [b, c], rownames=["a"], colnames=["b", "c"], margins=True, dropna=False
        )
        m = pd.MultiIndex.from_arrays(
            [
                ["one", "one", "two", "two", np.nan, np.nan, "All"],
                ["dull", "shiny", "dull", "shiny", "dull", "shiny", ""],
            ],
            names=["b", "c"],
        )
        expected = pd.DataFrame(
            [[1, 0, 1, 0, 0, 0, 2], [2, 0, 1, 1, 0, 1, 5], [3, 0, 2, 1, 0, 1, 7]],
            columns=m,
        )
        expected.index = pd.Index(["bar", "foo", "All"], name="a")
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab(
            [a, b], c, rownames=["a", "b"], colnames=["c"], margins=True, dropna=False
        )
        m = pd.MultiIndex.from_arrays(
            [
                ["bar", "bar", "bar", "foo", "foo", "foo", "All"],
                ["one", "two", np.nan, "one", "two", np.nan, ""],
            ],
            names=["a", "b"],
        )
        expected = pd.DataFrame(
            [
                [1, 0, 1.0],
                [1, 0, 1.0],
                [0, 0, np.nan],
                [2, 0, 2.0],
                [1, 1, 2.0],
                [0, 1, 1.0],
                [5, 2, 7.0],
            ],
            index=m,
        )
        expected.columns = pd.Index(["dull", "shiny", "All"], name="c")
        tm.assert_frame_equal(actual, expected)

        actual = pd.crosstab(
            [a, b], c, rownames=["a", "b"], colnames=["c"], margins=True, dropna=True
        )
        m = pd.MultiIndex.from_arrays(
            [["bar", "bar", "foo", "foo", "All"], ["one", "two", "one", "two", ""]],
            names=["a", "b"],
        )
        expected = pd.DataFrame(
            [[1, 0, 1], [1, 0, 1], [2, 0, 2], [1, 1, 2], [5, 1, 6]], index=m
        )
        expected.columns = pd.Index(["dull", "shiny", "All"], name="c")
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

    def test_crosstab_normalize_arrays(self):
        # GH#12578
        df = pd.DataFrame(
            {"a": [1, 2, 2, 2, 2], "b": [3, 3, 4, 4, 4], "c": [1, 1, np.nan, 1, 1]}
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

    def test_crosstab_with_categorical_columns(self):
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
                "D": np.random.default_rng(2).standard_normal(24),
                "E": np.random.default_rng(2).standard_normal(24),
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
        expected_column = pd.Index(["bar", "foo", "All"], name="C")
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
        # aggfunc is np.size, resulting in integers
        expected["All"] = expected["All"].astype("int64")
        tm.assert_frame_equal(result, expected)

    def test_crosstab_duplicate_names(self):
        # GH 13279 / 22529

        s1 = pd.Series(range(3), name="foo")
        s2_foo = pd.Series(range(1, 4), name="foo")
        s2_bar = pd.Series(range(1, 4), name="bar")
        s3 = pd.Series(range(3), name="waldo")

        # check result computed with duplicate labels against
        # result computed with unique labels, then relabelled
        mapper = {"bar": "foo"}

        # duplicate row, column labels
        result = pd.crosstab(s1, s2_foo)
        expected = pd.crosstab(s1, s2_bar).rename_axis(columns=mapper, axis=1)
        tm.assert_frame_equal(result, expected)

        # duplicate row, unique column labels
        result = pd.crosstab([s1, s2_foo], s3)
        expected = pd.crosstab([s1, s2_bar], s3).rename_axis(index=mapper, axis=0)
        tm.assert_frame_equal(result, expected)

        # unique row, duplicate column labels
        result = pd.crosstab(s3, [s1, s2_foo])
        expected = pd.crosstab(s3, [s1, s2_bar]).rename_axis(columns=mapper, axis=1)

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
        result = pd.crosstab(s1, s2)
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

    def test_crosstab_normalize_multiple_columns(self):
        # GH 15150
        df = pd.DataFrame(
            {
                "A": ["one", "one", "two", "three"] * 6,
                "B": ["A", "B", "C"] * 8,
                "C": ["foo", "foo", "foo", "bar", "bar", "bar"] * 4,
                "D": [0] * 24,
                "E": [0] * 24,
            }
        )

        result = pd.crosstab(
            [df.A, df.B],
            df.C,
            values=df.D,
            aggfunc=np.sum,
            normalize=True,
            margins=True,
        )
        expected = pd.DataFrame(
            np.array([0] * 29 + [1], dtype=float).reshape(10, 3),
            columns=pd.Index(["bar", "foo", "All"], name="C"),
            index=pd.MultiIndex.from_tuples(
                [
                    ("one", "A"),
                    ("one", "B"),
                    ("one", "C"),
                    ("three", "A"),
                    ("three", "B"),
                    ("three", "C"),
                    ("two", "A"),
                    ("two", "B"),
                    ("two", "C"),
                    ("All", ""),
                ],
                names=["A", "B"],
            ),
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
        expected.index = pd.MultiIndex(
            levels=[["Sub-Total", "bar", "foo"], ["", "one", "two"]],
            codes=[[1, 1, 2, 2, 0], [1, 2, 1, 2, 0]],
            names=["A", "B"],
        )
        expected.columns = pd.Index(["large", "small"], name="C")
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
        expected.columns = pd.Index(["large", "small", "Sub-Total"], name="C")
        expected.index = pd.MultiIndex(
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
        expected.columns = pd.Index(["large", "small", "Sub-Total"], name="C")
        expected.index = pd.MultiIndex(
            levels=[["Sub-Total", "bar", "foo"], ["", "one", "two"]],
            codes=[[1, 1, 2, 2, 0], [1, 2, 1, 2, 0]],
            names=["A", "B"],
        )
        tm.assert_frame_equal(result, expected)

    def test_margin_normalize_multiple_columns(self):
        # GH 35144
        # use multiple columns with margins and normalization
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
        result = pd.crosstab(
            index=df.C,
            columns=[df.A, df.B],
            margins=True,
            margins_name="margin",
            normalize=True,
        )
        expected = pd.DataFrame(
            [
                [0.111111, 0.111111, 0.222222, 0.000000, 0.444444],
                [0.111111, 0.111111, 0.111111, 0.222222, 0.555556],
                [0.222222, 0.222222, 0.333333, 0.222222, 1.0],
            ],
            index=["large", "small", "margin"],
        )
        expected.columns = pd.MultiIndex(
            levels=[["bar", "foo", "margin"], ["", "one", "two"]],
            codes=[[0, 0, 1, 1, 2], [1, 2, 1, 2, 0]],
            names=["A", "B"],
        )
        expected.index.name = "C"
        tm.assert_frame_equal(result, expected)

    def test_margin_support_Float(self):
        # GH 50313
        # use Float64 formats and function aggfunc with margins
        df = pd.DataFrame(
            {"A": [1, 2, 2, 1], "B": [3, 3, 4, 5], "C": [-1.0, 10.0, 1.0, 10.0]},
            dtype="Float64",
        )
        result = pd.crosstab(
            df["A"],
            df["B"],
            values=df["C"],
            aggfunc="sum",
            margins=True,
        )
        expected = pd.DataFrame(
            [
                [-1.0, pd.NA, 10.0, 9.0],
                [10.0, 1.0, pd.NA, 11.0],
                [9.0, 1.0, 10.0, 20.0],
            ],
            index=pd.Index([1.0, 2.0, "All"], dtype="object", name="A"),
            columns=pd.Index([3.0, 4.0, 5.0, "All"], dtype="object", name="B"),
            dtype="Float64",
        )
        tm.assert_frame_equal(result, expected)

    def test_margin_with_ordered_categorical_column(self):
        # GH 25278
        df = pd.DataFrame(
            {
                "First": ["B", "B", "C", "A", "B", "C"],
                "Second": ["C", "B", "B", "B", "C", "A"],
            }
        )
        df["First"] = df["First"].astype(pd.CategoricalDtype(ordered=True))
        customized_categories_order = ["C", "A", "B"]
        df["First"] = df["First"].cat.reorder_categories(customized_categories_order)
        result = pd.crosstab(df["First"], df["Second"], margins=True)

        expected_index = pd.Index(["C", "A", "B", "All"], name="First")
        expected_columns = pd.Index(["A", "B", "C", "All"], name="Second")
        expected_data = [[1, 1, 0, 2], [0, 1, 0, 1], [0, 1, 2, 3], [1, 3, 2, 6]]
        expected = pd.DataFrame(
            expected_data, index=expected_index, columns=expected_columns
        )
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("a_dtype", ["category", "int64"])
@pytest.mark.parametrize("b_dtype", ["category", "int64"])
def test_categoricals(a_dtype, b_dtype):
    # https://github.com/pandas-dev/pandas/issues/37465
    g = np.random.default_rng(2)
    a = pd.Series(g.integers(0, 3, size=100)).astype(a_dtype)
    b = pd.Series(g.integers(0, 2, size=100)).astype(b_dtype)
    result = pd.crosstab(a, b, margins=True, dropna=False)
    columns = pd.Index([0, 1, "All"], dtype="object", name="col_0")
    index = pd.Index([0, 1, 2, "All"], dtype="object", name="row_0")
    values = [[10, 18, 28], [23, 16, 39], [17, 16, 33], [50, 50, 100]]
    expected = pd.DataFrame(values, index, columns)
    tm.assert_frame_equal(result, expected)

    # Verify when categorical does not have all values present
    a.loc[a == 1] = 2
    a_is_cat = isinstance(a.dtype, pd.CategoricalDtype)
    assert not a_is_cat or a.value_counts().loc[1] == 0
    result = pd.crosstab(a, b, margins=True, dropna=False)
    values = [[10, 18, 28], [0, 0, 0], [40, 32, 72], [50, 50, 100]]
    expected = pd.DataFrame(values, index, columns)
    if not a_is_cat:
        expected = expected.loc[[0, 2, "All"]]
        expected["All"] = expected["All"].astype("int64")
    tm.assert_frame_equal(result, expected)
