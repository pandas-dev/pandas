import numpy as np
import pytest

import pandas._libs.index as libindex
from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.boolean import BooleanDtype


class TestMultiIndexBasic:
    def test_multiindex_perf_warn(self, performance_warning):
        df = pd.DataFrame(
            {
                "jim": [0, 0, 1, 1],
                "joe": ["x", "x", "z", "y"],
                "jolie": np.random.default_rng(2).random(4),
            }
        ).set_index(["jim", "joe"])

        # Full-length key uses the engine directly; no PerformanceWarning
        df.loc[(1, "z")]

        # Partial key still goes through lexsort path
        df = df.iloc[[2, 1, 3, 0]]
        msg = "indexing past lexsort depth may impact performance"
        with tm.assert_produces_warning(performance_warning, match=msg):
            df.loc[(0,)]

    @pytest.mark.parametrize("offset", [-5, 5])
    def test_indexing_over_hashtable_size_cutoff(self, monkeypatch, offset):
        size_cutoff = 20
        n = size_cutoff + offset

        with monkeypatch.context():
            monkeypatch.setattr(libindex, "_SIZE_CUTOFF", size_cutoff)
            s = pd.Series(
                np.arange(n), pd.MultiIndex.from_arrays((["a"] * n, np.arange(n)))
            )

            # hai it works!
            assert s[("a", 5)] == 5
            assert s[("a", 6)] == 6
            assert s[("a", 7)] == 7

    def test_multi_nan_indexing(self):
        # GH 3588
        df = pd.DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            }
        )
        result = df.set_index(["a", "b"], drop=False)
        expected = pd.DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            },
            index=[
                pd.Index(["R1", "R2", np.nan, "R4"], name="a"),
                pd.Index(["C1", "C2", "C3", "C4"], name="b"),
            ],
        )
        tm.assert_frame_equal(result, expected)

    def test_exclusive_nat_column_indexing(self):
        # GH 38025
        # test multi indexing when one column exclusively contains NaT values
        df = pd.DataFrame(
            {
                "a": [pd.NaT, pd.NaT, pd.NaT, pd.NaT],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            }
        )
        df = df.set_index(["a", "b"])
        expected = pd.DataFrame(
            {
                "c": [10, 15, np.nan, 20],
            },
            index=[
                pd.Index([pd.NaT, pd.NaT, pd.NaT, pd.NaT], name="a"),
                pd.Index(["C1", "C2", "C3", "C4"], name="b"),
            ],
        )
        tm.assert_frame_equal(df, expected)

    def test_nested_tuples_duplicates(self):
        # GH#30892

        dti = pd.to_datetime(["20190101", "20190101", "20190102"])
        idx = pd.Index(["a", "a", "c"])
        mi = pd.MultiIndex.from_arrays([dti, idx], names=["index1", "index2"])

        df = pd.DataFrame({"c1": [1, 2, 3], "c2": [np.nan, np.nan, np.nan]}, index=mi)

        expected = pd.DataFrame({"c1": df["c1"], "c2": [1.0, 1.0, np.nan]}, index=mi)

        df2 = df.copy(deep=True)
        df2.loc[(dti[0], "a"), "c2"] = 1.0
        tm.assert_frame_equal(df2, expected)

        df3 = df.copy(deep=True)
        df3.loc[[(dti[0], "a")], "c2"] = 1.0
        tm.assert_frame_equal(df3, expected)

    def test_multiindex_with_datatime_level_preserves_freq(self):
        # https://github.com/pandas-dev/pandas/issues/35563
        idx = pd.Index(range(2), name="A")
        dti = pd.date_range("2020-01-01", periods=7, freq="D", name="B")
        mi = pd.MultiIndex.from_product([idx, dti])
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((14, 2)), index=mi)
        result = df.loc[0].index
        tm.assert_index_equal(result, dti)
        assert result.freq == dti.freq

    def test_multiindex_complex(self):
        # GH#42145
        complex_data = [1 + 2j, 4 - 3j, 10 - 1j]
        non_complex_data = [3, 4, 5]
        result = pd.DataFrame(
            {
                "x": complex_data,
                "y": non_complex_data,
                "z": non_complex_data,
            }
        )
        result.set_index(["x", "y"], inplace=True)
        expected = pd.DataFrame(
            {"z": non_complex_data},
            index=pd.MultiIndex.from_arrays(
                [complex_data, non_complex_data],
                names=("x", "y"),
            ),
        )
        tm.assert_frame_equal(result, expected)

    def test_rename_multiindex_with_duplicates(self):
        # GH 38015
        mi = pd.MultiIndex.from_tuples([("A", "cat"), ("B", "cat"), ("B", "cat")])
        df = pd.DataFrame(index=mi)
        df = df.rename(index={"A": "Apple"}, level=0)

        mi2 = pd.MultiIndex.from_tuples([("Apple", "cat"), ("B", "cat"), ("B", "cat")])
        expected = pd.DataFrame(index=mi2)
        tm.assert_frame_equal(df, expected)

    def test_series_align_multiindex_with_nan_overlap_only(self):
        # GH 38439
        mi1 = pd.MultiIndex.from_arrays([[81.0, np.nan], [np.nan, np.nan]])
        mi2 = pd.MultiIndex.from_arrays([[np.nan, 82.0], [np.nan, np.nan]])
        ser1 = pd.Series([1, 2], index=mi1)
        ser2 = pd.Series([1, 2], index=mi2)
        result1, result2 = ser1.align(ser2)

        mi = pd.MultiIndex.from_arrays([[81.0, 82.0, np.nan], [np.nan, np.nan, np.nan]])
        expected1 = pd.Series([1.0, np.nan, 2.0], index=mi)
        expected2 = pd.Series([np.nan, 2.0, 1.0], index=mi)

        tm.assert_series_equal(result1, expected1)
        tm.assert_series_equal(result2, expected2)

    def test_series_align_multiindex_with_nan(self):
        # GH 38439
        mi1 = pd.MultiIndex.from_arrays([[81.0, np.nan], [np.nan, np.nan]])
        mi2 = pd.MultiIndex.from_arrays([[np.nan, 81.0], [np.nan, np.nan]])
        ser1 = pd.Series([1, 2], index=mi1)
        ser2 = pd.Series([1, 2], index=mi2)
        result1, result2 = ser1.align(ser2)

        mi = pd.MultiIndex.from_arrays([[81.0, np.nan], [np.nan, np.nan]])
        expected1 = pd.Series([1, 2], index=mi)
        expected2 = pd.Series([2, 1], index=mi)

        tm.assert_series_equal(result1, expected1)
        tm.assert_series_equal(result2, expected2)

    def test_nunique_smoke(self):
        # GH 34019
        n = pd.DataFrame([[1, 2], [1, 2]]).set_index([0, 1]).index.nunique()
        assert n == 1

    def test_multiindex_repeated_keys(self):
        # GH19414
        tm.assert_series_equal(
            pd.Series([1, 2], pd.MultiIndex.from_arrays([["a", "b"]])).loc[
                ["a", "a", "b", "b"]
            ],
            pd.Series([1, 1, 2, 2], pd.MultiIndex.from_arrays([["a", "a", "b", "b"]])),
        )

    def test_multiindex_with_na_missing_key(self):
        # GH46173
        df = pd.DataFrame.from_dict(
            {
                ("foo",): [1, 2, 3],
                ("bar",): [5, 6, 7],
                (None,): [8, 9, 0],
            }
        )
        with pytest.raises(KeyError, match="missing_key"):
            df[[("missing_key",)]]

    def test_multiindex_dtype_preservation(self):
        # GH51261
        columns = pd.MultiIndex.from_tuples([("A", "B")], names=["lvl1", "lvl2"])
        df = pd.DataFrame(["value"], columns=columns).astype("category")
        df_no_multiindex = df["A"]
        assert isinstance(df_no_multiindex["B"].dtype, pd.CategoricalDtype)

        # geopandas 1763 analogue
        df = pd.DataFrame(
            [[1, 0], [0, 1]],
            columns=[
                ["foo", "foo"],
                ["location", "location"],
                ["x", "y"],
            ],
        )
        with tm.assert_produces_warning(Pandas4Warning, match="Setting a new column"):
            df = df.assign(bools=pd.Series([True, False], dtype="boolean"))
        assert isinstance(df["bools"].dtype, BooleanDtype)

    def test_multiindex_from_tuples_with_nan(self):
        # GH#23578
        result = pd.MultiIndex.from_tuples([("a", "b", "c"), np.nan, ("d", "", "")])
        expected = pd.MultiIndex.from_tuples(
            [("a", "b", "c"), (np.nan, np.nan, np.nan), ("d", "", "")]
        )
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("operation", ["div", "mul", "add", "sub"])
    def test_groupyby_rename_categories_operation_with_multiindex(self, operation):
        # GH#51500
        data = pd.DataFrame(
            [["C", "B", "B"], ["B", "A", "A"], ["B", "A", "B"]], columns=["0", "1", "2"]
        )
        data["0"] = data["0"].astype("category")
        data["0"] = data["0"].cat.rename_categories({"C": "B", "B": "C"})

        a = data.groupby(by=["0", "1"])["2"].value_counts()
        b = data.groupby(by=["0", "1"]).size()

        result = getattr(a, operation)(b)
        expected = getattr(a, operation)(b.sort_index(ascending=False))

        tm.assert_series_equal(result, expected)

    def test_multiindex_assign_aligns_as_implicit_tuple(self):
        # GH 61841
        cols = pd.MultiIndex.from_tuples([("A", "B")])
        df1 = pd.DataFrame([[i] for i in range(3)], columns=cols)
        df2 = df1.copy()
        df3 = df1.copy()
        s1 = df1["A"].rolling(2).mean()
        s2 = s1.copy()
        s3 = s1.copy()

        msg = "Setting a new column"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df2["C"] = s2
        df3[("C", "")] = s3
        tm.assert_frame_equal(df2, df3)

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df1["C"] = s1
        tm.assert_frame_equal(df1, df2)
        tm.assert_frame_equal(df1, df3)

        # setting again, now an existing column, does not warn
        with tm.assert_produces_warning(None):
            df1["C"] = s1
        tm.assert_frame_equal(df1, df2)
        tm.assert_frame_equal(df1, df3)

        df1["C"] = s1 * 2
        tm.assert_series_equal(df1["C"], df2["C"] * 2)

    def test_multiindex_assign_alignment_with_non_string_dtype(self):
        # GH 62518
        columns = pd.MultiIndex.from_arrays(
            [["a", "a", "z", "z"], pd.Categorical([1, 2, 1, 2])]
        )

        meta = pd.DataFrame(columns=columns, dtype=object)
        meta["z"] = meta["z"].astype("int64")

        result = pd.DataFrame(
            data={
                ("a", 1): pd.Series([], dtype=object),
                ("a", 2): pd.Series([], dtype=object),
                ("z", 1): pd.Series([], dtype="int64"),
                ("z", 2): pd.Series([], dtype="int64"),
            },
            columns=columns,
        )

        tm.assert_frame_equal(meta, result)

    def test_multiindex_assign_alignment_with_object_dtype(self):
        # https://github.com/pandas-dev/pandas/issues/65118
        # second level of the multiindex is object dtype, and the value of
        # that level is 0 for the single column we are setting

        columns = pd.MultiIndex.from_tuples(
            [("A", "M"), ("A", 0), ("B", 1), ("B", 2), ("C", 0)]
        )
        df = pd.DataFrame(np.arange(20, dtype=float).reshape(4, 5), columns=columns)

        df["A"] = df["A"] / 100
        df["B"] = df["B"] / 100
        # this case specifically was the buggy one
        df["C"] = df["C"] / 100

        expected = pd.DataFrame(
            np.arange(20, dtype=float).reshape(4, 5) / 100, columns=columns
        )
        tm.assert_frame_equal(df, expected)
