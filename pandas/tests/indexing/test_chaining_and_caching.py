from string import ascii_letters

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    Series,
    Timestamp,
    date_range,
    option_context,
)
import pandas._testing as tm


class TestCaching:
    @pytest.mark.parametrize("do_ref", [True, False])
    def test_setitem_cache_updating(self, do_ref):
        # GH 5424
        cont = ["one", "two", "three", "four", "five", "six", "seven"]

        df = DataFrame({"a": cont, "b": cont[3:] + cont[:3], "c": np.arange(7)})

        # ref the cache
        if do_ref:
            df.loc[0, "c"]

        # set it
        df.loc[7, "c"] = 1

        assert df.loc[0, "c"] == 0.0
        assert df.loc[7, "c"] == 1.0

    def test_setitem_cache_updating_slices(self):
        # GH 7084
        # not updating cache on series setting with slices
        expected = DataFrame(
            {"A": [600, 600, 600]}, index=date_range("5/7/2014", "5/9/2014")
        )
        out = DataFrame({"A": [0, 0, 0]}, index=date_range("5/7/2014", "5/9/2014"))
        df = DataFrame({"C": ["A", "A", "A"], "D": [100, 200, 300]})

        # loop through df to update out
        six = Timestamp("5/7/2014")
        eix = Timestamp("5/9/2014")
        for ix, row in df.iterrows():
            out.loc[six:eix, row["C"]] = out.loc[six:eix, row["C"]] + row["D"]

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out["A"], expected["A"])

        # try via a chain indexing
        # this actually works
        out = DataFrame({"A": [0, 0, 0]}, index=date_range("5/7/2014", "5/9/2014"))
        out_original = out.copy()
        for ix, row in df.iterrows():
            v = out[row["C"]][six:eix] + row["D"]
            with tm.raises_chained_assignment_error():
                out[row["C"]][six:eix] = v

        tm.assert_frame_equal(out, out_original)
        tm.assert_series_equal(out["A"], out_original["A"])

        out = DataFrame({"A": [0, 0, 0]}, index=date_range("5/7/2014", "5/9/2014"))
        for ix, row in df.iterrows():
            out.loc[six:eix, row["C"]] += row["D"]

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out["A"], expected["A"])


class TestChaining:
    def test_setitem_chained_setfault(self):
        # GH6026
        data = ["right", "left", "left", "left", "right", "left", "timeout"]

        df = DataFrame({"response": np.array(data)})
        mask = df.response == "timeout"
        with tm.raises_chained_assignment_error():
            df.response[mask] = "none"
        tm.assert_frame_equal(df, DataFrame({"response": data}))

        recarray = np.rec.fromarrays([data], names=["response"])
        df = DataFrame(recarray)
        mask = df.response == "timeout"
        with tm.raises_chained_assignment_error():
            df.response[mask] = "none"
        tm.assert_frame_equal(df, DataFrame({"response": data}))

        df = DataFrame({"response": data, "response1": data})
        df_original = df.copy()
        mask = df.response == "timeout"
        with tm.raises_chained_assignment_error():
            df.response[mask] = "none"
        tm.assert_frame_equal(df, df_original)

        # GH 6056
        expected = DataFrame({"A": [np.nan, "bar", "bah", "foo", "bar"]})
        df = DataFrame({"A": np.array(["foo", "bar", "bah", "foo", "bar"])})
        with tm.raises_chained_assignment_error():
            df["A"].iloc[0] = np.nan
        expected = DataFrame({"A": ["foo", "bar", "bah", "foo", "bar"]})
        result = df.head()
        tm.assert_frame_equal(result, expected)

        df = DataFrame({"A": np.array(["foo", "bar", "bah", "foo", "bar"])})
        with tm.raises_chained_assignment_error():
            df.A.iloc[0] = np.nan
        result = df.head()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment(self):
        with option_context("chained_assignment", "raise"):
            # work with the chain
            df = DataFrame(
                np.arange(4).reshape(2, 2), columns=list("AB"), dtype="int64"
            )
            df_original = df.copy()

            with tm.raises_chained_assignment_error():
                df["A"][0] = -5
            with tm.raises_chained_assignment_error():
                df["A"][1] = -6
            tm.assert_frame_equal(df, df_original)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_raises(self):
        # test with the chaining
        df = DataFrame(
            {
                "A": Series(range(2), dtype="int64"),
                "B": np.array(np.arange(2, 4), dtype=np.float64),
            }
        )
        df_original = df.copy()
        with tm.raises_chained_assignment_error():
            df["A"][0] = -5
        with tm.raises_chained_assignment_error():
            df["A"][1] = -6
        tm.assert_frame_equal(df, df_original)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_fails(self):
        # Using a copy (the chain), fails
        df = DataFrame(
            {
                "A": Series(range(2), dtype="int64"),
                "B": np.array(np.arange(2, 4), dtype=np.float64),
            }
        )

        with tm.raises_chained_assignment_error():
            df.loc[0]["A"] = -5

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_doc_example(self):
        # Doc example
        df = DataFrame(
            {
                "a": ["one", "one", "two", "three", "two", "one", "six"],
                "c": Series(range(7), dtype="int64"),
            }
        )

        indexer = df.a.str.startswith("o")
        with tm.raises_chained_assignment_error():
            df[indexer]["c"] = 42

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_object_dtype(self):
        df = DataFrame(
            {"A": Series(["aaa", "bbb", "ccc"], dtype=object), "B": [1, 2, 3]}
        )
        df_original = df.copy()

        with tm.raises_chained_assignment_error():
            df["A"][0] = 111
        tm.assert_frame_equal(df, df_original)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_is_copy_pickle(self, temp_file):
        # gh-5475: Make sure that is_copy is picked up reconstruction
        df = DataFrame({"A": [1, 2]})

        path = str(temp_file)
        df.to_pickle(path)
        df2 = pd.read_pickle(path)
        df2["B"] = df2["A"]
        df2["B"] = df2["A"]

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_str(self):
        idxs = np.random.default_rng(2).integers(len(ascii_letters), size=(100, 2))
        idxs.sort(axis=1)
        strings = [ascii_letters[x[0] : x[1]] for x in idxs]

        df = DataFrame(strings, columns=["letters"])
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df.loc[indexer, "letters"] = df.loc[indexer, "letters"].apply(str.lower)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_sorting(self):
        df = DataFrame(np.random.default_rng(2).standard_normal((10, 4)))
        ser = df.iloc[:, 0].sort_values()

        tm.assert_series_equal(ser, df.iloc[:, 0].sort_values())
        tm.assert_series_equal(ser, df[0].sort_values())

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_false_positives(self):
        # see gh-6025: false positives
        df = DataFrame({"column1": ["a", "a", "a"], "column2": [4, 8, 9]})
        str(df)

        df["column1"] = df["column1"] + "b"
        str(df)

        df = df[df["column2"] != 8]
        str(df)

        df["column1"] = df["column1"] + "c"
        str(df)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_undefined_column(self):
        # from SO:
        # https://stackoverflow.com/questions/24054495/potential-bug-setting-value-for-undefined-column-using-iloc
        df = DataFrame(np.arange(0, 9), columns=["count"])
        df["group"] = "b"
        df_original = df.copy()
        with tm.raises_chained_assignment_error():
            df.iloc[0:5]["group"] = "a"
        tm.assert_frame_equal(df, df_original)

    @pytest.mark.arm_slow
    def test_detect_chained_assignment_changing_dtype(self):
        # Mixed type setting but same dtype & changing dtype
        df = DataFrame(
            {
                "A": date_range("20130101", periods=5),
                "B": np.random.default_rng(2).standard_normal(5),
                "C": np.arange(5, dtype="int64"),
                "D": ["a", "b", "c", "d", "e"],
            }
        )
        df_original = df.copy()

        with tm.raises_chained_assignment_error():
            df.loc[2]["D"] = "foo"
        with tm.raises_chained_assignment_error():
            df.loc[2]["C"] = "foo"
        tm.assert_frame_equal(df, df_original)
        # TODO: Use tm.raises_chained_assignment_error() when PDEP-6 is enforced
        with pytest.raises(TypeError, match="Invalid value"):
            with tm.raises_chained_assignment_error():
                df["C"][2] = "foo"

    def test_setting_with_copy_bug(self):
        # operating on a copy
        df = DataFrame(
            {"a": list(range(4)), "b": list("ab.."), "c": ["a", "b", np.nan, "d"]}
        )
        df_original = df.copy()
        mask = pd.isna(df.c)
        with tm.raises_chained_assignment_error():
            df[["c"]][mask] = df[["b"]][mask]
        tm.assert_frame_equal(df, df_original)

    def test_setting_with_copy_bug_no_warning(self):
        # invalid warning as we are returning a new object
        # GH 8730
        df1 = DataFrame({"x": Series(["a", "b", "c"]), "y": Series(["d", "e", "f"])})
        df2 = df1[["x"]]

        # this should not raise
        df2["y"] = ["g", "h", "i"]

    def test_detect_chained_assignment_warnings_errors(self):
        df = DataFrame({"A": ["aaa", "bbb", "ccc"], "B": [1, 2, 3]})
        with tm.raises_chained_assignment_error():
            df.loc[0]["A"] = 111

    @pytest.mark.parametrize("rhs", [3, DataFrame({0: [1, 2, 3, 4]})])
    def test_detect_chained_assignment_warning_stacklevel(self, rhs):
        # GH#42570
        df = DataFrame(np.arange(25).reshape(5, 5))
        df_original = df.copy()
        chained = df.loc[:3]
        chained[2] = rhs
        tm.assert_frame_equal(df, df_original)

    def test_chained_getitem_with_lists(self):
        # GH6394
        # Regression in chained getitem indexing with embedded list-like from
        # 0.12

        df = DataFrame({"A": 5 * [np.zeros(3)], "B": 5 * [np.ones(3)]})
        expected = df["A"].iloc[2]
        result = df.loc[2, "A"]
        tm.assert_numpy_array_equal(result, expected)
        result2 = df.iloc[2]["A"]
        tm.assert_numpy_array_equal(result2, expected)
        result3 = df["A"].loc[2]
        tm.assert_numpy_array_equal(result3, expected)
        result4 = df["A"].iloc[2]
        tm.assert_numpy_array_equal(result4, expected)

    def test_cache_updating(self):
        # GH 4939, make sure to update the cache on setitem

        df = DataFrame(
            np.zeros((10, 4)),
            columns=Index(list("ABCD"), dtype=object),
        )
        df["A"]  # cache series
        df.loc["Hello Friend"] = df.iloc[0]
        assert "Hello Friend" in df["A"].index
        assert "Hello Friend" in df["B"].index

    def test_cache_updating2(self):
        # 10264
        df = DataFrame(
            np.zeros((5, 5), dtype="int64"),
            columns=["a", "b", "c", "d", "e"],
            index=range(5),
        )
        df["f"] = 0
        df_orig = df.copy()
        with pytest.raises(ValueError, match="read-only"):
            df.f.values[3] = 1
        tm.assert_frame_equal(df, df_orig)

    def test_iloc_setitem_chained_assignment(self):
        # GH#3970
        with option_context("chained_assignment", None):
            df = DataFrame({"aa": range(5), "bb": [2.2] * 5})
            df["cc"] = 0.0

            ck = [True] * len(df)

            with tm.raises_chained_assignment_error():
                df["bb"].iloc[0] = 0.13

            # GH#3970 this lookup used to break the chained setting to 0.15
            df.iloc[ck]

            with tm.raises_chained_assignment_error():
                df["bb"].iloc[0] = 0.15

            assert df["bb"].iloc[0] == 2.2

    def test_getitem_loc_assignment_slice_state(self):
        # GH 13569
        df = DataFrame({"a": [10, 20, 30]})
        with tm.raises_chained_assignment_error():
            df["a"].loc[4] = 40
        tm.assert_frame_equal(df, DataFrame({"a": [10, 20, 30]}))
        tm.assert_series_equal(df["a"], Series([10, 20, 30], name="a"))
