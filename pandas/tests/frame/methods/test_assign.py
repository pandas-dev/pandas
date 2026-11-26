import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestAssign:
    def test_assign(self):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        original = df.copy()
        result = df.assign(C=df.B / df.A)
        expected = df.copy()
        expected["C"] = [4, 2.5, 2]
        tm.assert_frame_equal(result, expected)

        # lambda syntax
        result = df.assign(C=lambda x: x.B / x.A)
        tm.assert_frame_equal(result, expected)

        # original is unmodified
        tm.assert_frame_equal(df, original)

        # Non-Series array-like
        result = df.assign(C=[4, 2.5, 2])
        tm.assert_frame_equal(result, expected)
        # original is unmodified
        tm.assert_frame_equal(df, original)

        result = df.assign(B=df.B / df.A)
        expected = expected.drop("B", axis=1).rename(columns={"C": "B"})
        tm.assert_frame_equal(result, expected)

        # overwrite
        result = df.assign(A=df.A + df.B)
        expected = df.copy()
        expected["A"] = [5, 7, 9]
        tm.assert_frame_equal(result, expected)

        # lambda
        result = df.assign(A=lambda x: x.A + x.B)
        tm.assert_frame_equal(result, expected)

    def test_assign_multiple(self):
        df = pd.DataFrame([[1, 4], [2, 5], [3, 6]], columns=["A", "B"])
        result = df.assign(C=[7, 8, 9], D=df.A, E=lambda x: x.B)
        expected = pd.DataFrame(
            [[1, 4, 7, 1, 4], [2, 5, 8, 2, 5], [3, 6, 9, 3, 6]], columns=list("ABCDE")
        )
        tm.assert_frame_equal(result, expected)

    def test_assign_order(self):
        # GH 9818
        df = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        result = df.assign(D=df.A + df.B, C=df.A - df.B)

        expected = pd.DataFrame([[1, 2, 3, -1], [3, 4, 7, -1]], columns=list("ABDC"))
        tm.assert_frame_equal(result, expected)
        result = df.assign(C=df.A - df.B, D=df.A + df.B)

        expected = pd.DataFrame([[1, 2, -1, 3], [3, 4, -1, 7]], columns=list("ABCD"))

        tm.assert_frame_equal(result, expected)

    def test_assign_bad(self):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

        # non-keyword argument
        msg = r"assign\(\) takes 1 positional argument but 2 were given"
        with pytest.raises(TypeError, match=msg):
            df.assign(lambda x: x.A)
        msg = "'DataFrame' object has no attribute 'C'"
        with pytest.raises(AttributeError, match=msg):
            df.assign(C=df.A, D=df.A + df.C)

    def test_assign_dependent(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})

        result = df.assign(C=df.A, D=lambda x: x["A"] + x["C"])
        expected = pd.DataFrame([[1, 3, 1, 2], [2, 4, 2, 4]], columns=list("ABCD"))
        tm.assert_frame_equal(result, expected)

        result = df.assign(C=lambda df: df.A, D=lambda df: df["A"] + df["C"])
        expected = pd.DataFrame([[1, 3, 1, 2], [2, 4, 2, 4]], columns=list("ABCD"))
        tm.assert_frame_equal(result, expected)

    def test_assign_many_columns_not_fragmented(self):
        # GH#60075
        df = pd.DataFrame({"A": [1.0, 2.0]})
        values = {f"c{i}": np.array([i, i + 0.5]) for i in range(150)}
        with tm.assert_produces_warning(None):
            result = df.assign(**values)
        assert result._mgr.nblocks == 2
        expected = pd.DataFrame({"A": [1.0, 2.0], **values})
        tm.assert_frame_equal(result, expected)

    def test_assign_runs_between_callable_and_existing(self):
        # GH#60075 a callable sees the columns before it, an existing column
        #  keeps its position, and object column labels stay object
        cols = pd.Index(["A"], dtype=object)
        df = pd.DataFrame({"A": [1, 2]}, columns=cols)
        result = df.assign(
            B=[3, 4], C=[0, 0], D=lambda x: x["B"] * 2, A=[5, 6], E=7, F=8.5
        )
        expected = pd.DataFrame(
            {"A": [5, 6], "B": [3, 4], "C": [0, 0], "D": [6, 8], "E": 7, "F": 8.5},
            columns=pd.Index(list("ABCDEF"), dtype=object),
        )
        tm.assert_frame_equal(result, expected)

    def test_assign_bytes_columns_keep_dtype(self):
        # GH#60075 adding several columns at once must not store bytes as object
        df = pd.DataFrame({"A": [1, 2]})
        arr = np.array([b"x", b"y"])
        result = df.assign(B=arr, C=arr)
        expected = pd.Series(
            [np.dtype(np.int64), arr.dtype, arr.dtype], index=list("ABC")
        )
        tm.assert_series_equal(result.dtypes, expected)

    @pytest.mark.parametrize(
        "columns",
        [
            pd.RangeIndex(1),
            pd.Index(["A"], dtype=object, name="n"),
            pd.Index(["A"], dtype=pd.StringDtype(na_value=np.nan), name="n"),
            pd.Index(["A"], dtype="string[python]"),
            pd.CategoricalIndex(["A"], categories=["A", "B", "C", "D"]),
        ],
    )
    def test_assign_many_columns_labels_match_setitem(self, columns):
        # GH#60075
        df = pd.DataFrame([[1.0], [2.0]], columns=columns)
        result = df.assign(B=[3.0, 4.0], C=5.0, D=np.array([6.0, 7.0]))
        expected = df.copy()
        expected["B"] = [3.0, 4.0]
        expected["C"] = 5.0
        expected["D"] = np.array([6.0, 7.0])
        tm.assert_frame_equal(result, expected, check_column_type=True)
        assert result.columns.name == expected.columns.name

    @pytest.mark.parametrize(
        "columns, keys",
        [
            (pd.DatetimeIndex(["2020-01-01"]), ["2020-01-02", "2020-01-02 00:00:00"]),
            (pd.TimedeltaIndex(["1 day"]), ["2 days", "48h"]),
        ],
    )
    def test_assign_keys_matching_same_label(self, columns, keys):
        # GH#60075 the second key overwrites the column the first one added
        df = pd.DataFrame([[1.0]], columns=columns)
        result = df.assign(**{keys[0]: 2.0, keys[1]: 3.0})
        expected = df.copy()
        expected[keys[0]] = 2.0
        expected[keys[1]] = 3.0
        assert result.shape == (1, 2)
        tm.assert_frame_equal(result, expected)

    def test_assign_numeric_str_key_arrow_labels(self):
        # GH#60075 a numeric str key can keep Arrow numeric labels' dtype
        pytest.importorskip("pyarrow")
        df = pd.DataFrame([[1.0]], columns=pd.Index([1], dtype="int64[pyarrow]"))
        result = df.assign(**{"0": 2.0, "x": 3.0})
        expected = df.copy()
        expected["0"] = 2.0
        expected["x"] = 3.0
        tm.assert_frame_equal(result, expected)

    @pytest.mark.filterwarnings("ignore::pandas.errors.Pandas4Warning")
    def test_assign_multiindex_columns_pads_each_new_key(self):
        # GH#60075 MultiIndex columns keep adding one padded key at a time
        mi = pd.MultiIndex.from_tuples([("a", "b")])
        result = pd.DataFrame([[1], [2]], columns=mi).assign(x=[1, 2], y=[3, 4])
        assert result.columns.tolist() == [("a", "b"), ("x", ""), ("y", "")]

    def test_assign_empty_frame_takes_index_from_values(self):
        # GH#60075
        result = pd.DataFrame().assign(a=1, b=[1, 2])
        expected = pd.DataFrame()
        expected["a"] = 1
        expected["b"] = [1, 2]
        tm.assert_frame_equal(result, expected)

    def test_assign_uses_subclass_setitem(self):
        # GH#60075 a subclass overriding __setitem__ still has it called per key
        class UpperKeys(pd.DataFrame):
            @property
            def _constructor(self):
                return UpperKeys

            def __setitem__(self, key, value):
                super().__setitem__(str(key).upper(), value)

        result = UpperKeys({"a": [1, 2]}).assign(x=[1, 2], y=[3, 4])
        assert result.columns.tolist() == ["a", "X", "Y"]
