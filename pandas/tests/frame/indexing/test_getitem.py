import numpy as np
import pytest

from pandas import (
    Categorical,
    CategoricalDtype,
    CategoricalIndex,
    DataFrame,
    MultiIndex,
    Series,
    Timestamp,
    get_dummies,
    period_range,
)
import pandas._testing as tm
from pandas.core.arrays import SparseArray


class TestGetitem:
    def test_getitem_unused_level_raises(self):
        # GH#20410
        mi = MultiIndex(
            levels=[["a_lot", "onlyone", "notevenone"], [1970, ""]],
            codes=[[1, 0], [1, 0]],
        )
        df = DataFrame(-1, index=range(3), columns=mi)

        with pytest.raises(KeyError, match="notevenone"):
            df["notevenone"]

    def test_getitem_periodindex(self):
        rng = period_range("1/1/2000", periods=5)
        df = DataFrame(np.random.randn(10, 5), columns=rng)

        ts = df[rng[0]]
        tm.assert_series_equal(ts, df.iloc[:, 0])

        # GH#1211; smoketest unrelated to the rest of this test
        repr(df)

        ts = df["1/1/2000"]
        tm.assert_series_equal(ts, df.iloc[:, 0])

    def test_getitem_list_of_labels_categoricalindex_cols(self):
        # GH#16115
        cats = Categorical([Timestamp("12-31-1999"), Timestamp("12-31-2000")])

        expected = DataFrame(
            [[1, 0], [0, 1]], dtype="uint8", index=[0, 1], columns=cats
        )
        dummies = get_dummies(cats)
        result = dummies[list(dummies.columns)]
        tm.assert_frame_equal(result, expected)

    def test_getitem_sparse_column_return_type_and_dtype(self):
        # https://github.com/pandas-dev/pandas/issues/23559
        data = SparseArray([0, 1])
        df = DataFrame({"A": data})
        expected = Series(data, name="A")
        result = df["A"]
        tm.assert_series_equal(result, expected)

        # Also check iloc and loc while we're here
        result = df.iloc[:, 0]
        tm.assert_series_equal(result, expected)

        result = df.loc[:, "A"]
        tm.assert_series_equal(result, expected)


class TestGetitemListLike:
    def test_getitem_list_missing_key(self):
        # GH#13822, incorrect error string with non-unique columns when missing
        # column is accessed
        df = DataFrame({"x": [1.0], "y": [2.0], "z": [3.0]})
        df.columns = ["x", "x", "z"]

        # Check that we get the correct value in the KeyError
        with pytest.raises(KeyError, match=r"\['y'\] not in index"):
            df[["x", "y", "z"]]


class TestGetitemCallable:
    def test_getitem_callable(self, float_frame):
        # GH#12533
        result = float_frame[lambda x: "A"]
        expected = float_frame.loc[:, "A"]
        tm.assert_series_equal(result, expected)

        result = float_frame[lambda x: ["A", "B"]]
        expected = float_frame.loc[:, ["A", "B"]]
        tm.assert_frame_equal(result, float_frame.loc[:, ["A", "B"]])

        df = float_frame[:3]
        result = df[lambda x: [True, False, True]]
        expected = float_frame.iloc[[0, 2], :]
        tm.assert_frame_equal(result, expected)

    def test_loc_multiindex_columns_one_level(self):
        # GH#29749
        df = DataFrame([[1, 2]], columns=[["a", "b"]])
        expected = DataFrame([1], columns=[["a"]])

        result = df["a"]
        tm.assert_frame_equal(result, expected)

        result = df.loc[:, "a"]
        tm.assert_frame_equal(result, expected)


class TestGetitemBooleanMask:
    def test_getitem_bool_mask_categorical_index(self):

        df3 = DataFrame(
            {
                "A": np.arange(6, dtype="int64"),
            },
            index=CategoricalIndex(
                [1, 1, 2, 1, 3, 2],
                dtype=CategoricalDtype([3, 2, 1], ordered=True),
                name="B",
            ),
        )
        df4 = DataFrame(
            {
                "A": np.arange(6, dtype="int64"),
            },
            index=CategoricalIndex(
                [1, 1, 2, 1, 3, 2],
                dtype=CategoricalDtype([3, 2, 1], ordered=False),
                name="B",
            ),
        )

        result = df3[df3.index == "a"]
        expected = df3.iloc[[]]
        tm.assert_frame_equal(result, expected)

        result = df4[df4.index == "a"]
        expected = df4.iloc[[]]
        tm.assert_frame_equal(result, expected)

        result = df3[df3.index == 1]
        expected = df3.iloc[[0, 1, 3]]
        tm.assert_frame_equal(result, expected)

        result = df4[df4.index == 1]
        expected = df4.iloc[[0, 1, 3]]
        tm.assert_frame_equal(result, expected)

        # since we have an ordered categorical

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=True,
        #         name='B')
        result = df3[df3.index < 2]
        expected = df3.iloc[[4]]
        tm.assert_frame_equal(result, expected)

        result = df3[df3.index > 1]
        expected = df3.iloc[[]]
        tm.assert_frame_equal(result, expected)

        # unordered
        # cannot be compared

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=False,
        #         name='B')
        msg = "Unordered Categoricals can only compare equality or not"
        with pytest.raises(TypeError, match=msg):
            df4[df4.index < 2]
        with pytest.raises(TypeError, match=msg):
            df4[df4.index > 1]
