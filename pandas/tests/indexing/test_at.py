from datetime import (
    datetime,
    timezone,
)
from decimal import Decimal

import numpy as np
import pytest

from pandas.errors import (
    InvalidIndexError,
    Pandas4Warning,
)

from pandas import (
    CategoricalDtype,
    CategoricalIndex,
    DataFrame,
    DateOffset,
    DatetimeIndex,
    Index,
    MultiIndex,
    Series,
    Timestamp,
)
import pandas._testing as tm


def test_at_dateoffset_columns():
    # GH#20948 - .at with DateOffset columns
    offsets = Series(data=[-15, -10, -5, 0, 5, 10, 15], dtype=float).map(DateOffset)
    df = DataFrame(index=[0, 1], columns=Index(offsets))

    # read access
    result = df.at[0, offsets[0]]
    assert result is np.nan

    # write access
    df.at[0, offsets[0]] = 1
    assert df.at[0, offsets[0]] == 1


def test_at_multiindex_partial_date_string():
    # GH#43395 - .at with partial date string on MultiIndex with DatetimeIndex level
    timestamps = DatetimeIndex(
        ["2021-08-01", "2021-08-01 12:00", "2021-08-02", "2021-08-02 12:00"]
    )
    index = MultiIndex.from_product(
        [["A", "B"], timestamps], names=["ticker", "timestamp"]
    )
    df = DataFrame({"col": range(len(index))}, index=index)

    # DataFrame.at with partial date string should not raise TypeError
    result = df.at[("A", "2021-08-02"), "col"]
    expected = np.array([2, 3], dtype=np.int64)
    tm.assert_numpy_array_equal(result, expected)

    # Series.at with partial date string
    ser = df["col"]
    result = ser.at[("A", "2021-08-02")]
    expected = Series(
        [2, 3],
        index=MultiIndex.from_arrays(
            [
                ["A", "A"],
                DatetimeIndex(
                    ["2021-08-02", "2021-08-02 12:00"],
                    dtype="datetime64[us]",
                ),
            ],
            names=["ticker", "timestamp"],
        ),
        name="col",
    )
    tm.assert_series_equal(result, expected)


def test_at_incompatible_type_decimal():
    # GH#22740 - .at should not silently discard incompatible type
    df = DataFrame({"A": [1, 2, 3]})
    with pytest.raises(TypeError, match="Invalid value"):
        df.at[0, "A"] = Decimal("1")


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    result = DataFrame({"foo": [datetime(2000, 1, 1)]})
    with pytest.raises(TypeError, match="Invalid value"):
        result.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=timezone.utc)


def test_selection_methods_of_assigned_col():
    # GH 29282
    df = DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
    df2 = DataFrame(data={"c": [7, 8, 9]}, index=[2, 1, 0])
    df["c"] = df2["c"]
    df.at[1, "c"] = 11
    result = df
    expected = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [9, 11, 7]})
    tm.assert_frame_equal(result, expected)
    result = df.at[1, "c"]
    assert result == 11

    result = df["c"]
    expected = Series([9, 11, 7], name="c")
    tm.assert_series_equal(result, expected)

    result = df[["c"]]
    expected = DataFrame({"c": [9, 11, 7]})
    tm.assert_frame_equal(result, expected)


class TestAtSetItem:
    def test_at_setitem_mixed_index_assignment(self):
        # GH#19860
        ser = Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        ser.at["a"] = 11
        assert ser.iat[0] == 11
        ser.at[1] = 22
        assert ser.iat[3] == 22

    def test_at_setitem_categorical_missing(self):
        df = DataFrame(
            index=range(3), columns=range(3), dtype=CategoricalDtype(["foo", "bar"])
        )
        df.at[1, 1] = "foo"

        expected = DataFrame(
            [
                [np.nan, np.nan, np.nan],
                [np.nan, "foo", np.nan],
                [np.nan, np.nan, np.nan],
            ],
            dtype=CategoricalDtype(["foo", "bar"]),
        )

        tm.assert_frame_equal(df, expected)

    def test_at_setitem_multiindex(self):
        df = DataFrame(
            np.zeros((3, 2), dtype="int64"),
            columns=MultiIndex.from_tuples([("a", 0), ("a", 1)]),
        )
        df.at[0, "a"] = 10
        expected = DataFrame(
            [[10, 10], [0, 0], [0, 0]],
            columns=MultiIndex.from_tuples([("a", 0), ("a", 1)]),
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("row", (Timestamp("2019-01-01"), "2019-01-01"))
    def test_at_datetime_index(self, row):
        # Set float64 dtype to avoid upcast when setting .5
        df = DataFrame(
            data=[[1] * 2], index=DatetimeIndex(data=["2019-01-01", "2019-01-02"])
        ).astype({0: "float64"})
        expected = DataFrame(
            data=[[0.5, 1], [1.0, 1]],
            index=DatetimeIndex(data=["2019-01-01", "2019-01-02"]),
        )

        df.at[row, 0] = 0.5
        tm.assert_frame_equal(df, expected)

    def test_at_setitem_list_like_object_dtype_no_warning(self):
        # GH#61223 - list-like in object-dtype column succeeds without warning
        df = DataFrame({"A": Series([[1, 2], [3, 4], [5, 6]])})
        assert df["A"].dtype == np.dtype("object")
        with tm.assert_produces_warning(None):
            df.at[0, "A"] = [10, 20, 30]
        assert df.at[0, "A"] == [10, 20, 30]
        assert df["A"].dtype == np.dtype("object")


class TestAtSetItemWithExpansion:
    def test_at_setitem_expansion_series_dt64tz_value(self, tz_naive_fixture):
        # GH#25506
        ts = (
            Timestamp("2017-08-05 00:00:00+0100", tz=tz_naive_fixture)
            if tz_naive_fixture is not None
            else Timestamp("2017-08-05 00:00:00+0100")
        )
        result = Series(ts)
        with tm.assert_produces_warning(Pandas4Warning, match="does not exist"):
            result.at[1] = ts
        expected = Series([ts, ts])
        tm.assert_series_equal(result, expected)

    def test_at_setitem_expansion_deprecated_dataframe(self):
        # GH#48323
        df = DataFrame({"a": [1, 2]})
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[5, "a"] = 6
        expected = DataFrame({"a": [1, 2, 6]}, index=[0, 1, 5])
        tm.assert_frame_equal(df, expected)

    def test_at_setitem_expansion_deprecated_series(self):
        # GH#48323
        ser = Series([1, 2], index=[0, 1])
        msg = "Setting a value on a Series via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            ser.at[5] = 6
        expected = Series([1, 2, 6], index=[0, 1, 5])
        tm.assert_series_equal(ser, expected)

    def test_at_setitem_expansion_deprecated_new_column(self):
        # GH#48323 - new column (existing row) also expands
        df = DataFrame({"a": [1, 2]})
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[0, "b"] = 99
        assert "b" in df.columns

    def test_at_setitem_expansion_deprecated_new_column_multiindex_rows(self):
        # GH#48323 - new column on a frame with MultiIndex rows
        mi = MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 1)])
        df = DataFrame({"x": [1, 2, 3]}, index=mi)
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[("a", 1), "y"] = 99
        assert "y" in df.columns

    def test_at_setitem_expansion_deprecated_new_column_multiindex_cols(self):
        # GH#48323 - new tuple column on a frame with MultiIndex columns
        df = DataFrame(
            [[1, 2], [3, 4]],
            columns=MultiIndex.from_tuples([("a", 1), ("a", 2)]),
        )
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[0, ("b", 1)] = 99
        assert ("b", 1) in df.columns

    def test_at_setitem_no_warning_existing_key(self):
        # GH#48323 - no warning for existing keys
        df = DataFrame({"a": [1, 2]})
        with tm.assert_produces_warning(None):
            df.at[0, "a"] = 99
        assert df.at[0, "a"] == 99

        ser = Series([1, 2])
        with tm.assert_produces_warning(None):
            ser.at[0] = 99
        assert ser.at[0] == 99

    def test_at_setitem_list_like_missing_row_object_dtype(self):
        # GH#61223 - missing row falls back to .loc with Pandas4Warning
        df = DataFrame({"A": Series([[1], [2]], dtype=object)}, index=["a", "b"])
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at["c", "A"] = [3, 4]
        assert df.at["c", "A"] == [3, 4]
        assert df["A"].dtype == np.dtype("object")

    def test_at_setitem_list_like_new_column(self):
        # GH#61223 - list-like stored as single cell in new object-dtype column
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        with tm.assert_produces_warning(
            Pandas4Warning, match="column key 'B' that does not exist"
        ):
            df.at["a", "B"] = [10, 20, 30]
        assert df.at["a", "B"] == [10, 20, 30]
        assert df["B"].dtype == np.dtype("object")
        assert df.at["b", "B"] is np.nan

    def test_at_setitem_ndarray_new_column(self):
        # GH#61223 - np.ndarray stored as single cell in new object-dtype column
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        arr = np.array([10, 20, 30])
        with tm.assert_produces_warning(
            Pandas4Warning, match="column key 'B' that does not exist"
        ):
            df.at["a", "B"] = arr
        stored = df.at["a", "B"]
        assert isinstance(stored, np.ndarray)
        tm.assert_numpy_array_equal(stored, arr)
        assert df["B"].dtype == np.dtype("object")
        assert df.at["b", "B"] is np.nan

    def test_at_setitem_series_value_new_column(self):
        # GH#61223 - Series stored as single cell in new object-dtype column
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        val = Series([10, 20], index=["x", "y"])
        with tm.assert_produces_warning(
            Pandas4Warning, match="column key 'B' that does not exist"
        ):
            df.at["a", "B"] = val
        result = df.at["a", "B"]
        tm.assert_series_equal(result, val)
        assert df["B"].dtype == np.dtype("object")


class TestAtWithDuplicates:
    def test_at_with_duplicate_axes_requires_scalar_lookup(self):
        # GH#33041 check that falling back to loc doesn't allow non-scalar
        #  args to slip in

        arr = np.random.default_rng(2).standard_normal(6).reshape(3, 2)
        df = DataFrame(arr, columns=["A", "A"])

        msg = "Invalid call for scalar access"
        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]]
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]]
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"]

        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"] = 1

    def test_at_setitem_list_like_nonunique_columns_raises(self):
        # GH#61223 - non-unique column key with list-like raises KeyError
        df = DataFrame([[1, 2]], columns=["A", "A"])
        with pytest.raises(KeyError, match="A"):
            df.at[0, "A"] = [10, 20]

    def test_at_setitem_list_like_nonunique_columns_object_dtype_raises(self):
        # GH#61223 - non-unique column with object dtype also raises KeyError
        df = DataFrame([[[1], [2]]], columns=["A", "A"])
        with pytest.raises(KeyError, match="A"):
            df.at[0, "A"] = [10, 20]


class TestAtErrors:
    # TODO: De-duplicate/parametrize
    #  test_at_series_raises_key_error2, test_at_frame_raises_key_error2

    def test_at_series_raises_key_error(self, indexer_al):
        # GH#31724 .at should match .loc

        ser = Series([1, 2, 3], index=[3, 2, 1])
        result = indexer_al(ser)[1]
        assert result == 3

        with pytest.raises(KeyError, match="a"):
            indexer_al(ser)["a"]

    def test_at_frame_raises_key_error(self, indexer_al):
        # GH#31724 .at should match .loc

        df = DataFrame({0: [1, 2, 3]}, index=[3, 2, 1])

        result = indexer_al(df)[1, 0]
        assert result == 3

        with pytest.raises(KeyError, match="a"):
            indexer_al(df)["a", 0]

        with pytest.raises(KeyError, match="a"):
            indexer_al(df)[1, "a"]

    def test_at_series_raises_key_error2(self, indexer_al):
        # at should not fallback
        # GH#7814
        # GH#31724 .at should match .loc
        ser = Series([1, 2, 3], index=list("abc"))
        result = indexer_al(ser)["a"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            indexer_al(ser)[0]

    def test_at_frame_raises_key_error2(self, indexer_al):
        # GH#31724 .at should match .loc
        df = DataFrame({"A": [1, 2, 3]}, index=list("abc"))
        result = indexer_al(df)["a", "A"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            indexer_al(df)["a", 0]

    def test_at_frame_multiple_columns(self):
        # GH#48296 - at shouldn't modify multiple columns
        df = DataFrame({"a": [1, 2], "b": [3, 4]})
        new_row = [6, 7]
        with pytest.raises(
            InvalidIndexError,
            match=f"You can only assign a scalar value not a \\{type(new_row)}",
        ):
            df.at[5] = new_row

    def test_at_getitem_mixed_index_no_fallback(self):
        # GH#19860
        ser = Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        with pytest.raises(KeyError, match="^0$"):
            ser.at[0]
        with pytest.raises(KeyError, match="^4$"):
            ser.at[4]

    def test_at_categorical_integers(self):
        # CategoricalIndex with integer categories that don't happen to match
        #  the Categorical's codes
        ci = CategoricalIndex([3, 4])

        arr = np.arange(4).reshape(2, 2)
        frame = DataFrame(arr, index=ci)

        for df in [frame, frame.T]:
            for key in [0, 1]:
                with pytest.raises(KeyError, match=str(key)):
                    df.at[key, key]

    def test_at_applied_for_rows(self):
        # GH#48729 .at should raise InvalidIndexError when assigning rows
        df = DataFrame(index=["a"], columns=["col1", "col2"])
        new_row = [123, 15]
        with pytest.raises(
            InvalidIndexError,
            match=f"You can only assign a scalar value not a \\{type(new_row)}",
        ):
            df.at["a"] = new_row

    def test_at_setitem_list_like_wrong_dtype(self):
        # GH#61223 - list-like in non-object dtype column raises ValueError
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at["b", "A"] = [100, 200]

    def test_iat_setitem_list_like(self):
        # GH#61223 - .iat raises ValueError for non-object dtype
        df = DataFrame({"A": [1, 2, 3]})
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.iat[0, 0] = [99, 88, 77]

    def test_at_setitem_tuple_value(self):
        # GH#61223 - tuple raises ValueError for non-object dtype
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at["a", "A"] = (10, 20, 30)

    def test_at_setitem_list_like_multiindex_row(self):
        # GH#61223 - .at with MultiIndex row raises ValueError for non-object dtype
        idx = MultiIndex.from_tuples([("r1", 0), ("r2", 1)])
        df = DataFrame({"A": [1, 2]}, index=idx)
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at[("r1", 0), "A"] = [10, 20]

    def test_at_setitem_list_like_multiindex_partial_key(self):
        # GH#61223 - partial MultiIndex key raises KeyError (.at requires unique row)
        idx = MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 3)])
        df = DataFrame({"A": [[1], [2], [3]]}, index=idx)
        with pytest.raises(KeyError, match="a"):
            df.at["a", "A"] = [10, 20, 30]

    def test_at_setitem_list_like_partial_multiindex_key_matching_len(self):
        # GH#61223 - partial MultiIndex key; new column fires Pandas4Warning
        idx = MultiIndex.from_tuples([("r1", 0), ("r1", 1), ("r2", 0)])
        df = DataFrame({"A": [1, 2, 3]}, index=idx)
        msg_warn = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg_warn):
            with pytest.raises(KeyError, match="r1"):
                df.at["r1", "B"] = [10, 20]

    def test_at_setitem_ndarray_existing_column(self):
        # GH#61223 - np.ndarray in non-object dtype column raises ValueError
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at["b", "A"] = np.array([100, 200])

    def test_at_setitem_dict_value_not_affected(self):
        # GH#61223 - dict is not list-like; .at raises ValueError as before
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        with pytest.raises(ValueError, match="Incompatible indexer"):
            df.at["a", "A"] = {"key": "val"}

    def test_at_setitem_series_value_single_cell(self):
        # GH#61223 - Series raises ValueError for non-object dtype
        df = DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
        val = Series([10, 20], index=["x", "y"])
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at["a", "A"] = val

    def test_at_setitem_list_like_float_column(self):
        # GH#61223 - float64 column raises ValueError like int64
        df = DataFrame({"A": [1.0, 2.0, 3.0]}, index=["a", "b", "c"])
        with pytest.raises(ValueError, match="Cannot store a list-like value"):
            df.at["b", "A"] = [10.5, 20.5]

    def test_at_setitem_list_like_multiindex_row_object_dtype(self):
        # GH#61223 - full tuple key on unique MultiIndex stores list-like as single cell
        idx = MultiIndex.from_tuples([("r1", 0), ("r2", 1)])
        df = DataFrame({"A": Series([[1, 2], [3, 4]], index=idx)})
        assert df["A"].dtype == np.dtype("object")
        df.at[("r1", 0), "A"] = [10, 20, 30]
        assert df.at[("r1", 0), "A"] == [10, 20, 30]
        assert df["A"].dtype == np.dtype("object")

    def test_iat_setitem_list_like_object_dtype(self):
        # GH#61223 - .iat on object-dtype column stores list-like as single cell
        df = DataFrame({"A": Series([[1, 2], [3, 4], [5, 6]])})
        assert df["A"].dtype == np.dtype("object")
        df.iat[0, 0] = [10, 20, 30]
        assert df.iat[0, 0] == [10, 20, 30]
        assert df["A"].dtype == np.dtype("object")

    def test_at_setitem_tuple_value_object_dtype(self):
        # GH#61223 - tuple value in object-dtype column stores as single cell
        df = DataFrame({"A": Series([(1, 2), (3, 4)])})
        assert df["A"].dtype == np.dtype("object")
        df.at[0, "A"] = (10, 20, 30)
        assert df.at[0, "A"] == (10, 20, 30)

    def test_at_setitem_series_value_object_dtype_column(self):
        # GH#61223 - Series value in object-dtype column stores as single cell
        df = DataFrame({"A": Series([[1, 2], [3, 4]])})
        assert df["A"].dtype == np.dtype("object")
        val = Series([10, 20], index=["x", "y"])
        df.at[0, "A"] = val
        result = df.at[0, "A"]
        tm.assert_series_equal(result, val)
