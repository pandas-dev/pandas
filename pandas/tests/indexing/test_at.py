from datetime import (
    UTC,
    datetime,
)
from decimal import Decimal

import numpy as np
import pytest

from pandas.errors import (
    InvalidIndexError,
    Pandas4Warning,
)

import pandas as pd
import pandas._testing as tm


def test_at_dateoffset_columns():
    # GH#20948 - .at with DateOffset columns
    offsets = pd.Series(data=[-15, -10, -5, 0, 5, 10, 15], dtype=float).map(
        pd.DateOffset
    )
    df = pd.DataFrame(index=[0, 1], columns=pd.Index(offsets))

    # read access
    result = df.at[0, offsets[0]]
    assert result is np.nan

    # write access
    df.at[0, offsets[0]] = 1
    assert df.at[0, offsets[0]] == 1


def test_at_multiindex_partial_date_string():
    # GH#43395 - .at with partial date string on MultiIndex with DatetimeIndex level
    timestamps = pd.DatetimeIndex(
        ["2021-08-01", "2021-08-01 12:00", "2021-08-02", "2021-08-02 12:00"]
    )
    index = pd.MultiIndex.from_product(
        [["A", "B"], timestamps], names=["ticker", "timestamp"]
    )
    df = pd.DataFrame({"col": range(len(index))}, index=index)

    # DataFrame.at with partial date string should not raise TypeError
    result = df.at[("A", "2021-08-02"), "col"]
    expected = np.array([2, 3], dtype=np.int64)
    tm.assert_numpy_array_equal(result, expected)

    # Series.at with partial date string
    ser = df["col"]
    result = ser.at[("A", "2021-08-02")]
    expected = pd.Series(
        [2, 3],
        index=pd.MultiIndex.from_arrays(
            [
                ["A", "A"],
                pd.DatetimeIndex(
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
    df = pd.DataFrame({"A": [1, 2, 3]})
    with pytest.raises(TypeError, match="Invalid value"):
        df.at[0, "A"] = Decimal("1")


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    result = pd.DataFrame({"foo": [datetime(2000, 1, 1)]})
    with pytest.raises(TypeError, match="Invalid value"):
        result.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=UTC)


def test_selection_methods_of_assigned_col():
    # GH 29282
    df = pd.DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
    df2 = pd.DataFrame(data={"c": [7, 8, 9]}, index=[2, 1, 0])
    df["c"] = df2["c"]
    df.at[1, "c"] = 11
    result = df
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [9, 11, 7]})
    tm.assert_frame_equal(result, expected)
    result = df.at[1, "c"]
    assert result == 11

    result = df["c"]
    expected = pd.Series([9, 11, 7], name="c")
    tm.assert_series_equal(result, expected)

    result = df[["c"]]
    expected = pd.DataFrame({"c": [9, 11, 7]})
    tm.assert_frame_equal(result, expected)


class TestAtSetItem:
    def test_at_setitem_mixed_index_assignment(self):
        # GH#19860
        ser = pd.Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        ser.at["a"] = 11
        assert ser.iat[0] == 11
        ser.at[1] = 22
        assert ser.iat[3] == 22

    def test_at_setitem_categorical_missing(self):
        df = pd.DataFrame(
            index=range(3), columns=range(3), dtype=pd.CategoricalDtype(["foo", "bar"])
        )
        df.at[1, 1] = "foo"

        expected = pd.DataFrame(
            [
                [np.nan, np.nan, np.nan],
                [np.nan, "foo", np.nan],
                [np.nan, np.nan, np.nan],
            ],
            dtype=pd.CategoricalDtype(["foo", "bar"]),
        )

        tm.assert_frame_equal(df, expected)

    def test_at_setitem_multiindex(self):
        df = pd.DataFrame(
            np.zeros((3, 2), dtype="int64"),
            columns=pd.MultiIndex.from_tuples([("a", 0), ("a", 1)]),
        )
        df.at[0, "a"] = 10
        expected = pd.DataFrame(
            [[10, 10], [0, 0], [0, 0]],
            columns=pd.MultiIndex.from_tuples([("a", 0), ("a", 1)]),
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("row", (pd.Timestamp("2019-01-01"), "2019-01-01"))
    def test_at_datetime_index(self, row):
        # Set float64 dtype to avoid upcast when setting .5
        df = pd.DataFrame(
            data=[[1] * 2], index=pd.DatetimeIndex(data=["2019-01-01", "2019-01-02"])
        ).astype({0: "float64"})
        expected = pd.DataFrame(
            data=[[0.5, 1], [1.0, 1]],
            index=pd.DatetimeIndex(data=["2019-01-01", "2019-01-02"]),
        )

        df.at[row, 0] = 0.5
        tm.assert_frame_equal(df, expected)


class TestAtSetItemWithExpansion:
    def test_at_setitem_expansion_series_dt64tz_value(self, tz_naive_fixture):
        # GH#25506
        ts = (
            pd.Timestamp("2017-08-05 00:00:00+0100", tz=tz_naive_fixture)
            if tz_naive_fixture is not None
            else pd.Timestamp("2017-08-05 00:00:00+0100")
        )
        result = pd.Series(ts)
        with tm.assert_produces_warning(Pandas4Warning, match="does not exist"):
            result.at[1] = ts
        expected = pd.Series([ts, ts])
        tm.assert_series_equal(result, expected)

    def test_at_setitem_expansion_deprecated_dataframe(self):
        # GH#48323
        df = pd.DataFrame({"a": [1, 2]})
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[5, "a"] = 6
        expected = pd.DataFrame({"a": [1, 2, 6]}, index=[0, 1, 5])
        tm.assert_frame_equal(df, expected)

    def test_at_setitem_expansion_deprecated_series(self):
        # GH#48323
        ser = pd.Series([1, 2], index=[0, 1])
        msg = "Setting a value on a Series via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            ser.at[5] = 6
        expected = pd.Series([1, 2, 6], index=[0, 1, 5])
        tm.assert_series_equal(ser, expected)

    def test_at_setitem_expansion_deprecated_new_column(self):
        # GH#48323 - new column (existing row) also expands
        df = pd.DataFrame({"a": [1, 2]})
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[0, "b"] = 99
        assert "b" in df.columns

    def test_at_setitem_expansion_deprecated_new_column_multiindex_rows(self):
        # GH#48323 - new column on a frame with MultiIndex rows
        mi = pd.MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 1)])
        df = pd.DataFrame({"x": [1, 2, 3]}, index=mi)
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[("a", 1), "y"] = 99
        assert "y" in df.columns

    def test_at_setitem_expansion_deprecated_new_column_multiindex_cols(self):
        # GH#48323 - new tuple column on a frame with MultiIndex columns
        df = pd.DataFrame(
            [[1, 2], [3, 4]],
            columns=pd.MultiIndex.from_tuples([("a", 1), ("a", 2)]),
        )
        msg = "Setting a value on a DataFrame via .at with a key"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.at[0, ("b", 1)] = 99
        assert ("b", 1) in df.columns

    def test_at_setitem_no_warning_existing_key(self):
        # GH#48323 - no warning for existing keys
        df = pd.DataFrame({"a": [1, 2]})
        with tm.assert_produces_warning(None):
            df.at[0, "a"] = 99
        assert df.at[0, "a"] == 99

        ser = pd.Series([1, 2])
        with tm.assert_produces_warning(None):
            ser.at[0] = 99
        assert ser.at[0] == 99


class TestAtWithDuplicates:
    def test_at_with_duplicate_axes_requires_scalar_lookup(self):
        # GH#33041 check that falling back to loc doesn't allow non-scalar
        #  args to slip in

        arr = np.random.default_rng(2).standard_normal(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=["A", "A"])

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


class TestAtErrors:
    # TODO: De-duplicate/parametrize
    #  test_at_series_raises_key_error2, test_at_frame_raises_key_error2

    def test_at_series_raises_key_error(self, indexer_al):
        # GH#31724 .at should match .loc

        ser = pd.Series([1, 2, 3], index=[3, 2, 1])
        result = indexer_al(ser)[1]
        assert result == 3

        with pytest.raises(KeyError, match="a"):
            indexer_al(ser)["a"]

    def test_at_frame_raises_key_error(self, indexer_al):
        # GH#31724 .at should match .loc

        df = pd.DataFrame({0: [1, 2, 3]}, index=[3, 2, 1])

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
        ser = pd.Series([1, 2, 3], index=list("abc"))
        result = indexer_al(ser)["a"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            indexer_al(ser)[0]

    def test_at_frame_raises_key_error2(self, indexer_al):
        # GH#31724 .at should match .loc
        df = pd.DataFrame({"A": [1, 2, 3]}, index=list("abc"))
        result = indexer_al(df)["a", "A"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            indexer_al(df)["a", 0]

    def test_at_frame_multiple_columns(self):
        # GH#48296 - at shouldn't modify multiple columns
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        new_row = [6, 7]
        with pytest.raises(
            InvalidIndexError,
            match=".at-based indexing can only have scalar indexers; use .loc instead",
        ):
            df.at[5] = new_row

    def test_at_getitem_mixed_index_no_fallback(self):
        # GH#19860
        ser = pd.Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        with pytest.raises(KeyError, match="^0$"):
            ser.at[0]
        with pytest.raises(KeyError, match="^4$"):
            ser.at[4]

    def test_at_categorical_integers(self):
        # CategoricalIndex with integer categories that don't happen to match
        #  the Categorical's codes
        ci = pd.CategoricalIndex([3, 4])

        arr = np.arange(4).reshape(2, 2)
        frame = pd.DataFrame(arr, index=ci)

        for df in [frame, frame.T]:
            for key in [0, 1]:
                with pytest.raises(KeyError, match=str(key)):
                    df.at[key, key]

    def test_at_applied_for_rows(self):
        # GH#48729 .at should raise InvalidIndexError when assigning rows
        df = pd.DataFrame(index=["a"], columns=["col1", "col2"])
        new_row = [123, 15]
        with pytest.raises(
            InvalidIndexError,
            match=".at-based indexing can only have scalar indexers; use .loc instead",
        ):
            df.at["a"] = new_row

    @pytest.mark.parametrize("key", [lambda df: df["b"] == 6, [0, 1], np.array([0, 1])])
    def test_at_setitem_nonscalar_indexer_raises(self, key):
        # GH#51866 - .at setitem with a mask/list/array indexer should raise a
        #  clear error directing the user to .loc, not the misleading
        #  "scalar value" message
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        indexer = key(df) if callable(key) else key
        with pytest.raises(
            InvalidIndexError,
            match=".at-based indexing can only have scalar indexers; use .loc instead",
        ):
            df.at[indexer, "b"] = 9.0

    @pytest.mark.parametrize("key", [lambda ser: ser > 2, [0, 1], np.array([0, 1])])
    def test_at_setitem_series_nonscalar_indexer_raises(self, key):
        # GH#51866 - Series .at setitem with a non-scalar indexer should raise a
        #  clear error rather than dumping the raw indexer
        ser = pd.Series([1, 2, 3, 4])
        indexer = key(ser) if callable(key) else key
        with pytest.raises(
            InvalidIndexError,
            match=".at-based indexing can only have scalar indexers; use .loc instead",
        ):
            ser.at[indexer] = 9.0
