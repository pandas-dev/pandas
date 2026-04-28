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
