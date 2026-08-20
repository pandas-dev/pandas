import numpy as np
import pytest

from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    HDFStore,
    Index,
    Series,
    _testing as tm,
    concat,
    read_hdf,
)

pytestmark = [pytest.mark.single_cpu]


def test_categorical(temp_hdfstore):
    # Basic
    s = Series(
        Categorical(
            ["a", "b", "b", "a", "a", "c"],
            categories=["a", "b", "c", "d"],
            ordered=False,
        )
    )
    temp_hdfstore.append("s", s, format="table")
    result = temp_hdfstore.select("s")
    tm.assert_series_equal(s, result)

    s = Series(
        Categorical(
            ["a", "b", "b", "a", "a", "c"],
            categories=["a", "b", "c", "d"],
            ordered=True,
        )
    )
    temp_hdfstore.append("s_ordered", s, format="table")
    result = temp_hdfstore.select("s_ordered")
    tm.assert_series_equal(s, result)

    df = DataFrame({"s": s, "vals": [1, 2, 3, 4, 5, 6]})
    temp_hdfstore.append("df", df, format="table")
    result = temp_hdfstore.select("df")
    tm.assert_frame_equal(result, df)

    # Dtypes
    s = Series([1, 1, 2, 2, 3, 4, 5]).astype("category")
    temp_hdfstore.append("si", s)
    result = temp_hdfstore.select("si")
    tm.assert_series_equal(result, s)

    s = Series([1, 1, np.nan, 2, 3, 4, 5]).astype("category")
    temp_hdfstore.append("si2", s)
    result = temp_hdfstore.select("si2")
    tm.assert_series_equal(result, s)

    # Multiple
    df2 = df.copy()
    df2["s2"] = Series(list("abcdefg")).astype("category")
    temp_hdfstore.append("df2", df2)
    result = temp_hdfstore.select("df2")
    tm.assert_frame_equal(result, df2)

    # Make sure the metadata is OK
    info = temp_hdfstore.info()
    assert "/df2   " in info
    # df2._mgr.blocks[0] and df2._mgr.blocks[2] are Categorical
    assert "/df2/meta/values_block_0/meta" in info
    assert "/df2/meta/values_block_2/meta" in info

    # unordered
    s = Series(
        Categorical(
            ["a", "b", "b", "a", "a", "c"],
            categories=["a", "b", "c", "d"],
            ordered=False,
        )
    )
    temp_hdfstore.append("s2", s, format="table")
    result = temp_hdfstore.select("s2")
    tm.assert_series_equal(result, s)

    # Query
    temp_hdfstore.append("df3", df, data_columns=["s"])
    expected = df[df.s.isin(["b", "c"])]
    result = temp_hdfstore.select("df3", where=['s in ["b","c"]'])
    tm.assert_frame_equal(result, expected)

    expected = df[df.s.isin(["b", "c"])]
    result = temp_hdfstore.select("df3", where=['s = ["b","c"]'])
    tm.assert_frame_equal(result, expected)

    expected = df[df.s.isin(["d"])]
    result = temp_hdfstore.select("df3", where=['s in ["d"]'])
    tm.assert_frame_equal(result, expected)

    expected = df[df.s.isin(["f"])]
    result = temp_hdfstore.select("df3", where=['s in ["f"]'])
    tm.assert_frame_equal(result, expected)

    # Appending with same categories is ok
    temp_hdfstore.append("df3", df)

    df = concat([df, df])
    expected = df[df.s.isin(["b", "c"])]
    result = temp_hdfstore.select("df3", where=['s in ["b","c"]'])
    tm.assert_frame_equal(result, expected)

    # Appending must have the same categories
    df3 = df.copy()
    df3["s"] = df3["s"].cat.remove_unused_categories()

    msg = "cannot append a categorical with different categories to the existing"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df3", df3)

    # Remove, and make sure meta data is removed (its a recursive
    # removal so should be).
    result = temp_hdfstore.select("df3/meta/s/meta")
    assert result is not None
    temp_hdfstore.remove("df3")

    with pytest.raises(KeyError, match="'No object named df3/meta/s/meta in the file'"):
        temp_hdfstore.select("df3/meta/s/meta")


def test_categorical_conversion(temp_h5_path):
    # GH13322
    # Check that read_hdf with categorical columns doesn't return rows if
    # where criteria isn't met.
    obsids = ["ESP_012345_6789", "ESP_987654_3210"]
    imgids = ["APF00006np", "APF0001imm"]
    data = [4.3, 9.8]

    # Test without categories
    df = DataFrame({"obsids": obsids, "imgids": imgids, "data": data})

    # We are expecting an empty DataFrame matching types of df
    expected = df.iloc[[], :]
    df.to_hdf(temp_h5_path, key="df", format="table", data_columns=True)
    result = read_hdf(temp_h5_path, "df", where="obsids=B")
    tm.assert_frame_equal(result, expected)

    # Test with categories
    df.obsids = df.obsids.astype("category")
    df.imgids = df.imgids.astype("category")

    # We are expecting an empty DataFrame matching types of df
    expected = df.iloc[[], :]
    df.to_hdf(temp_h5_path, key="df", format="table", data_columns=True)
    result = read_hdf(temp_h5_path, "df", where="obsids=B")
    tm.assert_frame_equal(result, expected)


def test_categorical_nan_only_columns(temp_h5_path):
    # GH18413
    # Check that read_hdf with categorical columns with NaN-only values can
    # be read back.
    df = DataFrame(
        {
            "a": ["a", "b", "c", np.nan],
            "b": [np.nan, np.nan, np.nan, np.nan],
            "c": [1, 2, 3, 4],
            "d": Series([None] * 4, dtype=object),
        }
    )
    df["a"] = df.a.astype("category")
    df["b"] = df.b.astype("category")
    df["d"] = df.b.astype("category")
    expected = df
    df.to_hdf(temp_h5_path, key="df", format="table", data_columns=True)
    result = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(result, expected)


def test_categorical_nan_rep_collision(temp_h5_path):
    # GH#21741 - a category equal to the default nan_rep ("nan") used to
    # raise a broadcasting ValueError on read. It should instead be read
    # back as NaN, with the surviving codes renumbered correctly.
    df = DataFrame(
        {"A": Series(["aaa", "nan", "bbb", "aaa", "zzz", "bbb"]).astype("category")}
    )
    df.to_hdf(temp_h5_path, key="df", format="table")
    result = read_hdf(temp_h5_path, key="df")

    expected = DataFrame(
        {"A": Series(["aaa", np.nan, "bbb", "aaa", "zzz", "bbb"]).astype("category")}
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("where, expected", [["q", []], ["a", ["a"]]])
def test_convert_value(temp_h5_path, where: str, expected):
    # GH39420
    # Check that read_hdf with categorical columns can filter by where condition.
    df = DataFrame({"col": ["a", "b", "s"]})
    df.col = df.col.astype("category")
    max_widths = {"col": 1}
    categorical_values = sorted(df.col.unique())
    expected = DataFrame({"col": expected})
    expected.col = expected.col.astype("category")
    expected.col = expected.col.cat.set_categories(categorical_values)

    df.to_hdf(temp_h5_path, key="df", format="table", min_itemsize=max_widths)
    result = read_hdf(temp_h5_path, where=f'col=="{where}"')
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_roundtrip(fmt, temp_h5_path):
    # GH#33909 (fixed), GH#16118 (table)
    df = DataFrame(
        {"x": [1, 2, 3]}, index=Categorical(["A", "B", "A"], categories=["A", "B"])
    )
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_preserves_metadata(fmt, temp_h5_path):
    # GH#33909, GH#16118 - ordered flag, custom name, and a non-sorted
    # category order all need to round-trip.
    ci = CategoricalIndex(
        ["b", "a", "c"], categories=["c", "b", "a"], ordered=True, name="ix"
    )
    df = DataFrame({"v": [1.0, 2.0, 3.0]}, index=ci)
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(result, df)
    assert result.index.ordered is True
    assert list(result.index.categories) == ["c", "b", "a"]
    assert result.index.name == "ix"


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_with_nan(fmt, temp_h5_path):
    # GH#33909, GH#16118 - missing entries (code -1) round-trip
    df = DataFrame(
        {"x": [1, 2, 3]},
        index=CategoricalIndex(["a", None, "b"], categories=["a", "b"]),
    )
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_category_equal_nan_rep(fmt, temp_h5_path):
    # GH#65576 - a genuine "nan" category equals the default nan_rep string.
    # Reading such a file used to raise "Categorical categories cannot be
    # null" and permanently brick the key. fixed format stores categories
    # verbatim and keeps the category; table format encodes categories with
    # nan_rep, so "nan" decodes to NaN and is dropped, matching how a
    # categorical *column* behaves (GH#21741).
    df = DataFrame({"v": [1, 2, 3, 4]}, index=CategoricalIndex(["a", "b", "nan", "a"]))
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")

    if fmt == "fixed":
        expected = df
    else:
        expected = DataFrame(
            {"v": [1, 2, 3, 4]},
            index=CategoricalIndex(["a", "b", np.nan, "a"], categories=["a", "b"]),
        )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_all_nan(fmt, temp_h5_path):
    # GH#65576 - an all-NaN CategoricalIndex has zero categories, which cannot
    # be written as a metadata array. It must still round-trip as categorical
    # rather than as raw integer codes (table format used to read back
    # Index([-1, -1])).
    df = DataFrame({"v": [1, 2]}, index=CategoricalIndex([np.nan, np.nan]))
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")

    expected = DataFrame(
        {"v": [1, 2]},
        index=CategoricalIndex(
            Categorical.from_codes([-1, -1], categories=Index([], dtype="float64"))
        ),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_categorical_index_numeric_categories(fmt, temp_h5_path):
    # GH#33909, GH#16118 - non-string categories
    df = DataFrame(
        {"x": [1, 2, 3]},
        index=CategoricalIndex([10, 20, 30], categories=[10, 20, 30, 40]),
    )
    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize("fmt", ["fixed", "table"])
@pytest.mark.parametrize(
    "dtype, expected_dtype", [("Int64", "int64"), ("string", "str")]
)
def test_categorical_index_extension_categories(
    fmt, dtype, expected_dtype, temp_h5_path
):
    # GH#33909, GH#16118 - extension-dtype categories. Matches the
    # longstanding behavior for categorical columns: table format casts
    # the categories to the numpy equivalent; fixed raises for numeric EAs.
    values = [10, 20, 30] if dtype == "Int64" else ["a", "b", "c"]
    categories = Index(values, dtype=dtype)
    df = DataFrame(
        {"x": [1, 2, 3]}, index=CategoricalIndex(categories, categories=categories)
    )

    if fmt == "fixed" and dtype == "Int64":
        with pytest.raises(
            NotImplementedError, match="Cannot store an Index with dtype Int64"
        ):
            df.to_hdf(temp_h5_path, key="df", format=fmt)
        return

    df.to_hdf(temp_h5_path, key="df", format=fmt)
    result = read_hdf(temp_h5_path, key="df")
    expected_categories = categories.astype(expected_dtype)
    expected = DataFrame(
        {"x": [1, 2, 3]},
        index=CategoricalIndex(expected_categories, categories=expected_categories),
    )
    tm.assert_frame_equal(result, expected)


def test_categorical_columns_axis_fixed_format(temp_h5_path):
    # GH#33909 - the columns axis is also a CategoricalIndex; fixed format
    df = DataFrame(
        [[1, 2], [3, 4]],
        index=["r1", "r2"],
        columns=CategoricalIndex(["a", "b"]),
    )
    df.to_hdf(temp_h5_path, key="df")
    result = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(result, df)


def test_categorical_index_series_fixed_format(temp_h5_path):
    # GH#33909 - Series with a CategoricalIndex
    ser = Series(
        [10, 20, 30],
        index=CategoricalIndex(["x", "y", "z"], categories=["x", "y", "z"]),
        name="vals",
    )
    ser.to_hdf(temp_h5_path, key="ser")
    result = read_hdf(temp_h5_path, key="ser")
    tm.assert_series_equal(result, ser)


def test_categorical_index_table_append(temp_h5_path):
    # GH#16118 - appending rows with the same CategoricalIndex categories
    ci = CategoricalIndex(["a", "b", "c"], categories=["a", "b", "c"], name="ix")
    df = DataFrame({"v": [1.0, 2.0, 3.0]}, index=ci)
    df.to_hdf(temp_h5_path, key="df", format="table")
    df.to_hdf(temp_h5_path, key="df", format="table", append=True)
    result = read_hdf(temp_h5_path, key="df")
    expected = concat([df, df])
    tm.assert_frame_equal(result, expected)


def test_categorical_select_non_sorted_categories(temp_h5_path):
    # GH#38131 - select() returned wrong rows for non-sorted category order
    df = DataFrame(
        {"col1": Categorical(["a", "b", "c", "a"], categories=["c", "b", "a"])}
    )

    with HDFStore(temp_h5_path) as store:
        store.append("df", df, data_columns=df.columns)
        result = store.select("df", "col1='a'")

    expected = df[df["col1"] == "a"]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "where, mask",
    [
        ('col == "z"', lambda col: col == "z"),
        ('col in ["z", "y"]', lambda col: col.isin(["z", "y"])),
        ('col in ["a", "z"]', lambda col: col.isin(["a", "z"])),
        ('col != "z"', lambda col: col != "z"),
    ],
)
def test_categorical_where_non_category_with_nan(temp_h5_path, where, mask):
    # GH#22977 - querying for a value that is not one of the categories must
    # not match NaN rows (whose code is also -1)
    df = DataFrame({"col": Categorical(["a", "b", np.nan, "a"])})
    df.to_hdf(temp_h5_path, key="df", format="table", data_columns=["col"])

    result = read_hdf(temp_h5_path, "df", where=where)
    expected = df[mask(df["col"])]
    tm.assert_frame_equal(result, expected)
