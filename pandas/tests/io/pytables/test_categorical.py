import numpy as np
import pytest

from pandas import (
    Categorical,
    DataFrame,
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
