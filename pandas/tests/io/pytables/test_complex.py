import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm

from pandas.io.pytables import read_hdf


def test_complex_fixed(temp_h5_path):
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)).astype(np.complex64),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    df.to_hdf(temp_h5_path, key="df")
    reread = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, reread)

    df = DataFrame(
        np.random.default_rng(2).random((4, 5)).astype(np.complex128),
        index=list("abcd"),
        columns=list("ABCDE"),
    )
    df.to_hdf(temp_h5_path, key="df")
    reread = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, reread)


def test_complex_table(temp_h5_path):
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)).astype(np.complex64),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    df.to_hdf(temp_h5_path, key="df", format="table")
    reread = read_hdf(temp_h5_path, key="df")
    tm.assert_frame_equal(df, reread)

    df = DataFrame(
        np.random.default_rng(2).random((4, 5)).astype(np.complex128),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    df.to_hdf(temp_h5_path, key="df", format="table", mode="w")
    reread = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, reread)


def test_complex_mixed_fixed(temp_h5_path):
    complex64 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64
    )
    complex128 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex128
    )
    df = DataFrame(
        {
            "A": [1, 2, 3, 4],
            "B": ["a", "b", "c", "d"],
            "C": complex64,
            "D": complex128,
            "E": [1.0, 2.0, 3.0, 4.0],
        },
        index=list("abcd"),
    )
    df.to_hdf(temp_h5_path, key="df")
    reread = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, reread)


@pytest.mark.filterwarnings(
    "ignore:`alltrue` is deprecated as of NumPy 1.25.0:DeprecationWarning"
)
def test_complex_mixed_table_store_select(temp_hdfstore):
    complex64 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64
    )
    complex128 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex128
    )
    df = DataFrame(
        {
            "A": [1, 2, 3, 4],
            "B": ["a", "b", "c", "d"],
            "C": complex64,
            "D": complex128,
            "E": [1.0, 2.0, 3.0, 4.0],
        },
        index=list("abcd"),
    )

    temp_hdfstore.append("df", df, data_columns=["A", "B"])
    result = temp_hdfstore.select("df", where="A>2")
    tm.assert_frame_equal(df.loc[df.A > 2], result)


@pytest.mark.filterwarnings(
    "ignore:`alltrue` is deprecated as of NumPy 1.25.0:DeprecationWarning"
)
def test_complex_mixed_table_store_to_path(temp_h5_path):
    complex64 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64
    )
    complex128 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex128
    )
    df = DataFrame(
        {
            "A": [1, 2, 3, 4],
            "B": ["a", "b", "c", "d"],
            "C": complex64,
            "D": complex128,
            "E": [1.0, 2.0, 3.0, 4.0],
        },
        index=list("abcd"),
    )

    df.to_hdf(temp_h5_path, key="df", format="table")
    reread = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(df, reread)


def test_complex_across_dimensions_fixed(temp_h5_path):
    complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
    s = Series(complex128, index=list("abcd"))
    df = DataFrame({"A": s, "B": s})

    objs = [s, df]
    comps = [tm.assert_series_equal, tm.assert_frame_equal]
    for obj, comp in zip(objs, comps):
        obj.to_hdf(temp_h5_path, key="obj", format="fixed")
        reread = read_hdf(temp_h5_path, "obj")
        comp(obj, reread)


def test_complex_across_dimensions(temp_h5_path):
    complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
    s = Series(complex128, index=list("abcd"))
    df = DataFrame({"A": s, "B": s})

    df.to_hdf(temp_h5_path, key="obj", format="table")
    reread = read_hdf(temp_h5_path, "obj")
    tm.assert_frame_equal(df, reread)


def test_complex_indexing_error(temp_hdfstore):
    complex128 = np.array(
        [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex128
    )
    df = DataFrame(
        {"A": [1, 2, 3, 4], "B": ["a", "b", "c", "d"], "C": complex128},
        index=list("abcd"),
    )

    msg = (
        "Columns containing complex values can be stored "
        "but cannot be indexed when using table format. "
        "Either use fixed format, set index=False, "
        "or do not include the columns containing complex "
        "values to data_columns when initializing the table."
    )

    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.append("df", df, data_columns=["C"])


def test_complex_series_error(temp_h5_path):
    complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
    s = Series(complex128, index=list("abcd"))

    msg = (
        "Columns containing complex values can be stored "
        "but cannot be indexed when using table format. "
        "Either use fixed format, set index=False, "
        "or do not include the columns containing complex "
        "values to data_columns when initializing the table."
    )

    with pytest.raises(TypeError, match=msg):
        s.to_hdf(temp_h5_path, key="obj", format="t")

    s.to_hdf(temp_h5_path, key="obj", format="t", index=False)
    reread = read_hdf(temp_h5_path, "obj")
    tm.assert_series_equal(s, reread)


def test_complex_append(temp_hdfstore):
    df = DataFrame(
        {
            "a": np.random.default_rng(2).standard_normal(100).astype(np.complex128),
            "b": np.random.default_rng(2).standard_normal(100),
        }
    )

    temp_hdfstore.append("df", df, data_columns=["b"])
    temp_hdfstore.append("df", df)
    result = temp_hdfstore.select("df")
    tm.assert_frame_equal(pd.concat([df, df], axis=0), result)
