import numpy as np

from pandas import DataFrame
import pandas._testing as tm


def test_copy(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy()

    # the deep copy doesn't share memory
    assert not np.shares_memory(df_copy["a"].values, df["a"].values)
    if using_copy_on_write:
        assert df_copy._mgr.refs is None  # type: ignore[union-attr]

    # mutating copy doesn't mutate original
    df_copy.iloc[0, 0] = 0
    assert df.iloc[0, 0] == 1


def test_copy_shallow(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy(deep=False)

    # the shallow copy still shares memory
    assert np.shares_memory(df_copy["a"].values, df["a"].values)
    if using_copy_on_write:
        assert df_copy._mgr.refs is not None  # type: ignore[union-attr]

    if using_copy_on_write:
        # mutating shallow copy doesn't mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 1
        # mutating triggered a copy-on-write -> no longer shares memory
        assert not np.shares_memory(df_copy["a"].values, df["a"].values)
        # but still shares memory for the other columns/blocks
        assert np.shares_memory(df_copy["c"].values, df["c"].values)
    else:
        # mutating shallow copy does mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 0
        # and still shares memory
        assert np.shares_memory(df_copy["a"].values, df["a"].values)


# -----------------------------------------------------------------------------
# DataFrame methods returning new DataFrame using shallow copy


def test_reset_index(using_copy_on_write):
    # Case: resetting the index (i.e. adding a new column) + mutating the
    # resulting dataframe
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]}, index=[10, 11, 12]
    )
    df_orig = df.copy()
    df2 = df.reset_index()
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.shares_memory(df2["b"].values, df["b"].values)
        assert np.shares_memory(df2["c"].values, df["c"].values)
    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 2] = 0
    assert not np.shares_memory(df2["b"].values, df["b"].values)
    if using_copy_on_write:
        assert np.shares_memory(df2["c"].values, df["c"].values)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns(using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.rename(columns=str.upper)

    if using_copy_on_write:
        assert np.shares_memory(df2["A"].values, df["a"].values)
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(df2["A"].values, df["a"].values)
    if using_copy_on_write:
        assert np.shares_memory(df2["C"].values, df["c"].values)
    expected = DataFrame({"A": [0, 2, 3], "B": [4, 5, 6], "C": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns_modify_parent(using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the original (parent) dataframe
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df2 = df.rename(columns=str.upper)
    df2_orig = df2.copy()

    if using_copy_on_write:
        assert np.shares_memory(df2["A"].values, df["a"].values)
    else:
        assert not np.shares_memory(df2["A"].values, df["a"].values)
    df.iloc[0, 0] = 0
    assert not np.shares_memory(df2["A"].values, df["a"].values)
    if using_copy_on_write:
        assert np.shares_memory(df2["C"].values, df["c"].values)
    expected = DataFrame({"a": [0, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df, expected)
    tm.assert_frame_equal(df2, df2_orig)


def test_reindex_columns(using_copy_on_write):
    # Case: reindexing the column returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.reindex(columns=["a", "c"])

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.shares_memory(df2["a"].values, df["a"].values)
    else:
        assert not np.shares_memory(df2["a"].values, df["a"].values)
    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(df2["a"].values, df["a"].values)
    if using_copy_on_write:
        assert np.shares_memory(df2["c"].values, df["c"].values)
    tm.assert_frame_equal(df, df_orig)
