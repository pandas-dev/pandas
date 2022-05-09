import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
import pandas.core.common as com


def test_copy(using_copy_on_write):
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy()

    # the deep copy doesn't share memory
    assert not np.may_share_memory(df_copy["a"].values, df["a"].values)
    if using_copy_on_write:
        assert df_copy._mgr.refs is None

    # mutating copy doesn't mutate original
    df_copy.iloc[0, 0] = 0
    assert df.iloc[0, 0] == 1


def test_copy_shallow(using_copy_on_write):
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy(deep=False)

    # the shallow copy still shares memory
    assert np.may_share_memory(df_copy["a"].values, df["a"].values)
    if using_copy_on_write:
        assert df_copy._mgr.refs is not None

    if using_copy_on_write:
        # mutating shallow copy doesn't mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 1
        # mutating triggered a copy-on-write -> no longer shares memory
        assert not np.may_share_memory(df_copy["a"].values, df["a"].values)
        # but still shares memory for the other columns/blocks
        assert np.may_share_memory(df_copy["c"].values, df["c"].values)
    else:
        # mutating shallow copy does mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 0
        # and still shares memory
        assert np.may_share_memory(df_copy["a"].values, df["a"].values)


# -----------------------------------------------------------------------------
# DataFrame methods returning new DataFrame using shallow copy


def test_reset_index(using_copy_on_write):
    # Case: resetting the index (i.e. adding a new column) + mutating the
    # resulting dataframe
    df = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]}, index=[10, 11, 12]
    )
    df_orig = df.copy()
    df2 = df.reset_index()
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.may_share_memory(df2["b"].values, df["b"].values)
        assert np.may_share_memory(df2["c"].values, df["c"].values)
    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 2] = 0
    assert not np.may_share_memory(df2["b"].values, df["b"].values)
    if using_copy_on_write:
        assert np.may_share_memory(df2["c"].values, df["c"].values)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns(using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the result
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.rename(columns=str.upper)

    if using_copy_on_write:
        assert np.may_share_memory(df2["A"].values, df["a"].values)
    df2.iloc[0, 0] = 0
    assert not np.may_share_memory(df2["A"].values, df["a"].values)
    if using_copy_on_write:
        assert np.may_share_memory(df2["C"].values, df["c"].values)
    expected = pd.DataFrame({"A": [0, 2, 3], "B": [4, 5, 6], "C": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns_modify_parent(using_array_manager, using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the original (parent) dataframe
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df2 = df.rename(columns=str.upper)
    df2_orig = df2.copy()

    if using_copy_on_write:
        assert np.may_share_memory(df2["A"].values, df["a"].values)
    else:
        assert not np.may_share_memory(df2["A"].values, df["a"].values)
    df.iloc[0, 0] = 0
    assert not np.may_share_memory(df2["A"].values, df["a"].values)
    if using_array_manager:
        assert np.may_share_memory(df2["C"].values, df["c"].values)
    expected = pd.DataFrame({"a": [0, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df, expected)
    tm.assert_frame_equal(df2, df2_orig)


def test_reindex_columns(using_copy_on_write):
    # Case: reindexing the column returns a new dataframe
    # + afterwards modifying the result
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.reindex(columns=["a", "c"])

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.may_share_memory(df2["a"].values, df["a"].values)
    else:
        assert not np.may_share_memory(df2["a"].values, df["a"].values)
    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    assert not np.may_share_memory(df2["a"].values, df["a"].values)
    if using_copy_on_write:
        assert np.may_share_memory(df2["c"].values, df["c"].values)
    tm.assert_frame_equal(df, df_orig)


# # -----------------------------------------------------------------------------
# # Indexing operations taking subset + modifying the subset/parent


def test_subset_column_selection(using_copy_on_write):
    # Case: taking a subset of the columns of a DataFrame
    # + afterwards modifying the subset
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    subset = df[["a", "c"]]

    if using_copy_on_write:
        # the subset shares memory ...
        assert np.may_share_memory(subset["a"].values, df["a"].values)
        # ... but uses CoW when being modified
        subset.iloc[0, 0] = 0
    else:
        assert not np.may_share_memory(subset["a"].values, df["a"].values)
        # INFO this no longer raise warning since pandas 1.4
        # with pd.option_context("chained_assignment", "warn"):
        #     with tm.assert_produces_warning(com.SettingWithCopyWarning):
        subset.iloc[0, 0] = 0

    assert not np.may_share_memory(subset["a"].values, df["a"].values)

    expected = pd.DataFrame({"a": [0, 2, 3], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


def test_subset_column_selection_modify_parent(using_copy_on_write):
    # Case: taking a subset of the columns of a DataFrame
    # + afterwards modifying the parent
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})

    subset = df[["a", "c"]]
    if using_copy_on_write:
        # the subset shares memory ...
        assert np.may_share_memory(subset["a"].values, df["a"].values)
        # ... but parent uses CoW parent when it is modified
    df.iloc[0, 0] = 0

    assert not np.may_share_memory(subset["a"].values, df["a"].values)
    if using_copy_on_write:
        # different column/block still shares memory
        assert np.may_share_memory(subset["c"].values, df["c"].values)

    expected = pd.DataFrame({"a": [1, 2, 3], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(subset, expected)


def test_subset_row_slice(using_copy_on_write):
    # Case: taking a subset of the rows of a DataFrame using a slice
    # + afterwards modifying the subset
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    subset = df[1:3]
    subset._mgr._verify_integrity()

    assert np.may_share_memory(subset["a"].values, df["a"].values)

    if using_copy_on_write:
        subset.iloc[0, 0] = 0
        assert not np.may_share_memory(subset["a"].values, df["a"].values)

    else:
        # INFO this no longer raise warning since pandas 1.4
        # with pd.option_context("chained_assignment", "warn"):
        #     with tm.assert_produces_warning(com.SettingWithCopyWarning):
        subset.iloc[0, 0] = 0

    subset._mgr._verify_integrity()

    expected = pd.DataFrame(
        {"a": [0, 3], "b": [5, 6], "c": [0.2, 0.3]}, index=range(1, 3)
    )
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        # original parent dataframe is not modified (CoW)
        tm.assert_frame_equal(df, df_orig)
    else:
        # original parent dataframe is actually updated
        df_orig.iloc[1, 0] = 0
        tm.assert_frame_equal(df, df_orig)


def test_subset_column_slice(using_copy_on_write):
    # Case: taking a subset of the columns of a DataFrame using a slice
    # + afterwards modifying the subset
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    subset = df.iloc[:, 1:]
    subset._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.may_share_memory(subset["b"].values, df["b"].values)

        subset.iloc[0, 0] = 0
        assert not np.may_share_memory(subset["b"].values, df["b"].values)

    else:
        subset.iloc[0, 0] = 0

    expected = pd.DataFrame({"b": [0, 5, 6], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (also not for BlockManager case)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "indexer",
    [slice(0, 2), np.array([True, True, False]), np.array([0, 1])],
    ids=["slice", "mask", "array"],
)
def test_subset_set_with_row_indexer(indexer_si, indexer, using_copy_on_write):
    # Case: setting values with a row indexer on a viewing subset
    # subset[indexer] = value and subset.iloc[indexer] = value
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [4, 5, 6, 7], "c": [0.1, 0.2, 0.3, 0.4]})
    df_orig = df.copy()
    subset = df[1:4]

    if (
        indexer_si is tm.setitem
        and isinstance(indexer, np.ndarray)
        and indexer.dtype == "int"
    ):
        pytest.skip("setitem with labels selects on columns")

    if using_copy_on_write:
        indexer_si(subset)[indexer] = 0
    else:
        # INFO iloc no longer raises warning since pandas 1.4
        warn = com.SettingWithCopyWarning if indexer_si is tm.setitem else None
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(warn):
                indexer_si(subset)[indexer] = 0

    expected = pd.DataFrame(
        {"a": [0, 0, 4], "b": [0, 0, 7], "c": [0.0, 0.0, 0.4]}, index=range(1, 4)
    )
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        # original parent dataframe is not modified (CoW)
        tm.assert_frame_equal(df, df_orig)
    else:
        # original parent dataframe is actually updated
        df_orig[1:3] = 0
        tm.assert_frame_equal(df, df_orig)


def test_subset_set_with_mask(using_copy_on_write):
    # Case: setting values with a mask on a viewing subset: subset[mask] = value
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [4, 5, 6, 7], "c": [0.1, 0.2, 0.3, 0.4]})
    df_orig = df.copy()
    subset = df[1:4]

    mask = subset > 3

    if using_copy_on_write:
        subset[mask] = 0
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset[mask] = 0

    expected = pd.DataFrame(
        {"a": [2, 3, 0], "b": [0, 0, 0], "c": [0.20, 0.3, 0.4]}, index=range(1, 4)
    )
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        # original parent dataframe is not modified (CoW)
        tm.assert_frame_equal(df, df_orig)
    else:
        # original parent dataframe is actually updated
        df_orig.loc[3, "a"] = 0
        df_orig.loc[1:3, "b"] = 0
        tm.assert_frame_equal(df, df_orig)


def test_subset_set_column(using_copy_on_write):
    # Case: setting a single column on a viewing subset -> subset[col] = value
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    subset = df[1:3]

    if using_copy_on_write:
        subset["a"] = np.array([10, 11], dtype="int64")
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset["a"] = np.array([10, 11], dtype="int64")

    subset._mgr._verify_integrity()
    expected = pd.DataFrame(
        {"a": [10, 11], "b": [5, 6], "c": [0.2, 0.3]}, index=range(1, 3)
    )
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_set_column_with_loc(using_copy_on_write, dtype):
    # Case: setting a single column with loc on a viewing subset
    # -> subset.loc[:, col] = value
    df = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()
    subset = df[1:3]

    if using_copy_on_write:
        subset.loc[:, "a"] = np.array([10, 11], dtype="int64")
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset.loc[:, "a"] = np.array([10, 11], dtype="int64")

    subset._mgr._verify_integrity()
    expected = pd.DataFrame(
        {"a": [10, 11], "b": [5, 6], "c": np.array([8, 9], dtype=dtype)},
        index=range(1, 3),
    )
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        # original parent dataframe is not modified (CoW)
        tm.assert_frame_equal(df, df_orig)
    else:
        # original parent dataframe is actually updated
        df_orig.loc[1:3, "a"] = np.array([10, 11], dtype="int64")
        tm.assert_frame_equal(df, df_orig)


def test_subset_set_column_with_loc2(using_copy_on_write):
    # Case: setting a single column with loc on a viewing subset
    # -> subset.loc[:, col] = value
    # separate test for case of DataFrame of a single column -> takes a separate
    # code path
    df = pd.DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    subset = df[1:3]

    if using_copy_on_write:
        subset.loc[:, "a"] = 0
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset.loc[:, "a"] = 0

    subset._mgr._verify_integrity()
    expected = pd.DataFrame({"a": [0, 0]}, index=range(1, 3))
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        # original parent dataframe is not modified (CoW)
        tm.assert_frame_equal(df, df_orig)
    else:
        # original parent dataframe is actually updated
        df_orig.loc[1:3, "a"] = 0
        tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_set_columns(using_copy_on_write, dtype):
    # Case: setting multiple columns on a viewing subset
    # -> subset[[col1, col2]] = value
    df = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()
    subset = df[1:3]

    if using_copy_on_write:
        subset[["a", "c"]] = 0
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset[["a", "c"]] = 0

    subset._mgr._verify_integrity()
    if using_copy_on_write:
        # first and third column should certainly have no references anymore
        assert all(subset._mgr._has_no_reference(i) for i in [0, 2])
    expected = pd.DataFrame({"a": [0, 0], "b": [5, 6], "c": [0, 0]}, index=range(1, 3))
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "indexer",
    [slice("a", "b"), np.array([True, True, False]), ["a", "b"]],
    ids=["slice", "mask", "array"],
)
def test_subset_set_with_column_indexer(indexer, using_copy_on_write):
    # Case: setting multiple columns with a column indexer on a viewing subset
    # -> subset.loc[:, [col1, col2]] = value
    df = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "c": [4, 5, 6]})
    df_orig = df.copy()
    subset = df[1:3]

    if using_copy_on_write:
        subset.loc[:, indexer] = 0
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                subset.loc[:, indexer] = 0

    subset._mgr._verify_integrity()
    expected = pd.DataFrame(
        {"a": [0, 0], "b": [0.0, 0.0], "c": [5, 6]}, index=range(1, 3)
    )
    # TODO full row slice .loc[:, idx] update inplace instead of overwrite?
    expected["b"] = expected["b"].astype("int64")
    tm.assert_frame_equal(subset, expected)
    if using_copy_on_write:
        tm.assert_frame_equal(df, df_orig)
    else:
        # In the mixed case with BlockManager, only one of the two columns is
        # mutated in the parent frame ..
        df_orig.loc[1:2, ["a"]] = 0
        tm.assert_frame_equal(df, df_orig)


# TODO add more tests modifying the parent


# -----------------------------------------------------------------------------
# Series -- Indexing operations taking subset + modifying the subset/parent


def test_series_getitem_slice(using_copy_on_write):
    # Case: taking a slice of a Series + afterwards modifying the subset
    s = pd.Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()

    subset = s[:]
    assert np.may_share_memory(subset.values, s.values)

    subset.iloc[0] = 0

    if using_copy_on_write:
        assert not np.may_share_memory(subset.values, s.values)

    expected = pd.Series([0, 2, 3], index=["a", "b", "c"])
    tm.assert_series_equal(subset, expected)

    if using_copy_on_write:
        # original parent series is not modified (CoW)
        tm.assert_series_equal(s, s_orig)
    else:
        # original parent series is actually updated
        assert s.iloc[0] == 0


@pytest.mark.parametrize(
    "indexer",
    [slice(0, 2), np.array([True, True, False]), np.array([0, 1])],
    ids=["slice", "mask", "array"],
)
def test_series_subset_set_with_indexer(indexer_si, indexer, using_copy_on_write):
    # Case: setting values in a viewing Series with an indexer
    s = pd.Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()
    subset = s[:]

    indexer_si(subset)[indexer] = 0
    expected = pd.Series([0, 0, 3], index=["a", "b", "c"])
    tm.assert_series_equal(subset, expected)

    if using_copy_on_write:
        tm.assert_series_equal(s, s_orig)
    else:
        tm.assert_series_equal(s, expected)


# -----------------------------------------------------------------------------
# del operator


def test_del_frame(using_copy_on_write):
    # Case: deleting a column with `del` on a viewing child dataframe should
    # not modify parent + update the references
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df[:]

    assert np.may_share_memory(df["a"].values, df2["a"].values)

    del df2["b"]

    assert np.may_share_memory(df["a"].values, df2["a"].values)
    tm.assert_frame_equal(df, df_orig)
    tm.assert_frame_equal(df2, df_orig[["a", "c"]])
    df2._mgr._verify_integrity()

    # TODO in theory modifying column "b" of the parent wouldn't need a CoW
    # but the weakref is still alive and so we still perform CoW

    df2.loc[0, "a"] = 100
    if using_copy_on_write:
        # modifying child after deleting a column still doesn't update parent
        tm.assert_frame_equal(df, df_orig)
    else:
        assert df.loc[0, "a"] == 100


def test_del_series():
    s = pd.Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()
    s2 = s[:]

    assert np.may_share_memory(s.values, s2.values)

    del s2["a"]

    assert not np.may_share_memory(s.values, s2.values)
    tm.assert_series_equal(s, s_orig)
    tm.assert_series_equal(s2, s_orig[["b", "c"]])

    # modifying s2 doesn't need copy on write (due to `del`, s2 is backed by new array)
    values = s2.values
    s2.loc["b"] = 100
    assert values[0] == 100


# -----------------------------------------------------------------------------
# Accessing column as Series


def test_column_as_series(using_copy_on_write):
    # Case: selecting a single column now also uses Copy-on-Write
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    s = df["a"]

    assert np.may_share_memory(s.values, df["a"].values)

    if using_copy_on_write:
        s[0] = 0
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                s[0] = 0

    expected = pd.Series([0, 2, 3], name="a")
    tm.assert_series_equal(s, expected)
    if using_copy_on_write:
        # assert not np.may_share_memory(s.values, df["a"].values)
        tm.assert_frame_equal(df, df_orig)
        # ensure cached series on getitem is not the changed series
        tm.assert_series_equal(df["a"], df_orig["a"])
    else:
        df_orig.iloc[0, 0] = 0
        tm.assert_frame_equal(df, df_orig)


def test_column_as_series_set_with_upcast(using_copy_on_write):
    # Case: selecting a single column now also uses Copy-on-Write -> when
    # setting a value causes an upcast, we don't need to update the parent
    # DataFrame through the cache mechanism
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    s = df["a"]
    if using_copy_on_write:
        s[0] = "foo"
    else:
        with pd.option_context("chained_assignment", "warn"):
            with tm.assert_produces_warning(com.SettingWithCopyWarning):
                s[0] = "foo"

    expected = pd.Series(["foo", 2, 3], dtype=object, name="a")
    tm.assert_series_equal(s, expected)
    if using_copy_on_write:
        tm.assert_frame_equal(df, df_orig)
        # ensure cached series on getitem is not the changed series
        tm.assert_series_equal(df["a"], df_orig["a"])
    else:
        df_orig["a"] = expected
        tm.assert_frame_equal(df, df_orig)


# TODO add tests for other indexing methods on the Series


def test_dataframe_add_column_from_series():
    # Case: adding a new column to a DataFrame from an existing column/series
    # -> always already takes a copy on assignment
    # (no change in behaviour here)
    # TODO can we achieve the same behaviour with Copy-on-Write?
    df = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})

    s = pd.Series([10, 11, 12])
    df["new"] = s
    assert not np.may_share_memory(df["new"].values, s.values)

    # editing series -> doesn't modify column in frame
    s[0] = 0
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "new": [10, 11, 12]})
    tm.assert_frame_equal(df, expected)

    # editing column in frame -> doesn't modify series
    df.loc[2, "new"] = 100
    expected = pd.Series([0, 11, 12])
    tm.assert_series_equal(s, expected)


# TODO add tests for constructors
