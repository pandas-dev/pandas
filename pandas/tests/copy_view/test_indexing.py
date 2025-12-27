import numpy as np
import pytest

from pandas.core.dtypes.common import is_float_dtype

import pandas as pd
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.fixture(params=["numpy", "nullable"])
def backend(request):
    if request.param == "numpy":

        def make_dataframe(*args, **kwargs):
            return DataFrame(*args, **kwargs)

        def make_series(*args, **kwargs):
            return Series(*args, **kwargs)

    elif request.param == "nullable":

        def make_dataframe(*args, **kwargs):
            df = DataFrame(*args, **kwargs)
            df_nullable = df.convert_dtypes()
            # convert_dtypes will try to cast float to int if there is no loss in
            # precision -> undo that change
            for col in df.columns:
                if is_float_dtype(df[col].dtype) and not is_float_dtype(
                    df_nullable[col].dtype
                ):
                    df_nullable[col] = df_nullable[col].astype("Float64")
            # copy final result to ensure we start with a fully self-owning DataFrame
            return df_nullable.copy()

        def make_series(*args, **kwargs):
            ser = Series(*args, **kwargs)
            return ser.convert_dtypes().copy()

    return request.param, make_dataframe, make_series


# -----------------------------------------------------------------------------
# Indexing operations taking subset + modifying the subset/parent


def test_subset_column_selection(backend):
    # Case: taking a subset of the columns of a DataFrame
    # + afterwards modifying the subset
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    subset = df[["a", "c"]]

    assert subset.index is not df.index

    # the subset shares memory ...
    assert np.shares_memory(get_array(subset, "a"), get_array(df, "a"))
    # ... but uses CoW when being modified
    subset.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(subset, "a"), get_array(df, "a"))

    expected = DataFrame({"a": [0, 2, 3], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


def test_subset_column_selection_modify_parent(backend):
    # Case: taking a subset of the columns of a DataFrame
    # + afterwards modifying the parent
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})

    subset = df[["a", "c"]]

    # the subset shares memory ...
    assert np.shares_memory(get_array(subset, "a"), get_array(df, "a"))
    # ... but parent uses CoW parent when it is modified
    df.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(subset, "a"), get_array(df, "a"))
    # different column/block still shares memory
    assert np.shares_memory(get_array(subset, "c"), get_array(df, "c"))

    expected = DataFrame({"a": [1, 2, 3], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(subset, expected)


def test_subset_row_slice(backend):
    # Case: taking a subset of the rows of a DataFrame using a slice
    # + afterwards modifying the subset
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    subset = df[1:3]
    subset._mgr._verify_integrity()

    assert subset.columns is not df.columns
    assert np.shares_memory(get_array(subset, "a"), get_array(df, "a"))

    subset.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(subset, "a"), get_array(df, "a"))

    subset._mgr._verify_integrity()

    expected = DataFrame({"a": [0, 3], "b": [5, 6], "c": [0.2, 0.3]}, index=range(1, 3))
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (CoW)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_column_slice(backend, dtype):
    # Case: taking a subset of the columns of a DataFrame using a slice
    # + afterwards modifying the subset
    dtype_backend, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    subset = df.iloc[:, 1:]
    subset._mgr._verify_integrity()

    assert subset.index is not df.index
    assert np.shares_memory(get_array(subset, "b"), get_array(df, "b"))

    subset.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(subset, "b"), get_array(df, "b"))

    expected = DataFrame({"b": [0, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)})
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (also not for BlockManager case,
    # except for single block)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
@pytest.mark.parametrize(
    "row_indexer",
    [slice(1, 2), np.array([False, True, True]), np.array([1, 2])],
    ids=["slice", "mask", "array"],
)
@pytest.mark.parametrize(
    "column_indexer",
    [slice("b", "c"), np.array([False, True, True]), ["b", "c"]],
    ids=["slice", "mask", "array"],
)
def test_subset_loc_rows_columns(
    backend,
    dtype,
    row_indexer,
    column_indexer,
):
    # Case: taking a subset of the rows+columns of a DataFrame using .loc
    # + afterwards modifying the subset
    # Generic test for several combinations of row/column indexers, not all
    # of those could actually return a view / need CoW (so this test is not
    # checking memory sharing, only ensuring subsequent mutation doesn't
    # affect the parent dataframe)
    dtype_backend, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    subset = df.loc[row_indexer, column_indexer]

    assert subset.index is not df.index
    assert subset.columns is not df.columns

    # modifying the subset never modifies the parent
    subset.iloc[0, 0] = 0

    expected = DataFrame(
        {"b": [0, 6], "c": np.array([8, 9], dtype=dtype)}, index=range(1, 3)
    )
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
@pytest.mark.parametrize(
    "row_indexer",
    [slice(1, 3), np.array([False, True, True]), np.array([1, 2])],
    ids=["slice", "mask", "array"],
)
@pytest.mark.parametrize(
    "column_indexer",
    [slice(1, 3), np.array([False, True, True]), [1, 2]],
    ids=["slice", "mask", "array"],
)
def test_subset_iloc_rows_columns(
    backend,
    dtype,
    row_indexer,
    column_indexer,
):
    # Case: taking a subset of the rows+columns of a DataFrame using .iloc
    # + afterwards modifying the subset
    # Generic test for several combinations of row/column indexers, not all
    # of those could actually return a view / need CoW (so this test is not
    # checking memory sharing, only ensuring subsequent mutation doesn't
    # affect the parent dataframe)
    dtype_backend, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    subset = df.iloc[row_indexer, column_indexer]

    assert subset.index is not df.index
    assert subset.columns is not df.columns

    # modifying the subset never modifies the parent
    subset.iloc[0, 0] = 0

    expected = DataFrame(
        {"b": [0, 6], "c": np.array([8, 9], dtype=dtype)}, index=range(1, 3)
    )
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "indexer",
    [slice(0, 2), np.array([True, True, False]), np.array([0, 1])],
    ids=["slice", "mask", "array"],
)
def test_subset_set_with_row_indexer(backend, indexer_si, indexer):
    # Case: setting values with a row indexer on a viewing subset
    # subset[indexer] = value and subset.iloc[indexer] = value
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3, 4], "b": [4, 5, 6, 7], "c": [0.1, 0.2, 0.3, 0.4]})
    df_orig = df.copy()
    subset = df[1:4]

    if (
        indexer_si is tm.setitem
        and isinstance(indexer, np.ndarray)
        and indexer.dtype == "int"
    ):
        pytest.skip("setitem with labels selects on columns")

    indexer_si(subset)[indexer] = 0

    expected = DataFrame(
        {"a": [0, 0, 4], "b": [0, 0, 7], "c": [0.0, 0.0, 0.4]}, index=range(1, 4)
    )
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (CoW)
    tm.assert_frame_equal(df, df_orig)


def test_subset_set_with_mask(backend):
    # Case: setting values with a mask on a viewing subset: subset[mask] = value
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3, 4], "b": [4, 5, 6, 7], "c": [0.1, 0.2, 0.3, 0.4]})
    df_orig = df.copy()
    subset = df[1:4]

    mask = subset > 3

    subset[mask] = 0

    expected = DataFrame(
        {"a": [2, 3, 0], "b": [0, 0, 0], "c": [0.20, 0.3, 0.4]}, index=range(1, 4)
    )
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


def test_subset_set_column(backend):
    # Case: setting a single column on a viewing subset -> subset[col] = value
    dtype_backend, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    subset = df[1:3]

    if dtype_backend == "numpy":
        arr = np.array([10, 11], dtype="int64")
    else:
        arr = pd.array([10, 11], dtype="Int64")

    subset["a"] = arr
    subset._mgr._verify_integrity()
    expected = DataFrame(
        {"a": [10, 11], "b": [5, 6], "c": [0.2, 0.3]}, index=range(1, 3)
    )
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_set_column_with_loc(backend, dtype):
    # Case: setting a single column with loc on a viewing subset
    # -> subset.loc[:, col] = value
    _, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()
    subset = df[1:3]

    subset.loc[:, "a"] = np.array([10, 11], dtype="int64")

    subset._mgr._verify_integrity()
    expected = DataFrame(
        {"a": [10, 11], "b": [5, 6], "c": np.array([8, 9], dtype=dtype)},
        index=range(1, 3),
    )
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (CoW)
    tm.assert_frame_equal(df, df_orig)


def test_subset_set_column_with_loc2(backend):
    # Case: setting a single column with loc on a viewing subset
    # -> subset.loc[:, col] = value
    # separate test for case of DataFrame of a single column -> takes a separate
    # code path
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    subset = df[1:3]

    subset.loc[:, "a"] = 0

    subset._mgr._verify_integrity()
    expected = DataFrame({"a": [0, 0]}, index=range(1, 3))
    tm.assert_frame_equal(subset, expected)
    # original parent dataframe is not modified (CoW)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_set_columns(backend, dtype):
    # Case: setting multiple columns on a viewing subset
    # -> subset[[col1, col2]] = value
    dtype_backend, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()
    subset = df[1:3]

    subset[["a", "c"]] = 0

    subset._mgr._verify_integrity()
    # first and third column should certainly have no references anymore
    assert all(subset._mgr._has_no_reference(i) for i in [0, 2])
    expected = DataFrame({"a": [0, 0], "b": [5, 6], "c": [0, 0]}, index=range(1, 3))
    if dtype_backend == "nullable":
        # there is not yet a global option, so overriding a column by setting a scalar
        # defaults to numpy dtype even if original column was nullable
        expected["a"] = expected["a"].astype("int64")
        expected["c"] = expected["c"].astype("int64")

    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "indexer",
    [slice("a", "b"), np.array([True, True, False]), ["a", "b"]],
    ids=["slice", "mask", "array"],
)
def test_subset_set_with_column_indexer(backend, indexer):
    # Case: setting multiple columns with a column indexer on a viewing subset
    # -> subset.loc[:, [col1, col2]] = value
    _, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "c": [4, 5, 6]})
    df_orig = df.copy()
    subset = df[1:3]

    subset.loc[:, indexer] = 0

    subset._mgr._verify_integrity()
    expected = DataFrame({"a": [0, 0], "b": [0.0, 0.0], "c": [5, 6]}, index=range(1, 3))
    tm.assert_frame_equal(subset, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df[["a", "b"]][0:2],
        lambda df: df[0:2][["a", "b"]],
        lambda df: df[["a", "b"]].iloc[0:2],
        lambda df: df[["a", "b"]].loc[0:1],
        lambda df: df[0:2].iloc[:, 0:2],
        lambda df: df[0:2].loc[:, "a":"b"],  # type: ignore[misc]
    ],
    ids=[
        "row-getitem-slice",
        "column-getitem",
        "row-iloc-slice",
        "row-loc-slice",
        "column-iloc-slice",
        "column-loc-slice",
    ],
)
@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_chained_getitem(
    request,
    backend,
    method,
    dtype,
):
    # Case: creating a subset using multiple, chained getitem calls using views
    # still needs to guarantee proper CoW behaviour
    _, DataFrame, _ = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    # modify subset -> don't modify parent
    subset = method(df)

    subset.iloc[0, 0] = 0
    tm.assert_frame_equal(df, df_orig)

    # modify parent -> don't modify subset
    subset = method(df)
    df.iloc[0, 0] = 0
    expected = DataFrame({"a": [1, 2], "b": [4, 5]})
    tm.assert_frame_equal(subset, expected)


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_subset_chained_getitem_column(backend, dtype):
    # Case: creating a subset using multiple, chained getitem calls using views
    # still needs to guarantee proper CoW behaviour
    dtype_backend, DataFrame, Series = backend
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    # modify subset -> don't modify parent
    subset = df[:]["a"][0:2]
    subset.iloc[0] = 0
    tm.assert_frame_equal(df, df_orig)

    # modify parent -> don't modify subset
    subset = df[:]["a"][0:2]
    df.iloc[0, 0] = 0
    expected = Series([1, 2], name="a")
    tm.assert_series_equal(subset, expected)


@pytest.mark.parametrize(
    "method",
    [
        lambda s: s["a":"c"]["a":"b"],  # type: ignore[misc]
        lambda s: s.iloc[0:3].iloc[0:2],
        lambda s: s.loc["a":"c"].loc["a":"b"],  # type: ignore[misc]
        lambda s: s.loc["a":"c"]  # type: ignore[misc]
        .iloc[0:3]
        .iloc[0:2]
        .loc["a":"b"]  # type: ignore[misc]
        .iloc[0:1],
    ],
    ids=["getitem", "iloc", "loc", "long-chain"],
)
def test_subset_chained_getitem_series(backend, method):
    # Case: creating a subset using multiple, chained getitem calls using views
    # still needs to guarantee proper CoW behaviour
    _, _, Series = backend
    s = Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()

    # modify subset -> don't modify parent
    subset = method(s)
    subset.iloc[0] = 0
    tm.assert_series_equal(s, s_orig)

    # modify parent -> don't modify subset
    subset = s.iloc[0:3].iloc[0:2]
    s.iloc[0] = 0
    expected = Series([1, 2], index=["a", "b"])
    tm.assert_series_equal(subset, expected)


def test_subset_chained_single_block_row():
    # not parametrizing this for dtype backend, since this explicitly tests single block
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    df_orig = df.copy()

    # modify subset -> don't modify parent
    subset = df[:].iloc[0].iloc[0:2]
    subset.iloc[0] = 0
    tm.assert_frame_equal(df, df_orig)

    # modify parent -> don't modify subset
    subset = df[:].iloc[0].iloc[0:2]
    df.iloc[0, 0] = 0
    expected = Series([1, 4], index=["a", "b"], name=0)
    tm.assert_series_equal(subset, expected)


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df[:],
        lambda df: df.loc[:, :],
        lambda df: df.loc[:],
        lambda df: df.iloc[:, :],
        lambda df: df.iloc[:],
    ],
    ids=["getitem", "loc", "loc-rows", "iloc", "iloc-rows"],
)
def test_null_slice(backend, method):
    # Case: also all variants of indexing with a null slice (:) should return
    # new objects to ensure we correctly use CoW for the results
    dtype_backend, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    df_orig = df.copy()

    df2 = method(df)

    # we always return new objects (shallow copy), regardless of CoW or not
    assert df2 is not df
    assert df2.index is not df.index
    assert df2.columns is not df.columns

    # and those trigger CoW when mutated
    df2.iloc[0, 0] = 0
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "method",
    [
        lambda s: s[:],
        lambda s: s.loc[:],
        lambda s: s.iloc[:],
    ],
    ids=["getitem", "loc", "iloc"],
)
def test_null_slice_series(backend, method):
    _, _, Series = backend
    s = Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()

    s2 = method(s)

    # we always return new objects, regardless of CoW or not
    assert s2 is not s
    assert s2.index is not s.index

    # and those trigger CoW when mutated
    s2.iloc[0] = 0
    tm.assert_series_equal(s, s_orig)


# TODO add more tests modifying the parent


# -----------------------------------------------------------------------------
# Series -- Indexing operations taking subset + modifying the subset/parent


def test_series_getitem_slice(backend):
    # Case: taking a slice of a Series + afterwards modifying the subset
    _, _, Series = backend
    s = Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()

    subset = s[:]
    assert np.shares_memory(get_array(subset), get_array(s))
    assert subset.index is not s.index

    subset.iloc[0] = 0

    assert not np.shares_memory(get_array(subset), get_array(s))

    expected = Series([0, 2, 3], index=["a", "b", "c"])
    tm.assert_series_equal(subset, expected)

    # original parent series is not modified (CoW)
    tm.assert_series_equal(s, s_orig)


def test_series_getitem_ellipsis():
    # Case: taking a view of a Series using Ellipsis + afterwards modifying the subset
    s = Series([1, 2, 3])
    s_orig = s.copy()

    subset = s[...]
    assert np.shares_memory(get_array(subset), get_array(s))
    assert subset.index is not s.index

    subset.iloc[0] = 0

    assert not np.shares_memory(get_array(subset), get_array(s))

    expected = Series([0, 2, 3])
    tm.assert_series_equal(subset, expected)

    # original parent series is not modified (CoW)
    tm.assert_series_equal(s, s_orig)


@pytest.mark.parametrize(
    "indexer",
    [slice(0, 2), np.array([True, True, False]), np.array([0, 1])],
    ids=["slice", "mask", "array"],
)
def test_series_subset_set_with_indexer(backend, indexer_si, indexer):
    # Case: setting values in a viewing Series with an indexer
    _, _, Series = backend
    s = Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()
    subset = s[:]

    if (
        indexer_si is tm.setitem
        and isinstance(indexer, np.ndarray)
        and indexer.dtype.kind == "i"
    ):
        # In 3.0 we treat integers as always-labels
        with pytest.raises(KeyError):
            indexer_si(subset)[indexer] = 0
        return

    indexer_si(subset)[indexer] = 0
    expected = Series([0, 0, 3], index=["a", "b", "c"])
    tm.assert_series_equal(subset, expected)

    tm.assert_series_equal(s, s_orig)


# -----------------------------------------------------------------------------
# del operator


def test_del_frame(backend):
    # Case: deleting a column with `del` on a viewing child dataframe should
    # not modify parent + update the references
    dtype_backend, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df[:]

    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    del df2["b"]

    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    tm.assert_frame_equal(df, df_orig)
    tm.assert_frame_equal(df2, df_orig[["a", "c"]])
    df2._mgr._verify_integrity()

    df.loc[0, "b"] = 200
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    df_orig = df.copy()

    df2.loc[0, "a"] = 100
    # modifying child after deleting a column still doesn't update parent
    tm.assert_frame_equal(df, df_orig)


def test_del_series(backend):
    _, _, Series = backend
    s = Series([1, 2, 3], index=["a", "b", "c"])
    s_orig = s.copy()
    s2 = s[:]

    assert np.shares_memory(get_array(s), get_array(s2))

    del s2["a"]

    assert not np.shares_memory(get_array(s), get_array(s2))
    tm.assert_series_equal(s, s_orig)
    tm.assert_series_equal(s2, s_orig[["b", "c"]])

    # modifying s2 doesn't need copy on write (due to `del`, s2 is backed by new array)
    values = s2.values
    s2.loc["b"] = 100
    assert values[0] == 100


# -----------------------------------------------------------------------------
# Accessing column as Series


def test_column_as_series(backend):
    # Case: selecting a single column now also uses Copy-on-Write
    dtype_backend, DataFrame, Series = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    s = df["a"]

    assert s.index is not df.index
    assert np.shares_memory(get_array(s, "a"), get_array(df, "a"))

    s[0] = 0

    expected = Series([0, 2, 3], name="a")
    tm.assert_series_equal(s, expected)
    # assert not np.shares_memory(s.values, get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)
    # ensure cached series on getitem is not the changed series
    tm.assert_series_equal(df["a"], df_orig["a"])


def test_column_as_series_set_with_upcast(backend):
    # Case: selecting a single column now also uses Copy-on-Write -> when
    # setting a value causes an upcast, we don't need to update the parent
    # DataFrame through the cache mechanism
    dtype_backend, DataFrame, Series = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    s = df["a"]
    if dtype_backend == "nullable":
        with pytest.raises(TypeError, match="Invalid value"):
            s[0] = "foo"
        expected = Series([1, 2, 3], name="a")
        tm.assert_series_equal(s, expected)
        tm.assert_frame_equal(df, df_orig)
        # ensure cached series on getitem is not the changed series
        tm.assert_series_equal(df["a"], df_orig["a"])
    else:
        with pytest.raises(TypeError, match="Invalid value"):
            s[0] = "foo"


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df["a"],
        lambda df: df.loc[:, "a"],
        lambda df: df.iloc[:, 0],
    ],
    ids=["getitem", "loc", "iloc"],
)
def test_column_as_series_no_item_cache(request, backend, method):
    # Case: selecting a single column (which now also uses Copy-on-Write to protect
    # the view) should always give a new object (i.e. not make use of a cache)
    dtype_backend, DataFrame, _ = backend
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    s1 = method(df)
    s2 = method(df)

    assert s1 is not s2
    assert s1.index is not df.index
    assert s1.index is not s2.index

    s1.iloc[0] = 0

    tm.assert_series_equal(s2, df_orig["a"])
    tm.assert_frame_equal(df, df_orig)


# TODO add tests for other indexing methods on the Series


def test_dataframe_add_column_from_series(backend):
    # Case: adding a new column to a DataFrame from an existing column/series
    # -> delays copy under CoW
    _, DataFrame, Series = backend
    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})

    s = Series([10, 11, 12])
    df["new"] = s
    assert np.shares_memory(get_array(df, "new"), get_array(s))

    # editing series -> doesn't modify column in frame
    s[0] = 0
    expected = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "new": [10, 11, 12]})
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize("val", [100, "a"])
@pytest.mark.parametrize(
    "indexer_func, indexer",
    [
        (tm.loc, (0, "a")),
        (tm.iloc, (0, 0)),
        (tm.loc, ([0], "a")),
        (tm.iloc, ([0], 0)),
        (tm.loc, (slice(None), "a")),
        (tm.iloc, (slice(None), 0)),
    ],
)
@pytest.mark.parametrize(
    "col", [[0.1, 0.2, 0.3], [7, 8, 9]], ids=["mixed-block", "single-block"]
)
def test_set_value_copy_only_necessary_column(indexer_func, indexer, val, col):
    # When setting inplace, only copy column that is modified instead of the whole
    # block (by splitting the block)
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": col})
    df_orig = df.copy()
    view = df[:]

    if val == "a":
        with pytest.raises(TypeError, match="Invalid value"):
            indexer_func(df)[indexer] = val
    else:
        indexer_func(df)[indexer] = val

        assert np.shares_memory(get_array(df, "b"), get_array(view, "b"))
        assert not np.shares_memory(get_array(df, "a"), get_array(view, "a"))
        tm.assert_frame_equal(view, df_orig)


def test_series_midx_slice():
    ser = Series([1, 2, 3], index=pd.MultiIndex.from_arrays([[1, 1, 2], [3, 4, 5]]))
    ser_orig = ser.copy()
    result = ser[1]
    assert np.shares_memory(get_array(ser), get_array(result))
    result.iloc[0] = 100
    tm.assert_series_equal(ser, ser_orig)


def test_getitem_midx_slice():
    df = DataFrame({("a", "x"): [1, 2], ("a", "y"): 1, ("b", "x"): 2})
    df_orig = df.copy()
    new_df = df[("a",)]

    assert not new_df._mgr._has_no_reference(0)

    assert np.shares_memory(get_array(df, ("a", "x")), get_array(new_df, "x"))
    new_df.iloc[0, 0] = 100
    tm.assert_frame_equal(df_orig, df)


def test_series_midx_tuples_slice():
    ser = Series(
        [1, 2, 3],
        index=pd.MultiIndex.from_tuples([((1, 2), 3), ((1, 2), 4), ((2, 3), 4)]),
    )
    result = ser[(1, 2)]
    assert np.shares_memory(get_array(ser), get_array(result))
    result.iloc[0] = 100
    expected = Series(
        [1, 2, 3],
        index=pd.MultiIndex.from_tuples([((1, 2), 3), ((1, 2), 4), ((2, 3), 4)]),
    )
    tm.assert_series_equal(ser, expected)


def test_midx_read_only_bool_indexer():
    # GH#56635
    def mklbl(prefix, n):
        return [f"{prefix}{i}" for i in range(n)]

    idx = pd.MultiIndex.from_product(
        [mklbl("A", 4), mklbl("B", 2), mklbl("C", 4), mklbl("D", 2)]
    )
    cols = pd.MultiIndex.from_tuples(
        [("a", "foo"), ("a", "bar"), ("b", "foo"), ("b", "bah")], names=["lvl0", "lvl1"]
    )
    df = DataFrame(1, index=idx, columns=cols).sort_index().sort_index(axis=1)

    mask = df[("a", "foo")] == 1
    expected_mask = mask.copy()
    result = df.loc[pd.IndexSlice[mask, :, ["C1", "C3"]], :]
    expected = df.loc[pd.IndexSlice[:, :, ["C1", "C3"]], :]
    tm.assert_frame_equal(result, expected)
    tm.assert_series_equal(mask, expected_mask)


def test_loc_enlarging_with_dataframe():
    df = DataFrame({"a": [1, 2, 3]})
    rhs = DataFrame({"b": [1, 2, 3], "c": [4, 5, 6]})
    rhs_orig = rhs.copy()
    df.loc[:, ["b", "c"]] = rhs
    assert np.shares_memory(get_array(df, "b"), get_array(rhs, "b"))
    assert np.shares_memory(get_array(df, "c"), get_array(rhs, "c"))
    assert not df._mgr._has_no_reference(1)

    df.iloc[0, 1] = 100
    tm.assert_frame_equal(rhs, rhs_orig)
