import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_consolidate():
    # create unconsolidated DataFrame
    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})
    df["c"] = [4, 5, 6]

    # take a viewing subset
    subset = df[:]

    # each block of subset references a block of df
    assert all(blk.refs.has_reference() for blk in subset._mgr.blocks)

    # consolidate the two int64 blocks
    subset._consolidate_inplace()

    # the float64 block still references the parent one because it still a view
    assert subset._mgr.blocks[0].refs.has_reference()
    # equivalent of assert np.shares_memory(df["b"].values, subset["b"].values)
    # but avoids caching df["b"]
    assert np.shares_memory(get_array(df, "b"), get_array(subset, "b"))

    # the new consolidated int64 block does not reference another
    assert not subset._mgr.blocks[1].refs.has_reference()

    # the parent dataframe now also only is linked for the float column
    assert not df._mgr.blocks[0].refs.has_reference()
    assert df._mgr.blocks[1].refs.has_reference()
    assert not df._mgr.blocks[2].refs.has_reference()

    # and modifying subset still doesn't modify parent
    subset.iloc[0, 1] = 0.0
    assert not df._mgr.blocks[1].refs.has_reference()
    assert df.loc[0, "b"] == 0.1


@pytest.mark.parametrize("dtype", [np.intp, np.int8])
@pytest.mark.parametrize(
    "locs, arr",
    [
        ([0], np.array([-1, -2, -3])),
        ([1], np.array([-1, -2, -3])),
        ([5], np.array([-1, -2, -3])),
        ([0, 1], np.array([[-1, -2, -3], [-4, -5, -6]]).T),
        ([0, 2], np.array([[-1, -2, -3], [-4, -5, -6]]).T),
        ([0, 1, 2], np.array([[-1, -2, -3], [-4, -5, -6], [-4, -5, -6]]).T),
        ([1, 2], np.array([[-1, -2, -3], [-4, -5, -6]]).T),
        ([1, 3], np.array([[-1, -2, -3], [-4, -5, -6]]).T),
    ],
)
def test_iset_splits_blocks_inplace(locs, arr, dtype):
    # Nothing currently calls iset with
    # more than 1 loc with inplace=True (only happens with inplace=False)
    # but ensure that it works
    df = DataFrame(
        {
            "a": [1, 2, 3],
            "b": [4, 5, 6],
            "c": [7, 8, 9],
            "d": [10, 11, 12],
            "e": [13, 14, 15],
            "f": Series(["a", "b", "c"], dtype=object),
        },
    )
    arr = arr.astype(dtype)
    df_orig = df.copy()
    df2 = df.copy(deep=False)  # Trigger a CoW (if enabled, otherwise makes copy)
    df2._mgr.iset(locs, arr, inplace=True)

    tm.assert_frame_equal(df, df_orig)
    for i, col in enumerate(df.columns):
        if i not in locs:
            assert np.shares_memory(get_array(df, col), get_array(df2, col))


def test_exponential_backoff():
    # GH#55518
    df = DataFrame({"a": [1, 2, 3]})
    for i in range(490):
        df.copy(deep=False)

    assert len(df._mgr.blocks[0].refs.referenced_blocks) == 491

    df = DataFrame({"a": [1, 2, 3]})
    dfs = [df.copy(deep=False) for i in range(510)]

    for i in range(20):
        df.copy(deep=False)
    assert len(df._mgr.blocks[0].refs.referenced_blocks) == 531
    assert df._mgr.blocks[0].refs.clear_counter == 1000

    for i in range(500):
        df.copy(deep=False)

    # Don't reduce since we still have over 500 objects alive
    assert df._mgr.blocks[0].refs.clear_counter == 1000

    dfs = dfs[:300]
    for i in range(500):
        df.copy(deep=False)

    # Reduce since there are less than 500 objects alive
    assert df._mgr.blocks[0].refs.clear_counter == 500
