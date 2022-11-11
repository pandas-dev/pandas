"""
Testing interaction between the different managers (BlockManager, ArrayManager)
"""
import pytest

from pandas.core.dtypes.missing import array_equivalent

import pandas as pd
import pandas._testing as tm
from pandas.core.internals import (
    ArrayManager,
    BlockManager,
    SingleArrayManager,
    SingleBlockManager,
)


def test_equality_comparison_different_dataframe_managers():
    with pd.option_context("mode.data_manager", "array"):
        df_array = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
        df_array1 = df_array.copy()
    with pd.option_context("mode.data_manager", "block"):
        df_block = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})

    with pytest.raises(
        ValueError, match="BlockManager -> ArrayManager operators don't work."
    ):
        # Cannot compare two dataframes of different manager types using ==
        df_block == df_array
    with pytest.raises(
        ValueError, match="BlockManager -> ArrayManager operators don't work."
    ):
        df_block.equals(df_array)

    assert df_array.equals(df_block)
    assert df_array.equals(df_array1)

    # Test that most comparison operators work
    assert (df_array == df_block).all().all()
    assert (df_array >= df_block).all().all()
    assert not (df_array > df_block).any().any()
    assert not (df_array < df_block).any().any()

    # TODO: this direction doesn't work right now
    # assert df_block.equals(df_array)


def test_dataframe_creation():

    with pd.option_context("mode.data_manager", "block"):
        df_block = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "c": [4, 5, 6]})
    assert isinstance(df_block._mgr, BlockManager)

    with pd.option_context("mode.data_manager", "array"):
        df_array = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "c": [4, 5, 6]})
    assert isinstance(df_array._mgr, ArrayManager)

    # also ensure both are seen as equal
    tm.assert_frame_equal(df_block, df_array)

    # conversion from one manager to the other
    result = df_block._as_manager("block")
    assert isinstance(result._mgr, BlockManager)
    result = df_block._as_manager("array")
    assert isinstance(result._mgr, ArrayManager)
    tm.assert_frame_equal(result, df_block)
    assert all(
        array_equivalent(left, right)
        for left, right in zip(result._mgr.arrays, df_array._mgr.arrays)
    )

    result = df_array._as_manager("array")
    assert isinstance(result._mgr, ArrayManager)
    result = df_array._as_manager("block")
    assert isinstance(result._mgr, BlockManager)
    tm.assert_frame_equal(result, df_array)
    assert len(result._mgr.blocks) == 2


def test_series_creation():

    with pd.option_context("mode.data_manager", "block"):
        s_block = pd.Series([1, 2, 3], name="A", index=["a", "b", "c"])
    assert isinstance(s_block._mgr, SingleBlockManager)

    with pd.option_context("mode.data_manager", "array"):
        s_array = pd.Series([1, 2, 3], name="A", index=["a", "b", "c"])
    assert isinstance(s_array._mgr, SingleArrayManager)

    # also ensure both are seen as equal
    tm.assert_series_equal(s_block, s_array)

    # conversion from one manager to the other
    result = s_block._as_manager("block")
    assert isinstance(result._mgr, SingleBlockManager)
    result = s_block._as_manager("array")
    assert isinstance(result._mgr, SingleArrayManager)
    tm.assert_series_equal(result, s_block)

    result = s_array._as_manager("array")
    assert isinstance(result._mgr, SingleArrayManager)
    result = s_array._as_manager("block")
    assert isinstance(result._mgr, SingleBlockManager)
    tm.assert_series_equal(result, s_array)
