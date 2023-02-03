import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_assigning_to_same_variable_removes_references(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3]})
    df = df.reset_index()
    arr = get_array(df, "a")
    df.iloc[0, 1] = 100  # Write into a

    assert np.shares_memory(arr, get_array(df, "a"))


def test_setitem_dont_track_unnecessary_references(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 1})

    df["b"] = 100
    arr = get_array(df, "a")
    df.iloc[0, 0] = 100  # Check that we correctly removed the reference
    assert np.shares_memory(arr, get_array(df, "a"))


def test_setitem_with_view_copies(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 1})
    view = df[:]
    expected = df.copy()

    df["b"] = 100
    arr = get_array(df, "a")
    df.iloc[0, 0] = 100  # Check that we correctly track reference
    if using_copy_on_write:
        assert not np.shares_memory(arr, get_array(df, "a"))
        tm.assert_frame_equal(view, expected)


def test_setitem_with_view_invalidated_does_not_copy(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 1})
    view = df[:]

    df["b"] = 100
    arr = get_array(df, "a")
    view = None  # noqa
    df.iloc[0, 0] = 100  # Check that we correctly release the reference
    if using_copy_on_write:
        pytest.mark.xfail("blk.delete does not track references correctly")
        assert np.shares_memory(arr, get_array(df, "a"))


def test_out_of_scope(using_copy_on_write):
    def func():
        df = DataFrame({"a": [1, 2], "b": 1.5, "c": 1})
        # create some subset
        result = df[["a", "b"]]
        return result

    result = func()
    if using_copy_on_write:
        assert not result._mgr.blocks[0].refs.has_reference()
        assert not result._mgr.blocks[1].refs.has_reference()


def test_delete(using_copy_on_write):
    df = DataFrame(np.random.randn(4, 3), columns=["a", "b", "c"])
    del df["b"]
    if using_copy_on_write:
        # TODO: This should not have references, delete makes a shallow copy
        # but keeps the blocks alive
        assert df._mgr.blocks[0].refs.has_reference()
        assert df._mgr.blocks[1].refs.has_reference()

    df = df[["a"]]
    if using_copy_on_write:
        assert not df._mgr.blocks[0].refs.has_reference()
