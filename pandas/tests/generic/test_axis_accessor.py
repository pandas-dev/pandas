import pytest

import pandas as pd
import pandas._testing as tm
from pandas.api.extensions import register_index_accessor
from pandas.tests.test_register_accessor import ensure_removed


@pytest.fixture
def multiindex():
    return pd.MultiIndex.from_tuples([("x", 1), ("y", 2)], names=["k", "n"])


def test_ax_to_flat_index(frame_or_series, multiindex):
    obj = frame_or_series([1, 2], index=multiindex)
    result = obj.ax.to_flat_index()
    expected = frame_or_series(
        [1, 2], index=pd.Index([("x", 1), ("y", 2)], tupleize_cols=False)
    )
    tm.assert_equal(result, expected)
    # the original is unchanged
    assert obj.index is multiindex


@pytest.mark.parametrize(
    "method, args",
    [
        ("droplevel", ("n",)),
        ("swaplevel", ()),
        ("reorder_levels", ([1, 0],)),
    ],
)
@pytest.mark.parametrize("axis", [0, 1, "index", "columns"])
def test_ax_matches_ndframe_method(multiindex, method, args, axis):
    df = pd.DataFrame([[1, 2], [3, 4]], index=multiindex, columns=multiindex)
    result = getattr(df.ax, method)(*args, axis=axis)
    expected = getattr(df, method)(*args, axis=axis)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["tz_localize", "to_period"])
def test_ax_matches_ndframe_method_datetime(frame_or_series, method):
    obj = frame_or_series([1, 2], index=pd.date_range("2020", periods=2))
    arg = "UTC" if method == "tz_localize" else "D"
    result = getattr(obj.ax, method)(arg)
    expected = getattr(obj, method)(arg)
    tm.assert_equal(result, expected)


def test_ax_str_accessor_columns():
    df = pd.DataFrame({"a": [1], "b": [2]})
    result = df.ax.str.upper(axis=1)
    expected = pd.DataFrame({"A": [1], "B": [2]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method, args", [("delete", (0,)), ("get_loc", (("x", 1),))])
def test_ax_rejects_non_length_preserving(frame_or_series, multiindex, method, args):
    obj = frame_or_series([1, 2], index=multiindex)
    msg = f"MultiIndex.{method} did not return an Index of the same length"
    with pytest.raises(TypeError, match=msg):
        getattr(obj.ax, method)(*args)


def test_ax_rejects_non_length_preserving_columns(multiindex):
    # the message names the class of the axis that was used
    df = pd.DataFrame([[1, 2]], columns=multiindex)
    with pytest.raises(TypeError, match="MultiIndex.delete did not return"):
        df.ax.delete(0, axis=1)
    with pytest.raises(TypeError, match="Index.str.cat did not return"):
        pd.DataFrame({"a": [1]}).ax.str.cat(axis=1)


@pytest.mark.parametrize("name", ["names", "nonexistent", "swaplevel"])
def test_ax_rejects_non_methods(frame_or_series, name):
    # swaplevel is a MultiIndex method, but neither axis is a MultiIndex
    obj = frame_or_series([1, 2], index=["x", "y"])
    msg = f"no {type(obj).__name__} axis has an Index method '{name}'"
    with pytest.raises(AttributeError, match=msg):
        getattr(obj.ax, name)
    with pytest.raises(AttributeError, match=f"Index method 'str.{name}'"):
        getattr(obj.ax.str, name)


def test_ax_series_rejects_axis_1():
    ser = pd.Series([1, 2], index=pd.Index(["a", "b"]))
    with pytest.raises(ValueError, match="No axis named 1"):
        ser.ax.str.upper(axis=1)


def test_ax_dir(multiindex):
    df = pd.DataFrame({"a": [1, 2]}, index=multiindex)
    names = dir(df.ax)
    # MultiIndex-only methods are listed when either axis is a MultiIndex
    assert "swaplevel" in names
    assert "str" in names
    assert "names" not in names
    assert "upper" in dir(df.ax.str)
    assert df.ax.droplevel.__doc__ == pd.Index.droplevel.__doc__


def test_ax_rejects_inplace(frame_or_series):
    obj = frame_or_series([1, 2], index=pd.Index(["x", "y"], name="k"))
    with pytest.raises(ValueError, match="ax does not support inplace=True"):
        obj.ax.set_names("z", inplace=True)
    assert obj.index.name == "k"
    assert obj.ax.set_names("z", inplace=False).index.name == "z"


def test_ax_method_on_other_axis_only(multiindex):
    # swaplevel exists only on the columns' type
    df = pd.DataFrame([[1, 2]], columns=multiindex)
    result = df.ax.swaplevel(axis=1)
    tm.assert_frame_equal(result, df.swaplevel(axis=1))
    with pytest.raises(AttributeError, match="has no attribute 'swaplevel'"):
        df.ax.swaplevel()


def test_ax_registered_index_accessor():
    class Upper:
        def __init__(self, index) -> None:
            self._index = index

        def apply(self):
            return self._index.str.upper()

    with ensure_removed(pd.Index, "axtest"):
        register_index_accessor("axtest")(Upper)
        df = pd.DataFrame({"a": [1]}, index=["x"])
        result = df.ax.axtest.apply()
        expected = pd.DataFrame({"a": [1]}, index=["X"])
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "name",
    [
        "sort_values",
        "take",
        "union",
        "intersection",
        "difference",
        "symmetric_difference",
        "join",
    ],
)
def test_ax_rejects_reordering_methods(frame_or_series, name):
    # e.g. index.union(["a"]) sorts ["b", "a"] without moving the data
    obj = frame_or_series([10, 20], index=["b", "a"])
    with pytest.raises(AttributeError, match=f"does not support '{name}'"):
        getattr(obj.ax, name)
    assert not hasattr(obj.ax, name)
    assert name not in dir(obj.ax)


def test_ax_accessor_not_valid_on_any_axis():
    ser = pd.Series([1, 2], index=pd.date_range("2020", periods=2))
    assert "str" not in dir(ser.ax)
    with pytest.raises(AttributeError, match="no Series axis supports the 'str'"):
        ser.ax.str
