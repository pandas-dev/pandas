import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm


def test_none_to_nan():
    a = pd.arrays.DictArray._from_sequence([{"a": 1}, None, {"a": 1}])
    assert a[1] is not None
    assert np.isnan(a[1])


def test_setitem_validates():
    a = pd.arrays.DictArray._from_sequence([{"a": 1}, {"a": 2}])
    with pytest.raises(TypeError, match="Cannot set non-dict value"):
        a[0] = "not_a_dict"

    with pytest.raises(TypeError, match="Cannot set non-dict value"):
        a[:] = np.array([1, 2])


def test_constructor_raises():
    with pytest.raises(ValueError, match="sequence of dicts"):
        pd.arrays.DictArray(np.array(["a", "b"], dtype="S1"))

    with pytest.raises(ValueError, match="sequence of dicts"):
        pd.arrays.DictArray(np.array([]))


def test_reduce_raises():
    arr = pd.Series([{"a": 1}, {"a": 2}, {"a": 3}], dtype="dict")

    with pytest.raises(NotImplementedError):
        arr.sum()


def test_getitem():
    arr = pd.Series([{"a": 1}, {"a": 2}, {"a": 3}], dtype="dict")
    result = arr.dict["a"]
    expected = pd.Series([1, 2, 3])

    tm.assert_series_equal(result, expected)


def test_get():
    arr = pd.Series([{"a": 1}, {"b": 2}, {"a": 3}], dtype="dict")
    result = arr.dict.get("a")
    expected = pd.Series([1.0, np.nan, 3.0])

    tm.assert_series_equal(result, expected)
