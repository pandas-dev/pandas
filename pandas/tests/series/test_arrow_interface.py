import ctypes

import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow", minversion="16.0")


def test_series_arrow_interface():
    s = pd.Series([1, 4, 2])

    capsule = s.__arrow_c_stream__()
    assert (
        ctypes.pythonapi.PyCapsule_IsValid(
            ctypes.py_object(capsule), b"arrow_array_stream"
        )
        == 1
    )

    ca = pa.chunked_array(s)
    expected = pa.chunked_array([[1, 4, 2]])
    assert ca.equals(expected)
    ca = pa.chunked_array(s, type=pa.int32())
    expected = pa.chunked_array([[1, 4, 2]], type=pa.int32())
    assert ca.equals(expected)


def test_series_arrow_interface_arrow_dtypes():
    s = pd.Series([1, 4, 2], dtype="Int64[pyarrow]")

    capsule = s.__arrow_c_stream__()
    assert (
        ctypes.pythonapi.PyCapsule_IsValid(
            ctypes.py_object(capsule), b"arrow_array_stream"
        )
        == 1
    )

    ca = pa.chunked_array(s)
    expected = pa.chunked_array([[1, 4, 2]])
    assert ca.equals(expected)
    ca = pa.chunked_array(s, type=pa.int32())
    expected = pa.chunked_array([[1, 4, 2]], type=pa.int32())
    assert ca.equals(expected)


def test_series_arrow_interface_stringdtype():
    s = pd.Series(["foo", "bar"], dtype="string[pyarrow]")

    capsule = s.__arrow_c_stream__()
    assert (
        ctypes.pythonapi.PyCapsule_IsValid(
            ctypes.py_object(capsule), b"arrow_array_stream"
        )
        == 1
    )

    ca = pa.chunked_array(s)
    expected = pa.chunked_array([["foo", "bar"]], type=pa.large_string())
    assert ca.equals(expected)


class ArrowArrayWrapper:
    def __init__(self, array):
        self.array = array

    def __arrow_c_array__(self, requested_schema=None):
        return self.array.__arrow_c_array__(requested_schema)


class ArrowStreamWrapper:
    def __init__(self, chunked_array):
        self.stream = chunked_array

    def __arrow_c_stream__(self, requested_schema=None):
        return self.stream.__arrow_c_stream__(requested_schema)


@td.skip_if_no("pyarrow", min_version="14.0")
def test_dataframe_from_arrow():
    # objects with __arrow_c_stream__
    arr = pa.chunked_array([[1, 2, 3], [4, 5]])

    result = pd.Series.from_arrow(arr)
    expected = pd.Series([1, 2, 3, 4, 5])
    tm.assert_series_equal(result, expected)

    # not only pyarrow object are supported
    result = pd.Series.from_arrow(ArrowStreamWrapper(arr))
    tm.assert_series_equal(result, expected)

    # table works as well, but will be seen as a StructArray
    table = pa.table({"a": [1, 2, 3], "b": ["a", "b", "c"]})

    result = pd.Series.from_arrow(table)
    expected = pd.Series([{"a": 1, "b": "a"}, {"a": 2, "b": "b"}, {"a": 3, "b": "c"}])
    tm.assert_series_equal(result, expected)

    # objects with __arrow_c_array__
    arr = pa.array([1, 2, 3])

    expected = pd.Series([1, 2, 3])
    result = pd.Series.from_arrow(arr)
    tm.assert_series_equal(result, expected)

    result = pd.Series.from_arrow(ArrowArrayWrapper(arr))
    tm.assert_series_equal(result, expected)

    # only accept actual Arrow objects
    with pytest.raises(
        TypeError, match="Expected an Arrow-compatible array-like object"
    ):
        pd.Series.from_arrow([1, 2, 3])
