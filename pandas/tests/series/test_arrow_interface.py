import ctypes

import pytest

import pandas as pd

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
