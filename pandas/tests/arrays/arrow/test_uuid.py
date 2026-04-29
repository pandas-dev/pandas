import uuid

import pytest

from pandas.compat.pyarrow import pa_version_under18p0

import pandas as pd
from pandas.core.arrays import ArrowExtensionArray

pa = pytest.importorskip("pyarrow")

pytestmark = pytest.mark.skipif(
    pa_version_under18p0, reason="pa.uuid() not available before pyarrow 18.0.0"
)


def test_uuid_keeps_dtype_verify() -> None:
    #GH#63511: Verify that Arrow-backed UUIDs maintain their type during construction.
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    arr = pa.array([u.bytes, None], type=pa.uuid())
    s = pd.Series(arr)

    assert s.dtype != object

    ser_array = s.array
    assert isinstance(ser_array, ArrowExtensionArray)
    assert ser_array._pa_array.type == pa.uuid()

    d = pd.DataFrame({"id": arr})
    df_array = d["id"].array
    assert isinstance(df_array, ArrowExtensionArray)
    assert df_array._pa_array.type == pa.uuid()


def test_uuid_isna() -> None:
    #GH#63511: Verify that Arrow-back UUIDs array handles NULL values.
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    arr = pa.array([u.bytes, None], type=pa.uuid())
    s = pd.Series(arr)

    result = s.isna()
    assert not result[0]
    assert result[1]


def test_uuid_comparison_eq() -> None:
    #GH#63511: Verify that Arrow-back UUIDs array handles comparisons correctly. 
    u1 = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    u2 = uuid.UUID("550e8400-e29b-41d4-a716-446655440001")
    u3 = uuid.UUID("550e8400-e29b-41d4-a716-446655440002")
    arr1 = pa.array([u1.bytes, u2.bytes, u3.bytes, None], type=pa.uuid())
    arr2 = pa.array([u1.bytes, u1.bytes, None, None], type=pa.uuid())
    s1 = pd.Series(arr1)
    s2 = pd.Series(arr2)

    # Cast uuid.UUID to PyArrow type with length 16 bytes
    pa_binary16 = pd.ArrowDtype(pa.binary(16))

    out0 = s1[0] == s2[0]
    out1 = s1[1] == s2[1]
    out2 = s1[2] == s2[2]
    out3 = s1[3] == s2[3]
    out4 = (s1.astype(pa_binary16) == s2.astype(pa_binary16)).tolist()

    assert out0
    assert not out1
    assert not out2
    assert pd.isna(out3)
    assert out4 == [True, False, pd.NA, pd.NA]


def test_uuid_getitem_scalar() -> None:
    #GH#63511: Verify that iloc[i] works correctly with an Arrow-back UUID array.
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    arr = pa.array([u.bytes, None], type=pa.uuid())
    s = pd.Series(arr)

    # Extract individual elements
    v0 = s.iloc[0]
    v1 = s.iloc[1]

    assert isinstance(v0, uuid.UUID)
    assert v1 is pd.NA


def test_uuid_contains_behavior() -> None:
    #GH#63511: Verify that "x in s" checks index, "x in s.array" checks value
    u = uuid.uuid4()
    arr = pa.array([u.bytes], type=pa.uuid())
    s = pd.Series(arr)

    # Checks index
    assert (s.iloc[0] in s) is False

    # Checks values
    assert (s.iloc[0] in s.array) is True
    assert (None in s.array) is False


def test_series_from_pyarrow_uuid_chunkedarray() -> None:
    #GH#63511: Verify that chunked arrays correctly implement comparison and null values
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    chunk1 = pa.array([u.bytes], type=pa.uuid())
    chunk2 = pa.array([None], type=pa.uuid())
    carr = pa.chunked_array([chunk1, chunk2])
    s = pd.Series(carr)

    assert s.dtype != object

    ser_array = s.array
    assert isinstance(ser_array, ArrowExtensionArray)
    assert ser_array._pa_array.type == pa.uuid()
    assert ser_array.isna().tolist() == [False, True]
