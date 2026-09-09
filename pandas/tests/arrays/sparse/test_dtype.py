import re

import numpy as np
import pytest

import pandas as pd


@pytest.mark.parametrize(
    "dtype, fill_value",
    [
        ("int", 0),
        ("float", np.nan),
        ("bool", False),
        ("object", np.nan),
        ("datetime64[ns]", np.datetime64("NaT", "ns")),
        ("timedelta64[ns]", np.timedelta64("NaT", "ns")),
    ],
)
def test_inferred_dtype(dtype, fill_value):
    sparse_dtype = pd.SparseDtype(dtype)
    result = sparse_dtype.fill_value
    if pd.isna(fill_value):
        assert pd.isna(result) and type(result) == type(fill_value)
    else:
        assert result == fill_value


def test_from_sparse_dtype():
    dtype = pd.SparseDtype("float", 0)
    result = pd.SparseDtype(dtype)
    assert result.fill_value == 0


def test_from_sparse_dtype_fill_value():
    dtype = pd.SparseDtype("int", 1)
    result = pd.SparseDtype(dtype, fill_value=2)
    expected = pd.SparseDtype("int", 2)
    assert result == expected


@pytest.mark.parametrize(
    "dtype, fill_value",
    [
        ("int", None),
        ("float", None),
        ("bool", None),
        ("object", None),
        ("datetime64[ns]", None),
        ("timedelta64[ns]", None),
        ("int", np.nan),
        ("float", 0),
    ],
)
def test_equal(dtype, fill_value):
    a = pd.SparseDtype(dtype, fill_value)
    b = pd.SparseDtype(dtype, fill_value)
    assert a == b
    assert b == a


def test_nans_equal():
    a = pd.SparseDtype(float, float("nan"))
    b = pd.SparseDtype(float, np.nan)
    assert a == b
    assert b == a


def test_nans_not_equal():
    # GH 54770
    a = pd.SparseDtype(float, 0)
    b = pd.SparseDtype(float, pd.NA)
    assert a != b
    assert b != a


tups = [
    (pd.SparseDtype("float64"), pd.SparseDtype("float32")),
    (pd.SparseDtype("float64"), pd.SparseDtype("float64", 0)),
    (pd.SparseDtype("float64"), pd.SparseDtype("datetime64[ns]", np.nan)),
    (pd.SparseDtype("float64"), np.dtype("float64")),
]


@pytest.mark.parametrize(
    "a, b",
    tups,
)
def test_not_equal(a, b):
    assert a != b


def test_construct_from_string_raises():
    with pytest.raises(
        TypeError, match="Cannot construct a 'SparseDtype' from 'not a dtype'"
    ):
        pd.SparseDtype.construct_from_string("not a dtype")


@pytest.mark.parametrize(
    "dtype, expected",
    [
        (int, True),
        (float, True),
        (bool, True),
        (object, False),
        (str, False),
    ],
)
def test_is_numeric(dtype, expected):
    assert pd.SparseDtype(dtype)._is_numeric is expected


def test_str_uses_object():
    result = pd.SparseDtype(str).subtype
    assert result == np.dtype("object")


@pytest.mark.parametrize(
    "string, expected",
    [
        ("Sparse[float64]", pd.SparseDtype(np.dtype("float64"))),
        ("Sparse[float32]", pd.SparseDtype(np.dtype("float32"))),
        ("Sparse[int]", pd.SparseDtype(np.dtype("int"))),
        ("Sparse[str]", pd.SparseDtype(np.dtype("str"))),
        ("Sparse[datetime64[ns]]", pd.SparseDtype(np.dtype("datetime64[ns]"))),
        ("Sparse", pd.SparseDtype(np.dtype("float"), np.nan)),
    ],
)
def test_construct_from_string(string, expected):
    result = pd.SparseDtype.construct_from_string(string)
    assert result == expected


@pytest.mark.parametrize(
    "a, b, expected",
    [
        (pd.SparseDtype(float, 0.0), pd.SparseDtype(np.dtype("float"), 0.0), True),
        (pd.SparseDtype(int, 0), pd.SparseDtype(int, 0), True),
        (pd.SparseDtype(float, float("nan")), pd.SparseDtype(float, np.nan), True),
        (pd.SparseDtype(float, 0), pd.SparseDtype(float, np.nan), False),
        (pd.SparseDtype(int, 0.0), pd.SparseDtype(float, 0.0), False),
    ],
)
def test_hash_equal(a, b, expected):
    result = a == b
    assert result is expected

    result = hash(a) == hash(b)
    assert result is expected


@pytest.mark.parametrize(
    "string, expected",
    [
        ("Sparse[int]", "int"),
        ("Sparse[int, 0]", "int"),
        ("Sparse[int64]", "int64"),
        ("Sparse[int64, 0]", "int64"),
        ("Sparse[datetime64[ns], 0]", "datetime64[ns]"),
    ],
)
def test_parse_subtype(string, expected):
    subtype, _ = pd.SparseDtype._parse_subtype(string)
    assert subtype == expected


@pytest.mark.parametrize(
    "string", ["Sparse[int, 1]", "Sparse[float, 0.0]", "Sparse[bool, True]"]
)
def test_construct_from_string_fill_value_raises(string):
    with pytest.raises(TypeError, match="fill_value in the string is not"):
        pd.SparseDtype.construct_from_string(string)


@pytest.mark.parametrize(
    "original, dtype, expected",
    [
        (pd.SparseDtype(int, 0), float, pd.SparseDtype(float, 0.0)),
        (pd.SparseDtype(int, 1), float, pd.SparseDtype(float, 1.0)),
        (pd.SparseDtype(int, 1), np.str_, pd.SparseDtype(object, "1")),
        (pd.SparseDtype(float, 1.5), int, pd.SparseDtype(int, 1)),
    ],
)
def test_update_dtype(original, dtype, expected):
    result = original.update_dtype(dtype)
    assert result == expected


@pytest.mark.parametrize(
    "original, dtype, expected_error_msg",
    [
        (
            pd.SparseDtype(float, np.nan),
            int,
            re.escape("Cannot convert non-finite values (NA or inf) to integer"),
        ),
        (
            pd.SparseDtype(str, "abc"),
            int,
            r"invalid literal for int\(\) with base 10: ('abc'|np\.str_\('abc'\))",
        ),
    ],
)
def test_update_dtype_raises(original, dtype, expected_error_msg):
    with pytest.raises(ValueError, match=expected_error_msg):
        original.update_dtype(dtype)


def test_repr():
    # GH-34352
    result = str(pd.SparseDtype("int64", fill_value=0))
    expected = "Sparse[int64, 0]"
    assert result == expected

    result = str(pd.SparseDtype(object, fill_value="0"))
    expected = "Sparse[object, '0']"
    assert result == expected


def test_sparse_dtype_subtype_must_be_numpy_dtype():
    # GH#53160
    msg = "SparseDtype subtype must be a numpy dtype"
    with pytest.raises(TypeError, match=msg):
        pd.SparseDtype("category", fill_value="c")
