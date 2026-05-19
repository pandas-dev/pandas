import re

import pytest

from pandas import (
    ArrowDtype,
    Series,
)
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.mark.parametrize(
    "list_dtype",
    (
        pa.list_(pa.int64()),
        pa.list_(pa.int64(), list_size=3),
        pa.large_list(pa.int64()),
    ),
)
def test_list_getitem(list_dtype):
    ser = Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=ArrowDtype(list_dtype),
        name="a",
    )
    actual = ser.list[1]
    expected = Series([2, None, None], dtype="int64[pyarrow]", name="a")
    tm.assert_series_equal(actual, expected)


def test_list_getitem_index():
    # GH 58425
    ser = Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    actual = ser.list[1]
    expected = Series(
        [2, None, None],
        dtype="int64[pyarrow]",
        index=[1, 3, 7],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_getitem_slice():
    ser = Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    actual = ser.list[1:None:None]
    expected = Series(
        [[2, 3], [None, 5], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_len():
    ser = Series(
        [[1, 2, 3], [4, None], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
        name="a",
    )
    actual = ser.list.len()
    expected = Series([3, 2, None], dtype=ArrowDtype(pa.int32()), name="a")
    tm.assert_series_equal(actual, expected)


def test_list_flatten():
    ser = Series(
        [[1, 2, 3], None, [4, None], [], [7, 8]],
        dtype=ArrowDtype(pa.list_(pa.int64())),
        name="a",
    )
    actual = ser.list.flatten()
    expected = Series(
        [1, 2, 3, 4, None, 7, 8],
        dtype=ArrowDtype(pa.int64()),
        index=[0, 0, 0, 2, 2, 4, 4],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_getitem_slice_invalid():
    ser = Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
    )
    with tm.external_error_raised(pa.ArrowInvalid):
        ser.list[1:None:0]


def test_list_accessor_non_list_dtype():
    ser = Series(
        [1, 2, 4],
        dtype=ArrowDtype(pa.int64()),
    )
    with pytest.raises(
        AttributeError,
        match=re.escape(
            "Can only use the '.list' accessor with 'list[pyarrow]' dtype, "
            "not int64[pyarrow]."
        ),
    ):
        ser.list[1:None:0]


@pytest.mark.parametrize(
    "list_dtype",
    (
        pa.list_(pa.int64()),
        pa.list_(pa.int64(), list_size=3),
        pa.large_list(pa.int64()),
    ),
)
def test_list_getitem_invalid_index(list_dtype):
    ser = Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=ArrowDtype(list_dtype),
    )
    with tm.external_error_raised(pa.ArrowInvalid):
        ser.list[-1]
    with tm.external_error_raised(pa.ArrowInvalid):
        ser.list[5]
    with pytest.raises(ValueError, match="key must be an int or slice, got str"):
        ser.list["abc"]


def test_list_accessor_not_iterable():
    ser = Series(
        [[1, 2, 3], [4, None], None],
        dtype=ArrowDtype(pa.list_(pa.int64())),
    )
    with pytest.raises(TypeError, match="'ListAccessor' object is not iterable"):
        iter(ser.list)
