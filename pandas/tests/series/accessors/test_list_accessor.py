import re

import pytest

import pandas as pd
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
    ser = pd.Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=pd.ArrowDtype(list_dtype),
        name="a",
    )
    actual = ser.list[1]
    expected = pd.Series([2, None, None], dtype="int64[pyarrow]", name="a")
    tm.assert_series_equal(actual, expected)


def test_list_getitem_index():
    # GH 58425
    ser = pd.Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    actual = ser.list[1]
    expected = pd.Series(
        [2, None, None],
        dtype="int64[pyarrow]",
        index=[1, 3, 7],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_getitem_slice():
    ser = pd.Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    actual = ser.list[1:None:None]
    expected = pd.Series(
        [[2, 3], [None, 5], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
        index=[1, 3, 7],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_len():
    ser = pd.Series(
        [[1, 2, 3], [4, None], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
        name="a",
    )
    actual = ser.list.len()
    expected = pd.Series([3, 2, None], dtype=pd.ArrowDtype(pa.int32()), name="a")
    tm.assert_series_equal(actual, expected)


def test_list_flatten():
    ser = pd.Series(
        [[1, 2, 3], None, [4, None], [], [7, 8]],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
        name="a",
    )
    actual = ser.list.flatten()
    expected = pd.Series(
        [1, 2, 3, 4, None, 7, 8],
        dtype=pd.ArrowDtype(pa.int64()),
        index=[0, 0, 0, 2, 2, 4, 4],
        name="a",
    )
    tm.assert_series_equal(actual, expected)


def test_list_getitem_slice_invalid():
    ser = pd.Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
    )
    with tm.external_error_raised(pa.ArrowInvalid):
        ser.list[1:None:0]


def test_list_accessor_non_list_dtype():
    ser = pd.Series(
        [1, 2, 4],
        dtype=pd.ArrowDtype(pa.int64()),
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
    ser = pd.Series(
        [[1, 2, 3], [4, None, 5], None],
        dtype=pd.ArrowDtype(list_dtype),
    )
    with tm.external_error_raised(pa.ArrowInvalid):
        ser.list[5]
    with pytest.raises(IndexError, match="list index -4 out of range"):
        ser.list[-4]
    with pytest.raises(ValueError, match="key must be an int or slice, got str"):
        ser.list["abc"]
    result = ser.list[-1]
    expected = pd.Series(
        [3, 5, None],
        dtype=pd.ArrowDtype(pa.int64()),
    )
    tm.assert_series_equal(result, expected)


def test_list_accessor_not_iterable():
    ser = pd.Series(
        [[1, 2, 3], [4, None], None],
        dtype=pd.ArrowDtype(pa.list_(pa.int64())),
    )
    with pytest.raises(TypeError, match="'ListAccessor' object is not iterable"):
        iter(ser.list)

# GH#63221
def test_list_get_negative_index():
    ser = pd.Series(
        [["A", "B"], ["C", "D"]], dtype=pd.ArrowDtype(pa.list_(pa.string())), name="a"
    )
    result = ser.list[-1]
    expected = pd.Series(
        ["B", "D"],
        dtype=pd.ArrowDtype(pa.string()),  # item type, not list_dtype
        name="a",
    )
    tm.assert_series_equal(result, expected)


LIST_DTYPES = (
    pa.list_(pa.string()),
    pa.large_list(pa.string()),
)

@pytest.mark.parametrize("list_dtype", LIST_DTYPES)
@pytest.mark.parametrize("data", ([["A", "B"], ["C", "D"]], [["A", "B"], []]))
def test_list_getitem_negative_out_of_range(list_dtype, data):
    # GH#63221
    ser = pd.Series(data, dtype=pd.ArrowDtype(list_dtype))
    with pytest.raises(IndexError, match="list index -5 out of range"):
        ser.list[-5]


def test_list_getitem_negative_sliced_and_chunked():
    # GH#63221
    char_series = [["A", "B"], ["C", "D", "F"], None]
    ser = pd.Series(char_series, dtype=pd.ArrowDtype(pa.list_(pa.string())))
    result = ser.iloc[1:].list[-1]
    index = [1, 2]
    expected = pd.Series(["F", None], dtype=pd.ArrowDtype(pa.string()), index=index)
    tm.assert_series_equal(result, expected)

    chunked = pd.concat([ser, ser], ignore_index=True)
    result = chunked.list[-1]
    expected = pd.Series(
        ["B", "F", None, "B", "F", None], dtype=pd.ArrowDtype(pa.string())
    )
    tm.assert_series_equal(result, expected)
