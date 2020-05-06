import pytest

import pandas as pd
from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "data, index, drop_labels, axis, expected_data, expected_index",
    [
        # Unique Index
        ([1, 2], ["one", "two"], ["two"], 0, [1], ["one"]),
        ([1, 2], ["one", "two"], ["two"], "rows", [1], ["one"]),
        ([1, 1, 2], ["one", "two", "one"], ["two"], 0, [1, 2], ["one", "one"]),
        # GH 5248 Non-Unique Index
        ([1, 1, 2], ["one", "two", "one"], "two", 0, [1, 2], ["one", "one"]),
        ([1, 1, 2], ["one", "two", "one"], ["one"], 0, [1], ["two"]),
        ([1, 1, 2], ["one", "two", "one"], "one", 0, [1], ["two"]),
    ],
)
def test_drop_unique_and_non_unique_index(
    data, index, axis, drop_labels, expected_data, expected_index
):

    s = Series(data=data, index=index)
    result = s.drop(drop_labels, axis=axis)
    expected = Series(data=expected_data, index=expected_index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "data, index, drop_labels, axis, error_type, error_desc",
    [
        # single string/tuple-like
        (range(3), list("abc"), "bc", 0, KeyError, "not found in axis"),
        # bad axis
        (range(3), list("abc"), ("a",), 0, KeyError, "not found in axis"),
        (range(3), list("abc"), "one", "columns", ValueError, "No axis named columns"),
    ],
)
def test_drop_exception_raised(data, index, drop_labels, axis, error_type, error_desc):
    ser = Series(data, index=index)
    with pytest.raises(error_type, match=error_desc):
        ser.drop(drop_labels, axis=axis)


def test_drop_with_ignore_errors():
    # errors='ignore'
    s = Series(range(3), index=list("abc"))
    result = s.drop("bc", errors="ignore")
    tm.assert_series_equal(result, s)
    result = s.drop(["a", "d"], errors="ignore")
    expected = s.iloc[1:]
    tm.assert_series_equal(result, expected)

    # GH 8522
    s = Series([2, 3], index=[True, False])
    assert s.index.is_object()
    result = s.drop(True)
    expected = Series([3], index=[False])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("index", [[1, 2, 3], [1, 1, 3]])
@pytest.mark.parametrize("drop_labels", [[], [1], [3]])
def test_drop_empty_list(index, drop_labels):
    # GH 21494
    expected_index = [i for i in index if i not in drop_labels]
    series = pd.Series(index=index, dtype=object).drop(drop_labels)
    expected = pd.Series(index=expected_index, dtype=object)
    tm.assert_series_equal(series, expected)


@pytest.mark.parametrize(
    "data, index, drop_labels",
    [
        (None, [1, 2, 3], [1, 4]),
        (None, [1, 2, 2], [1, 4]),
        ([2, 3], [0, 1], [False, True]),
    ],
)
def test_drop_non_empty_list(data, index, drop_labels):
    # GH 21494 and GH 16877
    dtype = object if data is None else None
    ser = pd.Series(data=data, index=index, dtype=dtype)
    with pytest.raises(KeyError, match="not found in axis"):
        ser.drop(drop_labels)
