import pytest

from pandas import interval_range
import pandas._testing as tm


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(start=0, periods=4),
        dict(start=1, periods=5),
        dict(start=5, end=10, closed="left"),
    ],
)
def test_interval_array_equal(kwargs):
    arr = interval_range(**kwargs).values
    tm.assert_interval_array_equal(arr, arr)


def test_interval_array_equal_closed_mismatch():
    kwargs = dict(start=0, periods=5)
    arr1 = interval_range(closed="left", **kwargs).values
    arr2 = interval_range(closed="right", **kwargs).values

    msg = """\
IntervalArray are different

Attribute "closed" are different
\\[left\\]:  left
\\[right\\]: right"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_periods_mismatch():
    kwargs = dict(start=0)
    arr1 = interval_range(periods=5, **kwargs).values
    arr2 = interval_range(periods=6, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left length are different
\\[left\\]:  5, Int64Index\\(\\[0, 1, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: 6, Int64Index\\(\\[0, 1, 2, 3, 4, 5\\], dtype='int64'\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_end_mismatch():
    kwargs = dict(start=0, periods=5)
    arr1 = interval_range(end=10, **kwargs).values
    arr2 = interval_range(end=20, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(80.0 %\\)
\\[left\\]:  Int64Index\\(\\[0, 2, 4, 6, 8\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[0, 4, 8, 12, 16\\], dtype='int64'\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_start_mismatch():
    kwargs = dict(periods=4)
    arr1 = interval_range(start=0, **kwargs).values
    arr2 = interval_range(start=1, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(100.0 %\\)
\\[left\\]:  Int64Index\\(\\[0, 1, 2, 3\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)
