import numpy as np
import pytest

from pandas import interval_range
import pandas._testing as tm


@pytest.mark.parametrize(
    "kwargs",
    [
        {"start": 0.000001, "periods": 4},
        {"start": 1.000001, "periods": 5},
        {"start": np.datetime64("2003-04-08"), "periods": 3},
        {"start": np.timedelta64(1, "D"), "periods": 3},
        {"start": 5.000001, "end": 10, "inclusive": "left"},
    ],
)
def test_interval_array_almost_equal(kwargs):
    arr = interval_range(**kwargs).values
    tm.assert_interval_array_almost_equal(arr, arr)


def test_interval_array_almost_equal_within_tolerance():
    arr1 = interval_range(start=1.01, periods=5).values
    arr2 = interval_range(start=1.0, periods=5).values
    tm.assert_interval_array_almost_equal(arr1, arr2, atol=0, rtol=0.01)


def test_interval_array_almost_equal_outside_tolerance():
    arr1 = interval_range(start=1.01, periods=5).values
    arr2 = interval_range(start=1.0, periods=5).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(100.0 %\\)
\\[left\\]:  \\[1.01, 2.01, 3.01, 4.01, 5.01\\]
\\[right\\]: \\[1.0, 2.0, 3.0, 4.0, 5.0\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_almost_equal(arr1, arr2, atol=0, rtol=0.001)


def test_interval_array_almost_equal_closed_mismatch():
    kwargs = {"start": 0.000001, "periods": 5}
    arr1 = interval_range(inclusive="left", **kwargs).values
    arr2 = interval_range(inclusive="right", **kwargs).values

    msg = """\
IntervalArray are different

Attribute "inclusive" are different
\\[left\\]:  left
\\[right\\]: right"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_almost_equal(arr1, arr2)


def test_interval_array_almost_equal_periods_mismatch():
    kwargs = {"start": 0}
    arr1 = interval_range(periods=5, **kwargs).values
    arr2 = interval_range(periods=6, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left shapes are different
\\[left\\]:  \\(5,\\)
\\[right\\]: \\(6,\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_almost_equal(arr1, arr2)


def test_interval_array_almost_equal_end_mismatch():
    kwargs = {"start": 0, "periods": 5}
    arr1 = interval_range(end=10, **kwargs).values
    arr2 = interval_range(end=20, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(80.0 %\\)
\\[left\\]:  \\[0, 2, 4, 6, 8\\]
\\[right\\]: \\[0, 4, 8, 12, 16\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_almost_equal(arr1, arr2)


def test_interval_array_almost_equal_start_mismatch():
    kwargs = {"periods": 4}
    arr1 = interval_range(start=0, **kwargs).values
    arr2 = interval_range(start=1, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(100.0 %\\)
\\[left\\]:  \\[0, 1, 2, 3\\]
\\[right\\]: \\[1, 2, 3, 4\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_almost_equal(arr1, arr2)
