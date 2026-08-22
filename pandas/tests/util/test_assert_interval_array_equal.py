import pytest

from pandas import (
    Interval,
    Timestamp,
    interval_range,
)
import pandas._testing as tm
from pandas.arrays import IntervalArray


@pytest.mark.parametrize(
    "kwargs",
    [
        {"start": 0, "periods": 4},
        {"start": 1, "periods": 5},
        {"start": 5, "end": 10, "closed": "left"},
    ],
)
def test_interval_array_equal(kwargs):
    arr = interval_range(**kwargs).values
    tm.assert_interval_array_equal(arr, arr)


def test_interval_array_equal_closed_mismatch():
    kwargs = {"start": 0, "periods": 5}
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
    kwargs = {"start": 0}
    arr1 = interval_range(periods=5, **kwargs).values
    arr2 = interval_range(periods=6, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left shapes are different
\\[left\\]:  \\(5,\\)
\\[right\\]: \\(6,\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_end_mismatch():
    kwargs = {"start": 0, "periods": 5}
    arr1 = interval_range(end=10, **kwargs).values
    arr2 = interval_range(end=20, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(80.0 %\\)
\\[left\\]:  \\[0, 2, 4, 6, 8\\]
\\[right\\]: \\[0, 4, 8, 12, 16\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_start_mismatch():
    kwargs = {"periods": 4}
    arr1 = interval_range(start=0, **kwargs).values
    arr2 = interval_range(start=1, **kwargs).values

    msg = """\
IntervalArray.left are different

IntervalArray.left values are different \\(100.0 %\\)
\\[left\\]:  \\[0, 1, 2, 3\\]
\\[right\\]: \\[1, 2, 3, 4\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_end_mismatch_only():
    arr1 = IntervalArray([Interval(0, 1), Interval(0, 5)])
    arr2 = IntervalArray([Interval(0, 1), Interval(0, 6)])

    msg = """\
IntervalArray.right are different

IntervalArray.right values are different \\(50.0 %\\)
\\[left\\]:  \\[1, 5\\]
\\[right\\]: \\[1, 6\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_tolerance():
    # GH#43913 check_exact=False must apply rtol/atol to the endpoints
    arr1 = IntervalArray.from_tuples([(1.0, 2.0)])
    arr2 = IntervalArray.from_tuples([(1.0000000000001, 2.0)])

    tm.assert_interval_array_equal(arr1, arr2, check_exact=False, rtol=10.0, atol=10.0)


def test_interval_array_equal_tolerance_exceeded():
    # GH#43913 a difference larger than the tolerance must still be reported
    arr1 = IntervalArray.from_tuples([(1.0, 2.0)])
    arr2 = IntervalArray.from_tuples([(1.5, 2.0)])

    msg = r"IntervalArray\.left are different"

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(
            arr1, arr2, check_exact=False, rtol=0, atol=1e-14
        )


def test_interval_array_equal_check_exact_default():
    # GH#43913 the default stays exact, so a tiny difference is still reported
    arr1 = IntervalArray.from_tuples([(1.0, 2.0)])
    arr2 = IntervalArray.from_tuples([(1.0000000000001, 2.0)])

    msg = r"IntervalArray\.left are different"

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(arr1, arr2)


def test_interval_array_equal_tolerance_closed_mismatch():
    # GH#43913 closed is not numeric, a tolerance must not make it match
    arr1 = IntervalArray.from_tuples([(1.0, 2.0)], closed="left")
    arr2 = IntervalArray.from_tuples([(1.0, 2.0)], closed="right")

    msg = 'Attribute "closed" are different'

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(
            arr1, arr2, check_exact=False, rtol=10.0, atol=10.0
        )


def test_interval_array_equal_tolerance_datetimelike():
    # GH#43913 datetimelike endpoints keep the exact comparison
    arr1 = IntervalArray.from_tuples([(Timestamp("2020-01-01"), Timestamp("2021"))])
    arr2 = IntervalArray.from_tuples([(Timestamp("2020-01-02"), Timestamp("2021"))])

    msg = r"IntervalArray\.left\._ndarray are different"

    with pytest.raises(AssertionError, match=msg):
        tm.assert_interval_array_equal(
            arr1, arr2, check_exact=False, rtol=10.0, atol=10.0
        )
