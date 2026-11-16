import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import IntervalArray


@pytest.fixture(
    params=[
        [0, 1, 3, 6],
        [0.0, 1.5, 3.0, 4.5],
        pd.array([0, 1, 3, 6], dtype="Int64"),
        pd.date_range("2020-01-01", periods=4, tz="UTC"),
        pd.timedelta_range("1D", periods=4),
    ],
    ids=["int64", "float64", "Int64", "datetime64tz", "timedelta64"],
)
def ser(request):
    arr = IntervalArray.from_breaks(request.param)
    return pd.Series(arr, index=pd.Index(list("abc")), name="x")


@pytest.mark.parametrize("attr", ["left", "right", "mid", "length", "is_empty"])
def test_interval_accessor_properties(ser, attr):
    # GH#16401
    result = getattr(ser.interval, attr)
    expected = pd.Series(
        pd.Index(getattr(ser.array, attr)), index=ser.index, name=ser.name
    )
    tm.assert_series_equal(result, expected)


def test_interval_accessor_closed(ser):
    assert ser.interval.closed == "right"
    assert ser.interval.set_closed("both").interval.closed == "both"


def test_interval_accessor_methods(ser):
    arr = ser.array
    for method, arg in [
        ("contains", arr.mid[1]),
        ("overlaps", arr[1]),
        ("set_closed", "left"),
    ]:
        result = getattr(ser.interval, method)(arg)
        expected = pd.Series(
            getattr(arr, method)(arg), index=ser.index, name=ser.name, copy=True
        )
        tm.assert_series_equal(result, expected)
    assert ser.interval.contains(arr.mid[1]).tolist() == [False, True, False]


@pytest.mark.parametrize(
    "method", ["left", "right", "set_closed"], ids=["left", "right", "set_closed"]
)
def test_interval_accessor_result_does_not_share_data(ser, method):
    # mutating the result must not change the parent
    orig = ser.copy()
    if method == "set_closed":
        result = ser.interval.set_closed("both")
    else:
        result = getattr(ser.interval, method)
    result.iloc[0] = result.iloc[1]
    tm.assert_series_equal(ser, orig)


def test_interval_accessor_contains_series():
    ser = pd.Series(IntervalArray.from_breaks([0, 1, 3, 6]), index=[10, 11, 12])
    other = pd.Series([0.5, 2, 9], index=ser.index)
    result = ser.interval.contains(other)
    tm.assert_series_equal(result, pd.Series([True, True, False], index=ser.index))

    # a misaligned Series is not matched by position
    msg = "Can only compare identically-labeled Series objects"
    with pytest.raises(ValueError, match=msg):
        ser.interval.contains(other.iloc[::-1])
    with pytest.raises(ValueError, match=msg):
        ser.interval.contains(other=other.iloc[::-1])


def test_interval_accessor_missing():
    ser = pd.Series(IntervalArray([pd.Interval(0, 1), np.nan]))
    tm.assert_series_equal(ser.interval.left, pd.Series([0.0, np.nan]))
    tm.assert_series_equal(ser.interval.contains(0.5), pd.Series([True, False]))
    tm.assert_series_equal(
        ser.interval.overlaps(pd.Interval(0, 2)), pd.Series([True, False])
    )


@pytest.mark.parametrize(
    "data",
    [
        pd.Series([1, 2]),
        pd.Series(["a", "b"], dtype="category"),
        pd.Series([1, 2], dtype="category"),
        pd.Series([pd.Interval(0, 1), pd.Interval(1, 2)], dtype=object),
    ],
)
def test_interval_accessor_invalid(data):
    msg = "Can only use .interval accessor with an 'interval' dtype"
    with pytest.raises(AttributeError, match=msg):
        data.interval
    assert not hasattr(data, "interval")


def test_interval_accessor_categorical():
    # pd.cut returns a categorical with interval categories
    ser = pd.cut(pd.Series([0.5, 2, 4, np.nan], name="x"), [0, 1, 3, 6])
    dense = pd.Series(IntervalArray(ser), index=ser.index, name=ser.name)
    tm.assert_series_equal(ser.interval.mid, dense.interval.mid)
    tm.assert_series_equal(ser.interval.contains(2), dense.interval.contains(2))
    tm.assert_series_equal(
        ser.interval.set_closed("both"), dense.interval.set_closed("both")
    )
    assert ser.interval.mid.tolist()[:3] == [0.5, 2.0, 4.5]

    # without missing values the integer subtype is kept
    ser = pd.cut(pd.Series([0.5, 2]), [0, 1, 3])
    assert ser.interval.left.dtype == np.int64


def test_interval_accessor_label_attribute_access():
    # a label named "interval" stays reachable on a non-interval Series
    row = pd.DataFrame({"interval": [5, 10], "n": [1, 2]}).iloc[0]
    assert row.interval == 5
    row = pd.DataFrame({"interval": [5, 10]}).astype("category").iloc[0]
    assert row.interval == 5
