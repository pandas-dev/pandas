import numpy as np

from pandas import (
    Interval,
    IntervalIndex,
    Series,
    array,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_constructor_copy_input_interval_ea_default():
    # GH 63388
    arr = array([Interval(0, 1), Interval(1, 2)])
    idx = IntervalIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_intervalindex_readonly_data():
    # GH 63388
    arr = array([Interval(0, 1), Interval(1, 2)])
    arr._left.flags.writeable = False
    arr._right.flags.writeable = False
    ser = Series(IntervalIndex(arr))
    assert not np.shares_memory(arr._left, get_array(ser)._left)
    ser.iloc[0] = Interval(5, 6)
    expected = Series([Interval(5, 6), Interval(1, 2)], dtype="interval[int64, right]")
    tm.assert_series_equal(ser, expected)
