import numpy as np

import pandas as pd
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_constructor_copy_input_interval_ea_default():
    # GH 63388
    arr = pd.array([pd.Interval(0, 1), pd.Interval(1, 2)])
    idx = pd.IntervalIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_intervalindex_readonly_data():
    # GH 63388
    arr = pd.array([pd.Interval(0, 1), pd.Interval(1, 2)])
    arr._left.flags.writeable = False
    arr._right.flags.writeable = False
    ser = pd.Series(pd.IntervalIndex(arr))
    assert not np.shares_memory(arr._left, get_array(ser)._left)
    ser.iloc[0] = pd.Interval(5, 6)
    expected = pd.Series(
        [pd.Interval(5, 6), pd.Interval(1, 2)], dtype="interval[int64, right]"
    )
    tm.assert_series_equal(ser, expected)
