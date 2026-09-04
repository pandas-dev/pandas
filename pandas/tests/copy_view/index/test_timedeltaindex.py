import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize(
    "cons",
    [
        lambda x: pd.TimedeltaIndex(x),
        lambda x: pd.TimedeltaIndex(pd.TimedeltaIndex(x)),
    ],
)
def test_timedeltaindex(cons):
    dt = pd.timedelta_range("1 day", periods=3)
    ser = pd.Series(dt)
    idx = cons(ser)
    expected = idx.copy(deep=True)
    ser.iloc[0] = pd.Timedelta("5 days")
    tm.assert_index_equal(idx, expected)


def test_constructor_copy_input_timedelta_ndarray_default():
    # GH 63388
    arr = np.array([1, 2], dtype="timedelta64[ns]")
    idx = pd.TimedeltaIndex(arr)
    assert not np.shares_memory(arr, get_array(idx))


def test_constructor_copy_input_timedelta_ea_default():
    # GH 63388
    arr = pd.array([1, 2], dtype="timedelta64[ns]")
    idx = pd.TimedeltaIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_timedeltaindex_readonly_data():
    # GH 63388
    arr = np.array([1, 2], dtype="timedelta64[ns]")
    arr.flags.writeable = False
    ser = pd.Series(pd.TimedeltaIndex(arr))
    assert not np.shares_memory(arr, get_array(ser))
    ser.iloc[0] = pd.Timedelta(days=1)
    expected = pd.Series(
        [pd.Timedelta(days=1), pd.Timedelta(nanoseconds=2)], dtype="timedelta64[ns]"
    )
    tm.assert_series_equal(ser, expected)
