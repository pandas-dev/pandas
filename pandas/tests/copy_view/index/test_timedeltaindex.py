import numpy as np
import pytest

from pandas import (
    Series,
    Timedelta,
    TimedeltaIndex,
    array,
    timedelta_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

pytestmark = pytest.mark.filterwarnings(
    "ignore:Setting a value on a view:FutureWarning"
)


@pytest.mark.parametrize(
    "cons",
    [
        lambda x: TimedeltaIndex(x),
        lambda x: TimedeltaIndex(TimedeltaIndex(x)),
    ],
)
def test_timedeltaindex(cons):
    dt = timedelta_range("1 day", periods=3)
    ser = Series(dt)
    idx = cons(ser)
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timedelta("5 days")
    tm.assert_index_equal(idx, expected)


def test_constructor_copy_input_timedelta_ndarray_default():
    # GH 63388
    arr = np.array([1, 2], dtype="timedelta64[ns]")
    idx = TimedeltaIndex(arr)
    assert not np.shares_memory(arr, get_array(idx))


def test_constructor_copy_input_timedelta_ea_default():
    # GH 63388
    arr = array([1, 2], dtype="timedelta64[ns]")
    idx = TimedeltaIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_timedeltaindex_readonly_data():
    # GH 63388
    arr = np.array([1, 2], dtype="timedelta64[ns]")
    arr.flags.writeable = False
    ser = Series(TimedeltaIndex(arr))
    assert not np.shares_memory(arr, get_array(ser))
    ser.iloc[0] = Timedelta(days=1)
    expected = Series(
        [Timedelta(days=1), Timedelta(nanoseconds=2)], dtype="timedelta64[ns]"
    )
    tm.assert_series_equal(ser, expected)
