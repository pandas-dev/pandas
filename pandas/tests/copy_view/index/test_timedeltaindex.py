import pytest

from pandas import (
    Series,
    Timedelta,
    TimedeltaIndex,
    timedelta_range,
)
import pandas._testing as tm

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
