import pytest

from pandas import (
    Period,
    PeriodIndex,
    Series,
    period_range,
)
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Setting a value on a view:FutureWarning"
)


@pytest.mark.parametrize("box", [lambda x: x, PeriodIndex])
def test_periodindex(box):
    dt = period_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = box(PeriodIndex(ser))
    expected = idx.copy(deep=True)
    ser.iloc[0] = Period("2020-12-31")
    tm.assert_index_equal(idx, expected)
