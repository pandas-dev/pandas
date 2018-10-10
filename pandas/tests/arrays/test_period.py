import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas._libs.tslibs import iNaT
from pandas._libs.tslibs.period import IncompatibleFrequency
from pandas.core.arrays import PeriodArray, period_array


@pytest.mark.parametrize('key, value, expected', [
    ([0], pd.Period("2000", "D"), [10957, 1, 2]),
    ([0], None, [iNaT, 1, 2]),
    ([0, 1, 2], pd.Period("2000", "D"), [10957] * 3),
    ([0, 1, 2], [pd.Period("2000", "D"),
                 pd.Period("2001", "D"),
                 pd.Period("2002", "D")],
     [10957, 11323, 11688]),
])
def test_setitem(key, value, expected):
    arr = PeriodArray(np.arange(3), freq="D")
    expected = PeriodArray(expected, freq="D")
    arr[key] = value
    tm.assert_period_array_equal(arr, expected)


def test_setitem_raises():
    arr = PeriodArray(np.arange(3), freq="D")
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        arr[0] = pd.Period("2000", freq="A")

    with tm.assert_raises_regex(ValueError, "length"):
        arr[[0, 1]] = [pd.Period("2000", freq="D")]

    with tm.assert_raises_regex(TypeError, "int"):
        arr[0] = 1


@pytest.mark.parametrize("data, freq, expected", [
    ([pd.Period("2017", "D")], None, [17167]),
    ([pd.Period("2017", "D")], "D", [17167]),
    ([2017], "D", [17167]),
    (["2017"], "D", [17167]),
    ([pd.Period("2017", "D")], pd.tseries.offsets.Day(), [17167]),
    ([pd.Period("2017", "D"), None], None, [17167, iNaT]),
])
def test_to_period_ok(data, freq, expected):
    result = period_array(data, freq=freq).values
    expected = np.asarray(expected)
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("data, freq, msg", [
    ([pd.Period('2017', 'D'),
      pd.Period('2017', 'A')],
     None,
     "Input has different freq"),
    ([pd.Period('2017', 'D')],
     "A",
     "Input has different freq"),
])
def test_to_period_raises(data, freq, msg):
    with tm.assert_raises_regex(IncompatibleFrequency, msg):
        period_array(data, freq)
