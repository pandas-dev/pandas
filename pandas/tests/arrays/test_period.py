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
    ([0], np.nan, [iNaT, 1, 2]),
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


# period_array

@pytest.mark.parametrize("data, freq, expected", [
    ([pd.Period("2017", "D")], None, [17167]),
    ([pd.Period("2017", "D")], "D", [17167]),
    ([2017], "D", [17167]),
    (["2017"], "D", [17167]),
    ([pd.Period("2017", "D")], pd.tseries.offsets.Day(), [17167]),
    ([pd.Period("2017", "D"), None], None, [17167, iNaT]),
    (pd.Series(pd.date_range("2017", periods=3)), None,
     [17167, 17168, 17169]),
    (pd.date_range("2017", periods=3), None, [17167, 17168, 17169]),
])
def test_period_array_ok(data, freq, expected):
    result = period_array(data, freq=freq).values
    expected = np.asarray(expected, dtype=np.int64)
    tm.assert_numpy_array_equal(result, expected)


def test_from_datetime64_raises():
    arr = pd.date_range("2017", periods=3, freq="D")
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        PeriodArray._from_datetime64(arr, freq="M")


@pytest.mark.parametrize("data, freq, msg", [
    ([pd.Period('2017', 'D'),
      pd.Period('2017', 'A')],
     None,
     "Input has different freq"),
    ([pd.Period('2017', 'D')],
     "A",
     "Input has different freq"),
])
def test_period_array_raises(data, freq, msg):
    with tm.assert_raises_regex(IncompatibleFrequency, msg):
        period_array(data, freq)


def test_period_array_no_data():
    with tm.assert_raises_regex(ValueError, "one of"):
        period_array(None)


def test_asi8():
    result = period_array(['2000', '2001', None], freq='D').asi8
    expected = np.array([10957, 11323, iNaT])
    tm.assert_numpy_array_equal(result, expected)


def test_take_raises():
    arr = period_array(['2000', '2001'], freq='D')
    with tm.assert_raises_regex(IncompatibleFrequency, 'freq'):
        arr.take([0, -1], allow_fill=True,
                 fill_value=pd.Period('2000', freq='W'))
