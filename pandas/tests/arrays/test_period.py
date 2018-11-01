import numpy as np
import pytest

from pandas._libs.tslibs import iNaT
from pandas._libs.tslibs.period import IncompatibleFrequency

from pandas.core.dtypes.common import pandas_dtype
from pandas.core.dtypes.dtypes import PeriodDtype

import pandas as pd
from pandas.core.arrays import PeriodArray, period_array
import pandas.util.testing as tm

# ----------------------------------------------------------------------------
# Constructors

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
    result = period_array(data, freq=freq).asi8
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


def test_period_array_non_period_series_raies():
    ser = pd.Series([1, 2, 3])
    with tm.assert_raises_regex(TypeError, 'dtype'):
        PeriodArray(ser, freq='D')


def test_period_array_freq_mismatch():
    arr = period_array(['2000', '2001'], freq='D')
    with tm.assert_raises_regex(IncompatibleFrequency, 'freq'):
        PeriodArray(arr, freq='M')

    with tm.assert_raises_regex(IncompatibleFrequency, 'freq'):
        PeriodArray(arr, freq=pd.tseries.offsets.MonthEnd())


def test_asi8():
    result = period_array(['2000', '2001', None], freq='D').asi8
    expected = np.array([10957, 11323, iNaT])
    tm.assert_numpy_array_equal(result, expected)


def test_take_raises():
    arr = period_array(['2000', '2001'], freq='D')
    with tm.assert_raises_regex(IncompatibleFrequency, 'freq'):
        arr.take([0, -1], allow_fill=True,
                 fill_value=pd.Period('2000', freq='W'))

    with tm.assert_raises_regex(ValueError, 'foo'):
        arr.take([0, -1], allow_fill=True, fill_value='foo')


@pytest.mark.parametrize('dtype', [int, np.int32, np.int64])
def test_astype(dtype):
    # Need to ensure ordinals are astyped correctly for both
    # int32 and 64
    arr = period_array(['2000', '2001', None], freq='D')
    result = arr.astype(dtype)
    # need pandas_dtype to handle int32 vs. int64 correctly
    expected = pandas_dtype(dtype)
    assert result.dtype == expected


def test_astype_copies():
    arr = period_array(['2000', '2001', None], freq='D')
    result = arr.astype(np.int64, copy=False)
    assert result is arr._data

    result = arr.astype(np.int64, copy=True)
    assert result is not arr._data


def test_astype_categorical():
    arr = period_array(['2000', '2001', '2001', None], freq='D')
    result = arr.astype('category')
    categories = pd.PeriodIndex(['2000', '2001'], freq='D')
    expected = pd.Categorical.from_codes([0, 1, 1, -1], categories=categories)
    tm.assert_categorical_equal(result, expected)


def test_astype_period():
    arr = period_array(['2000', '2001', None], freq='D')
    result = arr.astype(PeriodDtype("M"))
    expected = period_array(['2000', '2001', None], freq='M')
    tm.assert_period_array_equal(result, expected)


@pytest.mark.parametrize('other', [
    'datetime64[ns]', 'timedelta64[ns]',
])
def test_astype_datetime(other):
    arr = period_array(['2000', '2001', None], freq='D')
    # slice off the [ns] so that the regex matches.
    with tm.assert_raises_regex(TypeError, other[:-4]):
        arr.astype(other)


def test_fillna_raises():
    arr = period_array(['2000', '2001', '2002'], freq='D')
    with tm.assert_raises_regex(ValueError, 'Length'):
        arr.fillna(arr[:2])


def test_fillna_copies():
    arr = period_array(['2000', '2001', '2002'], freq='D')
    result = arr.fillna(pd.Period("2000", "D"))
    assert result is not arr


# ----------------------------------------------------------------------------
# setitem

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


def test_setitem_raises_incompatible_freq():
    arr = PeriodArray(np.arange(3), freq="D")
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        arr[0] = pd.Period("2000", freq="A")

    other = period_array(['2000', '2001'], freq='A')
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        arr[[0, 1]] = other


def test_setitem_raises_length():
    arr = PeriodArray(np.arange(3), freq="D")
    with tm.assert_raises_regex(ValueError, "length"):
        arr[[0, 1]] = [pd.Period("2000", freq="D")]


def test_setitem_raises_type():
    arr = PeriodArray(np.arange(3), freq="D")
    with tm.assert_raises_regex(TypeError, "int"):
        arr[0] = 1


# ----------------------------------------------------------------------------
# Ops

def tet_sub_period():
    arr = period_array(['2000', '2001'], freq='D')
    other = pd.Period("2000", freq="M")
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        arr - other
