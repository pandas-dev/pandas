import numpy as np
import pytest

from pandas._libs.tslibs import iNaT
from pandas._libs.tslibs.period import IncompatibleFrequency

from pandas.core.dtypes.dtypes import PeriodDtype, registry

import pandas as pd
from pandas.core.arrays import PeriodArray, period_array
import pandas.util.testing as tm

# ----------------------------------------------------------------------------
# Dtype


def test_registered():
    assert PeriodDtype in registry.dtypes
    result = registry.find("Period[D]")
    expected = PeriodDtype("D")
    assert result == expected

# ----------------------------------------------------------------------------
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


def test_period_array_readonly_object():
    # https://github.com/pandas-dev/pandas/issues/25403
    pa = period_array([pd.Period('2019-01-01')])
    arr = np.asarray(pa, dtype='object')
    arr.setflags(write=False)

    result = period_array(arr)
    tm.assert_period_array_equal(result, pa)

    result = pd.Series(arr)
    tm.assert_series_equal(result, pd.Series(pa))

    result = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, pd.DataFrame({"A": pa}))


def test_from_datetime64_freq_changes():
    # https://github.com/pandas-dev/pandas/issues/23438
    arr = pd.date_range("2017", periods=3, freq="D")
    result = PeriodArray._from_datetime64(arr, freq="M")
    expected = period_array(['2017-01-01', '2017-01-01', '2017-01-01'],
                            freq="M")
    tm.assert_period_array_equal(result, expected)


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
    with pytest.raises(IncompatibleFrequency, match=msg):
        period_array(data, freq)


def test_period_array_non_period_series_raies():
    ser = pd.Series([1, 2, 3])
    with pytest.raises(TypeError, match='dtype'):
        PeriodArray(ser, freq='D')


def test_period_array_freq_mismatch():
    arr = period_array(['2000', '2001'], freq='D')
    with pytest.raises(IncompatibleFrequency, match='freq'):
        PeriodArray(arr, freq='M')

    with pytest.raises(IncompatibleFrequency, match='freq'):
        PeriodArray(arr, freq=pd.tseries.offsets.MonthEnd())


def test_asi8():
    result = period_array(['2000', '2001', None], freq='D').asi8
    expected = np.array([10957, 11323, iNaT])
    tm.assert_numpy_array_equal(result, expected)


def test_take_raises():
    arr = period_array(['2000', '2001'], freq='D')
    with pytest.raises(IncompatibleFrequency, match='freq'):
        arr.take([0, -1], allow_fill=True,
                 fill_value=pd.Period('2000', freq='W'))

    with pytest.raises(ValueError, match='foo'):
        arr.take([0, -1], allow_fill=True, fill_value='foo')


@pytest.mark.parametrize('dtype', [
    int, np.int32, np.int64, 'uint32', 'uint64',
])
def test_astype(dtype):
    # We choose to ignore the sign and size of integers for
    # Period/Datetime/Timedelta astype
    arr = period_array(['2000', '2001', None], freq='D')
    result = arr.astype(dtype)

    if np.dtype(dtype).kind == 'u':
        expected_dtype = np.dtype('uint64')
    else:
        expected_dtype = np.dtype('int64')
    expected = arr.astype(expected_dtype)

    assert result.dtype == expected_dtype
    tm.assert_numpy_array_equal(result, expected)


def test_astype_copies():
    arr = period_array(['2000', '2001', None], freq='D')
    result = arr.astype(np.int64, copy=False)
    # Add the `.base`, since we now use `.asi8` which returns a view.
    # We could maybe override it in PeriodArray to return ._data directly.
    assert result.base is arr._data

    result = arr.astype(np.int64, copy=True)
    assert result is not arr._data
    tm.assert_numpy_array_equal(result, arr._data.view('i8'))


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
    with pytest.raises(TypeError, match=other[:-4]):
        arr.astype(other)


def test_fillna_raises():
    arr = period_array(['2000', '2001', '2002'], freq='D')
    with pytest.raises(ValueError, match='Length'):
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
    with pytest.raises(IncompatibleFrequency, match="freq"):
        arr[0] = pd.Period("2000", freq="A")

    other = period_array(['2000', '2001'], freq='A')
    with pytest.raises(IncompatibleFrequency, match="freq"):
        arr[[0, 1]] = other


def test_setitem_raises_length():
    arr = PeriodArray(np.arange(3), freq="D")
    with pytest.raises(ValueError, match="length"):
        arr[[0, 1]] = [pd.Period("2000", freq="D")]


def test_setitem_raises_type():
    arr = PeriodArray(np.arange(3), freq="D")
    with pytest.raises(TypeError, match="int"):
        arr[0] = 1


# ----------------------------------------------------------------------------
# Ops

def test_sub_period():
    arr = period_array(['2000', '2001'], freq='D')
    other = pd.Period("2000", freq="M")
    with pytest.raises(IncompatibleFrequency, match="freq"):
        arr - other


# ----------------------------------------------------------------------------
# Methods

@pytest.mark.parametrize('other', [
    pd.Period('2000', freq='H'),
    period_array(['2000', '2001', '2000'], freq='H')
])
def test_where_different_freq_raises(other):
    ser = pd.Series(period_array(['2000', '2001', '2002'], freq='D'))
    cond = np.array([True, False, True])
    with pytest.raises(IncompatibleFrequency, match="freq"):
        ser.where(cond, other)


# ----------------------------------------------------------------------------
# Printing

def test_repr_small():
    arr = period_array(['2000', '2001'], freq='D')
    result = str(arr)
    expected = (
        "<PeriodArray>\n"
        "['2000-01-01', '2001-01-01']\n"
        "Length: 2, dtype: period[D]"
    )
    assert result == expected


def test_repr_large():
    arr = period_array(['2000', '2001'] * 500, freq='D')
    result = str(arr)
    expected = (
        "<PeriodArray>\n"
        "['2000-01-01', '2001-01-01', '2000-01-01', '2001-01-01', "
        "'2000-01-01',\n"
        " '2001-01-01', '2000-01-01', '2001-01-01', '2000-01-01', "
        "'2001-01-01',\n"
        " ...\n"
        " '2000-01-01', '2001-01-01', '2000-01-01', '2001-01-01', "
        "'2000-01-01',\n"
        " '2001-01-01', '2000-01-01', '2001-01-01', '2000-01-01', "
        "'2001-01-01']\n"
        "Length: 1000, dtype: period[D]"
    )
    assert result == expected


# ----------------------------------------------------------------------------
# Reductions

class TestReductions(object):

    def test_min_max(self):
        arr = period_array([
            '2000-01-03',
            '2000-01-03',
            'NaT',
            '2000-01-02',
            '2000-01-05',
            '2000-01-04',
        ], freq='D')

        result = arr.min()
        expected = pd.Period('2000-01-02', freq='D')
        assert result == expected

        result = arr.max()
        expected = pd.Period('2000-01-05', freq='D')
        assert result == expected

        result = arr.min(skipna=False)
        assert result is pd.NaT

        result = arr.max(skipna=False)
        assert result is pd.NaT

    @pytest.mark.parametrize('skipna', [True, False])
    def test_min_max_empty(self, skipna):
        arr = period_array([], freq='D')
        result = arr.min(skipna=skipna)
        assert result is pd.NaT

        result = arr.max(skipna=skipna)
        assert result is pd.NaT
