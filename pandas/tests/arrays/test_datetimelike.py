# -*- coding: utf-8 -*-
import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import (
    DatetimeArrayMixin as DatetimeArray, PeriodArray,
    TimedeltaArrayMixin as TimedeltaArray)
import pandas.util.testing as tm


# TODO: more freq variants
@pytest.fixture(params=['D', 'B', 'W', 'M', 'Q', 'Y'])
def period_index(request):
    """
    A fixture to provide PeriodIndex objects with different frequencies.

    Most PeriodArray behavior is already tested in PeriodIndex tests,
    so here we just test that the PeriodArray behavior matches
    the PeriodIndex behavior.
    """
    freqstr = request.param
    # TODO: non-monotone indexes; NaTs, different start dates
    pi = pd.period_range(start=pd.Timestamp('2000-01-01'),
                         periods=100,
                         freq=freqstr)
    return pi


@pytest.fixture(params=['D', 'B', 'W', 'M', 'Q', 'Y'])
def datetime_index(request):
    """
    A fixture to provide DatetimeIndex objects with different frequencies.

    Most DatetimeArray behavior is already tested in DatetimeIndex tests,
    so here we just test that the DatetimeArray behavior matches
    the DatetimeIndex behavior.
    """
    freqstr = request.param
    # TODO: non-monotone indexes; NaTs, different start dates, timezones
    pi = pd.date_range(start=pd.Timestamp('2000-01-01'),
                       periods=100,
                       freq=freqstr)
    return pi


@pytest.fixture
def timedelta_index(request):
    """
    A fixture to provide TimedeltaIndex objects with different frequencies.
     Most TimedeltaArray behavior is already tested in TimedeltaIndex tests,
    so here we just test that the TimedeltaArray behavior matches
    the TimedeltaIndex behavior.
    """
    # TODO: flesh this out
    return pd.TimedeltaIndex(['1 Day', '3 Hours', 'NaT'])


class SharedTests(object):
    index_cls = None

    def test_take(self):
        data = np.arange(100, dtype='i8')
        np.random.shuffle(data)

        idx = self.index_cls._simple_new(data, freq='D')
        arr = self.array_cls(idx)

        takers = [1, 4, 94]
        result = arr.take(takers)
        expected = idx.take(takers)

        tm.assert_index_equal(self.index_cls(result), expected)

        takers = np.array([1, 4, 94])
        result = arr.take(takers)
        expected = idx.take(takers)

        tm.assert_index_equal(self.index_cls(result), expected)

    def test_take_fill(self):
        data = np.arange(10, dtype='i8')

        idx = self.index_cls._simple_new(data, freq='D')
        arr = self.array_cls(idx)

        result = arr.take([-1, 1], allow_fill=True, fill_value=None)
        assert result[0] is pd.NaT

        result = arr.take([-1, 1], allow_fill=True, fill_value=np.nan)
        assert result[0] is pd.NaT

        result = arr.take([-1, 1], allow_fill=True, fill_value=pd.NaT)
        assert result[0] is pd.NaT

        with pytest.raises(ValueError):
            arr.take([0, 1], allow_fill=True, fill_value=2)

        with pytest.raises(ValueError):
            arr.take([0, 1], allow_fill=True, fill_value=2.0)

        with pytest.raises(ValueError):
            arr.take([0, 1], allow_fill=True,
                     fill_value=pd.Timestamp.now().time)

    def test_concat_same_type(self):
        data = np.arange(10, dtype='i8')

        idx = self.index_cls._simple_new(data, freq='D').insert(0, pd.NaT)
        arr = self.array_cls(idx)

        result = arr._concat_same_type([arr[:-1], arr[1:], arr])
        expected = idx._concat_same_dtype([idx[:-1], idx[1:], idx], None)

        tm.assert_index_equal(self.index_cls(result), expected)


class TestDatetimeArray(SharedTests):
    index_cls = pd.DatetimeIndex
    array_cls = DatetimeArray

    def test_array_object_dtype(self, tz_naive_fixture):
        # GH#23524
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArray(dti)

        expected = np.array(list(dti))

        result = np.array(arr, dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        # also test the DatetimeIndex method while we're at it
        result = np.array(dti, dtype=object)
        tm.assert_numpy_array_equal(result, expected)

    def test_array(self, tz_naive_fixture):
        # GH#23524
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArray(dti)

        expected = dti.asi8.view('M8[ns]')
        result = np.array(arr)
        tm.assert_numpy_array_equal(result, expected)

        # check that we are not making copies when setting copy=False
        result = np.array(arr, copy=False)
        assert result.base is expected.base
        assert result.base is not None

    def test_array_i8_dtype(self, tz_naive_fixture):
        # GH#23524
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArray(dti)

        expected = dti.asi8
        result = np.array(arr, dtype='i8')
        tm.assert_numpy_array_equal(result, expected)

        result = np.array(arr, dtype=np.int64)
        tm.assert_numpy_array_equal(result, expected)

        # check that we are not making copies when setting copy=False
        result = np.array(arr, dtype='i8', copy=False)
        assert result.base is expected.base
        assert result.base is not None

    def test_from_array_keeps_base(self):
        # Ensure that DatetimeArray._data.base isn't lost.
        arr = np.array(['2000-01-01', '2000-01-02'], dtype='M8[ns]')
        dta = DatetimeArray(arr)

        assert dta._data is arr
        dta = DatetimeArray(arr[:0])
        assert dta._data.base is arr

    def test_from_dti(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArray(dti)
        assert list(dti) == list(arr)

        # Check that Index.__new__ knows what to do with DatetimeArray
        dti2 = pd.Index(arr)
        assert isinstance(dti2, pd.DatetimeIndex)
        assert list(dti2) == list(arr)

    def test_astype_object(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArray(dti)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(dti)

    @pytest.mark.parametrize('freqstr', ['D', 'B', 'W', 'M', 'Q', 'Y'])
    def test_to_perioddelta(self, datetime_index, freqstr):
        # GH#23113
        dti = datetime_index
        arr = DatetimeArray(dti)

        expected = dti.to_perioddelta(freq=freqstr)
        result = arr.to_perioddelta(freq=freqstr)
        assert isinstance(result, TimedeltaArray)

        # placeholder until these become actual EA subclasses and we can use
        #  an EA-specific tm.assert_ function
        tm.assert_index_equal(pd.Index(result), pd.Index(expected))

    @pytest.mark.parametrize('freqstr', ['D', 'B', 'W', 'M', 'Q', 'Y'])
    def test_to_period(self, datetime_index, freqstr):
        dti = datetime_index
        arr = DatetimeArray(dti)

        expected = dti.to_period(freq=freqstr)
        result = arr.to_period(freq=freqstr)
        assert isinstance(result, PeriodArray)

        # placeholder until these become actual EA subclasses and we can use
        #  an EA-specific tm.assert_ function
        tm.assert_index_equal(pd.Index(result), pd.Index(expected))

    @pytest.mark.parametrize('propname',
                             pd.core.arrays.DatetimeArrayMixin._bool_ops)
    def test_bool_properties(self, datetime_index, propname):
        # in this case _bool_ops is just `is_leap_year`
        dti = datetime_index
        arr = DatetimeArray(dti)
        assert dti.freq == arr.freq

        result = getattr(arr, propname)
        expected = np.array(getattr(dti, propname), dtype=result.dtype)

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('propname',
                             pd.core.arrays.DatetimeArrayMixin._field_ops)
    def test_int_properties(self, datetime_index, propname):
        dti = datetime_index
        arr = DatetimeArray(dti)

        result = getattr(arr, propname)
        expected = np.array(getattr(dti, propname), dtype=result.dtype)

        tm.assert_numpy_array_equal(result, expected)

    def test_take_fill_valid(self, datetime_index, tz_naive_fixture):
        dti = datetime_index.tz_localize(tz_naive_fixture)
        arr = DatetimeArray(dti)

        now = pd.Timestamp.now().tz_localize(dti.tz)
        result = arr.take([-1, 1], allow_fill=True, fill_value=now)
        assert result[0] == now

        with pytest.raises(ValueError):
            # fill_value Timedelta invalid
            arr.take([-1, 1], allow_fill=True, fill_value=now - now)

        with pytest.raises(ValueError):
            # fill_value Period invalid
            arr.take([-1, 1], allow_fill=True, fill_value=pd.Period('2014Q1'))

        tz = None if dti.tz is not None else 'US/Eastern'
        now = pd.Timestamp.now().tz_localize(tz)
        with pytest.raises(TypeError):
            # Timestamp with mismatched tz-awareness
            arr.take([-1, 1], allow_fill=True, fill_value=now)

    def test_concat_same_type_invalid(self, datetime_index):
        # different timezones
        dti = datetime_index
        arr = DatetimeArray(dti)

        if arr.tz is None:
            other = arr.tz_localize('UTC')
        else:
            other = arr.tz_localize(None)

        with pytest.raises(AssertionError):
            arr._concat_same_type([arr, other])

    def test_concat_same_type_different_freq(self):
        # we *can* concatentate DTI with different freqs.
        a = DatetimeArray(pd.date_range('2000', periods=2, freq='D',
                                        tz='US/Central'))
        b = DatetimeArray(pd.date_range('2000', periods=2, freq='H',
                                        tz='US/Central'))
        result = DatetimeArray._concat_same_type([a, b])
        expected = DatetimeArray(pd.to_datetime([
            '2000-01-01 00:00:00', '2000-01-02 00:00:00',
            '2000-01-01 00:00:00', '2000-01-01 01:00:00',
        ]).tz_localize("US/Central"))

        tm.assert_datetime_array_equal(result, expected)


@pytest.mark.skip(reason="TODO")
class TestTimedeltaArray(SharedTests):
    index_cls = pd.TimedeltaIndex
    array_cls = TimedeltaArray

    def test_from_tdi(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArray(tdi)
        assert list(arr) == list(tdi)

        # Check that Index.__new__ knows what to do with TimedeltaArray
        tdi2 = pd.Index(arr)
        assert isinstance(tdi2, pd.TimedeltaIndex)
        assert list(tdi2) == list(arr)

    def test_astype_object(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArray(tdi)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(tdi)

    def test_to_pytimedelta(self, timedelta_index):
        tdi = timedelta_index
        arr = TimedeltaArray(tdi)

        expected = tdi.to_pytimedelta()
        result = arr.to_pytimedelta()

        tm.assert_numpy_array_equal(result, expected)

    def test_total_seconds(self, timedelta_index):
        tdi = timedelta_index
        arr = TimedeltaArray(tdi)

        expected = tdi.total_seconds()
        result = arr.total_seconds()

        tm.assert_numpy_array_equal(result, expected.values)

    @pytest.mark.parametrize('propname',
                             pd.core.arrays.TimedeltaArrayMixin._field_ops)
    def test_int_properties(self, timedelta_index, propname):
        tdi = timedelta_index
        arr = TimedeltaArray(tdi)

        result = getattr(arr, propname)
        expected = np.array(getattr(tdi, propname), dtype=result.dtype)

        tm.assert_numpy_array_equal(result, expected)

    def test_take_fill_valid(self, timedelta_index):
        tdi = timedelta_index
        arr = TimedeltaArray(tdi)

        td1 = pd.Timedelta(days=1)
        result = arr.take([-1, 1], allow_fill=True, fill_value=td1)
        assert result[0] == td1

        now = pd.Timestamp.now()
        with pytest.raises(ValueError):
            # fill_value Timestamp invalid
            arr.take([0, 1], allow_fill=True, fill_value=now)

        with pytest.raises(ValueError):
            # fill_value Period invalid
            arr.take([0, 1], allow_fill=True, fill_value=now.to_period('D'))


class TestPeriodArray(SharedTests):
    index_cls = pd.PeriodIndex
    array_cls = PeriodArray

    def test_from_pi(self, period_index):
        pi = period_index
        arr = PeriodArray(pi)
        assert list(arr) == list(pi)

        # Check that Index.__new__ knows what to do with PeriodArray
        pi2 = pd.Index(arr)
        assert isinstance(pi2, pd.PeriodIndex)
        assert list(pi2) == list(arr)

    def test_astype_object(self, period_index):
        pi = period_index
        arr = PeriodArray(pi)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(pi)

    @pytest.mark.parametrize('how', ['S', 'E'])
    def test_to_timestamp(self, how, period_index):
        pi = period_index
        arr = PeriodArray(pi)

        expected = DatetimeArray(pi.to_timestamp(how=how))
        result = arr.to_timestamp(how=how)
        assert isinstance(result, DatetimeArray)

        # placeholder until these become actual EA subclasses and we can use
        #  an EA-specific tm.assert_ function
        tm.assert_index_equal(pd.Index(result), pd.Index(expected))

    @pytest.mark.parametrize('propname', PeriodArray._bool_ops)
    def test_bool_properties(self, period_index, propname):
        # in this case _bool_ops is just `is_leap_year`
        pi = period_index
        arr = PeriodArray(pi)

        result = getattr(arr, propname)
        expected = np.array(getattr(pi, propname))

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('propname', PeriodArray._field_ops)
    def test_int_properties(self, period_index, propname):
        pi = period_index
        arr = PeriodArray(pi)

        result = getattr(arr, propname)
        expected = np.array(getattr(pi, propname))

        tm.assert_numpy_array_equal(result, expected)
