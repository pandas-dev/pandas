# -*- coding: utf-8 -*-
import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

from pandas.core.arrays import (DatetimeArrayMixin,
                                TimedeltaArrayMixin,
                                PeriodArray)


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
    so here we just test that the DatetimeIndex behavior matches
    the DatetimeIndex behavior.
    """
    freqstr = request.param
    # TODO: non-monotone indexes; NaTs, different start dates, timezones
    pi = pd.date_range(start=pd.Timestamp('2000-01-01'),
                       periods=100,
                       freq=freqstr)
    return pi


class TestDatetimeArray(object):

    def test_from_dti(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArrayMixin(dti)
        assert list(dti) == list(arr)

        # Check that Index.__new__ knows what to do with DatetimeArray
        dti2 = pd.Index(arr)
        assert isinstance(dti2, pd.DatetimeIndex)
        assert list(dti2) == list(arr)

    def test_astype_object(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArrayMixin(dti)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(dti)

    @pytest.mark.parametrize('freqstr', ['D', 'B', 'W', 'M', 'Q', 'Y'])
    def test_to_period(self, datetime_index, freqstr):
        dti = datetime_index
        arr = DatetimeArrayMixin(dti)

        expected = dti.to_period(freq=freqstr)
        result = arr.to_period(freq=freqstr)
        assert isinstance(result, PeriodArray)

        # placeholder until these become actual EA subclasses and we can use
        #  an EA-specific tm.assert_ function
        tm.assert_index_equal(pd.Index(result), pd.Index(expected))

    @pytest.mark.parametrize('propname', pd.DatetimeIndex._bool_ops)
    def test_bool_properties(self, datetime_index, propname):
        # in this case _bool_ops is just `is_leap_year`
        dti = datetime_index
        arr = DatetimeArrayMixin(dti)
        assert dti.freq == arr.freq

        result = getattr(arr, propname)
        expected = np.array(getattr(dti, propname), dtype=result.dtype)

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('propname', pd.DatetimeIndex._field_ops)
    def test_int_properties(self, datetime_index, propname):
        dti = datetime_index
        arr = DatetimeArrayMixin(dti)

        result = getattr(arr, propname)
        expected = np.array(getattr(dti, propname), dtype=result.dtype)

        tm.assert_numpy_array_equal(result, expected)


class TestTimedeltaArray(object):
    def test_from_tdi(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArrayMixin(tdi)
        assert list(arr) == list(tdi)

        # Check that Index.__new__ knows what to do with TimedeltaArray
        tdi2 = pd.Index(arr)
        assert isinstance(tdi2, pd.TimedeltaIndex)
        assert list(tdi2) == list(arr)

    def test_astype_object(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArrayMixin(tdi)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(tdi)


class TestPeriodArray(object):

    def test_from_pi(self, period_index):
        pi = period_index
        arr = PeriodArray(pi.values)
        assert list(arr) == list(pi)

        # Check that Index.__new__ knows what to do with PeriodArray
        pi2 = pd.Index(arr)
        assert isinstance(pi2, pd.PeriodIndex)
        assert list(pi2) == list(arr)

    def test_astype_object(self, period_index):
        pi = period_index
        arr = PeriodArray(pi.values)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(pi)

    @pytest.mark.parametrize('how', ['S', 'E'])
    def test_to_timestamp(self, how, period_index):
        pi = period_index
        arr = PeriodArray(pi.values)

        expected = DatetimeArrayMixin(pi.to_timestamp(how=how))
        result = arr.to_timestamp(how=how)
        assert isinstance(result, DatetimeArrayMixin)

        # placeholder until these become actual EA subclasses and we can use
        #  an EA-specific tm.assert_ function
        tm.assert_index_equal(pd.Index(result), pd.Index(expected))

    @pytest.mark.parametrize('propname', pd.core.arrays.PeriodArray._bool_ops)
    def test_bool_properties(self, period_index, propname):
        # in this case _bool_ops is just `is_leap_year`
        pi = period_index
        arr = PeriodArray(pi.values)

        result = getattr(arr, propname)
        expected = np.array(getattr(pi, propname))

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('propname',
                             pd.core.arrays.PeriodArray._field_ops)
    def test_int_properties(self, period_index, propname):
        pi = period_index
        arr = PeriodArray(pi.values)

        result = getattr(arr, propname)
        expected = np.array(getattr(pi, propname))

        tm.assert_numpy_array_equal(result, expected)
