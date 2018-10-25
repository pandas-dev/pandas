"""
Tests for pd.NaT behavior that depend on the implementations of Index/Series
"""
import pytest
import numpy as np

import pandas.util.testing as tm
from pandas import (
    Timestamp, Timedelta, Period, NaT,
    DatetimeIndex, TimedeltaIndex, PeriodIndex, Index, Series
)


@pytest.mark.parametrize('nat, idx', [(Timestamp('NaT'), DatetimeIndex),
                                      (Timedelta('NaT'), TimedeltaIndex),
                                      (Period('NaT', freq='M'), PeriodIndex)])
def test_nat_fields(nat, idx):
    for field in idx._field_ops:

        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == 'weekday':
            continue

        result = getattr(NaT, field)
        assert np.isnan(result)

        result = getattr(nat, field)
        assert np.isnan(result)

    for field in idx._bool_ops:

        result = getattr(NaT, field)
        assert result is False

        result = getattr(nat, field)
        assert result is False


def test_nat_vector_field_access():
    idx = DatetimeIndex(['1/1/2000', None, None, '1/4/2000'])

    for field in DatetimeIndex._field_ops:
        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == 'weekday':
            continue

        result = getattr(idx, field)
        expected = Index([getattr(x, field) for x in idx])
        tm.assert_index_equal(result, expected)

    ser = Series(idx)

    for field in DatetimeIndex._field_ops:

        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == 'weekday':
            continue

        result = getattr(ser.dt, field)
        expected = [getattr(x, field) for x in idx]
        tm.assert_series_equal(result, Series(expected))

    for field in DatetimeIndex._bool_ops:
        result = getattr(ser.dt, field)
        expected = [getattr(x, field) for x in idx]
        tm.assert_series_equal(result, Series(expected))


def test_nat_arithmetic_index():
    # GH#11718

    dti = DatetimeIndex(['2011-01-01', '2011-01-02'], name='x')
    exp = DatetimeIndex([NaT, NaT], name='x')
    tm.assert_index_equal(dti + NaT, exp)
    tm.assert_index_equal(NaT + dti, exp)

    dti_tz = DatetimeIndex(['2011-01-01', '2011-01-02'],
                           tz='US/Eastern', name='x')
    exp = DatetimeIndex([NaT, NaT], name='x', tz='US/Eastern')
    tm.assert_index_equal(dti_tz + NaT, exp)
    tm.assert_index_equal(NaT + dti_tz, exp)

    exp = TimedeltaIndex([NaT, NaT], name='x')
    for (left, right) in [(NaT, dti), (NaT, dti_tz)]:
        tm.assert_index_equal(left - right, exp)
        tm.assert_index_equal(right - left, exp)

    # timedelta # GH#19124
    tdi = TimedeltaIndex(['1 day', '2 day'], name='x')
    tdi_nat = TimedeltaIndex([NaT, NaT], name='x')

    tm.assert_index_equal(tdi + NaT, tdi_nat)
    tm.assert_index_equal(NaT + tdi, tdi_nat)
    tm.assert_index_equal(tdi - NaT, tdi_nat)
    tm.assert_index_equal(NaT - tdi, tdi_nat)


@pytest.mark.parametrize('box', [TimedeltaIndex, Series])
def test_nat_arithmetic_td64_vector(box):
    # GH#19124
    vec = box(['1 day', '2 day'], dtype='timedelta64[ns]')
    box_nat = box([NaT, NaT], dtype='timedelta64[ns]')

    tm.assert_equal(vec + NaT, box_nat)
    tm.assert_equal(NaT + vec, box_nat)
    tm.assert_equal(vec - NaT, box_nat)
    tm.assert_equal(NaT - vec, box_nat)
