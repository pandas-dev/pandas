# coding=utf-8
import pytest
import pytz

from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from pandas import (Series, NaT, date_range, timedelta_range,
                    Timestamp, Timedelta)

import pandas.util.testing as tm


class TestDatetimeSeriesArithmetic(object):
    @pytest.mark.parametrize(
        'box, assert_func',
        [(Series, tm.assert_series_equal),
         (pd.Index, tm.assert_index_equal)])
    def test_sub_datetime64_not_ns(self, box, assert_func):
        # GH#7996
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'

        obj = box(date_range('20130101', periods=3))
        res = obj - dt64
        expected = box([Timedelta(days=0), Timedelta(days=1),
                        Timedelta(days=2)])
        assert_func(res, expected)

        res = dt64 - obj
        assert_func(res, -expected)

    def test_operators_datetimelike(self):
        def run_ops(ops, get_ser, test_ser):

            # check that we are getting a TypeError
            # with 'operate' (from core/ops.py) for the ops that are not
            # defined
            for op_str in ops:
                op = getattr(get_ser, op_str, None)
                with tm.assert_raises_regex(TypeError, 'operate|cannot'):
                    op(test_ser)

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        ops = ['__add__', '__mul__', '__floordiv__', '__truediv__', '__div__',
               '__pow__', '__radd__', '__rmul__', '__rfloordiv__',
               '__rtruediv__', '__rdiv__', '__rpow__']
        run_ops(ops, dt1, dt2)
        dt1 - dt2
        dt2 - dt1

        # ## datetime64 with timetimedelta ###
        ops = ['__mul__', '__floordiv__', '__truediv__', '__div__', '__pow__',
               '__rmul__', '__rfloordiv__', '__rtruediv__', '__rdiv__',
               '__rpow__']
        run_ops(ops, dt1, td1)
        dt1 + td1
        td1 + dt1
        dt1 - td1
        # TODO: Decide if this ought to work.
        # td1 - dt1

        # ## timetimedelta with datetime64 ###
        ops = ['__sub__', '__mul__', '__floordiv__', '__truediv__', '__div__',
               '__pow__', '__rmul__', '__rfloordiv__', '__rtruediv__',
               '__rdiv__', '__rpow__']
        run_ops(ops, td1, dt1)
        td1 + dt1
        dt1 + td1

        # 8260, 10763
        # datetime64 with tz
        ops = ['__mul__', '__floordiv__', '__truediv__', '__div__', '__pow__',
               '__rmul__', '__rfloordiv__', '__rtruediv__', '__rdiv__',
               '__rpow__']

        tz = 'US/Eastern'
        dt1 = Series(date_range('2000-01-01 09:00:00', periods=5,
                                tz=tz), name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan
        td1 = Series(timedelta_range('1 days 1 min', periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan
        run_ops(ops, dt1, td1)

        result = dt1 + td1[0]
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 + td2[0]
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        # odd numpy behavior with scalar timedeltas
        result = td1[0] + dt1
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = td2[0] + dt2
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt1 - td1[0]
        exp = (dt1.dt.tz_localize(None) - td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)
        pytest.raises(TypeError, lambda: td1[0] - dt1)

        result = dt2 - td2[0]
        exp = (dt2.dt.tz_localize(None) - td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)
        pytest.raises(TypeError, lambda: td2[0] - dt2)

        result = dt1 + td1
        exp = (dt1.dt.tz_localize(None) + td1).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 + td2
        exp = (dt2.dt.tz_localize(None) + td2).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt1 - td1
        exp = (dt1.dt.tz_localize(None) - td1).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 - td2
        exp = (dt2.dt.tz_localize(None) - td2).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        pytest.raises(TypeError, lambda: td1 - dt1)
        pytest.raises(TypeError, lambda: td2 - dt2)

    def test_sub_single_tz(self):
        # GH12290
        s1 = Series([pd.Timestamp('2016-02-10', tz='America/Sao_Paulo')])
        s2 = Series([pd.Timestamp('2016-02-08', tz='America/Sao_Paulo')])
        result = s1 - s2
        expected = Series([Timedelta('2days')])
        tm.assert_series_equal(result, expected)
        result = s2 - s1
        expected = Series([Timedelta('-2days')])
        tm.assert_series_equal(result, expected)

    def test_dt64tz_series_sub_dtitz(self):
        # GH#19071 subtracting tzaware DatetimeIndex from tzaware Series
        # (with same tz) raises, fixed by #19024
        dti = pd.date_range('1999-09-30', periods=10, tz='US/Pacific')
        ser = pd.Series(dti)
        expected = pd.Series(pd.TimedeltaIndex(['0days'] * 10))

        res = dti - ser
        tm.assert_series_equal(res, expected)
        res = ser - dti
        tm.assert_series_equal(res, expected)

    def test_sub_datetime_compat(self):
        # see gh-14088
        s = Series([datetime(2016, 8, 23, 12, tzinfo=pytz.utc), pd.NaT])
        dt = datetime(2016, 8, 22, 12, tzinfo=pytz.utc)
        exp = Series([Timedelta('1 days'), pd.NaT])
        tm.assert_series_equal(s - dt, exp)
        tm.assert_series_equal(s - Timestamp(dt), exp)

    def test_dt64_series_with_timedelta(self):
        # scalar timedeltas/np.timedelta64 objects
        # operate with np.timedelta64 correctly
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + np.timedelta64(1, 's')
        result2 = np.timedelta64(1, 's') + s
        expected = Series([Timestamp('20130101 9:01:01'),
                           Timestamp('20130101 9:02:01')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        result = s + np.timedelta64(5, 'ms')
        result2 = np.timedelta64(5, 'ms') + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

    def test_dt64_series_add_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:01:05'),
                           Timestamp('20130101 9:02:05')])

        result = ser + pd.offsets.Second(5)
        tm.assert_series_equal(result, expected)

        result2 = pd.offsets.Second(5) + ser
        tm.assert_series_equal(result2, expected)

    def test_dt64_series_sub_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:00:55'),
                           Timestamp('20130101 9:01:55')])

        result = ser - pd.offsets.Second(5)
        tm.assert_series_equal(result, expected)

        result2 = -pd.offsets.Second(5) + ser
        tm.assert_series_equal(result2, expected)

        with pytest.raises(TypeError):
            pd.offsets.Second(5) - ser

    @pytest.mark.parametrize('cls_name', ['Day', 'Hour', 'Minute', 'Second',
                                          'Milli', 'Micro', 'Nano'])
    def test_dt64_series_with_tick_DateOffset_smoke(self, cls_name):
        # GH 4532
        # smoke tests for valid DateOffsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        offset_cls = getattr(pd.offsets, cls_name)
        ser + offset_cls(5)
        offset_cls(5) + ser

    def test_dt64_series_add_mixed_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + pd.offsets.Milli(5)
        result2 = pd.offsets.Milli(5) + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        result = s + pd.offsets.Minute(5) + pd.offsets.Milli(5)
        expected = Series([Timestamp('20130101 9:06:00.005'),
                           Timestamp('20130101 9:07:00.005')])
        tm.assert_series_equal(result, expected)

    def test_dt64_series_sub_NaT(self):
        # GH#18808
        dti = pd.DatetimeIndex([pd.NaT, pd.Timestamp('19900315')])
        ser = pd.Series(dti)
        res = ser - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

        dti_tz = dti.tz_localize('Asia/Tokyo')
        ser_tz = pd.Series(dti_tz)
        res = ser_tz - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

    def test_datetime64_ops_nat(self):
        # GH 11349
        datetime_series = Series([NaT, Timestamp('19900315')])
        nat_series_dtype_timestamp = Series([NaT, NaT], dtype='datetime64[ns]')
        single_nat_dtype_datetime = Series([NaT], dtype='datetime64[ns]')

        # subtraction
        tm.assert_series_equal(-NaT + datetime_series,
                               nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + datetime_series

        tm.assert_series_equal(-NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + nat_series_dtype_timestamp

        # addition
        tm.assert_series_equal(nat_series_dtype_timestamp + NaT,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

        tm.assert_series_equal(nat_series_dtype_timestamp + NaT,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

    @pytest.mark.parametrize('dt64_series', [
        Series([Timestamp('19900315'), Timestamp('19900315')]),
        Series([NaT, Timestamp('19900315')]),
        Series([NaT, NaT], dtype='datetime64[ns]')])
    @pytest.mark.parametrize('one', [1, 1.0, np.array(1)])
    def test_dt64_mul_div_numeric_invalid(self, one, dt64_series):
        # multiplication
        with pytest.raises(TypeError):
            dt64_series * one
        with pytest.raises(TypeError):
            one * dt64_series

        # division
        with pytest.raises(TypeError):
            dt64_series / one
        with pytest.raises(TypeError):
            one / dt64_series

    def test_dt64_series_arith_overflow(self):
        # GH#12534, fixed by #19024
        dt = pd.Timestamp('1700-01-31')
        td = pd.Timedelta('20000 Days')
        dti = pd.date_range('1949-09-30', freq='100Y', periods=4)
        ser = pd.Series(dti)
        with pytest.raises(OverflowError):
            ser - dt
        with pytest.raises(OverflowError):
            dt - ser
        with pytest.raises(OverflowError):
            ser + td
        with pytest.raises(OverflowError):
            td + ser

        ser.iloc[-1] = pd.NaT
        expected = pd.Series(['2004-10-03', '2104-10-04', '2204-10-04', 'NaT'],
                             dtype='datetime64[ns]')
        res = ser + td
        tm.assert_series_equal(res, expected)
        res = td + ser
        tm.assert_series_equal(res, expected)

        ser.iloc[1:] = pd.NaT
        expected = pd.Series(['91279 Days', 'NaT', 'NaT', 'NaT'],
                             dtype='timedelta64[ns]')
        res = ser - dt
        tm.assert_series_equal(res, expected)
        res = dt - ser
        tm.assert_series_equal(res, -expected)
