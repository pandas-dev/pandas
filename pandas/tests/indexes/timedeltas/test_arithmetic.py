# -*- coding: utf-8 -*-
import pytest
import numpy as np
from datetime import timedelta
from distutils.version import LooseVersion

import pandas as pd
import pandas.util.testing as tm
from pandas import (DatetimeIndex, TimedeltaIndex, Float64Index, Int64Index,
                    to_timedelta, timedelta_range, date_range,
                    Series,
                    Timestamp, Timedelta)


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=str)
def delta(request):
    # Several ways of representing two hours
    return request.param


@pytest.fixture(params=['B', 'D'])
def freq(request):
    return request.param


class TestTimedeltaIndexArithmetic(object):
    _holder = TimedeltaIndex

    @pytest.mark.xfail(reason='GH#18824 ufunc add cannot use operands...')
    def test_tdi_with_offset_array(self):
        # GH#18849
        tdi = pd.TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'])
        offs = np.array([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)])
        expected = pd.TimedeltaIndex(['1 days 01:00:00', '3 days 04:02:00'])

        res = tdi + offs
        tm.assert_index_equal(res, expected)

        res2 = offs + tdi
        tm.assert_index_equal(res2, expected)

        anchored = np.array([pd.offsets.QuarterEnd(),
                             pd.offsets.Week(weekday=2)])
        with pytest.raises(TypeError):
            tdi + anchored

    # TODO: Split by ops, better name
    def test_numeric_compat(self):
        idx = self._holder(np.arange(5, dtype='int64'))
        didx = self._holder(np.arange(5, dtype='int64') ** 2)
        result = idx * 1
        tm.assert_index_equal(result, idx)

        result = 1 * idx
        tm.assert_index_equal(result, idx)

        result = idx / 1
        tm.assert_index_equal(result, idx)

        result = idx // 1
        tm.assert_index_equal(result, idx)

        result = idx * np.array(5, dtype='int64')
        tm.assert_index_equal(result,
                              self._holder(np.arange(5, dtype='int64') * 5))

        result = idx * np.arange(5, dtype='int64')
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='int64'))
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='float64') + 0.1)
        tm.assert_index_equal(result, self._holder(np.arange(
            5, dtype='float64') * (np.arange(5, dtype='float64') + 0.1)))

        # invalid
        pytest.raises(TypeError, lambda: idx * idx)
        pytest.raises(ValueError, lambda: idx * self._holder(np.arange(3)))
        pytest.raises(ValueError, lambda: idx * np.array([1, 2]))

    def test_ufunc_coercions(self):
        # normal ops are also tested in tseries/test_timedeltas.py
        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                             freq='2H', name='x')

        for result in [idx * 2, np.multiply(idx, 2)]:
            assert isinstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['4H', '8H', '12H', '16H', '20H'],
                                 freq='4H', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '4H'

        for result in [idx / 2, np.divide(idx, 2)]:
            assert isinstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['1H', '2H', '3H', '4H', '5H'],
                                 freq='H', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == 'H'

        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                             freq='2H', name='x')
        for result in [-idx, np.negative(idx)]:
            assert isinstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['-2H', '-4H', '-6H', '-8H', '-10H'],
                                 freq='-2H', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '-2H'

        idx = TimedeltaIndex(['-2H', '-1H', '0H', '1H', '2H'],
                             freq='H', name='x')
        for result in [abs(idx), np.absolute(idx)]:
            assert isinstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['2H', '1H', '0H', '1H', '2H'],
                                 freq=None, name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq is None

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and integer

    def test_tdi_add_int(self, one):
        # Variants of `one` for #19012
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng + one
        expected = timedelta_range('1 days 10:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)

    def test_tdi_iadd_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        expected = timedelta_range('1 days 10:00:00', freq='H', periods=10)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_tdi_sub_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng - one
        expected = timedelta_range('1 days 08:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)

    def test_tdi_isub_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        expected = timedelta_range('1 days 08:00:00', freq='H', periods=10)
        rng -= one
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and timedelta-like

    def test_tdi_add_timedeltalike(self, delta):
        # only test adding/sub offsets as + is now numeric
        rng = timedelta_range('1 days', '10 days')
        result = rng + delta
        expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                   freq='D')
        tm.assert_index_equal(result, expected)

    def test_tdi_iadd_timedeltalike(self, delta):
        # only test adding/sub offsets as + is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                   freq='D')
        rng += delta
        tm.assert_index_equal(rng, expected)

    def test_tdi_sub_timedeltalike(self, delta):
        # only test adding/sub offsets as - is now numeric
        rng = timedelta_range('1 days', '10 days')
        result = rng - delta
        expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')
        tm.assert_index_equal(result, expected)

    def test_tdi_isub_timedeltalike(self, delta):
        # only test adding/sub offsets as - is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')
        rng -= delta
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and datetime-like

    def test_tdi_sub_timestamp_raises(self):
        idx = TimedeltaIndex(['1 day', '2 day'])
        msg = "cannot subtract a datelike from a TimedeltaIndex"
        with tm.assert_raises_regex(TypeError, msg):
            idx - Timestamp('2011-01-01')

    def test_tdi_add_timestamp(self):
        idx = TimedeltaIndex(['1 day', '2 day'])

        result = idx + Timestamp('2011-01-01')
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])
        tm.assert_index_equal(result, expected)

    def test_tdi_radd_timestamp(self):
        idx = TimedeltaIndex(['1 day', '2 day'])

        result = Timestamp('2011-01-01') + idx
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])
        tm.assert_index_equal(result, expected)

    # -------------------------------------------------------------

    # TODO: Split by operation, better name
    def test_ops_compat(self):

        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), Timedelta(hours=2)]

        rng = timedelta_range('1 days', '10 days', name='foo')

        # multiply
        for offset in offsets:
            pytest.raises(TypeError, lambda: rng * offset)

        # divide
        expected = Int64Index((np.arange(10) + 1) * 12, name='foo')
        for offset in offsets:
            result = rng / offset
            tm.assert_index_equal(result, expected, exact=False)

        # floor divide
        expected = Int64Index((np.arange(10) + 1) * 12, name='foo')
        for offset in offsets:
            result = rng // offset
            tm.assert_index_equal(result, expected, exact=False)

        # divide with nats
        rng = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        expected = Float64Index([12, np.nan, 24], name='foo')
        for offset in offsets:
            result = rng / offset
            tm.assert_index_equal(result, expected)

        # don't allow division by NaT (make could in the future)
        pytest.raises(TypeError, lambda: rng / pd.NaT)

    def test_subtraction_ops(self):
        # with datetimes/timedelta and tdi/dti
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')
        td = Timedelta('1 days')
        dt = Timestamp('20130101')

        pytest.raises(TypeError, lambda: tdi - dt)
        pytest.raises(TypeError, lambda: tdi - dti)
        pytest.raises(TypeError, lambda: td - dt)
        pytest.raises(TypeError, lambda: td - dti)

        result = dt - dti
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'], name='bar')
        tm.assert_index_equal(result, expected)

        result = dti - dt
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'], name='bar')
        tm.assert_index_equal(result, expected)

        result = tdi - td
        expected = TimedeltaIndex(['0 days', pd.NaT, '1 days'], name='foo')
        tm.assert_index_equal(result, expected, check_names=False)

        result = td - tdi
        expected = TimedeltaIndex(['0 days', pd.NaT, '-1 days'], name='foo')
        tm.assert_index_equal(result, expected, check_names=False)

        result = dti - td
        expected = DatetimeIndex(
            ['20121231', '20130101', '20130102'], name='bar')
        tm.assert_index_equal(result, expected, check_names=False)

        result = dt - tdi
        expected = DatetimeIndex(['20121231', pd.NaT, '20121230'], name='foo')
        tm.assert_index_equal(result, expected)

    def test_subtraction_ops_with_tz(self):

        # check that dt/dti subtraction ops with tz are validated
        dti = date_range('20130101', periods=3)
        ts = Timestamp('20130101')
        dt = ts.to_pydatetime()
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')
        ts_tz = Timestamp('20130101').tz_localize('US/Eastern')
        ts_tz2 = Timestamp('20130101').tz_localize('CET')
        dt_tz = ts_tz.to_pydatetime()
        td = Timedelta('1 days')

        def _check(result, expected):
            assert result == expected
            assert isinstance(result, Timedelta)

        # scalars
        result = ts - ts
        expected = Timedelta('0 days')
        _check(result, expected)

        result = dt_tz - ts_tz
        expected = Timedelta('0 days')
        _check(result, expected)

        result = ts_tz - dt_tz
        expected = Timedelta('0 days')
        _check(result, expected)

        # tz mismatches
        pytest.raises(TypeError, lambda: dt_tz - ts)
        pytest.raises(TypeError, lambda: dt_tz - dt)
        pytest.raises(TypeError, lambda: dt_tz - ts_tz2)
        pytest.raises(TypeError, lambda: dt - dt_tz)
        pytest.raises(TypeError, lambda: ts - dt_tz)
        pytest.raises(TypeError, lambda: ts_tz2 - ts)
        pytest.raises(TypeError, lambda: ts_tz2 - dt)
        pytest.raises(TypeError, lambda: ts_tz - ts_tz2)

        # with dti
        pytest.raises(TypeError, lambda: dti - ts_tz)
        pytest.raises(TypeError, lambda: dti_tz - ts)
        pytest.raises(TypeError, lambda: dti_tz - ts_tz2)

        result = dti_tz - dt_tz
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'])
        tm.assert_index_equal(result, expected)

        result = dt_tz - dti_tz
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'])
        tm.assert_index_equal(result, expected)

        result = dti_tz - ts_tz
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'])
        tm.assert_index_equal(result, expected)

        result = ts_tz - dti_tz
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'])
        tm.assert_index_equal(result, expected)

        result = td - td
        expected = Timedelta('0 days')
        _check(result, expected)

        result = dti_tz - td
        expected = DatetimeIndex(
            ['20121231', '20130101', '20130102'], tz='US/Eastern')
        tm.assert_index_equal(result, expected)

    def test_dti_tdi_numeric_ops(self):
        # These are normally union/diff set-like ops
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')

        # TODO(wesm): unused?
        # td = Timedelta('1 days')
        # dt = Timestamp('20130101')

        result = tdi - tdi
        expected = TimedeltaIndex(['0 days', pd.NaT, '0 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = tdi + tdi
        expected = TimedeltaIndex(['2 days', pd.NaT, '4 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = dti - tdi  # name will be reset
        expected = DatetimeIndex(['20121231', pd.NaT, '20130101'])
        tm.assert_index_equal(result, expected)

    def test_sub_period(self):
        # GH 13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        for freq in [None, 'H']:
            idx = pd.TimedeltaIndex(['1 hours', '2 hours'], freq=freq)

            with pytest.raises(TypeError):
                idx - p

            with pytest.raises(TypeError):
                p - idx

    def test_addition_ops(self):
        # with datetimes/timedelta and tdi/dti
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')
        td = Timedelta('1 days')
        dt = Timestamp('20130101')

        result = tdi + dt
        expected = DatetimeIndex(['20130102', pd.NaT, '20130103'], name='foo')
        tm.assert_index_equal(result, expected)

        result = dt + tdi
        expected = DatetimeIndex(['20130102', pd.NaT, '20130103'], name='foo')
        tm.assert_index_equal(result, expected)

        result = td + tdi
        expected = TimedeltaIndex(['2 days', pd.NaT, '3 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = tdi + td
        expected = TimedeltaIndex(['2 days', pd.NaT, '3 days'], name='foo')
        tm.assert_index_equal(result, expected)

        # unequal length
        pytest.raises(ValueError, lambda: tdi + dti[0:1])
        pytest.raises(ValueError, lambda: tdi[0:1] + dti)

        # random indexes
        pytest.raises(TypeError, lambda: tdi + Int64Index([1, 2, 3]))

        # this is a union!
        # pytest.raises(TypeError, lambda : Int64Index([1,2,3]) + tdi)

        result = tdi + dti  # name will be reset
        expected = DatetimeIndex(['20130102', pd.NaT, '20130105'])
        tm.assert_index_equal(result, expected)

        result = dti + tdi  # name will be reset
        expected = DatetimeIndex(['20130102', pd.NaT, '20130105'])
        tm.assert_index_equal(result, expected)

        result = dt + td
        expected = Timestamp('20130102')
        assert result == expected

        result = td + dt
        expected = Timestamp('20130102')
        assert result == expected

    def test_ops_ndarray(self):
        td = Timedelta('1 day')

        # timedelta, timedelta
        other = pd.to_timedelta(['1 day']).values
        expected = pd.to_timedelta(['2 days']).values
        tm.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= LooseVersion('1.8'):
            tm.assert_numpy_array_equal(other + td, expected)
        pytest.raises(TypeError, lambda: td + np.array([1]))
        pytest.raises(TypeError, lambda: np.array([1]) + td)

        expected = pd.to_timedelta(['0 days']).values
        tm.assert_numpy_array_equal(td - other, expected)
        if LooseVersion(np.__version__) >= LooseVersion('1.8'):
            tm.assert_numpy_array_equal(-other + td, expected)
        pytest.raises(TypeError, lambda: td - np.array([1]))
        pytest.raises(TypeError, lambda: np.array([1]) - td)

        expected = pd.to_timedelta(['2 days']).values
        tm.assert_numpy_array_equal(td * np.array([2]), expected)
        tm.assert_numpy_array_equal(np.array([2]) * td, expected)
        pytest.raises(TypeError, lambda: td * other)
        pytest.raises(TypeError, lambda: other * td)

        tm.assert_numpy_array_equal(td / other,
                                    np.array([1], dtype=np.float64))
        if LooseVersion(np.__version__) >= LooseVersion('1.8'):
            tm.assert_numpy_array_equal(other / td,
                                        np.array([1], dtype=np.float64))

        # timedelta, datetime
        other = pd.to_datetime(['2000-01-01']).values
        expected = pd.to_datetime(['2000-01-02']).values
        tm.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= LooseVersion('1.8'):
            tm.assert_numpy_array_equal(other + td, expected)

        expected = pd.to_datetime(['1999-12-31']).values
        tm.assert_numpy_array_equal(-td + other, expected)
        if LooseVersion(np.__version__) >= LooseVersion('1.8'):
            tm.assert_numpy_array_equal(other - td, expected)

    def test_ops_series(self):
        # regression test for GH8813
        td = Timedelta('1 day')
        other = pd.Series([1, 2])
        expected = pd.Series(pd.to_timedelta(['1 day', '2 days']))
        tm.assert_series_equal(expected, td * other)
        tm.assert_series_equal(expected, other * td)

    def test_ops_series_object(self):
        # GH 13043
        s = pd.Series([pd.Timestamp('2015-01-01', tz='US/Eastern'),
                       pd.Timestamp('2015-01-01', tz='Asia/Tokyo')],
                      name='xxx')
        assert s.dtype == object

        exp = pd.Series([pd.Timestamp('2015-01-02', tz='US/Eastern'),
                         pd.Timestamp('2015-01-02', tz='Asia/Tokyo')],
                        name='xxx')
        tm.assert_series_equal(s + pd.Timedelta('1 days'), exp)
        tm.assert_series_equal(pd.Timedelta('1 days') + s, exp)

        # object series & object series
        s2 = pd.Series([pd.Timestamp('2015-01-03', tz='US/Eastern'),
                        pd.Timestamp('2015-01-05', tz='Asia/Tokyo')],
                       name='xxx')
        assert s2.dtype == object
        exp = pd.Series([pd.Timedelta('2 days'), pd.Timedelta('4 days')],
                        name='xxx')
        tm.assert_series_equal(s2 - s, exp)
        tm.assert_series_equal(s - s2, -exp)

        s = pd.Series([pd.Timedelta('01:00:00'), pd.Timedelta('02:00:00')],
                      name='xxx', dtype=object)
        assert s.dtype == object

        exp = pd.Series([pd.Timedelta('01:30:00'), pd.Timedelta('02:30:00')],
                        name='xxx')
        tm.assert_series_equal(s + pd.Timedelta('00:30:00'), exp)
        tm.assert_series_equal(pd.Timedelta('00:30:00') + s, exp)

    def test_timedelta_ops_with_missing_values(self):
        # setup
        s1 = pd.to_timedelta(Series(['00:00:01']))
        s2 = pd.to_timedelta(Series(['00:00:02']))
        sn = pd.to_timedelta(Series([pd.NaT]))
        df1 = pd.DataFrame(['00:00:01']).apply(pd.to_timedelta)
        df2 = pd.DataFrame(['00:00:02']).apply(pd.to_timedelta)
        dfn = pd.DataFrame([pd.NaT]).apply(pd.to_timedelta)
        scalar1 = pd.to_timedelta('00:00:01')
        scalar2 = pd.to_timedelta('00:00:02')
        timedelta_NaT = pd.to_timedelta('NaT')
        NA = np.nan

        actual = scalar1 + scalar1
        assert actual == scalar2
        actual = scalar2 - scalar1
        assert actual == scalar1

        actual = s1 + s1
        tm.assert_series_equal(actual, s2)
        actual = s2 - s1
        tm.assert_series_equal(actual, s1)

        actual = s1 + scalar1
        tm.assert_series_equal(actual, s2)
        actual = scalar1 + s1
        tm.assert_series_equal(actual, s2)
        actual = s2 - scalar1
        tm.assert_series_equal(actual, s1)
        actual = -scalar1 + s2
        tm.assert_series_equal(actual, s1)

        actual = s1 + timedelta_NaT
        tm.assert_series_equal(actual, sn)
        actual = timedelta_NaT + s1
        tm.assert_series_equal(actual, sn)
        actual = s1 - timedelta_NaT
        tm.assert_series_equal(actual, sn)
        actual = -timedelta_NaT + s1
        tm.assert_series_equal(actual, sn)

        actual = s1 + NA
        tm.assert_series_equal(actual, sn)
        actual = NA + s1
        tm.assert_series_equal(actual, sn)
        actual = s1 - NA
        tm.assert_series_equal(actual, sn)
        actual = -NA + s1
        tm.assert_series_equal(actual, sn)

        actual = s1 + pd.NaT
        tm.assert_series_equal(actual, sn)
        actual = s2 - pd.NaT
        tm.assert_series_equal(actual, sn)

        actual = s1 + df1
        tm.assert_frame_equal(actual, df2)
        actual = s2 - df1
        tm.assert_frame_equal(actual, df1)
        actual = df1 + s1
        tm.assert_frame_equal(actual, df2)
        actual = df2 - s1
        tm.assert_frame_equal(actual, df1)

        actual = df1 + df1
        tm.assert_frame_equal(actual, df2)
        actual = df2 - df1
        tm.assert_frame_equal(actual, df1)

        actual = df1 + scalar1
        tm.assert_frame_equal(actual, df2)
        actual = df2 - scalar1
        tm.assert_frame_equal(actual, df1)

        actual = df1 + timedelta_NaT
        tm.assert_frame_equal(actual, dfn)
        actual = df1 - timedelta_NaT
        tm.assert_frame_equal(actual, dfn)

        actual = df1 + NA
        tm.assert_frame_equal(actual, dfn)
        actual = df1 - NA
        tm.assert_frame_equal(actual, dfn)

        actual = df1 + pd.NaT  # NaT is datetime, not timedelta
        tm.assert_frame_equal(actual, dfn)
        actual = df1 - pd.NaT
        tm.assert_frame_equal(actual, dfn)

    def test_add_overflow(self):
        # see gh-14068
        msg = "too (big|large) to convert"
        with tm.assert_raises_regex(OverflowError, msg):
            to_timedelta(106580, 'D') + Timestamp('2000')
        with tm.assert_raises_regex(OverflowError, msg):
            Timestamp('2000') + to_timedelta(106580, 'D')

        _NaT = int(pd.NaT) + 1
        msg = "Overflow in int64 addition"
        with tm.assert_raises_regex(OverflowError, msg):
            to_timedelta([106580], 'D') + Timestamp('2000')
        with tm.assert_raises_regex(OverflowError, msg):
            Timestamp('2000') + to_timedelta([106580], 'D')
        with tm.assert_raises_regex(OverflowError, msg):
            to_timedelta([_NaT]) - Timedelta('1 days')
        with tm.assert_raises_regex(OverflowError, msg):
            to_timedelta(['5 days', _NaT]) - Timedelta('1 days')
        with tm.assert_raises_regex(OverflowError, msg):
            (to_timedelta([_NaT, '5 days', '1 hours']) -
             to_timedelta(['7 seconds', _NaT, '4 hours']))

        # These should not overflow!
        exp = TimedeltaIndex([pd.NaT])
        result = to_timedelta([pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex(['4 days', pd.NaT])
        result = to_timedelta(['5 days', pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex([pd.NaT, pd.NaT, '5 hours'])
        result = (to_timedelta([pd.NaT, '5 days', '1 hours']) +
                  to_timedelta(['7 seconds', pd.NaT, '4 hours']))
        tm.assert_index_equal(result, exp)

    def test_timedeltaindex_add_timestamp_nat_masking(self):
        # GH17991 checking for overflow-masking with NaT
        tdinat = pd.to_timedelta(['24658 days 11:15:00', 'NaT'])

        tsneg = Timestamp('1950-01-01')
        ts_neg_variants = [tsneg,
                           tsneg.to_pydatetime(),
                           tsneg.to_datetime64().astype('datetime64[ns]'),
                           tsneg.to_datetime64().astype('datetime64[D]')]

        tspos = Timestamp('1980-01-01')
        ts_pos_variants = [tspos,
                           tspos.to_pydatetime(),
                           tspos.to_datetime64().astype('datetime64[ns]'),
                           tspos.to_datetime64().astype('datetime64[D]')]

        for variant in ts_neg_variants + ts_pos_variants:
            res = tdinat + variant
            assert res[1] is pd.NaT

    def test_tdi_ops_attributes(self):
        rng = timedelta_range('2 days', periods=5, freq='2D', name='x')

        result = rng + 1
        exp = timedelta_range('4 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '2D'

        result = rng - 2
        exp = timedelta_range('-2 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '2D'

        result = rng * 2
        exp = timedelta_range('4 days', periods=5, freq='4D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '4D'

        result = rng / 2
        exp = timedelta_range('1 days', periods=5, freq='D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == 'D'

        result = -rng
        exp = timedelta_range('-2 days', periods=5, freq='-2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '-2D'

        rng = pd.timedelta_range('-2 days', periods=5, freq='D', name='x')

        result = abs(rng)
        exp = TimedeltaIndex(['2 days', '1 days', '0 days', '1 days',
                              '2 days'], name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq is None

    # TODO: Needs more informative name, probably split up into
    # more targeted tests
    def test_timedelta(self, freq):
        index = date_range('1/1/2000', periods=50, freq=freq)

        shifted = index + timedelta(1)
        back = shifted + timedelta(-1)
        tm.assert_index_equal(index, back)

        if freq == 'D':
            expected = pd.tseries.offsets.Day(1)
            assert index.freq == expected
            assert shifted.freq == expected
            assert back.freq == expected
        else:  # freq == 'B'
            assert index.freq == pd.tseries.offsets.BusinessDay(1)
            assert shifted.freq is None
            assert back.freq == pd.tseries.offsets.BusinessDay(1)

        result = index - timedelta(1)
        expected = index + timedelta(-1)
        tm.assert_index_equal(result, expected)

        # GH4134, buggy with timedeltas
        rng = date_range('2013', '2014')
        s = Series(rng)
        result1 = rng - pd.offsets.Hour(1)
        result2 = DatetimeIndex(s - np.timedelta64(100000000))
        result3 = rng - np.timedelta64(100000000)
        result4 = DatetimeIndex(s - pd.offsets.Hour(1))
        tm.assert_index_equal(result1, result4)
        tm.assert_index_equal(result2, result3)
