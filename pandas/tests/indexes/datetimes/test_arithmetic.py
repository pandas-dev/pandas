# -*- coding: utf-8 -*-
from datetime import datetime, timedelta

import pytest
import pytz
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas.errors import PerformanceWarning, NullFrequencyError
from pandas import (Timestamp, Timedelta, Series,
                    DatetimeIndex,
                    date_range)


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=str)
def delta(request):
    # Several ways of representing two hours
    return request.param


class TestDatetimeIndexArithmetic(object):

    # -------------------------------------------------------------
    # DatetimeIndex.shift is used in integer addition

    def test_dti_shift_tzaware(self, tz_naive_fixture):
        # GH#9903
        tz = tz_naive_fixture
        idx = pd.DatetimeIndex([], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        tm.assert_index_equal(idx.shift(3, freq='H'), idx)

        idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01 11:00',
                                '2011-01-01 12:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        exp = pd.DatetimeIndex(['2011-01-01 13:00', '2011-01-01 14:00',
                                '2011-01-01 15:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(3, freq='H'), exp)
        exp = pd.DatetimeIndex(['2011-01-01 07:00', '2011-01-01 08:00',
                                '2011-01-01 09:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

    def test_dti_shift_freqs(self):
        # test shift for DatetimeIndex and non DatetimeIndex
        # GH#8083
        drange = pd.date_range('20130101', periods=5)
        result = drange.shift(1)
        expected = pd.DatetimeIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                     '2013-01-05',
                                     '2013-01-06'], freq='D')
        tm.assert_index_equal(result, expected)

        result = drange.shift(-1)
        expected = pd.DatetimeIndex(['2012-12-31', '2013-01-01', '2013-01-02',
                                     '2013-01-03', '2013-01-04'],
                                    freq='D')
        tm.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D')
        expected = pd.DatetimeIndex(['2013-01-07', '2013-01-08', '2013-01-09',
                                     '2013-01-10',
                                     '2013-01-11'], freq='D')
        tm.assert_index_equal(result, expected)

    def test_dti_shift_int(self):
        rng = date_range('1/1/2000', periods=20)

        result = rng + 5
        expected = rng.shift(5)
        tm.assert_index_equal(result, expected)

        result = rng - 5
        expected = rng.shift(-5)
        tm.assert_index_equal(result, expected)

    def test_dti_shift_no_freq(self):
        # GH#19147
        dti = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01'], freq=None)
        with pytest.raises(NullFrequencyError):
            dti.shift(2)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_shift_localized(self, tzstr):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize(tzstr)

        result = dr_tz.shift(1, '10T')
        assert result.tz == dr_tz.tz

    def test_dti_shift_across_dst(self):
        # GH 8616
        idx = date_range('2013-11-03', tz='America/Chicago',
                         periods=7, freq='H')
        s = Series(index=idx[:-1])
        result = s.shift(freq='H')
        expected = Series(index=idx[1:])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('shift, result_time', [
        [0, '2014-11-14 00:00:00'],
        [-1, '2014-11-13 23:00:00'],
        [1, '2014-11-14 01:00:00']])
    def test_dti_shift_near_midnight(self, shift, result_time):
        # GH 8616
        dt = datetime(2014, 11, 14, 0)
        dt_est = pytz.timezone('EST').localize(dt)
        s = Series(data=[1], index=[dt_est])
        result = s.shift(shift, freq='H')
        expected = Series(1, index=DatetimeIndex([result_time], tz='EST'))
        tm.assert_series_equal(result, expected)

    # -------------------------------------------------------------

    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '3D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == 'D'

    def test_datetimeindex_sub_timestamp_overflow(self):
        dtimax = pd.to_datetime(['now', pd.Timestamp.max])
        dtimin = pd.to_datetime(['now', pd.Timestamp.min])

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

        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        expected = pd.Timestamp.max.value - tspos.value
        for variant in ts_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        expected = pd.Timestamp.min.value - tsneg.value
        for variant in ts_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant

    @pytest.mark.parametrize('names', [('foo', None, None),
                                       ('baz', 'bar', None),
                                       ('bar', 'bar', 'bar')])
    @pytest.mark.parametrize('tz', [None, 'America/Chicago'])
    def test_dti_add_series(self, tz, names):
        # GH#13905
        index = DatetimeIndex(['2016-06-28 05:30', '2016-06-28 05:31'],
                              tz=tz, name=names[0])
        ser = Series([Timedelta(seconds=5)] * 2,
                     index=index, name=names[1])
        expected = Series(index + Timedelta(seconds=5),
                          index=index, name=names[2])

        # passing name arg isn't enough when names[2] is None
        expected.name = names[2]
        assert expected.dtype == index.dtype
        result = ser + index
        tm.assert_series_equal(result, expected)
        result2 = index + ser
        tm.assert_series_equal(result2, expected)

        expected = index + Timedelta(seconds=5)
        result3 = ser.values + index
        tm.assert_index_equal(result3, expected)
        result4 = index + ser.values
        tm.assert_index_equal(result4, expected)

    def test_dti_add_offset_array(self, tz_naive_fixture):
        # GH#18849
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_add_offset_index(self, tz_naive_fixture, names):
        # GH#18849, GH#19744
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    def test_dti_sub_offset_array(self, tz_naive_fixture):
        # GH#18824
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_sub_offset_index(self, tz_naive_fixture, names):
        # GH#18824, GH#19744
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_with_offset_series(self, tz_naive_fixture, names):
        # GH#18849
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = Series([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                       name=names[1])

        expected_add = Series([dti[n] + other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        tm.assert_series_equal(res, expected_add)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_series_equal(res2, expected_add)

        expected_sub = Series([dti[n] - other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res3 = dti - other
        tm.assert_series_equal(res3, expected_sub)

    def test_dti_add_offset_tzaware(self, tz_aware_fixture):
        timezone = tz_aware_fixture
        if timezone == 'US/Pacific':
            dates = date_range('2012-11-01', periods=3, tz=timezone)
            offset = dates + pd.offsets.Hour(5)
            assert dates[0] + pd.offsets.Hour(5) == offset[0]

        dates = date_range('2010-11-01 00:00',
                           periods=3, tz=timezone, freq='H')
        expected = DatetimeIndex(['2010-11-01 05:00', '2010-11-01 06:00',
                                  '2010-11-01 07:00'], freq='H', tz=timezone)

        offset = dates + pd.offsets.Hour(5)
        tm.assert_index_equal(offset, expected)
        offset = dates + np.timedelta64(5, 'h')
        tm.assert_index_equal(offset, expected)
        offset = dates + timedelta(hours=5)
        tm.assert_index_equal(offset, expected)
