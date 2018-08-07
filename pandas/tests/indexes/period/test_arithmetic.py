# -*- coding: utf-8 -*-

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import (Period, PeriodIndex,
                    _np_version_under1p10)
import pandas.core.indexes.period as period


class TestPeriodIndexArithmetic(object):
    # ---------------------------------------------------------------
    # PeriodIndex.shift is used by __add__ and __sub__

    def test_pi_shift_ndarray(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        result = idx.shift(np.array([1, 2, 3, 4]))
        expected = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)

        result = idx.shift(np.array([1, -2, 3, -4]))
        expected = PeriodIndex(['2011-02', '2010-12', 'NaT', '2010-12'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)

    def test_shift(self):
        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2002', end='12/1/2010')

        tm.assert_index_equal(pi1.shift(0), pi1)

        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2000', end='12/1/2008')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='2/1/2001', end='1/1/2010')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='12/1/2000', end='11/1/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='1/2/2001', end='12/2/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='12/31/2000', end='11/30/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

    def test_shift_corner_cases(self):
        # GH#9903
        idx = pd.PeriodIndex([], name='xxx', freq='H')

        with pytest.raises(TypeError):
            # period shift doesn't accept freq
            idx.shift(1, freq='H')

        tm.assert_index_equal(idx.shift(0), idx)
        tm.assert_index_equal(idx.shift(3), idx)

        idx = pd.PeriodIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                              '2011-01-01 12:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(0), idx)
        exp = pd.PeriodIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                              '2011-01-01 15:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(3), exp)
        exp = pd.PeriodIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                              '2011-01-01 09:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(-3), exp)

    def test_shift_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        result = idx.shift(1)
        expected = PeriodIndex(['2011-02', '2011-03', 'NaT', '2011-05'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

    def test_shift_gh8083(self):
        # test shift for PeriodIndex
        # GH#8083
        drange = pd.period_range('20130101', periods=5, freq='D')
        result = drange.shift(1)
        expected = PeriodIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                '2013-01-05', '2013-01-06'], freq='D')
        tm.assert_index_equal(result, expected)


class TestPeriodIndexSeriesMethods(object):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        if isinstance(expected, pd.Index):
            tm.assert_index_equal(result, expected)
        else:
            # comp op results in bool
            tm.assert_numpy_array_equal(result, expected)

        ser = pd.Series(values)
        result = func(ser)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_ops(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        expected = PeriodIndex(['2011-03', '2011-04', '2011-05', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        result = idx - Period('2011-01', freq='M')
        off = idx.freq
        exp = pd.Index([0 * off, 1 * off, 2 * off, 3 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = Period('2011-01', freq='M') - idx
        exp = pd.Index([0 * off, -1 * off, -2 * off, -3 * off], name='idx')
        tm.assert_index_equal(result, exp)

    @pytest.mark.parametrize('ng', ["str", 1.5])
    def test_pi_ops_errors(self, ng):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')
        ser = pd.Series(idx)

        msg = r"unsupported operand type\(s\)"

        for obj in [idx, ser]:
            with tm.assert_raises_regex(TypeError, msg):
                obj + ng

            with pytest.raises(TypeError):
                # error message differs between PY2 and 3
                ng + obj

            with tm.assert_raises_regex(TypeError, msg):
                obj - ng

            with pytest.raises(TypeError):
                np.add(obj, ng)

            if _np_version_under1p10:
                assert np.add(ng, obj) is NotImplemented
            else:
                with pytest.raises(TypeError):
                    np.add(ng, obj)

            with pytest.raises(TypeError):
                np.subtract(obj, ng)

            if _np_version_under1p10:
                assert np.subtract(ng, obj) is NotImplemented
            else:
                with pytest.raises(TypeError):
                    np.subtract(ng, obj)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        expected = PeriodIndex(['2011-03', '2011-04', 'NaT', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)
        self._check(idx, lambda x: np.add(x, 2), expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        self._check(idx + 2, lambda x: np.subtract(x, 2), idx)

        # freq with mult
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='2M', name='idx')
        expected = PeriodIndex(['2011-07', '2011-08', 'NaT', '2011-10'],
                               freq='2M', name='idx')
        self._check(idx, lambda x: x + 3, expected)
        self._check(idx, lambda x: 3 + x, expected)
        self._check(idx, lambda x: np.add(x, 3), expected)

        self._check(idx + 3, lambda x: x - 3, idx)
        self._check(idx + 3, lambda x: np.subtract(x, 3), idx)

    def test_pi_ops_array_int(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        f = lambda x: x + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.add(x, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT', '2011-06'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT', '2010-12'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.subtract(x, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT', '2011-06'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

    def test_pi_ops_offset(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        f = lambda x: x + pd.offsets.Day()
        exp = PeriodIndex(['2011-01-02', '2011-02-02', '2011-03-02',
                           '2011-04-02'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x + pd.offsets.Day(2)
        exp = PeriodIndex(['2011-01-03', '2011-02-03', '2011-03-03',
                           '2011-04-03'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - pd.offsets.Day(2)
        exp = PeriodIndex(['2010-12-30', '2011-01-30', '2011-02-27',
                           '2011-03-30'], freq='D', name='idx')
        self._check(idx, f, exp)

    def test_pi_offset_errors(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        ser = pd.Series(idx)

        # Series op is applied per Period instance, thus error is raised
        # from Period
        msg_idx = r"Input has different freq from PeriodIndex\(freq=D\)"
        msg_s = r"Input cannot be converted to Period\(freq=D\)"
        for obj, msg in [(idx, msg_idx), (ser, msg_s)]:
            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj + pd.offsets.Hour(2)

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                pd.offsets.Hour(2) + obj

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj - pd.offsets.Hour(2)

    def test_pi_sub_period(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        off = idx.freq
        exp = pd.Index([-12 * off, -11 * off, -10 * off, -9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(idx, pd.Period('2012-01', freq='M'))
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12 * off, 11 * off, 10 * off, 9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(pd.Period('2012-01', freq='M'), idx)
        if _np_version_under1p10:
            assert result is NotImplemented
        else:
            tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_sub_pdnat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        exp = pd.TimedeltaIndex([pd.NaT] * 4, name='idx')
        tm.assert_index_equal(pd.NaT - idx, exp)
        tm.assert_index_equal(idx - pd.NaT, exp)

    def test_pi_sub_period_nat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        off = idx.freq
        exp = pd.Index([-12 * off, pd.NaT, -10 * off, -9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12 * off, pd.NaT, 10 * off, 9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)
