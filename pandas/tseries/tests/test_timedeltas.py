# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import nose

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, Timestamp, isnull, notnull,
                    bdate_range, date_range, _np_version_under1p7)
import pandas.core.common as com
from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat, to_timedelta, tslib
from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type as ct
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_almost_equal,
                                 ensure_clean)
import pandas.util.testing as tm

def _skip_if_numpy_not_friendly():
    # not friendly for < 1.7
    if _np_version_under1p7:
        raise nose.SkipTest("numpy < 1.7")

class TestTimedeltas(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    def test_numeric_conversions(self):
        _skip_if_numpy_not_friendly()

        self.assertEqual(ct(0), np.timedelta64(0,'ns'))
        self.assertEqual(ct(10), np.timedelta64(10,'ns'))
        self.assertEqual(ct(10,unit='ns'), np.timedelta64(10,'ns').astype('m8[ns]'))

        self.assertEqual(ct(10,unit='us'), np.timedelta64(10,'us').astype('m8[ns]'))
        self.assertEqual(ct(10,unit='ms'), np.timedelta64(10,'ms').astype('m8[ns]'))
        self.assertEqual(ct(10,unit='s'), np.timedelta64(10,'s').astype('m8[ns]'))
        self.assertEqual(ct(10,unit='d'), np.timedelta64(10,'D').astype('m8[ns]'))

    def test_timedelta_conversions(self):
        _skip_if_numpy_not_friendly()

        self.assertEqual(ct(timedelta(seconds=1)), np.timedelta64(1,'s').astype('m8[ns]'))
        self.assertEqual(ct(timedelta(microseconds=1)), np.timedelta64(1,'us').astype('m8[ns]'))
        self.assertEqual(ct(timedelta(days=1)), np.timedelta64(1,'D').astype('m8[ns]'))

    def test_short_format_converters(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')

        self.assertEqual(ct('10'), np.timedelta64(10,'ns'))
        self.assertEqual(ct('10ns'), np.timedelta64(10,'ns'))
        self.assertEqual(ct('100'), np.timedelta64(100,'ns'))
        self.assertEqual(ct('100ns'), np.timedelta64(100,'ns'))

        self.assertEqual(ct('1000'), np.timedelta64(1000,'ns'))
        self.assertEqual(ct('1000ns'), np.timedelta64(1000,'ns'))
        self.assertEqual(ct('1000NS'), np.timedelta64(1000,'ns'))

        self.assertEqual(ct('10us'), np.timedelta64(10000,'ns'))
        self.assertEqual(ct('100us'), np.timedelta64(100000,'ns'))
        self.assertEqual(ct('1000us'), np.timedelta64(1000000,'ns'))
        self.assertEqual(ct('1000Us'), np.timedelta64(1000000,'ns'))
        self.assertEqual(ct('1000uS'), np.timedelta64(1000000,'ns'))

        self.assertEqual(ct('1ms'), np.timedelta64(1000000,'ns'))
        self.assertEqual(ct('10ms'), np.timedelta64(10000000,'ns'))
        self.assertEqual(ct('100ms'), np.timedelta64(100000000,'ns'))
        self.assertEqual(ct('1000ms'), np.timedelta64(1000000000,'ns'))

        self.assertEqual(ct('-1s'), -np.timedelta64(1000000000,'ns'))
        self.assertEqual(ct('1s'), np.timedelta64(1000000000,'ns'))
        self.assertEqual(ct('10s'), np.timedelta64(10000000000,'ns'))
        self.assertEqual(ct('100s'), np.timedelta64(100000000000,'ns'))
        self.assertEqual(ct('1000s'), np.timedelta64(1000000000000,'ns'))

        self.assertEqual(ct('1d'), conv(np.timedelta64(1,'D')))
        self.assertEqual(ct('-1d'), -conv(np.timedelta64(1,'D')))
        self.assertEqual(ct('1D'), conv(np.timedelta64(1,'D')))
        self.assertEqual(ct('10D'), conv(np.timedelta64(10,'D')))
        self.assertEqual(ct('100D'), conv(np.timedelta64(100,'D')))
        self.assertEqual(ct('1000D'), conv(np.timedelta64(1000,'D')))
        self.assertEqual(ct('10000D'), conv(np.timedelta64(10000,'D')))

        # space
        self.assertEqual(ct(' 10000D '), conv(np.timedelta64(10000,'D')))
        self.assertEqual(ct(' - 10000D '), -conv(np.timedelta64(10000,'D')))

        # invalid
        self.assertRaises(ValueError, ct, '1foo')
        self.assertRaises(ValueError, ct, 'foo')

    def test_full_format_converters(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')
        d1 = np.timedelta64(1,'D')

        self.assertEqual(ct('1days'), conv(d1))
        self.assertEqual(ct('1days,'), conv(d1))
        self.assertEqual(ct('- 1days,'), -conv(d1))

        self.assertEqual(ct('00:00:01'), conv(np.timedelta64(1,'s')))
        self.assertEqual(ct('06:00:01'), conv(np.timedelta64(6*3600+1,'s')))
        self.assertEqual(ct('06:00:01.0'), conv(np.timedelta64(6*3600+1,'s')))
        self.assertEqual(ct('06:00:01.01'), conv(np.timedelta64(1000*(6*3600+1)+10,'ms')))

        self.assertEqual(ct('- 1days, 00:00:01'), -conv(d1+np.timedelta64(1,'s')))
        self.assertEqual(ct('1days, 06:00:01'), conv(d1+np.timedelta64(6*3600+1,'s')))
        self.assertEqual(ct('1days, 06:00:01.01'), conv(d1+np.timedelta64(1000*(6*3600+1)+10,'ms')))

        # invalid
        self.assertRaises(ValueError, ct, '- 1days, 00')

    def test_nat_converters(self):
        _skip_if_numpy_not_friendly()

        self.assertEqual(to_timedelta('nat',box=False).astype('int64'), tslib.iNaT)
        self.assertEqual(to_timedelta('nan',box=False).astype('int64'), tslib.iNaT)

    def test_to_timedelta(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')
        d1 = np.timedelta64(1,'D')

        self.assertEqual(to_timedelta('1 days 06:05:01.00003',box=False), conv(d1+np.timedelta64(6*3600+5*60+1,'s')+np.timedelta64(30,'us')))
        self.assertEqual(to_timedelta('15.5us',box=False), conv(np.timedelta64(15500,'ns')))

        # empty string
        result = to_timedelta('',box=False)
        self.assertEqual(result.astype('int64'), tslib.iNaT)

        result = to_timedelta(['', ''])
        self.assertTrue(isnull(result).all())

        # pass thru
        result = to_timedelta(np.array([np.timedelta64(1,'s')]))
        expected = np.array([np.timedelta64(1,'s')])
        tm.assert_almost_equal(result,expected)

        # ints
        result = np.timedelta64(0,'ns')
        expected = to_timedelta(0,box=False)
        self.assertEqual(result, expected)

        # Series
        expected = Series([timedelta(days=1), timedelta(days=1, seconds=1)])
        result = to_timedelta(Series(['1d','1days 00:00:01']))
        tm.assert_series_equal(result, expected)

        # with units
        result = Series([ np.timedelta64(0,'ns'), np.timedelta64(10,'s').astype('m8[ns]') ],dtype='m8[ns]')
        expected = to_timedelta([0,10],unit='s')
        tm.assert_series_equal(result, expected)

        # single element conversion
        v = timedelta(seconds=1)
        result = to_timedelta(v,box=False)
        expected = np.timedelta64(timedelta(seconds=1))
        self.assertEqual(result, expected)

        v = np.timedelta64(timedelta(seconds=1))
        result = to_timedelta(v,box=False)
        expected = np.timedelta64(timedelta(seconds=1))
        self.assertEqual(result, expected)

        # arrays of various dtypes
        arr = np.array([1]*5,dtype='int64')
        result = to_timedelta(arr,unit='s')
        expected = Series([ np.timedelta64(1,'s') ]*5)
        tm.assert_series_equal(result, expected)

        arr = np.array([1]*5,dtype='int64')
        result = to_timedelta(arr,unit='m')
        expected = Series([ np.timedelta64(1,'m') ]*5)
        tm.assert_series_equal(result, expected)

        arr = np.array([1]*5,dtype='int64')
        result = to_timedelta(arr,unit='h')
        expected = Series([ np.timedelta64(1,'h') ]*5)
        tm.assert_series_equal(result, expected)

        arr = np.array([1]*5,dtype='timedelta64[s]')
        result = to_timedelta(arr)
        expected = Series([ np.timedelta64(1,'s') ]*5)
        tm.assert_series_equal(result, expected)

        arr = np.array([1]*5,dtype='timedelta64[D]')
        result = to_timedelta(arr)
        expected = Series([ np.timedelta64(1,'D') ]*5)
        tm.assert_series_equal(result, expected)

        # validate all units
        # GH 6855
        for unit in ['Y','M','W','D','y','w','d']:
            result = to_timedelta(np.arange(5),unit=unit)
            expected = Series([ np.timedelta64(i,unit.upper()) for i in np.arange(5).tolist() ])
            tm.assert_series_equal(result, expected)
        for unit in ['h','m','s','ms','us','ns','H','S','MS','US','NS']:
            result = to_timedelta(np.arange(5),unit=unit)
            expected = Series([ np.timedelta64(i,unit.lower()) for i in np.arange(5).tolist() ])
            tm.assert_series_equal(result, expected)

        # these will error
        self.assertRaises(ValueError, lambda : to_timedelta(['1h']))
        self.assertRaises(ValueError, lambda : to_timedelta(['1m']))

    def test_to_timedelta_via_apply(self):
        _skip_if_numpy_not_friendly()

        # GH 5458
        expected = Series([np.timedelta64(1,'s')])
        result = Series(['00:00:01']).apply(to_timedelta)
        tm.assert_series_equal(result, expected)

        result = Series([to_timedelta('00:00:01')])
        tm.assert_series_equal(result, expected)

    def test_timedelta_ops(self):
        _skip_if_numpy_not_friendly()

        # GH4984
        # make sure ops return timedeltas
        s = Series([Timestamp('20130101') + timedelta(seconds=i*i) for i in range(10) ])
        td = s.diff()

        result = td.mean()[0]
        # TODO This should have returned a scalar to begin with. Hack for now.
        expected = to_timedelta(timedelta(seconds=9))
        tm.assert_almost_equal(result, expected)

        result = td.quantile(.1)
        # This properly returned a scalar.
        expected = np.timedelta64(2599999999,'ns')
        tm.assert_almost_equal(result, expected)

        result = td.median()[0]
        # TODO This should have returned a scalar to begin with. Hack for now.
        expected = to_timedelta('00:00:08')
        tm.assert_almost_equal(result, expected)

        # GH 6462
        # consistency in returned values for sum
        result = td.sum()[0]
        expected = to_timedelta('00:01:21')
        tm.assert_almost_equal(result, expected)

    def test_timedelta_ops_scalar(self):
        _skip_if_numpy_not_friendly()

        # GH 6808
        base = pd.to_datetime('20130101 09:01:12.123456')
        expected_add = pd.to_datetime('20130101 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta(10,unit='s'),
                       timedelta(seconds=10),
                       np.timedelta64(10,'s'),
                       np.timedelta64(10000000000,'ns'),
                       pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

        base = pd.to_datetime('20130102 09:01:12.123456')
        expected_add = pd.to_datetime('20130103 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta('1 day, 00:00:10'),
                       pd.to_timedelta('1 days, 00:00:10'),
                       timedelta(days=1,seconds=10),
                       np.timedelta64(1,'D')+np.timedelta64(10,'s'),
                       pd.offsets.Day()+pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

    def test_to_timedelta_on_missing_values(self):
        _skip_if_numpy_not_friendly()

        # GH5438
        timedelta_NaT = np.timedelta64('NaT')

        actual = pd.to_timedelta(Series(['00:00:01', np.nan]))
        expected = Series([np.timedelta64(1000000000, 'ns'), timedelta_NaT], dtype='<m8[ns]')
        assert_series_equal(actual, expected)

        actual = pd.to_timedelta(Series(['00:00:01', pd.NaT]))
        assert_series_equal(actual, expected)

        actual = pd.to_timedelta(np.nan)
        self.assertEqual(actual.astype('int64'), timedelta_NaT.astype('int64'))

        actual = pd.to_timedelta(pd.NaT)
        self.assertEqual(actual.astype('int64'), timedelta_NaT.astype('int64'))

    def test_timedelta_ops_with_missing_values(self):
        _skip_if_numpy_not_friendly()

        # setup
        s1 = pd.to_timedelta(Series(['00:00:01']))
        s2 = pd.to_timedelta(Series(['00:00:02']))
        sn = pd.to_timedelta(Series([pd.NaT]))
        df1 = DataFrame(['00:00:01']).apply(pd.to_timedelta)
        df2 = DataFrame(['00:00:02']).apply(pd.to_timedelta)
        dfn = DataFrame([pd.NaT]).apply(pd.to_timedelta)
        scalar1 = pd.to_timedelta('00:00:01')
        scalar2 = pd.to_timedelta('00:00:02')
        timedelta_NaT = pd.to_timedelta('NaT')
        NA = np.nan

        actual = scalar1 + scalar1
        self.assertEqual(actual, scalar2)
        actual = scalar2 - scalar1
        self.assertEqual(actual, scalar1)

        actual = s1 + s1
        assert_series_equal(actual, s2)
        actual = s2 - s1
        assert_series_equal(actual, s1)

        actual = s1 + scalar1
        assert_series_equal(actual, s2)
        actual = s2 - scalar1
        assert_series_equal(actual, s1)

        actual = s1 + timedelta_NaT
        assert_series_equal(actual, sn)
        actual = s1 - timedelta_NaT
        assert_series_equal(actual, sn)

        actual = s1 + NA
        assert_series_equal(actual, sn)
        actual = s1 - NA
        assert_series_equal(actual, sn)

        actual = s1 + pd.NaT  # NaT is datetime, not timedelta
        assert_series_equal(actual, sn)
        actual = s2 - pd.NaT
        assert_series_equal(actual, sn)

        actual = s1 + df1
        assert_frame_equal(actual, df2)
        actual = s2 - df1
        assert_frame_equal(actual, df1)
        actual = df1 + s1
        assert_frame_equal(actual, df2)
        actual = df2 - s1
        assert_frame_equal(actual, df1)

        actual = df1 + df1
        assert_frame_equal(actual, df2)
        actual = df2 - df1
        assert_frame_equal(actual, df1)

        actual = df1 + scalar1
        assert_frame_equal(actual, df2)
        actual = df2 - scalar1
        assert_frame_equal(actual, df1)

        actual = df1 + timedelta_NaT
        assert_frame_equal(actual, dfn)
        actual = df1 - timedelta_NaT
        assert_frame_equal(actual, dfn)

        actual = df1 + NA
        assert_frame_equal(actual, dfn)
        actual = df1 - NA
        assert_frame_equal(actual, dfn)

        actual = df1 + pd.NaT  # NaT is datetime, not timedelta
        assert_frame_equal(actual, dfn)
        actual = df1 - pd.NaT
        assert_frame_equal(actual, dfn)

    def test_apply_to_timedelta(self):
        _skip_if_numpy_not_friendly()

        timedelta_NaT = pd.to_timedelta('NaT')

        list_of_valid_strings = ['00:00:01', '00:00:02']
        a = pd.to_timedelta(list_of_valid_strings)
        b = Series(list_of_valid_strings).apply(pd.to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

        list_of_strings = ['00:00:01', np.nan, pd.NaT, timedelta_NaT]
        a = pd.to_timedelta(list_of_strings)
        b = Series(list_of_strings).apply(pd.to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
