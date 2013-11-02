# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import nose
import unittest

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

class TestTimedeltas(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    def test_numeric_conversions(self):
        _skip_if_numpy_not_friendly()

        self.assert_(ct(0) == np.timedelta64(0,'ns'))
        self.assert_(ct(10) == np.timedelta64(10,'ns'))
        self.assert_(ct(10,unit='ns') == np.timedelta64(10,'ns').astype('m8[ns]'))

        self.assert_(ct(10,unit='us') == np.timedelta64(10,'us').astype('m8[ns]'))
        self.assert_(ct(10,unit='ms') == np.timedelta64(10,'ms').astype('m8[ns]'))
        self.assert_(ct(10,unit='s') == np.timedelta64(10,'s').astype('m8[ns]'))
        self.assert_(ct(10,unit='d') == np.timedelta64(10,'D').astype('m8[ns]'))

    def test_timedelta_conversions(self):
        _skip_if_numpy_not_friendly()

        self.assert_(ct(timedelta(seconds=1)) == np.timedelta64(1,'s').astype('m8[ns]'))
        self.assert_(ct(timedelta(microseconds=1)) == np.timedelta64(1,'us').astype('m8[ns]'))
        self.assert_(ct(timedelta(days=1)) == np.timedelta64(1,'D').astype('m8[ns]'))

    def test_short_format_converters(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')

        self.assert_(ct('10') == np.timedelta64(10,'ns'))
        self.assert_(ct('10ns') == np.timedelta64(10,'ns'))
        self.assert_(ct('100') == np.timedelta64(100,'ns'))
        self.assert_(ct('100ns') == np.timedelta64(100,'ns'))

        self.assert_(ct('1000') == np.timedelta64(1000,'ns'))
        self.assert_(ct('1000ns') == np.timedelta64(1000,'ns'))
        self.assert_(ct('1000NS') == np.timedelta64(1000,'ns'))

        self.assert_(ct('10us') == np.timedelta64(10000,'ns'))
        self.assert_(ct('100us') == np.timedelta64(100000,'ns'))
        self.assert_(ct('1000us') == np.timedelta64(1000000,'ns'))
        self.assert_(ct('1000Us') == np.timedelta64(1000000,'ns'))
        self.assert_(ct('1000uS') == np.timedelta64(1000000,'ns'))

        self.assert_(ct('1ms') == np.timedelta64(1000000,'ns'))
        self.assert_(ct('10ms') == np.timedelta64(10000000,'ns'))
        self.assert_(ct('100ms') == np.timedelta64(100000000,'ns'))
        self.assert_(ct('1000ms') == np.timedelta64(1000000000,'ns'))

        self.assert_(ct('-1s') == -np.timedelta64(1000000000,'ns'))
        self.assert_(ct('1s') == np.timedelta64(1000000000,'ns'))
        self.assert_(ct('10s') == np.timedelta64(10000000000,'ns'))
        self.assert_(ct('100s') == np.timedelta64(100000000000,'ns'))
        self.assert_(ct('1000s') == np.timedelta64(1000000000000,'ns'))

        self.assert_(ct('1d') == conv(np.timedelta64(1,'D')))
        self.assert_(ct('-1d') == -conv(np.timedelta64(1,'D')))
        self.assert_(ct('1D') == conv(np.timedelta64(1,'D')))
        self.assert_(ct('10D') == conv(np.timedelta64(10,'D')))
        self.assert_(ct('100D') == conv(np.timedelta64(100,'D')))
        self.assert_(ct('1000D') == conv(np.timedelta64(1000,'D')))
        self.assert_(ct('10000D') == conv(np.timedelta64(10000,'D')))

        # space
        self.assert_(ct(' 10000D ') == conv(np.timedelta64(10000,'D')))
        self.assert_(ct(' - 10000D ') == -conv(np.timedelta64(10000,'D')))

        # invalid
        self.assertRaises(ValueError, ct, '1foo')
        self.assertRaises(ValueError, ct, 'foo')

    def test_full_format_converters(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')
        d1 = np.timedelta64(1,'D')

        self.assert_(ct('1days') == conv(d1))
        self.assert_(ct('1days,') == conv(d1))
        self.assert_(ct('- 1days,') == -conv(d1))

        self.assert_(ct('00:00:01') == conv(np.timedelta64(1,'s')))
        self.assert_(ct('06:00:01') == conv(np.timedelta64(6*3600+1,'s')))
        self.assert_(ct('06:00:01.0') == conv(np.timedelta64(6*3600+1,'s')))
        self.assert_(ct('06:00:01.01') == conv(np.timedelta64(1000*(6*3600+1)+10,'ms')))

        self.assert_(ct('- 1days, 00:00:01') == -conv(d1+np.timedelta64(1,'s')))
        self.assert_(ct('1days, 06:00:01') == conv(d1+np.timedelta64(6*3600+1,'s')))
        self.assert_(ct('1days, 06:00:01.01') == conv(d1+np.timedelta64(1000*(6*3600+1)+10,'ms')))

        # invalid
        self.assertRaises(ValueError, ct, '- 1days, 00')

    def test_nat_converters(self):
        _skip_if_numpy_not_friendly()

        self.assert_(to_timedelta('nat',box=False) == tslib.iNaT)
        self.assert_(to_timedelta('nan',box=False) == tslib.iNaT)

    def test_to_timedelta(self):
        _skip_if_numpy_not_friendly()

        def conv(v):
            return v.astype('m8[ns]')
        d1 = np.timedelta64(1,'D')

        self.assert_(to_timedelta('1 days 06:05:01.00003',box=False) == conv(d1+np.timedelta64(6*3600+5*60+1,'s')+np.timedelta64(30,'us')))
        self.assert_(to_timedelta('15.5us',box=False) == conv(np.timedelta64(15500,'ns')))

        # empty string
        result = to_timedelta('',box=False)
        self.assert_(result == tslib.iNaT)

        result = to_timedelta(['', ''])
        self.assert_(isnull(result).all())

        # pass thru
        result = to_timedelta(np.array([np.timedelta64(1,'s')]))
        expected = np.array([np.timedelta64(1,'s')])
        tm.assert_almost_equal(result,expected)

        # ints
        result = np.timedelta64(0,'ns')
        expected = to_timedelta(0,box=False)
        self.assert_(result == expected)

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
        expected = to_timedelta([v])

        v = np.timedelta64(timedelta(seconds=1))
        result = to_timedelta(v,box=False)
        expected = to_timedelta([v])

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
        expected = to_timedelta('00:00:02.6')
        tm.assert_almost_equal(result, expected)

        result = td.median()[0]
        # TODO This should have returned a scalar to begin with. Hack for now.
        expected = to_timedelta('00:00:08')
        tm.assert_almost_equal(result, expected)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
