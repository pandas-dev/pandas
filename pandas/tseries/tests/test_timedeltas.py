# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import nose
import unittest

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, isnull, notnull,
                    bdate_range, date_range, _np_version_under1p7)
import pandas.core.common as com
from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat
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

        # ns not converted properly
        self.assert_(ct(0) == np.timedelta64(0,'ns'))
        self.assert_(ct(10) == np.timedelta64(0,'ns'))
        self.assert_(ct(10,unit='ns') == np.timedelta64(0,'ns').astype('m8[ns]'))

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

        # ns not converted properly
        self.assert_(ct('10') == np.timedelta64(0,'ns'))
        self.assert_(ct('10ns') == np.timedelta64(0,'ns'))
        self.assert_(ct('100') == np.timedelta64(0,'ns'))
        self.assert_(ct('100ns') == np.timedelta64(0,'ns'))

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


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
