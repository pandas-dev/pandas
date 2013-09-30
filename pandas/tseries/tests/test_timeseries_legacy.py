# pylint: disable-msg=E1101,W0612
from datetime import datetime, time, timedelta
import sys
import os
import unittest

import nose

import numpy as np
randn = np.random.randn

from pandas import (Index, Series, TimeSeries, DataFrame,
                    isnull, date_range, Timestamp, DatetimeIndex,
                    Int64Index, to_datetime, bdate_range)

from pandas.core.daterange import DateRange
import pandas.core.datetools as datetools
import pandas.tseries.offsets as offsets
import pandas.tseries.frequencies as fmod
import pandas as pd

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

from pandas.tslib import NaT, iNaT
import pandas.lib as lib
import pandas.tslib as tslib

import pandas.index as _index

from pandas.compat import(
    range, long, StringIO, lrange, lmap, map, zip, cPickle as pickle, product
)
from pandas import read_pickle
import pandas.core.datetools as dt
from numpy.random import rand
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal
import pandas.compat as compat
from pandas.core.datetools import BDay
import pandas.core.common as com
from pandas import concat

from numpy.testing.decorators import slow


def _skip_if_no_pytz():
    try:
        import pytz
    except ImportError:
        raise nose.SkipTest("pytz not installed")

# infortunately, too much has changed to handle these legacy pickles
# class TestLegacySupport(unittest.TestCase):
class LegacySupport(object):

    _multiprocess_can_split_ = True

    @classmethod
    def setUpClass(cls):
        if compat.PY3:
            raise nose.SkipTest("not compatible with Python >= 3")

        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'frame.pickle')

        with open(filepath, 'rb') as f:
            cls.frame = pickle.load(f)

        filepath = os.path.join(pth, 'data', 'series.pickle')
        with open(filepath, 'rb') as f:
            cls.series = pickle.load(f)

    def test_pass_offset_warn(self):
        buf = StringIO()

        sys.stderr = buf
        DatetimeIndex(start='1/1/2000', periods=10, offset='H')
        sys.stderr = sys.__stderr__

    def test_unpickle_legacy_frame(self):
        dtindex = DatetimeIndex(start='1/3/2005', end='1/14/2005',
                                freq=BDay(1))

        unpickled = self.frame

        self.assertEquals(type(unpickled.index), DatetimeIndex)
        self.assertEquals(len(unpickled), 10)
        self.assert_((unpickled.columns == Int64Index(np.arange(5))).all())
        self.assert_((unpickled.index == dtindex).all())
        self.assertEquals(unpickled.index.offset, BDay(1, normalize=True))

    def test_unpickle_legacy_series(self):
        from pandas.core.datetools import BDay

        unpickled = self.series

        dtindex = DatetimeIndex(start='1/3/2005', end='1/14/2005',
                                freq=BDay(1))

        self.assertEquals(type(unpickled.index), DatetimeIndex)
        self.assertEquals(len(unpickled), 10)
        self.assert_((unpickled.index == dtindex).all())
        self.assertEquals(unpickled.index.offset, BDay(1, normalize=True))

    def test_unpickle_legacy_len0_daterange(self):
        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'series_daterange0.pickle')

        result = pd.read_pickle(filepath)

        ex_index = DatetimeIndex([], freq='B')

        self.assert_(result.index.equals(ex_index))
        tm.assert_isinstance(result.index.freq, offsets.BDay)
        self.assert_(len(result) == 0)

    def test_arithmetic_interaction(self):
        index = self.frame.index
        obj_index = index.asobject

        dseries = Series(rand(len(index)), index=index)
        oseries = Series(dseries.values, index=obj_index)

        result = dseries + oseries
        expected = dseries * 2
        tm.assert_isinstance(result.index, DatetimeIndex)
        assert_series_equal(result, expected)

        result = dseries + oseries[:5]
        expected = dseries + dseries[:5]
        tm.assert_isinstance(result.index, DatetimeIndex)
        assert_series_equal(result, expected)

    def test_join_interaction(self):
        index = self.frame.index
        obj_index = index.asobject

        def _check_join(left, right, how='inner'):
            ra, rb, rc = left.join(right, how=how, return_indexers=True)
            ea, eb, ec = left.join(DatetimeIndex(right), how=how,
                                   return_indexers=True)

            tm.assert_isinstance(ra, DatetimeIndex)
            self.assert_(ra.equals(ea))

            assert_almost_equal(rb, eb)
            assert_almost_equal(rc, ec)

        _check_join(index[:15], obj_index[5:], how='inner')
        _check_join(index[:15], obj_index[5:], how='outer')
        _check_join(index[:15], obj_index[5:], how='right')
        _check_join(index[:15], obj_index[5:], how='left')

    def test_join_nonunique(self):
        idx1 = to_datetime(['2012-11-06 16:00:11.477563',
                            '2012-11-06 16:00:11.477563'])
        idx2 = to_datetime(['2012-11-06 15:11:09.006507',
                            '2012-11-06 15:11:09.006507'])
        rs = idx1.join(idx2, how='outer')
        self.assert_(rs.is_monotonic)

    def test_unpickle_daterange(self):
        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'daterange_073.pickle')

        rng = read_pickle(filepath)
        tm.assert_isinstance(rng[0], datetime)
        tm.assert_isinstance(rng.offset, offsets.BDay)
        self.assert_(rng.values.dtype == object)

    def test_setops(self):
        index = self.frame.index
        obj_index = index.asobject

        result = index[:5].union(obj_index[5:])
        expected = index
        tm.assert_isinstance(result, DatetimeIndex)
        self.assert_(result.equals(expected))

        result = index[:10].intersection(obj_index[5:])
        expected = index[5:10]
        tm.assert_isinstance(result, DatetimeIndex)
        self.assert_(result.equals(expected))

        result = index[:10] - obj_index[5:]
        expected = index[:5]
        tm.assert_isinstance(result, DatetimeIndex)
        self.assert_(result.equals(expected))

    def test_index_conversion(self):
        index = self.frame.index
        obj_index = index.asobject

        conv = DatetimeIndex(obj_index)
        self.assert_(conv.equals(index))

        self.assertRaises(ValueError, DatetimeIndex, ['a', 'b', 'c', 'd'])

    def test_tolist(self):
        rng = date_range('1/1/2000', periods=10)

        result = rng.tolist()
        tm.assert_isinstance(result[0], Timestamp)

    def test_object_convert_fail(self):
        idx = DatetimeIndex([NaT])
        self.assertRaises(ValueError, idx.astype, 'O')

    def test_setops_conversion_fail(self):
        index = self.frame.index

        right = Index(['a', 'b', 'c', 'd'])

        result = index.union(right)
        expected = Index(np.concatenate([index.asobject, right]))
        self.assert_(result.equals(expected))

        result = index.intersection(right)
        expected = Index([])
        self.assert_(result.equals(expected))

    def test_legacy_time_rules(self):
        rules = [('WEEKDAY', 'B'),
                 ('EOM', 'BM'),
                 ('W@MON', 'W-MON'), ('W@TUE', 'W-TUE'), ('W@WED', 'W-WED'),
                 ('W@THU', 'W-THU'), ('W@FRI', 'W-FRI'),
                 ('Q@JAN', 'BQ-JAN'), ('Q@FEB', 'BQ-FEB'), ('Q@MAR', 'BQ-MAR'),
                 ('A@JAN', 'BA-JAN'), ('A@FEB', 'BA-FEB'), ('A@MAR', 'BA-MAR'),
                 ('A@APR', 'BA-APR'), ('A@MAY', 'BA-MAY'), ('A@JUN', 'BA-JUN'),
                 ('A@JUL', 'BA-JUL'), ('A@AUG', 'BA-AUG'), ('A@SEP', 'BA-SEP'),
                 ('A@OCT', 'BA-OCT'), ('A@NOV', 'BA-NOV'), ('A@DEC', 'BA-DEC'),
                 ('WOM@1FRI', 'WOM-1FRI'), ('WOM@2FRI', 'WOM-2FRI'),
                 ('WOM@3FRI', 'WOM-3FRI'), ('WOM@4FRI', 'WOM-4FRI')]

        start, end = '1/1/2000', '1/1/2010'

        for old_freq, new_freq in rules:
            old_rng = date_range(start, end, freq=old_freq)
            new_rng = date_range(start, end, freq=new_freq)
            self.assert_(old_rng.equals(new_rng))

            # test get_legacy_offset_name
            offset = datetools.get_offset(new_freq)
            old_name = datetools.get_legacy_offset_name(offset)
            self.assertEquals(old_name, old_freq)

    def test_ms_vs_MS(self):
        left = datetools.get_offset('ms')
        right = datetools.get_offset('MS')
        self.assert_(left == datetools.Milli())
        self.assert_(right == datetools.MonthBegin())

    def test_rule_aliases(self):
        rule = datetools.to_offset('10us')
        self.assert_(rule == datetools.Micro(10))


class TestLegacyCompat(unittest.TestCase):

    def setUp(self):
        # suppress deprecation warnings
        sys.stderr = StringIO()

    def test_inferTimeRule(self):
        from pandas.tseries.frequencies import inferTimeRule

        index1 = [datetime(2010, 1, 29, 0, 0),
                  datetime(2010, 2, 26, 0, 0),
                  datetime(2010, 3, 31, 0, 0)]

        index2 = [datetime(2010, 3, 26, 0, 0),
                  datetime(2010, 3, 29, 0, 0),
                  datetime(2010, 3, 30, 0, 0)]

        index3 = [datetime(2010, 3, 26, 0, 0),
                  datetime(2010, 3, 27, 0, 0),
                  datetime(2010, 3, 29, 0, 0)]

        # LEGACY
        assert inferTimeRule(index1) == 'EOM'
        assert inferTimeRule(index2) == 'WEEKDAY'

        self.assertRaises(Exception, inferTimeRule, index1[:2])
        self.assertRaises(Exception, inferTimeRule, index3)

    def test_time_rule(self):
        result = DateRange('1/1/2000', '1/30/2000', time_rule='WEEKDAY')
        result2 = DateRange('1/1/2000', '1/30/2000', timeRule='WEEKDAY')
        expected = date_range('1/1/2000', '1/30/2000', freq='B')

        self.assert_(result.equals(expected))
        self.assert_(result2.equals(expected))

    def tearDown(self):
        sys.stderr = sys.__stderr__


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
