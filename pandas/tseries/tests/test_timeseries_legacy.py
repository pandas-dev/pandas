# pylint: disable-msg=E1101,W0612
from datetime import datetime
import sys
import os
import nose
import numpy as np

from pandas import (Index, Series, date_range, Timestamp,
                    DatetimeIndex, Int64Index, to_datetime)

from pandas.tseries.frequencies import get_offset, to_offset
from pandas.tseries.offsets import BDay, Micro, Milli, MonthBegin
import pandas as pd

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

from pandas.compat import StringIO, cPickle as pickle
from pandas import read_pickle
from numpy.random import rand
import pandas.compat as compat

randn = np.random.randn


# Unfortunately, too much has changed to handle these legacy pickles
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

        self.assertEqual(type(unpickled.index), DatetimeIndex)
        self.assertEqual(len(unpickled), 10)
        self.assertTrue((unpickled.columns == Int64Index(np.arange(5))).all())
        self.assertTrue((unpickled.index == dtindex).all())
        self.assertEqual(unpickled.index.offset, BDay(1, normalize=True))

    def test_unpickle_legacy_series(self):
        unpickled = self.series

        dtindex = DatetimeIndex(start='1/3/2005', end='1/14/2005',
                                freq=BDay(1))

        self.assertEqual(type(unpickled.index), DatetimeIndex)
        self.assertEqual(len(unpickled), 10)
        self.assertTrue((unpickled.index == dtindex).all())
        self.assertEqual(unpickled.index.offset, BDay(1, normalize=True))

    def test_unpickle_legacy_len0_daterange(self):
        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'series_daterange0.pickle')

        result = pd.read_pickle(filepath)

        ex_index = DatetimeIndex([], freq='B')

        self.assert_index_equal(result.index, ex_index)
        tm.assertIsInstance(result.index.freq, BDay)
        self.assertEqual(len(result), 0)

    def test_arithmetic_interaction(self):
        index = self.frame.index
        obj_index = index.asobject

        dseries = Series(rand(len(index)), index=index)
        oseries = Series(dseries.values, index=obj_index)

        result = dseries + oseries
        expected = dseries * 2
        tm.assertIsInstance(result.index, DatetimeIndex)
        assert_series_equal(result, expected)

        result = dseries + oseries[:5]
        expected = dseries + dseries[:5]
        tm.assertIsInstance(result.index, DatetimeIndex)
        assert_series_equal(result, expected)

    def test_join_interaction(self):
        index = self.frame.index
        obj_index = index.asobject

        def _check_join(left, right, how='inner'):
            ra, rb, rc = left.join(right, how=how, return_indexers=True)
            ea, eb, ec = left.join(DatetimeIndex(right), how=how,
                                   return_indexers=True)

            tm.assertIsInstance(ra, DatetimeIndex)
            self.assert_index_equal(ra, ea)

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
        self.assertTrue(rs.is_monotonic)

    def test_unpickle_daterange(self):
        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'daterange_073.pickle')

        rng = read_pickle(filepath)
        tm.assertIsInstance(rng[0], datetime)
        tm.assertIsInstance(rng.offset, BDay)
        self.assertEqual(rng.values.dtype, object)

    def test_setops(self):
        index = self.frame.index
        obj_index = index.asobject

        result = index[:5].union(obj_index[5:])
        expected = index
        tm.assertIsInstance(result, DatetimeIndex)
        self.assert_index_equal(result, expected)

        result = index[:10].intersection(obj_index[5:])
        expected = index[5:10]
        tm.assertIsInstance(result, DatetimeIndex)
        self.assert_index_equal(result, expected)

        result = index[:10] - obj_index[5:]
        expected = index[:5]
        tm.assertIsInstance(result, DatetimeIndex)
        self.assert_index_equal(result, expected)

    def test_index_conversion(self):
        index = self.frame.index
        obj_index = index.asobject

        conv = DatetimeIndex(obj_index)
        self.assert_index_equal(conv, index)

        self.assertRaises(ValueError, DatetimeIndex, ['a', 'b', 'c', 'd'])

    def test_tolist(self):
        rng = date_range('1/1/2000', periods=10)

        result = rng.tolist()
        tm.assertIsInstance(result[0], Timestamp)

    def test_object_convert_fail(self):
        idx = DatetimeIndex([np.NaT])
        self.assertRaises(ValueError, idx.astype, 'O')

    def test_setops_conversion_fail(self):
        index = self.frame.index

        right = Index(['a', 'b', 'c', 'd'])

        result = index.union(right)
        expected = Index(np.concatenate([index.asobject, right]))
        self.assert_index_equal(result, expected)

        result = index.intersection(right)
        expected = Index([])
        self.assert_index_equal(result, expected)

    def test_legacy_time_rules(self):
        rules = [('WEEKDAY', 'B'), ('EOM', 'BM'), ('W@MON', 'W-MON'),
                 ('W@TUE', 'W-TUE'), ('W@WED', 'W-WED'), ('W@THU', 'W-THU'),
                 ('W@FRI', 'W-FRI'), ('Q@JAN', 'BQ-JAN'), ('Q@FEB', 'BQ-FEB'),
                 ('Q@MAR', 'BQ-MAR'), ('A@JAN', 'BA-JAN'), ('A@FEB', 'BA-FEB'),
                 ('A@MAR', 'BA-MAR'), ('A@APR', 'BA-APR'), ('A@MAY', 'BA-MAY'),
                 ('A@JUN', 'BA-JUN'), ('A@JUL', 'BA-JUL'), ('A@AUG', 'BA-AUG'),
                 ('A@SEP', 'BA-SEP'), ('A@OCT', 'BA-OCT'), ('A@NOV', 'BA-NOV'),
                 ('A@DEC', 'BA-DEC'), ('WOM@1FRI', 'WOM-1FRI'),
                 ('WOM@2FRI', 'WOM-2FRI'), ('WOM@3FRI', 'WOM-3FRI'),
                 ('WOM@4FRI', 'WOM-4FRI')]

        start, end = '1/1/2000', '1/1/2010'

        for old_freq, new_freq in rules:
            old_rng = date_range(start, end, freq=old_freq)
            new_rng = date_range(start, end, freq=new_freq)
            self.assert_index_equal(old_rng, new_rng)

    def test_ms_vs_MS(self):
        left = get_offset('ms')
        right = get_offset('MS')
        self.assertEqual(left, Milli())
        self.assertEqual(right, MonthBegin())

    def test_rule_aliases(self):
        rule = to_offset('10us')
        self.assertEqual(rule, Micro(10))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
