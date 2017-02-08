import numpy as np

import pandas as pd
import pandas.util.testing as tm
import pandas.tseries.period as period
from pandas import period_range, PeriodIndex, Index, date_range


def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


class TestPeriodIndex(tm.TestCase):

    def setUp(self):
        pass

    def test_joins(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        for kind in ['inner', 'outer', 'left', 'right']:
            joined = index.join(index[:-5], how=kind)

            tm.assertIsInstance(joined, PeriodIndex)
            self.assertEqual(joined.freq, index.freq)

    def test_join_self(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        for kind in ['inner', 'outer', 'left', 'right']:
            res = index.join(index, how=kind)
            self.assertIs(index, res)

    def test_join_does_not_recur(self):
        df = tm.makeCustomDataframe(
            3, 2, data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type='p', r_idx_type='dt')
        s = df.iloc[:2, 0]

        res = s.index.join(df.columns, how='outer')
        expected = Index([s.index[0], s.index[1],
                          df.columns[0], df.columns[1]], object)
        tm.assert_index_equal(res, expected)

    def test_union(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].union(index[10:])
        tm.assert_index_equal(result, index)

        # not in order
        result = _permute(index[:-5]).union(_permute(index[10:]))
        tm.assert_index_equal(result, index)

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.union(index2)

        msg = 'can only call with other PeriodIndex-ed objects'
        with tm.assertRaisesRegexp(ValueError, msg):
            index.join(index.to_timestamp())

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.join(index3)

    def test_union_dataframe_index(self):
        rng1 = pd.period_range('1/1/1999', '1/1/2012', freq='M')
        s1 = pd.Series(np.random.randn(len(rng1)), rng1)

        rng2 = pd.period_range('1/1/1980', '12/1/2001', freq='M')
        s2 = pd.Series(np.random.randn(len(rng2)), rng2)
        df = pd.DataFrame({'s1': s1, 's2': s2})

        exp = pd.period_range('1/1/1980', '1/1/2012', freq='M')
        self.assert_index_equal(df.index, exp)

    def test_intersection(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].intersection(index[10:])
        tm.assert_index_equal(result, index[10:-5])

        # not in order
        left = _permute(index[:-5])
        right = _permute(index[10:])
        result = left.intersection(right).sort_values()
        tm.assert_index_equal(result, index[10:-5])

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.intersection(index2)

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.intersection(index3)

    def test_intersection_cases(self):
        base = period_range('6/1/2000', '6/30/2000', freq='D', name='idx')

        # if target has the same name, it is preserved
        rng2 = period_range('5/15/2000', '6/20/2000', freq='D', name='idx')
        expected2 = period_range('6/1/2000', '6/20/2000', freq='D',
                                 name='idx')

        # if target name is different, it will be reset
        rng3 = period_range('5/15/2000', '6/20/2000', freq='D', name='other')
        expected3 = period_range('6/1/2000', '6/20/2000', freq='D',
                                 name=None)

        rng4 = period_range('7/1/2000', '7/31/2000', freq='D', name='idx')
        expected4 = PeriodIndex([], name='idx', freq='D')

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        # non-monotonic
        base = PeriodIndex(['2011-01-05', '2011-01-04', '2011-01-02',
                            '2011-01-03'], freq='D', name='idx')

        rng2 = PeriodIndex(['2011-01-04', '2011-01-02',
                            '2011-02-02', '2011-02-03'],
                           freq='D', name='idx')
        expected2 = PeriodIndex(['2011-01-04', '2011-01-02'], freq='D',
                                name='idx')

        rng3 = PeriodIndex(['2011-01-04', '2011-01-02', '2011-02-02',
                            '2011-02-03'],
                           freq='D', name='other')
        expected3 = PeriodIndex(['2011-01-04', '2011-01-02'], freq='D',
                                name=None)

        rng4 = period_range('7/1/2000', '7/31/2000', freq='D', name='idx')
        expected4 = PeriodIndex([], freq='D', name='idx')

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, 'D')

        # empty same freq
        rng = date_range('6/1/2000', '6/15/2000', freq='T')
        result = rng[0:0].intersection(rng)
        self.assertEqual(len(result), 0)

        result = rng.intersection(rng[0:0])
        self.assertEqual(len(result), 0)
