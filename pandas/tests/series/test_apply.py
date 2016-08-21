# coding=utf-8
# pylint: disable-msg=E1101,W0612

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, isnull)
from pandas.compat import lrange
from pandas import compat
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesApply(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

    def test_apply(self):
        with np.errstate(all='ignore'):
            assert_series_equal(self.ts.apply(np.sqrt), np.sqrt(self.ts))

            # elementwise-apply
            import math
            assert_series_equal(self.ts.apply(math.exp), np.exp(self.ts))

            # how to handle Series result, #2316
            result = self.ts.apply(lambda x: Series(
                [x, x ** 2], index=['x', 'x^2']))
            expected = DataFrame({'x': self.ts, 'x^2': self.ts ** 2})
            tm.assert_frame_equal(result, expected)

        # empty series
        s = Series(dtype=object, name='foo', index=pd.Index([], name='bar'))
        rs = s.apply(lambda x: x)
        tm.assert_series_equal(s, rs)
        # check all metadata (GH 9322)
        self.assertIsNot(s, rs)
        self.assertIs(s.index, rs.index)
        self.assertEqual(s.dtype, rs.dtype)
        self.assertEqual(s.name, rs.name)

        # index but no data
        s = Series(index=[1, 2, 3])
        rs = s.apply(lambda x: x)
        tm.assert_series_equal(s, rs)

    def test_apply_same_length_inference_bug(self):
        s = Series([1, 2])
        f = lambda x: (x, x + 1)

        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

        s = Series([1, 2, 3])
        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

    def test_apply_dont_convert_dtype(self):
        s = Series(np.random.randn(10))

        f = lambda x: x if x > 0 else np.nan
        result = s.apply(f, convert_dtype=False)
        self.assertEqual(result.dtype, object)

    def test_apply_args(self):
        s = Series(['foo,bar'])

        result = s.apply(str.split, args=(',', ))
        self.assertEqual(result[0], ['foo', 'bar'])
        tm.assertIsInstance(result[0], list)

    def test_apply_box(self):
        # ufunc will not be boxed. Same test cases as the test_map_box
        vals = [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'datetime64[ns]')
        # boxed value must be Timestamp instance
        res = s.apply(lambda x: '{0}_{1}_{2}'.format(x.__class__.__name__,
                                                     x.day, x.tz))
        exp = pd.Series(['Timestamp_1_None', 'Timestamp_2_None'])
        tm.assert_series_equal(res, exp)

        vals = [pd.Timestamp('2011-01-01', tz='US/Eastern'),
                pd.Timestamp('2011-01-02', tz='US/Eastern')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'datetime64[ns, US/Eastern]')
        res = s.apply(lambda x: '{0}_{1}_{2}'.format(x.__class__.__name__,
                                                     x.day, x.tz))
        exp = pd.Series(['Timestamp_1_US/Eastern', 'Timestamp_2_US/Eastern'])
        tm.assert_series_equal(res, exp)

        # timedelta
        vals = [pd.Timedelta('1 days'), pd.Timedelta('2 days')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'timedelta64[ns]')
        res = s.apply(lambda x: '{0}_{1}'.format(x.__class__.__name__, x.days))
        exp = pd.Series(['Timedelta_1', 'Timedelta_2'])
        tm.assert_series_equal(res, exp)

        # period (object dtype, not boxed)
        vals = [pd.Period('2011-01-01', freq='M'),
                pd.Period('2011-01-02', freq='M')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'object')
        res = s.apply(lambda x: '{0}_{1}'.format(x.__class__.__name__,
                                                 x.freqstr))
        exp = pd.Series(['Period_M', 'Period_M'])
        tm.assert_series_equal(res, exp)

    def test_apply_datetimetz(self):
        values = pd.date_range('2011-01-01', '2011-01-02',
                               freq='H').tz_localize('Asia/Tokyo')
        s = pd.Series(values, name='XX')

        result = s.apply(lambda x: x + pd.offsets.Day())
        exp_values = pd.date_range('2011-01-02', '2011-01-03',
                                   freq='H').tz_localize('Asia/Tokyo')
        exp = pd.Series(exp_values, name='XX')
        tm.assert_series_equal(result, exp)

        # change dtype
        result = s.apply(lambda x: x.hour)
        exp = pd.Series(list(range(24)) + [0], name='XX', dtype=np.int32)
        tm.assert_series_equal(result, exp)

        # not vectorized
        def f(x):
            if not isinstance(x, pd.Timestamp):
                raise ValueError
            return str(x.tz)

        result = s.map(f)
        exp = pd.Series(['Asia/Tokyo'] * 25, name='XX')
        tm.assert_series_equal(result, exp)


class TestSeriesMap(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

    def test_map(self):
        index, data = tm.getMixedTypeDict()

        source = Series(data['B'], index=data['C'])
        target = Series(data['C'][:4], index=data['D'][:4])

        merged = target.map(source)

        for k, v in compat.iteritems(merged):
            self.assertEqual(v, source[target[k]])

        # input could be a dict
        merged = target.map(source.to_dict())

        for k, v in compat.iteritems(merged):
            self.assertEqual(v, source[target[k]])

        # function
        result = self.ts.map(lambda x: x * 2)
        self.assert_series_equal(result, self.ts * 2)

        # GH 10324
        a = Series([1, 2, 3, 4])
        b = Series(["even", "odd", "even", "odd"], dtype="category")
        c = Series(["even", "odd", "even", "odd"])

        exp = Series(["odd", "even", "odd", np.nan], dtype="category")
        self.assert_series_equal(a.map(b), exp)
        exp = Series(["odd", "even", "odd", np.nan])
        self.assert_series_equal(a.map(c), exp)

        a = Series(['a', 'b', 'c', 'd'])
        b = Series([1, 2, 3, 4],
                   index=pd.CategoricalIndex(['b', 'c', 'd', 'e']))
        c = Series([1, 2, 3, 4], index=Index(['b', 'c', 'd', 'e']))

        exp = Series([np.nan, 1, 2, 3])
        self.assert_series_equal(a.map(b), exp)
        exp = Series([np.nan, 1, 2, 3])
        self.assert_series_equal(a.map(c), exp)

        a = Series(['a', 'b', 'c', 'd'])
        b = Series(['B', 'C', 'D', 'E'], dtype='category',
                   index=pd.CategoricalIndex(['b', 'c', 'd', 'e']))
        c = Series(['B', 'C', 'D', 'E'], index=Index(['b', 'c', 'd', 'e']))

        exp = Series(pd.Categorical([np.nan, 'B', 'C', 'D'],
                                    categories=['B', 'C', 'D', 'E']))
        self.assert_series_equal(a.map(b), exp)
        exp = Series([np.nan, 'B', 'C', 'D'])
        self.assert_series_equal(a.map(c), exp)

    def test_map_compat(self):
        # related GH 8024
        s = Series([True, True, False], index=[1, 2, 3])
        result = s.map({True: 'foo', False: 'bar'})
        expected = Series(['foo', 'foo', 'bar'], index=[1, 2, 3])
        assert_series_equal(result, expected)

    def test_map_int(self):
        left = Series({'a': 1., 'b': 2., 'c': 3., 'd': 4})
        right = Series({1: 11, 2: 22, 3: 33})

        self.assertEqual(left.dtype, np.float_)
        self.assertTrue(issubclass(right.dtype.type, np.integer))

        merged = left.map(right)
        self.assertEqual(merged.dtype, np.float_)
        self.assertTrue(isnull(merged['d']))
        self.assertTrue(not isnull(merged['c']))

    def test_map_type_inference(self):
        s = Series(lrange(3))
        s2 = s.map(lambda x: np.where(x == 0, 0, 1))
        self.assertTrue(issubclass(s2.dtype.type, np.integer))

    def test_map_decimal(self):
        from decimal import Decimal

        result = self.series.map(lambda x: Decimal(str(x)))
        self.assertEqual(result.dtype, np.object_)
        tm.assertIsInstance(result[0], Decimal)

    def test_map_na_exclusion(self):
        s = Series([1.5, np.nan, 3, np.nan, 5])

        result = s.map(lambda x: x * 2, na_action='ignore')
        exp = s * 2
        assert_series_equal(result, exp)

    def test_map_dict_with_tuple_keys(self):
        """
        Due to new MultiIndex-ing behaviour in v0.14.0,
        dicts with tuple keys passed to map were being
        converted to a multi-index, preventing tuple values
        from being mapped properly.
        """
        df = pd.DataFrame({'a': [(1, ), (2, ), (3, 4), (5, 6)]})
        label_mappings = {(1, ): 'A', (2, ): 'B', (3, 4): 'A', (5, 6): 'B'}
        df['labels'] = df['a'].map(label_mappings)
        df['expected_labels'] = pd.Series(['A', 'B', 'A', 'B'], index=df.index)
        # All labels should be filled now
        tm.assert_series_equal(df['labels'], df['expected_labels'],
                               check_names=False)

    def test_map_box(self):
        vals = [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'datetime64[ns]')
        # boxed value must be Timestamp instance
        res = s.map(lambda x: '{0}_{1}_{2}'.format(x.__class__.__name__,
                                                   x.day, x.tz))
        exp = pd.Series(['Timestamp_1_None', 'Timestamp_2_None'])
        tm.assert_series_equal(res, exp)

        vals = [pd.Timestamp('2011-01-01', tz='US/Eastern'),
                pd.Timestamp('2011-01-02', tz='US/Eastern')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'datetime64[ns, US/Eastern]')
        res = s.map(lambda x: '{0}_{1}_{2}'.format(x.__class__.__name__,
                                                   x.day, x.tz))
        exp = pd.Series(['Timestamp_1_US/Eastern', 'Timestamp_2_US/Eastern'])
        tm.assert_series_equal(res, exp)

        # timedelta
        vals = [pd.Timedelta('1 days'), pd.Timedelta('2 days')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'timedelta64[ns]')
        res = s.map(lambda x: '{0}_{1}'.format(x.__class__.__name__, x.days))
        exp = pd.Series(['Timedelta_1', 'Timedelta_2'])
        tm.assert_series_equal(res, exp)

        # period (object dtype, not boxed)
        vals = [pd.Period('2011-01-01', freq='M'),
                pd.Period('2011-01-02', freq='M')]
        s = pd.Series(vals)
        self.assertEqual(s.dtype, 'object')
        res = s.map(lambda x: '{0}_{1}'.format(x.__class__.__name__,
                                               x.freqstr))
        exp = pd.Series(['Period_M', 'Period_M'])
        tm.assert_series_equal(res, exp)

    def test_map_categorical(self):
        values = pd.Categorical(list('ABBABCD'), categories=list('DCBA'),
                                ordered=True)
        s = pd.Series(values, name='XX', index=list('abcdefg'))

        result = s.map(lambda x: x.lower())
        exp_values = pd.Categorical(list('abbabcd'), categories=list('dcba'),
                                    ordered=True)
        exp = pd.Series(exp_values, name='XX', index=list('abcdefg'))
        tm.assert_series_equal(result, exp)
        tm.assert_categorical_equal(result.values, exp_values)

        result = s.map(lambda x: 'A')
        exp = pd.Series(['A'] * 7, name='XX', index=list('abcdefg'))
        tm.assert_series_equal(result, exp)
        self.assertEqual(result.dtype, np.object)

        with tm.assertRaises(NotImplementedError):
            s.map(lambda x: x, na_action='ignore')

    def test_map_datetimetz(self):
        values = pd.date_range('2011-01-01', '2011-01-02',
                               freq='H').tz_localize('Asia/Tokyo')
        s = pd.Series(values, name='XX')

        # keep tz
        result = s.map(lambda x: x + pd.offsets.Day())
        exp_values = pd.date_range('2011-01-02', '2011-01-03',
                                   freq='H').tz_localize('Asia/Tokyo')
        exp = pd.Series(exp_values, name='XX')
        tm.assert_series_equal(result, exp)

        # change dtype
        result = s.map(lambda x: x.hour)
        exp = pd.Series(list(range(24)) + [0], name='XX', dtype=np.int32)
        tm.assert_series_equal(result, exp)

        with tm.assertRaises(NotImplementedError):
            s.map(lambda x: x, na_action='ignore')

        # not vectorized
        def f(x):
            if not isinstance(x, pd.Timestamp):
                raise ValueError
            return str(x.tz)

        result = s.map(f)
        exp = pd.Series(['Asia/Tokyo'] * 25, name='XX')
        tm.assert_series_equal(result, exp)
