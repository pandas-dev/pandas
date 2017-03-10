""" test scalar indexing, including at and iat """

import numpy as np

from pandas import (Series, DataFrame, Timestamp,
                    Timedelta, date_range)
from pandas.util import testing as tm
from pandas.tests.indexing.common import Base


class TestScalar(Base, tm.TestCase):

    def test_at_and_iat_get(self):
        def _check(f, func, values=False):

            if f is not None:
                indicies = self.generate_indices(f, values)
                for i in indicies:
                    result = getattr(f, func)[i]
                    expected = self.get_value(f, i, values)
                    tm.assert_almost_equal(result, expected)

        for o in self._objs:

            d = getattr(self, o)

            # iat
            for f in [d['ints'], d['uints']]:
                _check(f, 'iat', values=True)

            for f in [d['labels'], d['ts'], d['floats']]:
                if f is not None:
                    self.assertRaises(ValueError, self.check_values, f, 'iat')

            # at
            for f in [d['ints'], d['uints'], d['labels'],
                      d['ts'], d['floats']]:
                _check(f, 'at')

    def test_at_and_iat_set(self):
        def _check(f, func, values=False):

            if f is not None:
                indicies = self.generate_indices(f, values)
                for i in indicies:
                    getattr(f, func)[i] = 1
                    expected = self.get_value(f, i, values)
                    tm.assert_almost_equal(expected, 1)

        for t in self._objs:

            d = getattr(self, t)

            # iat
            for f in [d['ints'], d['uints']]:
                _check(f, 'iat', values=True)

            for f in [d['labels'], d['ts'], d['floats']]:
                if f is not None:
                    self.assertRaises(ValueError, _check, f, 'iat')

            # at
            for f in [d['ints'], d['uints'], d['labels'],
                      d['ts'], d['floats']]:
                _check(f, 'at')

    def test_at_iat_coercion(self):

        # as timestamp is not a tuple!
        dates = date_range('1/1/2000', periods=8)
        df = DataFrame(np.random.randn(8, 4),
                       index=dates,
                       columns=['A', 'B', 'C', 'D'])
        s = df['A']

        result = s.at[dates[5]]
        xp = s.values[5]
        self.assertEqual(result, xp)

        # GH 7729
        # make sure we are boxing the returns
        s = Series(['2014-01-01', '2014-02-02'], dtype='datetime64[ns]')
        expected = Timestamp('2014-02-02')

        for r in [lambda: s.iat[1], lambda: s.iloc[1]]:
            result = r()
            self.assertEqual(result, expected)

        s = Series(['1 days', '2 days'], dtype='timedelta64[ns]')
        expected = Timedelta('2 days')

        for r in [lambda: s.iat[1], lambda: s.iloc[1]]:
            result = r()
            self.assertEqual(result, expected)

    def test_iat_invalid_args(self):
        pass

    def test_imethods_with_dups(self):

        # GH6493
        # iat/iloc with dups

        s = Series(range(5), index=[1, 1, 2, 2, 3], dtype='int64')
        result = s.iloc[2]
        self.assertEqual(result, 2)
        result = s.iat[2]
        self.assertEqual(result, 2)

        self.assertRaises(IndexError, lambda: s.iat[10])
        self.assertRaises(IndexError, lambda: s.iat[-10])

        result = s.iloc[[2, 3]]
        expected = Series([2, 3], [2, 2], dtype='int64')
        tm.assert_series_equal(result, expected)

        df = s.to_frame()
        result = df.iloc[2]
        expected = Series(2, index=[0], name=2)
        tm.assert_series_equal(result, expected)

        result = df.iat[2, 0]
        expected = 2
        self.assertEqual(result, 2)

    def test_at_to_fail(self):
        # at should not fallback
        # GH 7814
        s = Series([1, 2, 3], index=list('abc'))
        result = s.at['a']
        self.assertEqual(result, 1)
        self.assertRaises(ValueError, lambda: s.at[0])

        df = DataFrame({'A': [1, 2, 3]}, index=list('abc'))
        result = df.at['a', 'A']
        self.assertEqual(result, 1)
        self.assertRaises(ValueError, lambda: df.at['a', 0])

        s = Series([1, 2, 3], index=[3, 2, 1])
        result = s.at[1]
        self.assertEqual(result, 3)
        self.assertRaises(ValueError, lambda: s.at['a'])

        df = DataFrame({0: [1, 2, 3]}, index=[3, 2, 1])
        result = df.at[1, 0]
        self.assertEqual(result, 3)
        self.assertRaises(ValueError, lambda: df.at['a', 0])

        # GH 13822, incorrect error string with non-unique columns when missing
        # column is accessed
        df = DataFrame({'x': [1.], 'y': [2.], 'z': [3.]})
        df.columns = ['x', 'x', 'z']

        # Check that we get the correct value in the KeyError
        self.assertRaisesRegexp(KeyError, r"\['y'\] not in index",
                                lambda: df[['x', 'y', 'z']])
