# -*- coding: utf-8 -*-

from __future__ import print_function

from pandas.compat import range, lrange
import numpy as np

from pandas import DataFrame, Series

from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


# Column add, remove, delete.


class TestDataFrameMutateColumns(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_assign(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        original = df.copy()
        result = df.assign(C=df.B / df.A)
        expected = df.copy()
        expected['C'] = [4, 2.5, 2]
        assert_frame_equal(result, expected)

        # lambda syntax
        result = df.assign(C=lambda x: x.B / x.A)
        assert_frame_equal(result, expected)

        # original is unmodified
        assert_frame_equal(df, original)

        # Non-Series array-like
        result = df.assign(C=[4, 2.5, 2])
        assert_frame_equal(result, expected)
        # original is unmodified
        assert_frame_equal(df, original)

        result = df.assign(B=df.B / df.A)
        expected = expected.drop('B', axis=1).rename(columns={'C': 'B'})
        assert_frame_equal(result, expected)

        # overwrite
        result = df.assign(A=df.A + df.B)
        expected = df.copy()
        expected['A'] = [5, 7, 9]
        assert_frame_equal(result, expected)

        # lambda
        result = df.assign(A=lambda x: x.A + x.B)
        assert_frame_equal(result, expected)

    def test_assign_multiple(self):
        df = DataFrame([[1, 4], [2, 5], [3, 6]], columns=['A', 'B'])
        result = df.assign(C=[7, 8, 9], D=df.A, E=lambda x: x.B)
        expected = DataFrame([[1, 4, 7, 1, 4], [2, 5, 8, 2, 5],
                              [3, 6, 9, 3, 6]], columns=list('ABCDE'))
        assert_frame_equal(result, expected)

    def test_assign_alphabetical(self):
        # GH 9818
        df = DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        result = df.assign(D=df.A + df.B, C=df.A - df.B)
        expected = DataFrame([[1, 2, -1, 3], [3, 4, -1, 7]],
                             columns=list('ABCD'))
        assert_frame_equal(result, expected)
        result = df.assign(C=df.A - df.B, D=df.A + df.B)
        assert_frame_equal(result, expected)

    def test_assign_bad(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        # non-keyword argument
        with tm.assertRaises(TypeError):
            df.assign(lambda x: x.A)
        with tm.assertRaises(AttributeError):
            df.assign(C=df.A, D=df.A + df.C)
        with tm.assertRaises(KeyError):
            df.assign(C=lambda df: df.A, D=lambda df: df['A'] + df['C'])
        with tm.assertRaises(KeyError):
            df.assign(C=df.A, D=lambda x: x['A'] + x['C'])

    def test_insert_error_msmgs(self):

        # GH 7432
        df = DataFrame({'foo': ['a', 'b', 'c'], 'bar': [
                       1, 2, 3], 'baz': ['d', 'e', 'f']}).set_index('foo')
        s = DataFrame({'foo': ['a', 'b', 'c', 'a'], 'fiz': [
                      'g', 'h', 'i', 'j']}).set_index('foo')
        msg = 'cannot reindex from a duplicate axis'
        with assertRaisesRegexp(ValueError, msg):
            df['newcol'] = s

        # GH 4107, more descriptive error message
        df = DataFrame(np.random.randint(0, 2, (4, 4)),
                       columns=['a', 'b', 'c', 'd'])

        msg = 'incompatible index of inserted column with frame index'
        with assertRaisesRegexp(TypeError, msg):
            df['gr'] = df.groupby(['b', 'c']).count()

    def test_insert_benchmark(self):
        # from the vb_suite/frame_methods/frame_insert_columns
        N = 10
        K = 5
        df = DataFrame(index=lrange(N))
        new_col = np.random.randn(N)
        for i in range(K):
            df[i] = new_col
        expected = DataFrame(np.repeat(new_col, K).reshape(N, K),
                             index=lrange(N))
        assert_frame_equal(df, expected)

    def test_insert(self):
        df = DataFrame(np.random.randn(5, 3), index=np.arange(5),
                       columns=['c', 'b', 'a'])

        df.insert(0, 'foo', df['a'])
        self.assert_numpy_array_equal(df.columns, ['foo', 'c', 'b', 'a'])
        tm.assert_series_equal(df['a'], df['foo'], check_names=False)

        df.insert(2, 'bar', df['c'])
        self.assert_numpy_array_equal(df.columns,
                                      ['foo', 'c', 'bar', 'b', 'a'])
        tm.assert_almost_equal(df['c'], df['bar'], check_names=False)

        # diff dtype

        # new item
        df['x'] = df['a'].astype('float32')
        result = Series(dict(float64=5, float32=1))
        self.assertTrue((df.get_dtype_counts() == result).all())

        # replacing current (in different block)
        df['a'] = df['a'].astype('float32')
        result = Series(dict(float64=4, float32=2))
        self.assertTrue((df.get_dtype_counts() == result).all())

        df['y'] = df['a'].astype('int32')
        result = Series(dict(float64=4, float32=2, int32=1))
        self.assertTrue((df.get_dtype_counts() == result).all())

        with assertRaisesRegexp(ValueError, 'already exists'):
            df.insert(1, 'a', df['b'])
        self.assertRaises(ValueError, df.insert, 1, 'c', df['b'])

        df.columns.name = 'some_name'
        # preserve columns name field
        df.insert(0, 'baz', df['c'])
        self.assertEqual(df.columns.name, 'some_name')

    def test_delitem(self):
        del self.frame['A']
        self.assertNotIn('A', self.frame)

    def test_pop(self):
        self.frame.columns.name = 'baz'

        self.frame.pop('A')
        self.assertNotIn('A', self.frame)

        self.frame['foo'] = 'bar'
        self.frame.pop('foo')
        self.assertNotIn('foo', self.frame)
        # TODO self.assertEqual(self.frame.columns.name, 'baz')

        # 10912
        # inplace ops cause caching issue
        a = DataFrame([[1, 2, 3], [4, 5, 6]], columns=[
                      'A', 'B', 'C'], index=['X', 'Y'])
        b = a.pop('B')
        b += 1

        # original frame
        expected = DataFrame([[1, 3], [4, 6]], columns=[
                             'A', 'C'], index=['X', 'Y'])
        assert_frame_equal(a, expected)

        # result
        expected = Series([2, 5], index=['X', 'Y'], name='B') + 1
        assert_series_equal(b, expected)

    def test_pop_non_unique_cols(self):
        df = DataFrame({0: [0, 1], 1: [0, 1], 2: [4, 5]})
        df.columns = ["a", "b", "a"]

        res = df.pop("a")
        self.assertEqual(type(res), DataFrame)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(df.columns), 1)
        self.assertTrue("b" in df.columns)
        self.assertFalse("a" in df.columns)
        self.assertEqual(len(df.index), 2)

    def test_insert_column_bug_4032(self):

        # GH4032, inserting a column and renaming causing errors
        df = DataFrame({'b': [1.1, 2.2]})
        df = df.rename(columns={})
        df.insert(0, 'a', [1, 2])

        result = df.rename(columns={})
        str(result)
        expected = DataFrame([[1, 1.1], [2, 2.2]], columns=['a', 'b'])
        assert_frame_equal(result, expected)
        df.insert(0, 'c', [1.3, 2.3])

        result = df.rename(columns={})
        str(result)

        expected = DataFrame([[1.3, 1, 1.1], [2.3, 2, 2.2]],
                             columns=['c', 'a', 'b'])
        assert_frame_equal(result, expected)
