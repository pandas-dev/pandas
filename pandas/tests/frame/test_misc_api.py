# -*- coding: utf-8 -*-

from __future__ import print_function
# pylint: disable-msg=W0612,E1101
from copy import deepcopy
import sys
import nose
from distutils.version import LooseVersion

from pandas.compat import range, lrange
from pandas import compat

from numpy.random import randn
import numpy as np

from pandas import DataFrame, Series
import pandas as pd

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


class SharedWithSparse(object):

    _multiprocess_can_split_ = True

    def test_copy_index_name_checking(self):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy
        for attr in ('index', 'columns'):
            ind = getattr(self.frame, attr)
            ind.name = None
            cp = self.frame.copy()
            getattr(cp, attr).name = 'foo'
            self.assertIsNone(getattr(self.frame, attr).name)

    def test_getitem_pop_assign_name(self):
        s = self.frame['A']
        self.assertEqual(s.name, 'A')

        s = self.frame.pop('A')
        self.assertEqual(s.name, 'A')

        s = self.frame.ix[:, 'B']
        self.assertEqual(s.name, 'B')

        s2 = s.ix[:]
        self.assertEqual(s2.name, 'B')

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

    def test_join_index(self):
        # left / right

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2)
        self.assertTrue(f.index.equals(joined.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='left')
        self.assertTrue(joined.index.equals(f.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='right')
        self.assertTrue(joined.index.equals(f2.index))
        self.assertEqual(len(joined.columns), 4)

        # inner

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='inner')
        self.assertTrue(joined.index.equals(f.index.intersection(f2.index)))
        self.assertEqual(len(joined.columns), 4)

        # outer

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='outer')
        self.assertTrue(tm.equalContents(self.frame.index, joined.index))
        self.assertEqual(len(joined.columns), 4)

        assertRaisesRegexp(ValueError, 'join method', f.join, f2, how='foo')

        # corner case - overlapping columns
        for how in ('outer', 'left', 'inner'):
            with assertRaisesRegexp(ValueError, 'columns overlap but '
                                    'no suffix'):
                self.frame.join(self.frame, how=how)

    def test_join_index_more(self):
        af = self.frame.ix[:, ['A', 'B']]
        bf = self.frame.ix[::2, ['C', 'D']]

        expected = af.copy()
        expected['C'] = self.frame['C'][::2]
        expected['D'] = self.frame['D'][::2]

        result = af.join(bf)
        assert_frame_equal(result, expected)

        result = af.join(bf, how='right')
        assert_frame_equal(result, expected[::2])

        result = bf.join(af, how='right')
        assert_frame_equal(result, expected.ix[:, result.columns])

    def test_join_index_series(self):
        df = self.frame.copy()
        s = df.pop(self.frame.columns[-1])
        joined = df.join(s)

        # TODO should this check_names ?
        assert_frame_equal(joined, self.frame, check_names=False)

        s.name = None
        assertRaisesRegexp(ValueError, 'must have a name', df.join, s)

    def test_join_overlap(self):
        df1 = self.frame.ix[:, ['A', 'B', 'C']]
        df2 = self.frame.ix[:, ['B', 'C', 'D']]

        joined = df1.join(df2, lsuffix='_df1', rsuffix='_df2')
        df1_suf = df1.ix[:, ['B', 'C']].add_suffix('_df1')
        df2_suf = df2.ix[:, ['B', 'C']].add_suffix('_df2')

        no_overlap = self.frame.ix[:, ['A', 'D']]
        expected = df1_suf.join(df2_suf).join(no_overlap)

        # column order not necessarily sorted
        assert_frame_equal(joined, expected.ix[:, joined.columns])

    def test_add_prefix_suffix(self):
        with_prefix = self.frame.add_prefix('foo#')
        expected = ['foo#%s' % c for c in self.frame.columns]
        self.assert_numpy_array_equal(with_prefix.columns, expected)

        with_suffix = self.frame.add_suffix('#foo')
        expected = ['%s#foo' % c for c in self.frame.columns]
        self.assert_numpy_array_equal(with_suffix.columns, expected)


class TestDataFrameMisc(tm.TestCase, SharedWithSparse, TestData):

    klass = DataFrame

    _multiprocess_can_split_ = True

    def test_get_axis(self):
        f = self.frame
        self.assertEqual(f._get_axis_number(0), 0)
        self.assertEqual(f._get_axis_number(1), 1)
        self.assertEqual(f._get_axis_number('index'), 0)
        self.assertEqual(f._get_axis_number('rows'), 0)
        self.assertEqual(f._get_axis_number('columns'), 1)

        self.assertEqual(f._get_axis_name(0), 'index')
        self.assertEqual(f._get_axis_name(1), 'columns')
        self.assertEqual(f._get_axis_name('index'), 'index')
        self.assertEqual(f._get_axis_name('rows'), 'index')
        self.assertEqual(f._get_axis_name('columns'), 'columns')

        self.assertIs(f._get_axis(0), f.index)
        self.assertIs(f._get_axis(1), f.columns)

        assertRaisesRegexp(ValueError, 'No axis named', f._get_axis_number, 2)
        assertRaisesRegexp(ValueError, 'No axis.*foo', f._get_axis_name, 'foo')
        assertRaisesRegexp(ValueError, 'No axis.*None', f._get_axis_name, None)
        assertRaisesRegexp(ValueError, 'No axis named', f._get_axis_number,
                           None)

    def test_keys(self):
        getkeys = self.frame.keys
        self.assertIs(getkeys(), self.frame.columns)

    def test_column_contains_typeerror(self):
        try:
            self.frame.columns in self.frame
        except TypeError:
            pass

    def test_not_hashable(self):
        df = pd.DataFrame([1])
        self.assertRaises(TypeError, hash, df)
        self.assertRaises(TypeError, hash, self.empty)

    def test_new_empty_index(self):
        df1 = DataFrame(randn(0, 3))
        df2 = DataFrame(randn(0, 3))
        df1.index.name = 'foo'
        self.assertIsNone(df2.index.name)

    def test_array_interface(self):
        result = np.sqrt(self.frame)
        tm.assertIsInstance(result, type(self.frame))
        self.assertIs(result.index, self.frame.index)
        self.assertIs(result.columns, self.frame.columns)

        assert_frame_equal(result, self.frame.apply(np.sqrt))

    def test_get_agg_axis(self):
        cols = self.frame._get_agg_axis(0)
        self.assertIs(cols, self.frame.columns)

        idx = self.frame._get_agg_axis(1)
        self.assertIs(idx, self.frame.index)

        self.assertRaises(ValueError, self.frame._get_agg_axis, 2)

    def test_nonzero(self):
        self.assertTrue(self.empty.empty)

        self.assertFalse(self.frame.empty)
        self.assertFalse(self.mixed_frame.empty)

        # corner case
        df = DataFrame({'A': [1., 2., 3.],
                        'B': ['a', 'b', 'c']},
                       index=np.arange(3))
        del df['A']
        self.assertFalse(df.empty)

    def test_iteritems(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=['a', 'a', 'b'])
        for k, v in compat.iteritems(df):
            self.assertEqual(type(v), Series)

    def test_iter(self):
        self.assertTrue(tm.equalContents(list(self.frame), self.frame.columns))

    def test_iterrows(self):
        for i, (k, v) in enumerate(self.frame.iterrows()):
            exp = self.frame.xs(self.frame.index[i])
            assert_series_equal(v, exp)

        for i, (k, v) in enumerate(self.mixed_frame.iterrows()):
            exp = self.mixed_frame.xs(self.mixed_frame.index[i])
            assert_series_equal(v, exp)

    def test_itertuples(self):
        for i, tup in enumerate(self.frame.itertuples()):
            s = Series(tup[1:])
            s.name = tup[0]
            expected = self.frame.ix[i, :].reset_index(drop=True)
            assert_series_equal(s, expected)

        df = DataFrame({'floats': np.random.randn(5),
                        'ints': lrange(5)}, columns=['floats', 'ints'])

        for tup in df.itertuples(index=False):
            tm.assertIsInstance(tup[1], np.integer)

        df = DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
        dfaa = df[['a', 'a']]
        self.assertEqual(list(dfaa.itertuples()), [
                         (0, 1, 1), (1, 2, 2), (2, 3, 3)])

        self.assertEqual(repr(list(df.itertuples(name=None))),
                         '[(0, 1, 4), (1, 2, 5), (2, 3, 6)]')

        tup = next(df.itertuples(name='TestName'))

        # no support for field renaming in Python 2.6, regular tuples are
        # returned
        if sys.version >= LooseVersion('2.7'):
            self.assertEqual(tup._fields, ('Index', 'a', 'b'))
            self.assertEqual((tup.Index, tup.a, tup.b), tup)
            self.assertEqual(type(tup).__name__, 'TestName')

        df.columns = ['def', 'return']
        tup2 = next(df.itertuples(name='TestName'))
        self.assertEqual(tup2, (0, 1, 4))

        if sys.version >= LooseVersion('2.7'):
            self.assertEqual(tup2._fields, ('Index', '_1', '_2'))

        df3 = DataFrame(dict(('f' + str(i), [i]) for i in range(1024)))
        # will raise SyntaxError if trying to create namedtuple
        tup3 = next(df3.itertuples())
        self.assertFalse(hasattr(tup3, '_fields'))
        self.assertIsInstance(tup3, tuple)

    def test_len(self):
        self.assertEqual(len(self.frame), len(self.frame.index))

    def test_as_matrix(self):
        frame = self.frame
        mat = frame.as_matrix()

        frameCols = frame.columns
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if np.isnan(value):
                    self.assertTrue(np.isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])

        # mixed type
        mat = self.mixed_frame.as_matrix(['foo', 'A'])
        self.assertEqual(mat[0, 0], 'bar')

        df = DataFrame({'real': [1, 2, 3], 'complex': [1j, 2j, 3j]})
        mat = df.as_matrix()
        self.assertEqual(mat[0, 0], 1j)

        # single block corner case
        mat = self.frame.as_matrix(['A', 'B'])
        expected = self.frame.reindex(columns=['A', 'B']).values
        assert_almost_equal(mat, expected)

    def test_values(self):
        self.frame.values[:, 0] = 5.
        self.assertTrue((self.frame.values[:, 0] == 5).all())

    def test_deepcopy(self):
        cp = deepcopy(self.frame)
        series = cp['A']
        series[:] = 10
        for idx, value in compat.iteritems(series):
            self.assertNotEqual(self.frame['A'][idx], value)

    # ---------------------------------------------------------------------
    # Transposing

    def test_transpose(self):
        frame = self.frame
        dft = frame.T
        for idx, series in compat.iteritems(dft):
            for col, value in compat.iteritems(series):
                if np.isnan(value):
                    self.assertTrue(np.isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

        # mixed type
        index, data = tm.getMixedTypeDict()
        mixed = DataFrame(data, index=index)

        mixed_T = mixed.T
        for col, s in compat.iteritems(mixed_T):
            self.assertEqual(s.dtype, np.object_)

    def test_transpose_get_view(self):
        dft = self.frame.T
        dft.values[:, 5:10] = 5

        self.assertTrue((self.frame.values[5:10] == 5).all())

    def test_swapaxes(self):
        df = DataFrame(np.random.randn(10, 5))
        assert_frame_equal(df.T, df.swapaxes(0, 1))
        assert_frame_equal(df.T, df.swapaxes(1, 0))
        assert_frame_equal(df, df.swapaxes(0, 0))
        self.assertRaises(ValueError, df.swapaxes, 2, 5)

    def test_axis_aliases(self):
        f = self.frame

        # reg name
        expected = f.sum(axis=0)
        result = f.sum(axis='index')
        assert_series_equal(result, expected)

        expected = f.sum(axis=1)
        result = f.sum(axis='columns')
        assert_series_equal(result, expected)

    def test_more_asMatrix(self):
        values = self.mixed_frame.as_matrix()
        self.assertEqual(values.shape[1], len(self.mixed_frame.columns))

    def test_repr_with_mi_nat(self):
        df = DataFrame({'X': [1, 2]},
                       index=[[pd.NaT, pd.Timestamp('20130101')], ['a', 'b']])
        res = repr(df)
        exp = '              X\nNaT        a  1\n2013-01-01 b  2'
        nose.tools.assert_equal(res, exp)

    def test_iterkv_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            self.mixed_float.iterkv()

    def test_iterkv_names(self):
        for k, v in compat.iteritems(self.mixed_frame):
            self.assertEqual(v.name, k)

    def test_series_put_names(self):
        series = self.mixed_frame._series
        for k, v in compat.iteritems(series):
            self.assertEqual(v.name, k)

    def test_empty_nonzero(self):
        df = DataFrame([1, 2, 3])
        self.assertFalse(df.empty)
        df = DataFrame(index=['a', 'b'], columns=['c', 'd']).dropna()
        self.assertTrue(df.empty)
        self.assertTrue(df.T.empty)

    def test_inplace_return_self(self):
        # re #1893

        data = DataFrame({'a': ['foo', 'bar', 'baz', 'qux'],
                          'b': [0, 0, 1, 1],
                          'c': [1, 2, 3, 4]})

        def _check_f(base, f):
            result = f(base)
            self.assertTrue(result is None)

        # -----DataFrame-----

        # set_index
        f = lambda x: x.set_index('a', inplace=True)
        _check_f(data.copy(), f)

        # reset_index
        f = lambda x: x.reset_index(inplace=True)
        _check_f(data.set_index('a'), f)

        # drop_duplicates
        f = lambda x: x.drop_duplicates(inplace=True)
        _check_f(data.copy(), f)

        # sort
        f = lambda x: x.sort_values('b', inplace=True)
        _check_f(data.copy(), f)

        # sort_index
        f = lambda x: x.sort_index(inplace=True)
        _check_f(data.copy(), f)

        # sortlevel
        f = lambda x: x.sortlevel(0, inplace=True)
        _check_f(data.set_index(['a', 'b']), f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(data.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(data.copy(), f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(data.copy(), f)

        # -----Series-----
        d = data.copy()['c']

        # reset_index
        f = lambda x: x.reset_index(inplace=True, drop=True)
        _check_f(data.set_index('a')['c'], f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(d.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(d.copy(), f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(d.copy(), f)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
