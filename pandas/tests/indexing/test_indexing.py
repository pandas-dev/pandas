# -*- coding: utf-8 -*-
# pylint: disable-msg=W0612,E1101
import itertools
import warnings
from warnings import catch_warnings
from datetime import datetime

from pandas.types.common import (is_integer_dtype,
                                 is_float_dtype,
                                 is_scalar)
from pandas.compat import range, lrange, lzip, StringIO, lmap
from pandas.tslib import NaT
from numpy import nan
from numpy.random import randn
import numpy as np

import pandas as pd
from pandas import option_context
from pandas.core.indexing import _non_reducing_slice, _maybe_numeric_slice
from pandas.core.api import (DataFrame, Index, Series, Panel, isnull,
                             MultiIndex, Timestamp, Timedelta, UInt64Index)
from pandas.formats.printing import pprint_thing
from pandas import concat
from pandas.core.common import PerformanceWarning
from pandas.tests.indexing.common import _mklbl
import pandas.util.testing as tm
from pandas import date_range


_verbose = False

# ------------------------------------------------------------------------
# Indexing test cases


def _generate_indices(f, values=False):
    """ generate the indicies
          if values is True , use the axis values
                    is False, use the range
                    """

    axes = f.axes
    if values:
        axes = [lrange(len(a)) for a in axes]

    return itertools.product(*axes)


def _get_value(f, i, values=False):
    """ return the value for the location i """

    # check agains values
    if values:
        return f.values[i]

    # this is equiv of f[col][row].....
    # v = f
    # for a in reversed(i):
    #    v = v.__getitem__(a)
    # return v
    with catch_warnings(record=True):
        return f.ix[i]


def _get_result(obj, method, key, axis):
    """ return the result for this obj with this key and this axis """

    if isinstance(key, dict):
        key = key[axis]

    # use an artifical conversion to map the key as integers to the labels
    # so ix can work for comparisions
    if method == 'indexer':
        method = 'ix'
        key = obj._get_axis(axis)[key]

    # in case we actually want 0 index slicing
    try:
        xp = getattr(obj, method).__getitem__(_axify(obj, key, axis))
    except:
        xp = getattr(obj, method).__getitem__(key)

    return xp


def _axify(obj, key, axis):
    # create a tuple accessor
    axes = [slice(None)] * obj.ndim
    axes[axis] = key
    return tuple(axes)


class TestIndexing(tm.TestCase):

    _objs = set(['series', 'frame', 'panel'])
    _typs = set(['ints', 'uints', 'labels', 'mixed',
                 'ts', 'floats', 'empty', 'ts_rev'])

    def setUp(self):

        self.series_ints = Series(np.random.rand(4), index=lrange(0, 8, 2))
        self.frame_ints = DataFrame(np.random.randn(4, 4),
                                    index=lrange(0, 8, 2),
                                    columns=lrange(0, 12, 3))
        self.panel_ints = Panel(np.random.rand(4, 4, 4),
                                items=lrange(0, 8, 2),
                                major_axis=lrange(0, 12, 3),
                                minor_axis=lrange(0, 16, 4))

        self.series_uints = Series(np.random.rand(4),
                                   index=UInt64Index(lrange(0, 8, 2)))
        self.frame_uints = DataFrame(np.random.randn(4, 4),
                                     index=UInt64Index(lrange(0, 8, 2)),
                                     columns=UInt64Index(lrange(0, 12, 3)))
        self.panel_uints = Panel(np.random.rand(4, 4, 4),
                                 items=UInt64Index(lrange(0, 8, 2)),
                                 major_axis=UInt64Index(lrange(0, 12, 3)),
                                 minor_axis=UInt64Index(lrange(0, 16, 4)))

        self.series_labels = Series(np.random.randn(4), index=list('abcd'))
        self.frame_labels = DataFrame(np.random.randn(4, 4),
                                      index=list('abcd'), columns=list('ABCD'))
        self.panel_labels = Panel(np.random.randn(4, 4, 4),
                                  items=list('abcd'),
                                  major_axis=list('ABCD'),
                                  minor_axis=list('ZYXW'))

        self.series_mixed = Series(np.random.randn(4), index=[2, 4, 'null', 8])
        self.frame_mixed = DataFrame(np.random.randn(4, 4),
                                     index=[2, 4, 'null', 8])
        self.panel_mixed = Panel(np.random.randn(4, 4, 4),
                                 items=[2, 4, 'null', 8])

        self.series_ts = Series(np.random.randn(4),
                                index=date_range('20130101', periods=4))
        self.frame_ts = DataFrame(np.random.randn(4, 4),
                                  index=date_range('20130101', periods=4))
        self.panel_ts = Panel(np.random.randn(4, 4, 4),
                              items=date_range('20130101', periods=4))

        dates_rev = (date_range('20130101', periods=4)
                     .sort_values(ascending=False))
        self.series_ts_rev = Series(np.random.randn(4),
                                    index=dates_rev)
        self.frame_ts_rev = DataFrame(np.random.randn(4, 4),
                                      index=dates_rev)
        self.panel_ts_rev = Panel(np.random.randn(4, 4, 4),
                                  items=dates_rev)

        self.frame_empty = DataFrame({})
        self.series_empty = Series({})
        self.panel_empty = Panel({})

        # form agglomerates
        for o in self._objs:

            d = dict()
            for t in self._typs:
                d[t] = getattr(self, '%s_%s' % (o, t), None)

            setattr(self, o, d)

    def check_values(self, f, func, values=False):

        if f is None:
            return
        axes = f.axes
        indicies = itertools.product(*axes)

        for i in indicies:
            result = getattr(f, func)[i]

            # check agains values
            if values:
                expected = f.values[i]
            else:
                expected = f
                for a in reversed(i):
                    expected = expected.__getitem__(a)

            tm.assert_almost_equal(result, expected)

    def check_result(self, name, method1, key1, method2, key2, typs=None,
                     objs=None, axes=None, fails=None):
        def _eq(t, o, a, obj, k1, k2):
            """ compare equal for these 2 keys """

            if a is not None and a > obj.ndim - 1:
                return

            def _print(result, error=None):
                if error is not None:
                    error = str(error)
                v = ("%-16.16s [%-16.16s]: [typ->%-8.8s,obj->%-8.8s,"
                     "key1->(%-4.4s),key2->(%-4.4s),axis->%s] %s" %
                     (name, result, t, o, method1, method2, a, error or ''))
                if _verbose:
                    pprint_thing(v)

            try:
                rs = getattr(obj, method1).__getitem__(_axify(obj, k1, a))

                try:
                    xp = _get_result(obj, method2, k2, a)
                except:
                    result = 'no comp'
                    _print(result)
                    return

                detail = None

                try:
                    if is_scalar(rs) and is_scalar(xp):
                        self.assertEqual(rs, xp)
                    elif xp.ndim == 1:
                        tm.assert_series_equal(rs, xp)
                    elif xp.ndim == 2:
                        tm.assert_frame_equal(rs, xp)
                    elif xp.ndim == 3:
                        tm.assert_panel_equal(rs, xp)
                    result = 'ok'
                except AssertionError as e:
                    detail = str(e)
                    result = 'fail'

                # reverse the checks
                if fails is True:
                    if result == 'fail':
                        result = 'ok (fail)'

                _print(result)
                if not result.startswith('ok'):
                    raise AssertionError(detail)

            except AssertionError:
                raise
            except Exception as detail:

                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if isinstance(detail, fails):
                        result = 'ok (%s)' % type(detail).__name__
                        _print(result)
                        return

                result = type(detail).__name__
                raise AssertionError(_print(result, error=detail))

        if typs is None:
            typs = self._typs

        if objs is None:
            objs = self._objs

        if axes is not None:
            if not isinstance(axes, (tuple, list)):
                axes = [axes]
            else:
                axes = list(axes)
        else:
            axes = [0, 1, 2]

        # check
        for o in objs:
            if o not in self._objs:
                continue

            d = getattr(self, o)
            for a in axes:
                for t in typs:
                    if t not in self._typs:
                        continue

                    obj = d[t]
                    if obj is not None:
                        obj = obj.copy()

                        k2 = key2
                        _eq(t, o, a, obj, key1, k2)

    def test_ix_deprecation(self):
        # GH 15114

        df = DataFrame({'A': [1, 2, 3]})
        with tm.assert_produces_warning(DeprecationWarning,
                                        check_stacklevel=False):
            df.ix[1, 'A']

    def test_indexer_caching(self):
        # GH5727
        # make sure that indexers are in the _internal_names_set
        n = 1000001
        arrays = [lrange(n), lrange(n)]
        index = MultiIndex.from_tuples(lzip(*arrays))
        s = Series(np.zeros(n), index=index)
        str(s)

        # setitem
        expected = Series(np.ones(n), index=index)
        s = Series(np.zeros(n), index=index)
        s[s == 0] = 1
        tm.assert_series_equal(s, expected)

    def test_at_and_iat_get(self):
        def _check(f, func, values=False):

            if f is not None:
                indicies = _generate_indices(f, values)
                for i in indicies:
                    result = getattr(f, func)[i]
                    expected = _get_value(f, i, values)
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
                indicies = _generate_indices(f, values)
                for i in indicies:
                    getattr(f, func)[i] = 1
                    expected = _get_value(f, i, values)
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
        df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
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

    def test_repeated_getitem_dups(self):
        # GH 5678
        # repeated gettitems on a dup index returing a ndarray
        df = DataFrame(
            np.random.random_sample((20, 5)),
            index=['ABCDE' [x % 5] for x in range(20)])
        expected = df.loc['A', 0]
        result = df.loc[:, 0].loc['A']
        tm.assert_series_equal(result, expected)

    def test_iloc_exceeds_bounds(self):

        # GH6296
        # iloc should allow indexers that exceed the bounds
        df = DataFrame(np.random.random_sample((20, 5)), columns=list('ABCDE'))
        expected = df

        # lists of positions should raise IndexErrror!
        with tm.assertRaisesRegexp(IndexError,
                                   'positional indexers are out-of-bounds'):
            df.iloc[:, [0, 1, 2, 3, 4, 5]]
        self.assertRaises(IndexError, lambda: df.iloc[[1, 30]])
        self.assertRaises(IndexError, lambda: df.iloc[[1, -30]])
        self.assertRaises(IndexError, lambda: df.iloc[[100]])

        s = df['A']
        self.assertRaises(IndexError, lambda: s.iloc[[100]])
        self.assertRaises(IndexError, lambda: s.iloc[[-100]])

        # still raise on a single indexer
        msg = 'single positional indexer is out-of-bounds'
        with tm.assertRaisesRegexp(IndexError, msg):
            df.iloc[30]
        self.assertRaises(IndexError, lambda: df.iloc[-30])

        # GH10779
        # single positive/negative indexer exceeding Series bounds should raise
        # an IndexError
        with tm.assertRaisesRegexp(IndexError, msg):
            s.iloc[30]
        self.assertRaises(IndexError, lambda: s.iloc[-30])

        # slices are ok
        result = df.iloc[:, 4:10]  # 0 < start < len < stop
        expected = df.iloc[:, 4:]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -4:-10]  # stop < 0 < start < len
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:4:-1]  # 0 < stop < len < start (down)
        expected = df.iloc[:, :4:-1]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 4:-10:-1]  # stop < 0 < start < len (down)
        expected = df.iloc[:, 4::-1]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -10:4]  # start < 0 < stop < len
        expected = df.iloc[:, :4]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:4]  # 0 < stop < len < start
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -10:-11:-1]  # stop < start < 0 < len (down)
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:11]  # 0 < len < start < stop
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        # slice bounds exceeding is ok
        result = s.iloc[18:30]
        expected = s.iloc[18:]
        tm.assert_series_equal(result, expected)

        result = s.iloc[30:]
        expected = s.iloc[:0]
        tm.assert_series_equal(result, expected)

        result = s.iloc[30::-1]
        expected = s.iloc[::-1]
        tm.assert_series_equal(result, expected)

        # doc example
        def check(result, expected):
            str(result)
            result.dtypes
            tm.assert_frame_equal(result, expected)

        dfl = DataFrame(np.random.randn(5, 2), columns=list('AB'))
        check(dfl.iloc[:, 2:3], DataFrame(index=dfl.index))
        check(dfl.iloc[:, 1:3], dfl.iloc[:, [1]])
        check(dfl.iloc[4:6], dfl.iloc[[4]])

        self.assertRaises(IndexError, lambda: dfl.iloc[[4, 5, 6]])
        self.assertRaises(IndexError, lambda: dfl.iloc[:, 4])

    def test_iloc_getitem_int(self):

        # integer
        self.check_result('integer', 'iloc', 2, 'ix',
                          {0: 4, 1: 6, 2: 8}, typs=['ints', 'uints'])
        self.check_result('integer', 'iloc', 2, 'indexer', 2,
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

    def test_iloc_getitem_neg_int(self):

        # neg integer
        self.check_result('neg int', 'iloc', -1, 'ix',
                          {0: 6, 1: 9, 2: 12}, typs=['ints', 'uints'])
        self.check_result('neg int', 'iloc', -1, 'indexer', -1,
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

    def test_iloc_getitem_list_int(self):

        # list of ints
        self.check_result('list int', 'iloc', [0, 1, 2], 'ix',
                          {0: [0, 2, 4], 1: [0, 3, 6], 2: [0, 4, 8]},
                          typs=['ints', 'uints'])
        self.check_result('list int', 'iloc', [2], 'ix',
                          {0: [4], 1: [6], 2: [8]}, typs=['ints', 'uints'])
        self.check_result('list int', 'iloc', [0, 1, 2], 'indexer', [0, 1, 2],
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

        # array of ints (GH5006), make sure that a single indexer is returning
        # the correct type
        self.check_result('array int', 'iloc', np.array([0, 1, 2]), 'ix',
                          {0: [0, 2, 4],
                           1: [0, 3, 6],
                           2: [0, 4, 8]}, typs=['ints', 'uints'])
        self.check_result('array int', 'iloc', np.array([2]), 'ix',
                          {0: [4], 1: [6], 2: [8]}, typs=['ints', 'uints'])
        self.check_result('array int', 'iloc', np.array([0, 1, 2]), 'indexer',
                          [0, 1, 2],
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

    def test_iloc_getitem_neg_int_can_reach_first_index(self):
        # GH10547 and GH10779
        # negative integers should be able to reach index 0
        df = DataFrame({'A': [2, 3, 5], 'B': [7, 11, 13]})
        s = df['A']

        expected = df.iloc[0]
        result = df.iloc[-3]
        tm.assert_series_equal(result, expected)

        expected = df.iloc[[0]]
        result = df.iloc[[-3]]
        tm.assert_frame_equal(result, expected)

        expected = s.iloc[0]
        result = s.iloc[-3]
        self.assertEqual(result, expected)

        expected = s.iloc[[0]]
        result = s.iloc[[-3]]
        tm.assert_series_equal(result, expected)

        # check the length 1 Series case highlighted in GH10547
        expected = pd.Series(['a'], index=['A'])
        result = expected.iloc[[-1]]
        tm.assert_series_equal(result, expected)

    def test_iloc_getitem_dups(self):

        # no dups in panel (bug?)
        self.check_result('list int (dups)', 'iloc', [0, 1, 1, 3], 'ix',
                          {0: [0, 2, 2, 6], 1: [0, 3, 3, 9]},
                          objs=['series', 'frame'], typs=['ints', 'uints'])

        # GH 6766
        df1 = DataFrame([{'A': None, 'B': 1}, {'A': 2, 'B': 2}])
        df2 = DataFrame([{'A': 3, 'B': 3}, {'A': 4, 'B': 4}])
        df = concat([df1, df2], axis=1)

        # cross-sectional indexing
        result = df.iloc[0, 0]
        self.assertTrue(isnull(result))

        result = df.iloc[0, :]
        expected = Series([np.nan, 1, 3, 3], index=['A', 'B', 'A', 'B'],
                          name=0)
        tm.assert_series_equal(result, expected)

    def test_iloc_getitem_array(self):

        # array like
        s = Series(index=lrange(1, 4))
        self.check_result('array like', 'iloc', s.index, 'ix',
                          {0: [2, 4, 6], 1: [3, 6, 9], 2: [4, 8, 12]},
                          typs=['ints', 'uints'])

    def test_iloc_getitem_bool(self):

        # boolean indexers
        b = [True, False, True, False, ]
        self.check_result('bool', 'iloc', b, 'ix', b, typs=['ints', 'uints'])
        self.check_result('bool', 'iloc', b, 'ix', b,
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

    def test_iloc_getitem_slice(self):

        # slices
        self.check_result('slice', 'iloc', slice(1, 3), 'ix',
                          {0: [2, 4], 1: [3, 6], 2: [4, 8]},
                          typs=['ints', 'uints'])
        self.check_result('slice', 'iloc', slice(1, 3), 'indexer',
                          slice(1, 3),
                          typs=['labels', 'mixed', 'ts', 'floats', 'empty'],
                          fails=IndexError)

    def test_iloc_getitem_slice_dups(self):

        df1 = DataFrame(np.random.randn(10, 4), columns=['A', 'A', 'B', 'B'])
        df2 = DataFrame(np.random.randint(0, 10, size=20).reshape(10, 2),
                        columns=['A', 'C'])

        # axis=1
        df = concat([df1, df2], axis=1)
        tm.assert_frame_equal(df.iloc[:, :4], df1)
        tm.assert_frame_equal(df.iloc[:, 4:], df2)

        df = concat([df2, df1], axis=1)
        tm.assert_frame_equal(df.iloc[:, :2], df2)
        tm.assert_frame_equal(df.iloc[:, 2:], df1)

        exp = concat([df2, df1.iloc[:, [0]]], axis=1)
        tm.assert_frame_equal(df.iloc[:, 0:3], exp)

        # axis=0
        df = concat([df, df], axis=0)
        tm.assert_frame_equal(df.iloc[0:10, :2], df2)
        tm.assert_frame_equal(df.iloc[0:10, 2:], df1)
        tm.assert_frame_equal(df.iloc[10:, :2], df2)
        tm.assert_frame_equal(df.iloc[10:, 2:], df1)

    def test_iloc_setitem(self):
        df = self.frame_ints

        df.iloc[1, 1] = 1
        result = df.iloc[1, 1]
        self.assertEqual(result, 1)

        df.iloc[:, 2:3] = 0
        expected = df.iloc[:, 2:3]
        result = df.iloc[:, 2:3]
        tm.assert_frame_equal(result, expected)

        # GH5771
        s = Series(0, index=[4, 5, 6])
        s.iloc[1:2] += 1
        expected = Series([0, 1, 0], index=[4, 5, 6])
        tm.assert_series_equal(s, expected)

    def test_loc_setitem_slice(self):
        # GH10503

        # assigning the same type should not change the type
        df1 = DataFrame({'a': [0, 1, 1],
                         'b': Series([100, 200, 300], dtype='uint32')})
        ix = df1['a'] == 1
        newb1 = df1.loc[ix, 'b'] + 1
        df1.loc[ix, 'b'] = newb1
        expected = DataFrame({'a': [0, 1, 1],
                              'b': Series([100, 201, 301], dtype='uint32')})
        tm.assert_frame_equal(df1, expected)

        # assigning a new type should get the inferred type
        df2 = DataFrame({'a': [0, 1, 1], 'b': [100, 200, 300]},
                        dtype='uint64')
        ix = df1['a'] == 1
        newb2 = df2.loc[ix, 'b']
        df1.loc[ix, 'b'] = newb2
        expected = DataFrame({'a': [0, 1, 1], 'b': [100, 200, 300]},
                             dtype='uint64')
        tm.assert_frame_equal(df2, expected)

    def test_ix_loc_setitem_consistency(self):

        # GH 5771
        # loc with slice and series
        s = Series(0, index=[4, 5, 6])
        s.loc[4:5] += 1
        expected = Series([1, 1, 0], index=[4, 5, 6])
        tm.assert_series_equal(s, expected)

        # GH 5928
        # chained indexing assignment
        df = DataFrame({'a': [0, 1, 2]})
        expected = df.copy()
        with catch_warnings(record=True):
            expected.ix[[0, 1, 2], 'a'] = -expected.ix[[0, 1, 2], 'a']

        with catch_warnings(record=True):
            df['a'].ix[[0, 1, 2]] = -df['a'].ix[[0, 1, 2]]
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'a': [0, 1, 2], 'b': [0, 1, 2]})
        with catch_warnings(record=True):
            df['a'].ix[[0, 1, 2]] = -df['a'].ix[[0, 1, 2]].astype(
                'float64') + 0.5
        expected = DataFrame({'a': [0.5, -0.5, -1.5], 'b': [0, 1, 2]})
        tm.assert_frame_equal(df, expected)

        # GH 8607
        # ix setitem consistency
        df = DataFrame({'timestamp': [1413840976, 1413842580, 1413760580],
                        'delta': [1174, 904, 161],
                        'elapsed': [7673, 9277, 1470]})
        expected = DataFrame({'timestamp': pd.to_datetime(
            [1413840976, 1413842580, 1413760580], unit='s'),
            'delta': [1174, 904, 161],
            'elapsed': [7673, 9277, 1470]})

        df2 = df.copy()
        df2['timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

        df2 = df.copy()
        df2.loc[:, 'timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

        df2 = df.copy()
        with catch_warnings(record=True):
            df2.ix[:, 2] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

    def test_ix_loc_consistency(self):

        # GH 8613
        # some edge cases where ix/loc should return the same
        # this is not an exhaustive case

        def compare(result, expected):
            if is_scalar(expected):
                self.assertEqual(result, expected)
            else:
                self.assertTrue(expected.equals(result))

        # failure cases for .loc, but these work for .ix
        df = pd.DataFrame(np.random.randn(5, 4), columns=list('ABCD'))
        for key in [slice(1, 3), tuple([slice(0, 2), slice(0, 2)]),
                    tuple([slice(0, 2), df.columns[0:2]])]:

            for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                          tm.makeDateIndex, tm.makePeriodIndex,
                          tm.makeTimedeltaIndex]:
                df.index = index(len(df.index))
                with catch_warnings(record=True):
                    df.ix[key]

                self.assertRaises(TypeError, lambda: df.loc[key])

        df = pd.DataFrame(np.random.randn(5, 4), columns=list('ABCD'),
                          index=pd.date_range('2012-01-01', periods=5))

        for key in ['2012-01-03',
                    '2012-01-31',
                    slice('2012-01-03', '2012-01-03'),
                    slice('2012-01-03', '2012-01-04'),
                    slice('2012-01-03', '2012-01-06', 2),
                    slice('2012-01-03', '2012-01-31'),
                    tuple([[True, True, True, False, True]]), ]:

            # getitem

            # if the expected raises, then compare the exceptions
            try:
                with catch_warnings(record=True):
                    expected = df.ix[key]
            except KeyError:
                self.assertRaises(KeyError, lambda: df.loc[key])
                continue

            result = df.loc[key]
            compare(result, expected)

            # setitem
            df1 = df.copy()
            df2 = df.copy()

            with catch_warnings(record=True):
                df1.ix[key] = 10
            df2.loc[key] = 10
            compare(df2, df1)

        # edge cases
        s = Series([1, 2, 3, 4], index=list('abde'))

        result1 = s['a':'c']
        with catch_warnings(record=True):
            result2 = s.ix['a':'c']
        result3 = s.loc['a':'c']
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

        # now work rather than raising KeyError
        s = Series(range(5), [-2, -1, 1, 2, 3])

        with catch_warnings(record=True):
            result1 = s.ix[-10:3]
        result2 = s.loc[-10:3]
        tm.assert_series_equal(result1, result2)

        with catch_warnings(record=True):
            result1 = s.ix[0:3]
        result2 = s.loc[0:3]
        tm.assert_series_equal(result1, result2)

    def test_loc_setitem_dups(self):

        # GH 6541
        df_orig = DataFrame(
            {'me': list('rttti'),
             'foo': list('aaade'),
             'bar': np.arange(5, dtype='float64') * 1.34 + 2,
             'bar2': np.arange(5, dtype='float64') * -.34 + 2}).set_index('me')

        indexer = tuple(['r', ['bar', 'bar2']])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_series_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

        indexer = tuple(['r', 'bar'])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        self.assertEqual(df.loc[indexer], 2.0 * df_orig.loc[indexer])

        indexer = tuple(['t', ['bar', 'bar2']])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_frame_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

    def test_iloc_setitem_dups(self):

        # GH 6766
        # iloc with a mask aligning from another iloc
        df1 = DataFrame([{'A': None, 'B': 1}, {'A': 2, 'B': 2}])
        df2 = DataFrame([{'A': 3, 'B': 3}, {'A': 4, 'B': 4}])
        df = concat([df1, df2], axis=1)

        expected = df.fillna(3)
        expected['A'] = expected['A'].astype('float64')
        inds = np.isnan(df.iloc[:, 0])
        mask = inds[inds].index
        df.iloc[mask, 0] = df.iloc[mask, 2]
        tm.assert_frame_equal(df, expected)

        # del a dup column across blocks
        expected = DataFrame({0: [1, 2], 1: [3, 4]})
        expected.columns = ['B', 'B']
        del df['A']
        tm.assert_frame_equal(df, expected)

        # assign back to self
        df.iloc[[0, 1], [0, 1]] = df.iloc[[0, 1], [0, 1]]
        tm.assert_frame_equal(df, expected)

        # reversed x 2
        df.iloc[[1, 0], [0, 1]] = df.iloc[[1, 0], [0, 1]].reset_index(
            drop=True)
        df.iloc[[1, 0], [0, 1]] = df.iloc[[1, 0], [0, 1]].reset_index(
            drop=True)
        tm.assert_frame_equal(df, expected)

    def test_chained_getitem_with_lists(self):

        # GH6394
        # Regression in chained getitem indexing with embedded list-like from
        # 0.12
        def check(result, expected):
            tm.assert_numpy_array_equal(result, expected)
            tm.assertIsInstance(result, np.ndarray)

        df = DataFrame({'A': 5 * [np.zeros(3)], 'B': 5 * [np.ones(3)]})
        expected = df['A'].iloc[2]
        result = df.loc[2, 'A']
        check(result, expected)
        result2 = df.iloc[2]['A']
        check(result2, expected)
        result3 = df['A'].loc[2]
        check(result3, expected)
        result4 = df['A'].iloc[2]
        check(result4, expected)

    def test_loc_getitem_int(self):

        # int label
        self.check_result('int label', 'loc', 2, 'ix', 2,
                          typs=['ints', 'uints'], axes=0)
        self.check_result('int label', 'loc', 3, 'ix', 3,
                          typs=['ints', 'uints'], axes=1)
        self.check_result('int label', 'loc', 4, 'ix', 4,
                          typs=['ints', 'uints'], axes=2)
        self.check_result('int label', 'loc', 2, 'ix', 2,
                          typs=['label'], fails=KeyError)

    def test_loc_getitem_label(self):

        # label
        self.check_result('label', 'loc', 'c', 'ix', 'c', typs=['labels'],
                          axes=0)
        self.check_result('label', 'loc', 'null', 'ix', 'null', typs=['mixed'],
                          axes=0)
        self.check_result('label', 'loc', 8, 'ix', 8, typs=['mixed'], axes=0)
        self.check_result('label', 'loc', Timestamp('20130102'), 'ix', 1,
                          typs=['ts'], axes=0)
        self.check_result('label', 'loc', 'c', 'ix', 'c', typs=['empty'],
                          fails=KeyError)

    def test_loc_getitem_label_out_of_range(self):

        # out of range label
        self.check_result('label range', 'loc', 'f', 'ix', 'f',
                          typs=['ints', 'uints', 'labels', 'mixed', 'ts'],
                          fails=KeyError)
        self.check_result('label range', 'loc', 'f', 'ix', 'f',
                          typs=['floats'], fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20,
                          typs=['ints', 'uints', 'mixed'], fails=KeyError)
        self.check_result('label range', 'loc', 20, 'ix', 20,
                          typs=['labels'], fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs=['ts'],
                          axes=0, fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs=['floats'],
                          axes=0, fails=TypeError)

    def test_loc_getitem_label_list(self):

        # list of labels
        self.check_result('list lbl', 'loc', [0, 2, 4], 'ix', [0, 2, 4],
                          typs=['ints', 'uints'], axes=0)
        self.check_result('list lbl', 'loc', [3, 6, 9], 'ix', [3, 6, 9],
                          typs=['ints', 'uints'], axes=1)
        self.check_result('list lbl', 'loc', [4, 8, 12], 'ix', [4, 8, 12],
                          typs=['ints', 'uints'], axes=2)
        self.check_result('list lbl', 'loc', ['a', 'b', 'd'], 'ix',
                          ['a', 'b', 'd'], typs=['labels'], axes=0)
        self.check_result('list lbl', 'loc', ['A', 'B', 'C'], 'ix',
                          ['A', 'B', 'C'], typs=['labels'], axes=1)
        self.check_result('list lbl', 'loc', ['Z', 'Y', 'W'], 'ix',
                          ['Z', 'Y', 'W'], typs=['labels'], axes=2)
        self.check_result('list lbl', 'loc', [2, 8, 'null'], 'ix',
                          [2, 8, 'null'], typs=['mixed'], axes=0)
        self.check_result('list lbl', 'loc',
                          [Timestamp('20130102'), Timestamp('20130103')], 'ix',
                          [Timestamp('20130102'), Timestamp('20130103')],
                          typs=['ts'], axes=0)

        self.check_result('list lbl', 'loc', [0, 1, 2], 'indexer', [0, 1, 2],
                          typs=['empty'], fails=KeyError)
        self.check_result('list lbl', 'loc', [0, 2, 3], 'ix', [0, 2, 3],
                          typs=['ints', 'uints'], axes=0, fails=KeyError)
        self.check_result('list lbl', 'loc', [3, 6, 7], 'ix', [3, 6, 7],
                          typs=['ints', 'uints'], axes=1, fails=KeyError)
        self.check_result('list lbl', 'loc', [4, 8, 10], 'ix', [4, 8, 10],
                          typs=['ints', 'uints'], axes=2, fails=KeyError)

    def test_loc_getitem_label_list_fails(self):
        # fails
        self.check_result('list lbl', 'loc', [20, 30, 40], 'ix', [20, 30, 40],
                          typs=['ints', 'uints'], axes=1, fails=KeyError)
        self.check_result('list lbl', 'loc', [20, 30, 40], 'ix', [20, 30, 40],
                          typs=['ints', 'uints'], axes=2, fails=KeyError)

    def test_loc_getitem_label_array_like(self):
        # array like
        self.check_result('array like', 'loc', Series(index=[0, 2, 4]).index,
                          'ix', [0, 2, 4], typs=['ints', 'uints'], axes=0)
        self.check_result('array like', 'loc', Series(index=[3, 6, 9]).index,
                          'ix', [3, 6, 9], typs=['ints', 'uints'], axes=1)
        self.check_result('array like', 'loc', Series(index=[4, 8, 12]).index,
                          'ix', [4, 8, 12], typs=['ints', 'uints'], axes=2)

    def test_loc_getitem_bool(self):
        # boolean indexers
        b = [True, False, True, False]
        self.check_result('bool', 'loc', b, 'ix', b,
                          typs=['ints', 'uints', 'labels',
                                'mixed', 'ts', 'floats'])
        self.check_result('bool', 'loc', b, 'ix', b, typs=['empty'],
                          fails=KeyError)

    def test_loc_getitem_int_slice(self):

        # ok
        self.check_result('int slice2', 'loc', slice(2, 4), 'ix', [2, 4],
                          typs=['ints', 'uints'], axes=0)
        self.check_result('int slice2', 'loc', slice(3, 6), 'ix', [3, 6],
                          typs=['ints', 'uints'], axes=1)
        self.check_result('int slice2', 'loc', slice(4, 8), 'ix', [4, 8],
                          typs=['ints', 'uints'], axes=2)

        # GH 3053
        # loc should treat integer slices like label slices
        from itertools import product

        index = MultiIndex.from_tuples([t for t in product(
            [6, 7, 8], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[6:8, :]
        with catch_warnings(record=True):
            expected = df.ix[6:8, :]
        tm.assert_frame_equal(result, expected)

        index = MultiIndex.from_tuples([t
                                        for t in product(
                                            [10, 20, 30], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[20:30, :]
        with catch_warnings(record=True):
            expected = df.ix[20:30, :]
        tm.assert_frame_equal(result, expected)

        # doc examples
        result = df.loc[10, :]
        with catch_warnings(record=True):
            expected = df.ix[10, :]
        tm.assert_frame_equal(result, expected)

        result = df.loc[:, 10]
        # expected = df.ix[:,10] (this fails)
        expected = df[10]
        tm.assert_frame_equal(result, expected)

    def test_loc_to_fail(self):

        # GH3449
        df = DataFrame(np.random.random((3, 3)),
                       index=['a', 'b', 'c'],
                       columns=['e', 'f', 'g'])

        # raise a KeyError?
        self.assertRaises(KeyError, df.loc.__getitem__,
                          tuple([[1, 2], [1, 2]]))

        # GH  7496
        # loc should not fallback

        s = Series()
        s.loc[1] = 1
        s.loc['a'] = 2

        self.assertRaises(KeyError, lambda: s.loc[-1])
        self.assertRaises(KeyError, lambda: s.loc[[-1, -2]])

        self.assertRaises(KeyError, lambda: s.loc[['4']])

        s.loc[-1] = 3
        result = s.loc[[-1, -2]]
        expected = Series([3, np.nan], index=[-1, -2])
        tm.assert_series_equal(result, expected)

        s['a'] = 2
        self.assertRaises(KeyError, lambda: s.loc[[-2]])

        del s['a']

        def f():
            s.loc[[-2]] = 0

        self.assertRaises(KeyError, f)

        # inconsistency between .loc[values] and .loc[values,:]
        # GH 7999
        df = DataFrame([['a'], ['b']], index=[1, 2], columns=['value'])

        def f():
            df.loc[[3], :]

        self.assertRaises(KeyError, f)

        def f():
            df.loc[[3]]

        self.assertRaises(KeyError, f)

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

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)
        self.check_result('lab slice', 'loc', slice(1, 3),
                          'ix', slice(1, 3),
                          typs=['labels', 'mixed', 'empty', 'ts', 'floats'],
                          fails=TypeError)

        # real label slices
        self.check_result('lab slice', 'loc', slice('a', 'c'),
                          'ix', slice('a', 'c'), typs=['labels'], axes=0)
        self.check_result('lab slice', 'loc', slice('A', 'C'),
                          'ix', slice('A', 'C'), typs=['labels'], axes=1)
        self.check_result('lab slice', 'loc', slice('W', 'Z'),
                          'ix', slice('W', 'Z'), typs=['labels'], axes=2)

        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=0)
        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=1, fails=TypeError)
        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=2, fails=TypeError)

        # GH 14316
        self.check_result('ts slice rev', 'loc', slice('20130104', '20130102'),
                          'indexer', [0, 1, 2], typs=['ts_rev'], axes=0)

        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=0, fails=TypeError)
        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=1, fails=KeyError)
        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=2, fails=KeyError)

        self.check_result('mixed slice', 'loc', slice(2, 4, 2), 'ix', slice(
            2, 4, 2), typs=['mixed'], axes=0, fails=TypeError)

    def test_loc_general(self):

        df = DataFrame(
            np.random.rand(4, 4), columns=['A', 'B', 'C', 'D'],
            index=['A', 'B', 'C', 'D'])

        # want this to work
        result = df.loc[:, "A":"B"].iloc[0:2, :]
        self.assertTrue((result.columns == ['A', 'B']).all())
        self.assertTrue((result.index == ['A', 'B']).all())

        # mixed type
        result = DataFrame({'a': [Timestamp('20130101')], 'b': [1]}).iloc[0]
        expected = Series([Timestamp('20130101'), 1], index=['a', 'b'], name=0)
        tm.assert_series_equal(result, expected)
        self.assertEqual(result.dtype, object)

    def test_loc_setitem_consistency(self):
        # GH 6149
        # coerce similary for setitem and loc when rows have a null-slice
        expected = DataFrame({'date': Series(0, index=range(5),
                                             dtype=np.int64),
                              'val': Series(range(5), dtype=np.int64)})

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(
                            range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 0
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = np.array(0, dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = np.array([0, 0, 0, 0, 0], dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        expected = DataFrame({'date': Series('foo', index=range(5)),
                              'val': Series(range(5), dtype=np.int64)})
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 'foo'
        tm.assert_frame_equal(df, expected)

        expected = DataFrame({'date': Series(1.0, index=range(5)),
                              'val': Series(range(5), dtype=np.int64)})
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 1.0
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_empty(self):
        # empty (essentially noops)
        expected = DataFrame(columns=['x', 'y'])
        expected['x'] = expected['x'].astype(np.int64)
        df = DataFrame(columns=['x', 'y'])
        df.loc[:, 'x'] = 1
        tm.assert_frame_equal(df, expected)

        df = DataFrame(columns=['x', 'y'])
        df['x'] = 1
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_slice_column_len(self):
        # .loc[:,column] setting with slice == len of the column
        # GH10408
        data = """Level_0,,,Respondent,Respondent,Respondent,OtherCat,OtherCat
Level_1,,,Something,StartDate,EndDate,Yes/No,SomethingElse
Region,Site,RespondentID,,,,,
Region_1,Site_1,3987227376,A,5/25/2015 10:59,5/25/2015 11:22,Yes,
Region_1,Site_1,3980680971,A,5/21/2015 9:40,5/21/2015 9:52,Yes,Yes
Region_1,Site_2,3977723249,A,5/20/2015 8:27,5/20/2015 8:41,Yes,
Region_1,Site_2,3977723089,A,5/20/2015 8:33,5/20/2015 9:09,Yes,No"""

        df = pd.read_csv(StringIO(data), header=[0, 1], index_col=[0, 1, 2])
        df.loc[:, ('Respondent', 'StartDate')] = pd.to_datetime(df.loc[:, (
            'Respondent', 'StartDate')])
        df.loc[:, ('Respondent', 'EndDate')] = pd.to_datetime(df.loc[:, (
            'Respondent', 'EndDate')])
        df.loc[:, ('Respondent', 'Duration')] = df.loc[:, (
            'Respondent', 'EndDate')] - df.loc[:, ('Respondent', 'StartDate')]

        df.loc[:, ('Respondent', 'Duration')] = df.loc[:, (
            'Respondent', 'Duration')].astype('timedelta64[s]')
        expected = Series([1380, 720, 840, 2160.], index=df.index,
                          name=('Respondent', 'Duration'))
        tm.assert_series_equal(df[('Respondent', 'Duration')], expected)

    def test_loc_setitem_frame(self):
        df = self.frame_labels

        result = df.iloc[0, 0]

        df.loc['a', 'A'] = 1
        result = df.loc['a', 'A']
        self.assertEqual(result, 1)

        result = df.iloc[0, 0]
        self.assertEqual(result, 1)

        df.loc[:, 'B':'D'] = 0
        expected = df.loc[:, 'B':'D']
        with catch_warnings(record=True):
            result = df.ix[:, 1:]
        tm.assert_frame_equal(result, expected)

        # GH 6254
        # setting issue
        df = DataFrame(index=[3, 5, 4], columns=['A'])
        df.loc[[4, 3, 5], 'A'] = np.array([1, 2, 3], dtype='int64')
        expected = DataFrame(dict(A=Series(
            [1, 2, 3], index=[4, 3, 5]))).reindex(index=[3, 5, 4])
        tm.assert_frame_equal(df, expected)

        # GH 6252
        # setting with an empty frame
        keys1 = ['@' + str(i) for i in range(5)]
        val1 = np.arange(5, dtype='int64')

        keys2 = ['@' + str(i) for i in range(4)]
        val2 = np.arange(4, dtype='int64')

        index = list(set(keys1).union(keys2))
        df = DataFrame(index=index)
        df['A'] = nan
        df.loc[keys1, 'A'] = val1

        df['B'] = nan
        df.loc[keys2, 'B'] = val2

        expected = DataFrame(dict(A=Series(val1, index=keys1), B=Series(
            val2, index=keys2))).reindex(index=index)
        tm.assert_frame_equal(df, expected)

        # GH 8669
        # invalid coercion of nan -> int
        df = DataFrame({'A': [1, 2, 3], 'B': np.nan})
        df.loc[df.B > df.A, 'B'] = df.A
        expected = DataFrame({'A': [1, 2, 3], 'B': np.nan})
        tm.assert_frame_equal(df, expected)

        # GH 6546
        # setting with mixed labels
        df = DataFrame({1: [1, 2], 2: [3, 4], 'a': ['a', 'b']})

        result = df.loc[0, [1, 2]]
        expected = Series([1, 3], index=[1, 2], dtype=object, name=0)
        tm.assert_series_equal(result, expected)

        expected = DataFrame({1: [5, 2], 2: [6, 4], 'a': ['a', 'b']})
        df.loc[0, [1, 2]] = [5, 6]
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_frame_multiples(self):
        # multiple setting
        df = DataFrame({'A': ['foo', 'bar', 'baz'],
                        'B': Series(
                            range(3), dtype=np.int64)})
        rhs = df.loc[1:2]
        rhs.index = df.index[0:2]
        df.loc[0:1] = rhs
        expected = DataFrame({'A': ['bar', 'baz', 'baz'],
                              'B': Series(
                                  [1, 2, 2], dtype=np.int64)})
        tm.assert_frame_equal(df, expected)

        # multiple setting with frame on rhs (with M8)
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(
                            range(5), dtype=np.int64)})
        expected = DataFrame({'date': [Timestamp('20000101'), Timestamp(
            '20000102'), Timestamp('20000101'), Timestamp('20000102'),
            Timestamp('20000103')],
            'val': Series(
            [0, 1, 0, 1, 2], dtype=np.int64)})
        rhs = df.loc[0:2]
        rhs.index = df.index[2:5]
        df.loc[2:4] = rhs
        tm.assert_frame_equal(df, expected)

    def test_iloc_getitem_frame(self):
        df = DataFrame(np.random.randn(10, 4), index=lrange(0, 20, 2),
                       columns=lrange(0, 8, 2))

        result = df.iloc[2]
        with catch_warnings(record=True):
            exp = df.ix[4]
        tm.assert_series_equal(result, exp)

        result = df.iloc[2, 2]
        with catch_warnings(record=True):
            exp = df.ix[4, 4]
        self.assertEqual(result, exp)

        # slice
        result = df.iloc[4:8]
        with catch_warnings(record=True):
            expected = df.ix[8:14]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 2:3]
        with catch_warnings(record=True):
            expected = df.ix[:, 4:5]
        tm.assert_frame_equal(result, expected)

        # list of integers
        result = df.iloc[[0, 1, 3]]
        with catch_warnings(record=True):
            expected = df.ix[[0, 2, 6]]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[[0, 1, 3], [0, 1]]
        with catch_warnings(record=True):
            expected = df.ix[[0, 2, 6], [0, 2]]
        tm.assert_frame_equal(result, expected)

        # neg indicies
        result = df.iloc[[-1, 1, 3], [-1, 1]]
        with catch_warnings(record=True):
            expected = df.ix[[18, 2, 6], [6, 2]]
        tm.assert_frame_equal(result, expected)

        # dups indicies
        result = df.iloc[[-1, -1, 1, 3], [-1, 1]]
        with catch_warnings(record=True):
            expected = df.ix[[18, 18, 2, 6], [6, 2]]
        tm.assert_frame_equal(result, expected)

        # with index-like
        s = Series(index=lrange(1, 5))
        result = df.iloc[s.index]
        with catch_warnings(record=True):
            expected = df.ix[[2, 4, 6, 8]]
        tm.assert_frame_equal(result, expected)

    def test_iloc_getitem_labelled_frame(self):
        # try with labelled frame
        df = DataFrame(np.random.randn(10, 4),
                       index=list('abcdefghij'), columns=list('ABCD'))

        result = df.iloc[1, 1]
        exp = df.loc['b', 'B']
        self.assertEqual(result, exp)

        result = df.iloc[:, 2:3]
        expected = df.loc[:, ['C']]
        tm.assert_frame_equal(result, expected)

        # negative indexing
        result = df.iloc[-1, -1]
        exp = df.loc['j', 'D']
        self.assertEqual(result, exp)

        # out-of-bounds exception
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([10, 5]))

        # trying to use a label
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple(['j', 'D']))

    def test_iloc_getitem_doc_issue(self):

        # multi axis slicing issue with single block
        # surfaced in GH 6059

        arr = np.random.randn(6, 4)
        index = date_range('20130101', periods=6)
        columns = list('ABCD')
        df = DataFrame(arr, index=index, columns=columns)

        # defines ref_locs
        df.describe()

        result = df.iloc[3:5, 0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5, 0:2], index=index[3:5],
                             columns=columns[0:2])
        tm.assert_frame_equal(result, expected)

        # for dups
        df.columns = list('aaaa')
        result = df.iloc[3:5, 0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5, 0:2], index=index[3:5],
                             columns=list('aa'))
        tm.assert_frame_equal(result, expected)

        # related
        arr = np.random.randn(6, 4)
        index = list(range(0, 12, 2))
        columns = list(range(0, 8, 2))
        df = DataFrame(arr, index=index, columns=columns)

        df._data.blocks[0].mgr_locs
        result = df.iloc[1:5, 2:4]
        str(result)
        result.dtypes
        expected = DataFrame(arr[1:5, 2:4], index=index[1:5],
                             columns=columns[2:4])
        tm.assert_frame_equal(result, expected)

    def test_setitem_ndarray_1d(self):
        # GH5508

        # len of indexer vs length of the 1d ndarray
        df = DataFrame(index=Index(lrange(1, 11)))
        df['foo'] = np.zeros(10, dtype=np.float64)
        df['bar'] = np.zeros(10, dtype=np.complex)

        # invalid
        def f():
            with catch_warnings(record=True):
                df.ix[2:5, 'bar'] = np.array([2.33j, 1.23 + 0.1j, 2.2])

        self.assertRaises(ValueError, f)

        def f():
            df.loc[df.index[2:5], 'bar'] = np.array([2.33j, 1.23 + 0.1j,
                                                     2.2, 1.0])

        self.assertRaises(ValueError, f)

        # valid
        df.loc[df.index[2:6], 'bar'] = np.array([2.33j, 1.23 + 0.1j,
                                                 2.2, 1.0])

        result = df.loc[df.index[2:6], 'bar']
        expected = Series([2.33j, 1.23 + 0.1j, 2.2, 1.0], index=[3, 4, 5, 6],
                          name='bar')
        tm.assert_series_equal(result, expected)

        # dtype getting changed?
        df = DataFrame(index=Index(lrange(1, 11)))
        df['foo'] = np.zeros(10, dtype=np.float64)
        df['bar'] = np.zeros(10, dtype=np.complex)

        def f():
            df[2:5] = np.arange(1, 4) * 1j

        self.assertRaises(ValueError, f)

    def test_iloc_setitem_series(self):
        df = DataFrame(np.random.randn(10, 4), index=list('abcdefghij'),
                       columns=list('ABCD'))

        df.iloc[1, 1] = 1
        result = df.iloc[1, 1]
        self.assertEqual(result, 1)

        df.iloc[:, 2:3] = 0
        expected = df.iloc[:, 2:3]
        result = df.iloc[:, 2:3]
        tm.assert_frame_equal(result, expected)

        s = Series(np.random.randn(10), index=lrange(0, 20, 2))

        s.iloc[1] = 1
        result = s.iloc[1]
        self.assertEqual(result, 1)

        s.iloc[:4] = 0
        expected = s.iloc[:4]
        result = s.iloc[:4]
        tm.assert_series_equal(result, expected)

        s = Series([-1] * 6)
        s.iloc[0::2] = [0, 2, 4]
        s.iloc[1::2] = [1, 3, 5]
        result = s
        expected = Series([0, 1, 2, 3, 4, 5])
        tm.assert_series_equal(result, expected)

    def test_iloc_setitem_list_of_lists(self):

        # GH 7551
        # list-of-list is set incorrectly in mixed vs. single dtyped frames
        df = DataFrame(dict(A=np.arange(5, dtype='int64'),
                            B=np.arange(5, 10, dtype='int64')))
        df.iloc[2:4] = [[10, 11], [12, 13]]
        expected = DataFrame(dict(A=[0, 1, 10, 12, 4], B=[5, 6, 11, 13, 9]))
        tm.assert_frame_equal(df, expected)

        df = DataFrame(
            dict(A=list('abcde'), B=np.arange(5, 10, dtype='int64')))
        df.iloc[2:4] = [['x', 11], ['y', 13]]
        expected = DataFrame(dict(A=['a', 'b', 'x', 'y', 'e'],
                                  B=[5, 6, 11, 13, 9]))
        tm.assert_frame_equal(df, expected)

    def test_ix_general(self):

        # ix general issues

        # GH 2817
        data = {'amount': {0: 700, 1: 600, 2: 222, 3: 333, 4: 444},
                'col': {0: 3.5, 1: 3.5, 2: 4.0, 3: 4.0, 4: 4.0},
                'year': {0: 2012, 1: 2011, 2: 2012, 3: 2012, 4: 2012}}
        df = DataFrame(data).set_index(keys=['col', 'year'])
        key = 4.0, 2012

        # emits a PerformanceWarning, ok
        with self.assert_produces_warning(PerformanceWarning):
            tm.assert_frame_equal(df.loc[key], df.iloc[2:])

        # this is ok
        df.sort_index(inplace=True)
        res = df.loc[key]

        # col has float dtype, result should be Float64Index
        index = MultiIndex.from_arrays([[4.] * 3, [2012] * 3],
                                       names=['col', 'year'])
        expected = DataFrame({'amount': [222, 333, 444]}, index=index)
        tm.assert_frame_equal(res, expected)

    def test_ix_weird_slicing(self):
        # http://stackoverflow.com/q/17056560/1240268
        df = DataFrame({'one': [1, 2, 3, np.nan, np.nan],
                        'two': [1, 2, 3, 4, 5]})
        df.loc[df['one'] > 1, 'two'] = -df['two']

        expected = DataFrame({'one': {0: 1.0,
                                      1: 2.0,
                                      2: 3.0,
                                      3: nan,
                                      4: nan},
                              'two': {0: 1,
                                      1: -2,
                                      2: -3,
                                      3: 4,
                                      4: 5}})
        tm.assert_frame_equal(df, expected)

    def test_loc_coerceion(self):

        # 12411
        df = DataFrame({'date': [pd.Timestamp('20130101').tz_localize('UTC'),
                                 pd.NaT]})
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 12045
        import datetime
        df = DataFrame({'date': [datetime.datetime(2012, 1, 1),
                                 datetime.datetime(1012, 1, 2)]})
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 11594
        df = DataFrame({'text': ['some words'] + [None] * 9})
        expected = df.dtypes

        result = df.iloc[0:2]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[3:]
        tm.assert_series_equal(result.dtypes, expected)

    def test_setitem_dtype_upcast(self):

        # GH3216
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df['c'] = np.nan
        self.assertEqual(df['c'].dtype, np.float64)

        df.loc[0, 'c'] = 'foo'
        expected = DataFrame([{"a": 1, "c": 'foo'},
                              {"a": 3, "b": 2, "c": np.nan}])
        tm.assert_frame_equal(df, expected)

        # GH10280
        df = DataFrame(np.arange(6, dtype='int64').reshape(2, 3),
                       index=list('ab'),
                       columns=['foo', 'bar', 'baz'])

        for val in [3.14, 'wxyz']:
            left = df.copy()
            left.loc['a', 'bar'] = val
            right = DataFrame([[0, val, 2], [3, 4, 5]], index=list('ab'),
                              columns=['foo', 'bar', 'baz'])

            tm.assert_frame_equal(left, right)
            self.assertTrue(is_integer_dtype(left['foo']))
            self.assertTrue(is_integer_dtype(left['baz']))

        left = DataFrame(np.arange(6, dtype='int64').reshape(2, 3) / 10.0,
                         index=list('ab'),
                         columns=['foo', 'bar', 'baz'])
        left.loc['a', 'bar'] = 'wxyz'

        right = DataFrame([[0, 'wxyz', .2], [.3, .4, .5]], index=list('ab'),
                          columns=['foo', 'bar', 'baz'])

        tm.assert_frame_equal(left, right)
        self.assertTrue(is_float_dtype(left['foo']))
        self.assertTrue(is_float_dtype(left['baz']))

    def test_setitem_iloc(self):

        # setitem with an iloc list
        df = DataFrame(np.arange(9).reshape((3, 3)), index=["A", "B", "C"],
                       columns=["A", "B", "C"])
        df.iloc[[0, 1], [1, 2]]
        df.iloc[[0, 1], [1, 2]] += 100

        expected = DataFrame(
            np.array([0, 101, 102, 3, 104, 105, 6, 7, 8]).reshape((3, 3)),
            index=["A", "B", "C"], columns=["A", "B", "C"])
        tm.assert_frame_equal(df, expected)

    def test_dups_fancy_indexing(self):

        # GH 3455
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(10, 3)
        df.columns = ['a', 'a', 'b']
        result = df[['b', 'a']].columns
        expected = Index(['b', 'a', 'a'])
        self.assert_index_equal(result, expected)

        # across dtypes
        df = DataFrame([[1, 2, 1., 2., 3., 'foo', 'bar']],
                       columns=list('aaaaaaa'))
        df.head()
        str(df)
        result = DataFrame([[1, 2, 1., 2., 3., 'foo', 'bar']])
        result.columns = list('aaaaaaa')

        # TODO(wesm): unused?
        df_v = df.iloc[:, 4]  # noqa
        res_v = result.iloc[:, 4]  # noqa

        tm.assert_frame_equal(df, result)

        # GH 3561, dups not in selected order
        df = DataFrame(
            {'test': [5, 7, 9, 11],
             'test1': [4., 5, 6, 7],
             'other': list('abcd')}, index=['A', 'A', 'B', 'C'])
        rows = ['C', 'B']
        expected = DataFrame(
            {'test': [11, 9],
             'test1': [7., 6],
             'other': ['d', 'c']}, index=rows)
        result = df.loc[rows]
        tm.assert_frame_equal(result, expected)

        result = df.loc[Index(rows)]
        tm.assert_frame_equal(result, expected)

        rows = ['C', 'B', 'E']
        expected = DataFrame(
            {'test': [11, 9, np.nan],
             'test1': [7., 6, np.nan],
             'other': ['d', 'c', np.nan]}, index=rows)

        result = df.loc[rows]
        tm.assert_frame_equal(result, expected)

        # see GH5553, make sure we use the right indexer
        rows = ['F', 'G', 'H', 'C', 'B', 'E']
        expected = DataFrame({'test': [np.nan, np.nan, np.nan, 11, 9, np.nan],
                              'test1': [np.nan, np.nan, np.nan, 7., 6, np.nan],
                              'other': [np.nan, np.nan, np.nan,
                                        'd', 'c', np.nan]},
                             index=rows)
        result = df.loc[rows]
        tm.assert_frame_equal(result, expected)

        # inconsistent returns for unique/duplicate indices when values are
        # missing
        df = DataFrame(randn(4, 3), index=list('ABCD'))
        expected = df.ix[['E']]

        dfnu = DataFrame(randn(5, 3), index=list('AABCD'))
        result = dfnu.ix[['E']]
        tm.assert_frame_equal(result, expected)

        # ToDo: check_index_type can be True after GH 11497

        # GH 4619; duplicate indexer with missing label
        df = DataFrame({"A": [0, 1, 2]})
        result = df.ix[[0, 8, 0]]
        expected = DataFrame({"A": [0, np.nan, 0]}, index=[0, 8, 0])
        tm.assert_frame_equal(result, expected, check_index_type=False)

        df = DataFrame({"A": list('abc')})
        result = df.ix[[0, 8, 0]]
        expected = DataFrame({"A": ['a', np.nan, 'a']}, index=[0, 8, 0])
        tm.assert_frame_equal(result, expected, check_index_type=False)

        # non unique with non unique selector
        df = DataFrame({'test': [5, 7, 9, 11]}, index=['A', 'A', 'B', 'C'])
        expected = DataFrame(
            {'test': [5, 7, 5, 7, np.nan]}, index=['A', 'A', 'A', 'A', 'E'])
        result = df.ix[['A', 'A', 'E']]
        tm.assert_frame_equal(result, expected)

        # GH 5835
        # dups on index and missing values
        df = DataFrame(
            np.random.randn(5, 5), columns=['A', 'B', 'B', 'B', 'A'])

        expected = pd.concat(
            [df.ix[:, ['A', 'B']], DataFrame(np.nan, columns=['C'],
                                             index=df.index)], axis=1)
        result = df.ix[:, ['A', 'B', 'C']]
        tm.assert_frame_equal(result, expected)

        # GH 6504, multi-axis indexing
        df = DataFrame(np.random.randn(9, 2),
                       index=[1, 1, 1, 2, 2, 2, 3, 3, 3], columns=['a', 'b'])

        expected = df.iloc[0:6]
        result = df.loc[[1, 2]]
        tm.assert_frame_equal(result, expected)

        expected = df
        result = df.loc[:, ['a', 'b']]
        tm.assert_frame_equal(result, expected)

        expected = df.iloc[0:6, :]
        result = df.loc[[1, 2], ['a', 'b']]
        tm.assert_frame_equal(result, expected)

    def test_indexing_mixed_frame_bug(self):

        # GH3492
        df = DataFrame({'a': {1: 'aaa', 2: 'bbb', 3: 'ccc'},
                        'b': {1: 111, 2: 222, 3: 333}})

        # this works, new column is created correctly
        df['test'] = df['a'].apply(lambda x: '_' if x == 'aaa' else x)

        # this does not work, ie column test is not changed
        idx = df['test'] == '_'
        temp = df.ix[idx, 'a'].apply(lambda x: '-----' if x == 'aaa' else x)
        df.ix[idx, 'test'] = temp
        self.assertEqual(df.iloc[0, 2], '-----')

        # if I look at df, then element [0,2] equals '_'. If instead I type
        # df.ix[idx,'test'], I get '-----', finally by typing df.iloc[0,2] I
        # get '_'.

    def test_multitype_list_index_access(self):
        # GH 10610
        df = pd.DataFrame(np.random.random((10, 5)),
                          columns=["a"] + [20, 21, 22, 23])

        with self.assertRaises(KeyError):
            df[[22, 26, -8]]
        self.assertEqual(df[21].shape[0], df.shape[0])

    def test_set_index_nan(self):

        # GH 3586
        df = DataFrame({'PRuid': {17: 'nonQC',
                                  18: 'nonQC',
                                  19: 'nonQC',
                                  20: '10',
                                  21: '11',
                                  22: '12',
                                  23: '13',
                                  24: '24',
                                  25: '35',
                                  26: '46',
                                  27: '47',
                                  28: '48',
                                  29: '59',
                                  30: '10'},
                        'QC': {17: 0.0,
                               18: 0.0,
                               19: 0.0,
                               20: nan,
                               21: nan,
                               22: nan,
                               23: nan,
                               24: 1.0,
                               25: nan,
                               26: nan,
                               27: nan,
                               28: nan,
                               29: nan,
                               30: nan},
                        'data': {17: 7.9544899999999998,
                                 18: 8.0142609999999994,
                                 19: 7.8591520000000008,
                                 20: 0.86140349999999999,
                                 21: 0.87853110000000001,
                                 22: 0.8427041999999999,
                                 23: 0.78587700000000005,
                                 24: 0.73062459999999996,
                                 25: 0.81668560000000001,
                                 26: 0.81927080000000008,
                                 27: 0.80705009999999999,
                                 28: 0.81440240000000008,
                                 29: 0.80140849999999997,
                                 30: 0.81307740000000006},
                        'year': {17: 2006,
                                 18: 2007,
                                 19: 2008,
                                 20: 1985,
                                 21: 1985,
                                 22: 1985,
                                 23: 1985,
                                 24: 1985,
                                 25: 1985,
                                 26: 1985,
                                 27: 1985,
                                 28: 1985,
                                 29: 1985,
                                 30: 1986}}).reset_index()

        result = df.set_index(['year', 'PRuid', 'QC']).reset_index().reindex(
            columns=df.columns)
        tm.assert_frame_equal(result, df)

    def test_multi_nan_indexing(self):

        # GH 3588
        df = DataFrame({"a": ['R1', 'R2', np.nan, 'R4'],
                        'b': ["C1", "C2", "C3", "C4"],
                        "c": [10, 15, np.nan, 20]})
        result = df.set_index(['a', 'b'], drop=False)
        expected = DataFrame({"a": ['R1', 'R2', np.nan, 'R4'],
                              'b': ["C1", "C2", "C3", "C4"],
                              "c": [10, 15, np.nan, 20]},
                             index=[Index(['R1', 'R2', np.nan, 'R4'],
                                          name='a'),
                                    Index(['C1', 'C2', 'C3', 'C4'], name='b')])
        tm.assert_frame_equal(result, expected)

    def test_multi_assign(self):

        # GH 3626, an assignement of a sub-df to a df
        df = DataFrame({'FC': ['a', 'b', 'a', 'b', 'a', 'b'],
                        'PF': [0, 0, 0, 0, 1, 1],
                        'col1': lrange(6),
                        'col2': lrange(6, 12)})
        df.ix[1, 0] = np.nan
        df2 = df.copy()

        mask = ~df2.FC.isnull()
        cols = ['col1', 'col2']

        dft = df2 * 2
        dft.ix[3, 3] = np.nan

        expected = DataFrame({'FC': ['a', np.nan, 'a', 'b', 'a', 'b'],
                              'PF': [0, 0, 0, 0, 1, 1],
                              'col1': Series([0, 1, 4, 6, 8, 10]),
                              'col2': [12, 7, 16, np.nan, 20, 22]})

        # frame on rhs
        df2.ix[mask, cols] = dft.ix[mask, cols]
        tm.assert_frame_equal(df2, expected)

        df2.ix[mask, cols] = dft.ix[mask, cols]
        tm.assert_frame_equal(df2, expected)

        # with an ndarray on rhs
        df2 = df.copy()
        df2.ix[mask, cols] = dft.ix[mask, cols].values
        tm.assert_frame_equal(df2, expected)
        df2.ix[mask, cols] = dft.ix[mask, cols].values
        tm.assert_frame_equal(df2, expected)

        # broadcasting on the rhs is required
        df = DataFrame(dict(A=[1, 2, 0, 0, 0], B=[0, 0, 0, 10, 11], C=[
                       0, 0, 0, 10, 11], D=[3, 4, 5, 6, 7]))

        expected = df.copy()
        mask = expected['A'] == 0
        for col in ['A', 'B']:
            expected.loc[mask, col] = df['D']

        df.loc[df['A'] == 0, ['A', 'B']] = df['D']
        tm.assert_frame_equal(df, expected)

    def test_ix_assign_column_mixed(self):
        # GH #1142
        df = DataFrame(tm.getSeriesData())
        df['foo'] = 'bar'

        orig = df.ix[:, 'B'].copy()
        df.ix[:, 'B'] = df.ix[:, 'B'] + 1
        tm.assert_series_equal(df.B, orig + 1)

        # GH 3668, mixed frame with series value
        df = DataFrame({'x': lrange(10), 'y': lrange(10, 20), 'z': 'bar'})
        expected = df.copy()

        for i in range(5):
            indexer = i * 2
            v = 1000 + i * 200
            expected.ix[indexer, 'y'] = v
            self.assertEqual(expected.ix[indexer, 'y'], v)

        df.ix[df.x % 2 == 0, 'y'] = df.ix[df.x % 2 == 0, 'y'] * 100
        tm.assert_frame_equal(df, expected)

        # GH 4508, making sure consistency of assignments
        df = DataFrame({'a': [1, 2, 3], 'b': [0, 1, 2]})
        df.ix[[0, 2, ], 'b'] = [100, -100]
        expected = DataFrame({'a': [1, 2, 3], 'b': [100, 1, -100]})
        tm.assert_frame_equal(df, expected)

        df = pd.DataFrame({'a': lrange(4)})
        df['b'] = np.nan
        df.ix[[1, 3], 'b'] = [100, -100]
        expected = DataFrame({'a': [0, 1, 2, 3],
                              'b': [np.nan, 100, np.nan, -100]})
        tm.assert_frame_equal(df, expected)

        # ok, but chained assignments are dangerous
        # if we turn off chained assignement it will work
        with option_context('chained_assignment', None):
            df = pd.DataFrame({'a': lrange(4)})
            df['b'] = np.nan
            df['b'].ix[[1, 3]] = [100, -100]
            tm.assert_frame_equal(df, expected)

    def test_ix_get_set_consistency(self):

        # GH 4544
        # ix/loc get/set not consistent when
        # a mixed int/string index
        df = DataFrame(np.arange(16).reshape((4, 4)),
                       columns=['a', 'b', 8, 'c'],
                       index=['e', 7, 'f', 'g'])

        self.assertEqual(df.ix['e', 8], 2)
        self.assertEqual(df.loc['e', 8], 2)

        df.ix['e', 8] = 42
        self.assertEqual(df.ix['e', 8], 42)
        self.assertEqual(df.loc['e', 8], 42)

        df.loc['e', 8] = 45
        self.assertEqual(df.ix['e', 8], 45)
        self.assertEqual(df.loc['e', 8], 45)

    def test_setitem_list(self):

        # GH 6043
        # ix with a list
        df = DataFrame(index=[0, 1], columns=[0])
        df.ix[1, 0] = [1, 2, 3]
        df.ix[1, 0] = [1, 2]

        result = DataFrame(index=[0, 1], columns=[0])
        result.ix[1, 0] = [1, 2]

        tm.assert_frame_equal(result, df)

        # ix with an object
        class TO(object):

            def __init__(self, value):
                self.value = value

            def __str__(self):
                return "[{0}]".format(self.value)

            __repr__ = __str__

            def __eq__(self, other):
                return self.value == other.value

            def view(self):
                return self

        df = DataFrame(index=[0, 1], columns=[0])
        df.ix[1, 0] = TO(1)
        df.ix[1, 0] = TO(2)

        result = DataFrame(index=[0, 1], columns=[0])
        result.ix[1, 0] = TO(2)

        tm.assert_frame_equal(result, df)

        # remains object dtype even after setting it back
        df = DataFrame(index=[0, 1], columns=[0])
        df.ix[1, 0] = TO(1)
        df.ix[1, 0] = np.nan
        result = DataFrame(index=[0, 1], columns=[0])

        tm.assert_frame_equal(result, df)

    def test_iloc_mask(self):

        # GH 3631, iloc with a mask (of a series) should raise
        df = DataFrame(lrange(5), list('ABCDE'), columns=['a'])
        mask = (df.a % 2 == 0)
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple([mask]))
        mask.index = lrange(len(mask))
        self.assertRaises(NotImplementedError, df.iloc.__getitem__,
                          tuple([mask]))

        # ndarray ok
        result = df.iloc[np.array([True] * len(mask), dtype=bool)]
        tm.assert_frame_equal(result, df)

        # the possibilities
        locs = np.arange(4)
        nums = 2 ** locs
        reps = lmap(bin, nums)
        df = DataFrame({'locs': locs, 'nums': nums}, reps)

        expected = {
            (None, ''): '0b1100',
            (None, '.loc'): '0b1100',
            (None, '.iloc'): '0b1100',
            ('index', ''): '0b11',
            ('index', '.loc'): '0b11',
            ('index', '.iloc'): ('iLocation based boolean indexing '
                                 'cannot use an indexable as a mask'),
            ('locs', ''): 'Unalignable boolean Series provided as indexer '
                          '(index of the boolean Series and of the indexed '
                          'object do not match',
            ('locs', '.loc'): 'Unalignable boolean Series provided as indexer '
                              '(index of the boolean Series and of the '
                              'indexed object do not match',
            ('locs', '.iloc'): ('iLocation based boolean indexing on an '
                                'integer type is not available'),
        }

        # UserWarnings from reindex of a boolean mask
        with warnings.catch_warnings(record=True):
            result = dict()
            for idx in [None, 'index', 'locs']:
                mask = (df.nums > 2).values
                if idx:
                    mask = Series(mask, list(reversed(getattr(df, idx))))
                for method in ['', '.loc', '.iloc']:
                    try:
                        if method:
                            accessor = getattr(df, method[1:])
                        else:
                            accessor = df
                        ans = str(bin(accessor[mask]['nums'].sum()))
                    except Exception as e:
                        ans = str(e)

                    key = tuple([idx, method])
                    r = expected.get(key)
                    if r != ans:
                        raise AssertionError(
                            "[%s] does not match [%s], received [%s]"
                            % (key, ans, r))

    def test_ix_slicing_strings(self):
        # GH3836
        data = {'Classification':
                ['SA EQUITY CFD', 'bbb', 'SA EQUITY', 'SA SSF', 'aaa'],
                'Random': [1, 2, 3, 4, 5],
                'X': ['correct', 'wrong', 'correct', 'correct', 'wrong']}
        df = DataFrame(data)
        x = df[~df.Classification.isin(['SA EQUITY CFD', 'SA EQUITY', 'SA SSF'
                                        ])]
        df.ix[x.index, 'X'] = df['Classification']

        expected = DataFrame({'Classification': {0: 'SA EQUITY CFD',
                                                 1: 'bbb',
                                                 2: 'SA EQUITY',
                                                 3: 'SA SSF',
                                                 4: 'aaa'},
                              'Random': {0: 1,
                                         1: 2,
                                         2: 3,
                                         3: 4,
                                         4: 5},
                              'X': {0: 'correct',
                                    1: 'bbb',
                                    2: 'correct',
                                    3: 'correct',
                                    4: 'aaa'}})  # bug was 4: 'bbb'

        tm.assert_frame_equal(df, expected)

    def test_non_unique_loc(self):
        # GH3659
        # non-unique indexer with loc slice
        # https://groups.google.com/forum/?fromgroups#!topic/pydata/zTm2No0crYs

        # these are going to raise becuase the we are non monotonic
        df = DataFrame({'A': [1, 2, 3, 4, 5, 6],
                        'B': [3, 4, 5, 6, 7, 8]}, index=[0, 1, 0, 1, 2, 3])
        self.assertRaises(KeyError, df.loc.__getitem__,
                          tuple([slice(1, None)]))
        self.assertRaises(KeyError, df.loc.__getitem__,
                          tuple([slice(0, None)]))
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(1, 2)]))

        # monotonic are ok
        df = DataFrame({'A': [1, 2, 3, 4, 5, 6],
                        'B': [3, 4, 5, 6, 7, 8]},
                       index=[0, 1, 0, 1, 2, 3]).sort_index(axis=0)
        result = df.loc[1:]
        expected = DataFrame({'A': [2, 4, 5, 6], 'B': [4, 6, 7, 8]},
                             index=[1, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

        result = df.loc[0:]
        tm.assert_frame_equal(result, df)

        result = df.loc[1:2]
        expected = DataFrame({'A': [2, 4, 5], 'B': [4, 6, 7]},
                             index=[1, 1, 2])
        tm.assert_frame_equal(result, expected)

    def test_loc_name(self):
        # GH 3880
        df = DataFrame([[1, 1], [1, 1]])
        df.index.name = 'index_name'
        result = df.iloc[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

        result = df.ix[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

        result = df.loc[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

    def test_iloc_non_unique_indexing(self):

        # GH 4017, non-unique indexing (on the axis)
        df = DataFrame({'A': [0.1] * 3000, 'B': [1] * 3000})
        idx = np.array(lrange(30)) * 99
        expected = df.iloc[idx]

        df3 = pd.concat([df, 2 * df, 3 * df])
        result = df3.iloc[idx]

        tm.assert_frame_equal(result, expected)

        df2 = DataFrame({'A': [0.1] * 1000, 'B': [1] * 1000})
        df2 = pd.concat([df2, 2 * df2, 3 * df2])

        sidx = df2.index.to_series()
        expected = df2.iloc[idx[idx <= sidx.max()]]

        new_list = []
        for r, s in expected.iterrows():
            new_list.append(s)
            new_list.append(s * 2)
            new_list.append(s * 3)

        expected = DataFrame(new_list)
        expected = pd.concat([expected, DataFrame(index=idx[idx > sidx.max()])
                              ])
        result = df2.loc[idx]
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_string_slice(self):
        # GH 14424
        # string indexing against datetimelike with object
        # dtype should properly raises KeyError
        df = pd.DataFrame([1], pd.Index([pd.Timestamp('2011-01-01')],
                                        dtype=object))
        self.assertTrue(df.index.is_all_dates)
        with tm.assertRaises(KeyError):
            df['2011']

        with tm.assertRaises(KeyError):
            df.loc['2011', 0]

        df = pd.DataFrame()
        self.assertFalse(df.index.is_all_dates)
        with tm.assertRaises(KeyError):
            df['2011']

        with tm.assertRaises(KeyError):
            df.loc['2011', 0]

    def test_mi_access(self):

        # GH 4145
        data = """h1 main  h3 sub  h5
0  a    A   1  A1   1
1  b    B   2  B1   2
2  c    B   3  A1   3
3  d    A   4  B2   4
4  e    A   5  B2   5
5  f    B   6  A2   6
"""

        df = pd.read_csv(StringIO(data), sep=r'\s+', index_col=0)
        df2 = df.set_index(['main', 'sub']).T.sort_index(1)
        index = Index(['h1', 'h3', 'h5'])
        columns = MultiIndex.from_tuples([('A', 'A1')], names=['main', 'sub'])
        expected = DataFrame([['a', 1, 1]], index=columns, columns=index).T

        result = df2.loc[:, ('A', 'A1')]
        tm.assert_frame_equal(result, expected)

        result = df2[('A', 'A1')]
        tm.assert_frame_equal(result, expected)

        # GH 4146, not returning a block manager when selecting a unique index
        # from a duplicate index
        # as of 4879, this returns a Series (which is similar to what happens
        # with a non-unique)
        expected = Series(['a', 1, 1], index=['h1', 'h3', 'h5'], name='A1')
        result = df2['A']['A1']
        tm.assert_series_equal(result, expected)

        # selecting a non_unique from the 2nd level
        expected = DataFrame([['d', 4, 4], ['e', 5, 5]],
                             index=Index(['B2', 'B2'], name='sub'),
                             columns=['h1', 'h3', 'h5'], ).T
        result = df2['A']['B2']
        tm.assert_frame_equal(result, expected)

    def test_non_unique_loc_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list('ABCDEFG')

        def gen_test(l, l2):
            return pd.concat([DataFrame(randn(l, len(columns)),
                                        index=lrange(l), columns=columns),
                              DataFrame(np.ones((l2, len(columns))),
                                        index=[0] * l2, columns=columns)])

        def gen_expected(df, mask):
            l = len(mask)
            return pd.concat([df.take([0], convert=False),
                              DataFrame(np.ones((l, len(columns))),
                                        index=[0] * l,
                                        columns=columns),
                              df.take(mask[1:], convert=False)])

        df = gen_test(900, 100)
        self.assertFalse(df.index.is_unique)

        mask = np.arange(100)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

        df = gen_test(900000, 100000)
        self.assertFalse(df.index.is_unique)

        mask = np.arange(100000)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

    def test_astype_assignment(self):

        # GH4312 (iloc)
        df_orig = DataFrame([['1', '2', '3', '.4', 5, 6., 'foo']],
                            columns=list('ABCDEFG'))

        df = df_orig.copy()
        df.iloc[:, 0:2] = df.iloc[:, 0:2].astype(np.int64)
        expected = DataFrame([[1, 2, '3', '.4', 5, 6., 'foo']],
                             columns=list('ABCDEFG'))
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.iloc[:, 0:2] = df.iloc[:, 0:2]._convert(datetime=True, numeric=True)
        expected = DataFrame([[1, 2, '3', '.4', 5, 6., 'foo']],
                             columns=list('ABCDEFG'))
        tm.assert_frame_equal(df, expected)

        # GH5702 (loc)
        df = df_orig.copy()
        df.loc[:, 'A'] = df.loc[:, 'A'].astype(np.int64)
        expected = DataFrame([[1, '2', '3', '.4', 5, 6., 'foo']],
                             columns=list('ABCDEFG'))
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[:, ['B', 'C']] = df.loc[:, ['B', 'C']].astype(np.int64)
        expected = DataFrame([['1', 2, 3, '.4', 5, 6., 'foo']],
                             columns=list('ABCDEFG'))
        tm.assert_frame_equal(df, expected)

        # full replacements / no nans
        df = DataFrame({'A': [1., 2., 3., 4.]})
        df.iloc[:, 0] = df['A'].astype(np.int64)
        expected = DataFrame({'A': [1, 2, 3, 4]})
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'A': [1., 2., 3., 4.]})
        df.loc[:, 'A'] = df['A'].astype(np.int64)
        expected = DataFrame({'A': [1, 2, 3, 4]})
        tm.assert_frame_equal(df, expected)

    def test_astype_assignment_with_dups(self):

        # GH 4686
        # assignment with dups that has a dtype change
        cols = pd.MultiIndex.from_tuples([('A', '1'), ('B', '1'), ('A', '2')])
        df = DataFrame(np.arange(3).reshape((1, 3)),
                       columns=cols, dtype=object)
        index = df.index.copy()

        df['A'] = df['A'].astype(np.float64)
        self.assert_index_equal(df.index, index)

        # TODO(wesm): unused variables
        # result = df.get_dtype_counts().sort_index()
        # expected = Series({'float64': 2, 'object': 1}).sort_index()

    def test_dups_loc(self):

        # GH4726
        # dup indexing with iloc/loc
        df = DataFrame([[1, 2, 'foo', 'bar', Timestamp('20130101')]],
                       columns=['a', 'a', 'a', 'a', 'a'], index=[1])
        expected = Series([1, 2, 'foo', 'bar', Timestamp('20130101')],
                          index=['a', 'a', 'a', 'a', 'a'], name=1)

        result = df.iloc[0]
        tm.assert_series_equal(result, expected)

        result = df.loc[1]
        tm.assert_series_equal(result, expected)

    def test_partial_setting(self):

        # GH2578, allow ix and friends to partially set

        # series
        s_orig = Series([1, 2, 3])

        s = s_orig.copy()
        s[5] = 5
        expected = Series([1, 2, 3, 5], index=[0, 1, 2, 5])
        tm.assert_series_equal(s, expected)

        s = s_orig.copy()
        s.loc[5] = 5
        expected = Series([1, 2, 3, 5], index=[0, 1, 2, 5])
        tm.assert_series_equal(s, expected)

        s = s_orig.copy()
        s[5] = 5.
        expected = Series([1, 2, 3, 5.], index=[0, 1, 2, 5])
        tm.assert_series_equal(s, expected)

        s = s_orig.copy()
        s.loc[5] = 5.
        expected = Series([1, 2, 3, 5.], index=[0, 1, 2, 5])
        tm.assert_series_equal(s, expected)

        # iloc/iat raise
        s = s_orig.copy()

        def f():
            s.iloc[3] = 5.

        self.assertRaises(IndexError, f)

        def f():
            s.iat[3] = 5.

        self.assertRaises(IndexError, f)

        # ## frame ##

        df_orig = DataFrame(
            np.arange(6).reshape(3, 2), columns=['A', 'B'], dtype='int64')

        # iloc/iat raise
        df = df_orig.copy()

        def f():
            df.iloc[4, 2] = 5.

        self.assertRaises(IndexError, f)

        def f():
            df.iat[4, 2] = 5.

        self.assertRaises(IndexError, f)

        # row setting where it exists
        expected = DataFrame(dict({'A': [0, 4, 4], 'B': [1, 5, 5]}))
        df = df_orig.copy()
        df.iloc[1] = df.iloc[2]
        tm.assert_frame_equal(df, expected)

        expected = DataFrame(dict({'A': [0, 4, 4], 'B': [1, 5, 5]}))
        df = df_orig.copy()
        df.loc[1] = df.loc[2]
        tm.assert_frame_equal(df, expected)

        # like 2578, partial setting with dtype preservation
        expected = DataFrame(dict({'A': [0, 2, 4, 4], 'B': [1, 3, 5, 5]}))
        df = df_orig.copy()
        df.loc[3] = df.loc[2]
        tm.assert_frame_equal(df, expected)

        # single dtype frame, overwrite
        expected = DataFrame(dict({'A': [0, 2, 4], 'B': [0, 2, 4]}))
        df = df_orig.copy()
        df.ix[:, 'B'] = df.ix[:, 'A']
        tm.assert_frame_equal(df, expected)

        # mixed dtype frame, overwrite
        expected = DataFrame(dict({'A': [0, 2, 4], 'B': Series([0, 2, 4])}))
        df = df_orig.copy()
        df['B'] = df['B'].astype(np.float64)
        df.ix[:, 'B'] = df.ix[:, 'A']
        tm.assert_frame_equal(df, expected)

        # single dtype frame, partial setting
        expected = df_orig.copy()
        expected['C'] = df['A']
        df = df_orig.copy()
        df.ix[:, 'C'] = df.ix[:, 'A']
        tm.assert_frame_equal(df, expected)

        # mixed frame, partial setting
        expected = df_orig.copy()
        expected['C'] = df['A']
        df = df_orig.copy()
        df.ix[:, 'C'] = df.ix[:, 'A']
        tm.assert_frame_equal(df, expected)

        # ## panel ##
        p_orig = Panel(np.arange(16).reshape(2, 4, 2),
                       items=['Item1', 'Item2'],
                       major_axis=pd.date_range('2001/1/12', periods=4),
                       minor_axis=['A', 'B'], dtype='float64')

        # panel setting via item
        p_orig = Panel(np.arange(16).reshape(2, 4, 2),
                       items=['Item1', 'Item2'],
                       major_axis=pd.date_range('2001/1/12', periods=4),
                       minor_axis=['A', 'B'], dtype='float64')
        expected = p_orig.copy()
        expected['Item3'] = expected['Item1']
        p = p_orig.copy()
        p.loc['Item3'] = p['Item1']
        tm.assert_panel_equal(p, expected)

        # panel with aligned series
        expected = p_orig.copy()
        expected = expected.transpose(2, 1, 0)
        expected['C'] = DataFrame({'Item1': [30, 30, 30, 30],
                                   'Item2': [32, 32, 32, 32]},
                                  index=p_orig.major_axis)
        expected = expected.transpose(2, 1, 0)
        p = p_orig.copy()
        p.loc[:, :, 'C'] = Series([30, 32], index=p_orig.items)
        tm.assert_panel_equal(p, expected)

        # GH 8473
        dates = date_range('1/1/2000', periods=8)
        df_orig = DataFrame(np.random.randn(8, 4), index=dates,
                            columns=['A', 'B', 'C', 'D'])

        expected = pd.concat([df_orig, DataFrame(
            {'A': 7}, index=[dates[-1] + 1])])
        df = df_orig.copy()
        df.loc[dates[-1] + 1, 'A'] = 7
        tm.assert_frame_equal(df, expected)
        df = df_orig.copy()
        df.at[dates[-1] + 1, 'A'] = 7
        tm.assert_frame_equal(df, expected)

        exp_other = DataFrame({0: 7}, index=[dates[-1] + 1])
        expected = pd.concat([df_orig, exp_other], axis=1)

        df = df_orig.copy()
        df.loc[dates[-1] + 1, 0] = 7
        tm.assert_frame_equal(df, expected)
        df = df_orig.copy()
        df.at[dates[-1] + 1, 0] = 7
        tm.assert_frame_equal(df, expected)

    def test_partial_setting_mixed_dtype(self):

        # in a mixed dtype environment, try to preserve dtypes
        # by appending
        df = DataFrame([[True, 1], [False, 2]], columns=["female", "fitness"])

        s = df.loc[1].copy()
        s.name = 2
        expected = df.append(s)

        df.loc[2] = df.loc[1]
        tm.assert_frame_equal(df, expected)

        # columns will align
        df = DataFrame(columns=['A', 'B'])
        df.loc[0] = Series(1, index=range(4))
        tm.assert_frame_equal(df, DataFrame(columns=['A', 'B'], index=[0]))

        # columns will align
        df = DataFrame(columns=['A', 'B'])
        df.loc[0] = Series(1, index=['B'])

        exp = DataFrame([[np.nan, 1]], columns=['A', 'B'],
                        index=[0], dtype='float64')
        tm.assert_frame_equal(df, exp)

        # list-like must conform
        df = DataFrame(columns=['A', 'B'])

        def f():
            df.loc[0] = [1, 2, 3]

        self.assertRaises(ValueError, f)

        # these are coerced to float unavoidably (as its a list-like to begin)
        df = DataFrame(columns=['A', 'B'])
        df.loc[3] = [6, 7]

        exp = DataFrame([[6, 7]], index=[3], columns=['A', 'B'],
                        dtype='float64')
        tm.assert_frame_equal(df, exp)

    def test_series_partial_set(self):
        # partial set with new index
        # Regression from GH4825
        ser = Series([0.1, 0.2], index=[1, 2])

        # loc
        expected = Series([np.nan, 0.2, np.nan], index=[3, 2, 3])
        result = ser.loc[[3, 2, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([np.nan, 0.2, np.nan, np.nan], index=[3, 2, 3, 'x'])
        result = ser.loc[[3, 2, 3, 'x']]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([0.2, 0.2, 0.1], index=[2, 2, 1])
        result = ser.loc[[2, 2, 1]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([0.2, 0.2, np.nan, 0.1], index=[2, 2, 'x', 1])
        result = ser.loc[[2, 2, 'x', 1]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        # raises as nothing in in the index
        self.assertRaises(KeyError, lambda: ser.loc[[3, 3, 3]])

        expected = Series([0.2, 0.2, np.nan], index=[2, 2, 3])
        result = ser.loc[[2, 2, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([0.3, np.nan, np.nan], index=[3, 4, 4])
        result = Series([0.1, 0.2, 0.3], index=[1, 2, 3]).loc[[3, 4, 4]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([np.nan, 0.3, 0.3], index=[5, 3, 3])
        result = Series([0.1, 0.2, 0.3, 0.4],
                        index=[1, 2, 3, 4]).loc[[5, 3, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([np.nan, 0.4, 0.4], index=[5, 4, 4])
        result = Series([0.1, 0.2, 0.3, 0.4],
                        index=[1, 2, 3, 4]).loc[[5, 4, 4]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([0.4, np.nan, np.nan], index=[7, 2, 2])
        result = Series([0.1, 0.2, 0.3, 0.4],
                        index=[4, 5, 6, 7]).loc[[7, 2, 2]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        expected = Series([0.4, np.nan, np.nan], index=[4, 5, 5])
        result = Series([0.1, 0.2, 0.3, 0.4],
                        index=[1, 2, 3, 4]).loc[[4, 5, 5]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        # iloc
        expected = Series([0.2, 0.2, 0.1, 0.1], index=[2, 2, 1, 1])
        result = ser.iloc[[1, 1, 0, 0]]
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_series_partial_set_with_name(self):
        # GH 11497

        idx = Index([1, 2], dtype='int64', name='idx')
        ser = Series([0.1, 0.2], index=idx, name='s')

        # loc
        exp_idx = Index([3, 2, 3], dtype='int64', name='idx')
        expected = Series([np.nan, 0.2, np.nan], index=exp_idx, name='s')
        result = ser.loc[[3, 2, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([3, 2, 3, 'x'], dtype='object', name='idx')
        expected = Series([np.nan, 0.2, np.nan, np.nan], index=exp_idx,
                          name='s')
        result = ser.loc[[3, 2, 3, 'x']]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([2, 2, 1], dtype='int64', name='idx')
        expected = Series([0.2, 0.2, 0.1], index=exp_idx, name='s')
        result = ser.loc[[2, 2, 1]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([2, 2, 'x', 1], dtype='object', name='idx')
        expected = Series([0.2, 0.2, np.nan, 0.1], index=exp_idx, name='s')
        result = ser.loc[[2, 2, 'x', 1]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        # raises as nothing in in the index
        self.assertRaises(KeyError, lambda: ser.loc[[3, 3, 3]])

        exp_idx = Index([2, 2, 3], dtype='int64', name='idx')
        expected = Series([0.2, 0.2, np.nan], index=exp_idx, name='s')
        result = ser.loc[[2, 2, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([3, 4, 4], dtype='int64', name='idx')
        expected = Series([0.3, np.nan, np.nan], index=exp_idx, name='s')
        idx = Index([1, 2, 3], dtype='int64', name='idx')
        result = Series([0.1, 0.2, 0.3], index=idx, name='s').loc[[3, 4, 4]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([5, 3, 3], dtype='int64', name='idx')
        expected = Series([np.nan, 0.3, 0.3], index=exp_idx, name='s')
        idx = Index([1, 2, 3, 4], dtype='int64', name='idx')
        result = Series([0.1, 0.2, 0.3, 0.4], index=idx,
                        name='s').loc[[5, 3, 3]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([5, 4, 4], dtype='int64', name='idx')
        expected = Series([np.nan, 0.4, 0.4], index=exp_idx, name='s')
        idx = Index([1, 2, 3, 4], dtype='int64', name='idx')
        result = Series([0.1, 0.2, 0.3, 0.4], index=idx,
                        name='s').loc[[5, 4, 4]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([7, 2, 2], dtype='int64', name='idx')
        expected = Series([0.4, np.nan, np.nan], index=exp_idx, name='s')
        idx = Index([4, 5, 6, 7], dtype='int64', name='idx')
        result = Series([0.1, 0.2, 0.3, 0.4], index=idx,
                        name='s').loc[[7, 2, 2]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        exp_idx = Index([4, 5, 5], dtype='int64', name='idx')
        expected = Series([0.4, np.nan, np.nan], index=exp_idx, name='s')
        idx = Index([1, 2, 3, 4], dtype='int64', name='idx')
        result = Series([0.1, 0.2, 0.3, 0.4], index=idx,
                        name='s').loc[[4, 5, 5]]
        tm.assert_series_equal(result, expected, check_index_type=True)

        # iloc
        exp_idx = Index([2, 2, 1, 1], dtype='int64', name='idx')
        expected = Series([0.2, 0.2, 0.1, 0.1], index=exp_idx, name='s')
        result = ser.iloc[[1, 1, 0, 0]]
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_partial_set_invalid(self):

        # GH 4940
        # allow only setting of 'valid' values

        orig = tm.makeTimeDataFrame()
        df = orig.copy()

        # don't allow not string inserts
        def f():
            df.loc[100.0, :] = df.ix[0]

        self.assertRaises(TypeError, f)

        def f():
            df.loc[100, :] = df.ix[0]

        self.assertRaises(TypeError, f)

        def f():
            df.ix[100.0, :] = df.ix[0]

        self.assertRaises(TypeError, f)

        def f():
            df.ix[100, :] = df.ix[0]

        self.assertRaises(ValueError, f)

        # allow object conversion here
        df = orig.copy()
        df.loc['a', :] = df.ix[0]
        exp = orig.append(pd.Series(df.ix[0], name='a'))
        tm.assert_frame_equal(df, exp)
        tm.assert_index_equal(df.index,
                              pd.Index(orig.index.tolist() + ['a']))
        self.assertEqual(df.index.dtype, 'object')

    def test_partial_set_empty_series(self):

        # GH5226

        # partially set with an empty object series
        s = Series()
        s.loc[1] = 1
        tm.assert_series_equal(s, Series([1], index=[1]))
        s.loc[3] = 3
        tm.assert_series_equal(s, Series([1, 3], index=[1, 3]))

        s = Series()
        s.loc[1] = 1.
        tm.assert_series_equal(s, Series([1.], index=[1]))
        s.loc[3] = 3.
        tm.assert_series_equal(s, Series([1., 3.], index=[1, 3]))

        s = Series()
        s.loc['foo'] = 1
        tm.assert_series_equal(s, Series([1], index=['foo']))
        s.loc['bar'] = 3
        tm.assert_series_equal(s, Series([1, 3], index=['foo', 'bar']))
        s.loc[3] = 4
        tm.assert_series_equal(s, Series([1, 3, 4], index=['foo', 'bar', 3]))

    def test_partial_set_empty_frame(self):

        # partially set with an empty object
        # frame
        df = DataFrame()

        def f():
            df.loc[1] = 1

        self.assertRaises(ValueError, f)

        def f():
            df.loc[1] = Series([1], index=['foo'])

        self.assertRaises(ValueError, f)

        def f():
            df.loc[:, 1] = 1

        self.assertRaises(ValueError, f)

        # these work as they don't really change
        # anything but the index
        # GH5632
        expected = DataFrame(columns=['foo'], index=pd.Index(
            [], dtype='int64'))

        def f():
            df = DataFrame()
            df['foo'] = Series([], dtype='object')
            return df

        tm.assert_frame_equal(f(), expected)

        def f():
            df = DataFrame()
            df['foo'] = Series(df.index)
            return df

        tm.assert_frame_equal(f(), expected)

        def f():
            df = DataFrame()
            df['foo'] = df.index
            return df

        tm.assert_frame_equal(f(), expected)

        expected = DataFrame(columns=['foo'],
                             index=pd.Index([], dtype='int64'))
        expected['foo'] = expected['foo'].astype('float64')

        def f():
            df = DataFrame()
            df['foo'] = []
            return df

        tm.assert_frame_equal(f(), expected)

        def f():
            df = DataFrame()
            df['foo'] = Series(range(len(df)))
            return df

        tm.assert_frame_equal(f(), expected)

        def f():
            df = DataFrame()
            tm.assert_index_equal(df.index, pd.Index([], dtype='object'))
            df['foo'] = range(len(df))
            return df

        expected = DataFrame(columns=['foo'],
                             index=pd.Index([], dtype='int64'))
        expected['foo'] = expected['foo'].astype('float64')
        tm.assert_frame_equal(f(), expected)

        df = DataFrame()
        tm.assert_index_equal(df.columns, pd.Index([], dtype=object))
        df2 = DataFrame()
        df2[1] = Series([1], index=['foo'])
        df.loc[:, 1] = Series([1], index=['foo'])
        tm.assert_frame_equal(df, DataFrame([[1]], index=['foo'], columns=[1]))
        tm.assert_frame_equal(df, df2)

        # no index to start
        expected = DataFrame({0: Series(1, index=range(4))},
                             columns=['A', 'B', 0])

        df = DataFrame(columns=['A', 'B'])
        df[0] = Series(1, index=range(4))
        df.dtypes
        str(df)
        tm.assert_frame_equal(df, expected)

        df = DataFrame(columns=['A', 'B'])
        df.loc[:, 0] = Series(1, index=range(4))
        df.dtypes
        str(df)
        tm.assert_frame_equal(df, expected)

    def test_partial_set_empty_frame_row(self):
        # GH5720, GH5744
        # don't create rows when empty
        expected = DataFrame(columns=['A', 'B', 'New'],
                             index=pd.Index([], dtype='int64'))
        expected['A'] = expected['A'].astype('int64')
        expected['B'] = expected['B'].astype('float64')
        expected['New'] = expected['New'].astype('float64')

        df = DataFrame({"A": [1, 2, 3], "B": [1.2, 4.2, 5.2]})
        y = df[df.A > 5]
        y['New'] = np.nan
        tm.assert_frame_equal(y, expected)
        # tm.assert_frame_equal(y,expected)

        expected = DataFrame(columns=['a', 'b', 'c c', 'd'])
        expected['d'] = expected['d'].astype('int64')
        df = DataFrame(columns=['a', 'b', 'c c'])
        df['d'] = 3
        tm.assert_frame_equal(df, expected)
        tm.assert_series_equal(df['c c'], Series(name='c c', dtype=object))

        # reindex columns is ok
        df = DataFrame({"A": [1, 2, 3], "B": [1.2, 4.2, 5.2]})
        y = df[df.A > 5]
        result = y.reindex(columns=['A', 'B', 'C'])
        expected = DataFrame(columns=['A', 'B', 'C'],
                             index=pd.Index([], dtype='int64'))
        expected['A'] = expected['A'].astype('int64')
        expected['B'] = expected['B'].astype('float64')
        expected['C'] = expected['C'].astype('float64')
        tm.assert_frame_equal(result, expected)

    def test_partial_set_empty_frame_set_series(self):
        # GH 5756
        # setting with empty Series
        df = DataFrame(Series())
        tm.assert_frame_equal(df, DataFrame({0: Series()}))

        df = DataFrame(Series(name='foo'))
        tm.assert_frame_equal(df, DataFrame({'foo': Series()}))

    def test_partial_set_empty_frame_empty_copy_assignment(self):
        # GH 5932
        # copy on empty with assignment fails
        df = DataFrame(index=[0])
        df = df.copy()
        df['a'] = 0
        expected = DataFrame(0, index=[0], columns=['a'])
        tm.assert_frame_equal(df, expected)

    def test_partial_set_empty_frame_empty_consistencies(self):
        # GH 6171
        # consistency on empty frames
        df = DataFrame(columns=['x', 'y'])
        df['x'] = [1, 2]
        expected = DataFrame(dict(x=[1, 2], y=[np.nan, np.nan]))
        tm.assert_frame_equal(df, expected, check_dtype=False)

        df = DataFrame(columns=['x', 'y'])
        df['x'] = ['1', '2']
        expected = DataFrame(
            dict(x=['1', '2'], y=[np.nan, np.nan]), dtype=object)
        tm.assert_frame_equal(df, expected)

        df = DataFrame(columns=['x', 'y'])
        df.loc[0, 'x'] = 1
        expected = DataFrame(dict(x=[1], y=[np.nan]))
        tm.assert_frame_equal(df, expected, check_dtype=False)

    def test_cache_updating(self):
        # GH 4939, make sure to update the cache on setitem

        df = tm.makeDataFrame()
        df['A']  # cache series
        df.ix["Hello Friend"] = df.ix[0]
        self.assertIn("Hello Friend", df['A'].index)
        self.assertIn("Hello Friend", df['B'].index)

        panel = tm.makePanel()
        panel.ix[0]  # get first item into cache
        panel.ix[:, :, 'A+1'] = panel.ix[:, :, 'A'] + 1
        self.assertIn("A+1", panel.ix[0].columns)
        self.assertIn("A+1", panel.ix[1].columns)

        # 5216
        # make sure that we don't try to set a dead cache
        a = np.random.rand(10, 3)
        df = DataFrame(a, columns=['x', 'y', 'z'])
        tuples = [(i, j) for i in range(5) for j in range(2)]
        index = MultiIndex.from_tuples(tuples)
        df.index = index

        # setting via chained assignment
        # but actually works, since everything is a view
        df.loc[0]['z'].iloc[0] = 1.
        result = df.loc[(0, 0), 'z']
        self.assertEqual(result, 1)

        # correct setting
        df.loc[(0, 0), 'z'] = 2
        result = df.loc[(0, 0), 'z']
        self.assertEqual(result, 2)

        # 10264
        df = DataFrame(np.zeros((5, 5), dtype='int64'), columns=[
                       'a', 'b', 'c', 'd', 'e'], index=range(5))
        df['f'] = 0
        df.f.values[3] = 1

        # TODO(wesm): unused?
        # y = df.iloc[np.arange(2, len(df))]

        df.f.values[3] = 2
        expected = DataFrame(np.zeros((5, 6), dtype='int64'), columns=[
                             'a', 'b', 'c', 'd', 'e', 'f'], index=range(5))
        expected.at[3, 'f'] = 2
        tm.assert_frame_equal(df, expected)
        expected = Series([0, 0, 0, 2, 0], name='f')
        tm.assert_series_equal(df.f, expected)

    def test_set_ix_out_of_bounds_axis_0(self):
        df = pd.DataFrame(
            randn(2, 5), index=["row%s" % i for i in range(2)],
            columns=["col%s" % i for i in range(5)])
        self.assertRaises(ValueError, df.ix.__setitem__, (2, 0), 100)

    def test_set_ix_out_of_bounds_axis_1(self):
        df = pd.DataFrame(
            randn(5, 2), index=["row%s" % i for i in range(5)],
            columns=["col%s" % i for i in range(2)])
        self.assertRaises(ValueError, df.ix.__setitem__, (0, 2), 100)

    def test_iloc_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        tm.assert_frame_equal(df.iloc[:, []], df.iloc[:, :0],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.iloc[[], :], df.iloc[:0, :],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.iloc[[]], df.iloc[:0, :],
                              check_index_type=True,
                              check_column_type=True)

    def test_loc_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        tm.assert_frame_equal(df.loc[:, []], df.iloc[:, :0],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.loc[[], :], df.iloc[:0, :],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.loc[[]], df.iloc[:0, :],
                              check_index_type=True,
                              check_column_type=True)

    def test_ix_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        tm.assert_frame_equal(df.ix[:, []], df.iloc[:, :0],
                              check_index_type=True,
                              check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.ix[[], :], df.iloc[:0, :],
                              check_index_type=True,
                              check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.ix[[]], df.iloc[:0, :],
                              check_index_type=True,
                              check_column_type=True)

    def test_index_type_coercion(self):

        # GH 11836
        # if we have an index type and set it with something that looks
        # to numpy like the same, but is actually, not
        # (e.g. setting with a float or string '0')
        # then we need to coerce to object

        # integer indexes
        for s in [Series(range(5)),
                  Series(range(5), index=range(1, 6))]:

            self.assertTrue(s.index.is_integer())

            for indexer in [lambda x: x.ix,
                            lambda x: x.loc,
                            lambda x: x]:
                s2 = s.copy()
                indexer(s2)[0.1] = 0
                self.assertTrue(s2.index.is_floating())
                self.assertTrue(indexer(s2)[0.1] == 0)

                s2 = s.copy()
                indexer(s2)[0.0] = 0
                exp = s.index
                if 0 not in s:
                    exp = Index(s.index.tolist() + [0])
                tm.assert_index_equal(s2.index, exp)

                s2 = s.copy()
                indexer(s2)['0'] = 0
                self.assertTrue(s2.index.is_object())

        for s in [Series(range(5), index=np.arange(5.))]:

            self.assertTrue(s.index.is_floating())

            for idxr in [lambda x: x.ix,
                         lambda x: x.loc,
                         lambda x: x]:

                s2 = s.copy()
                idxr(s2)[0.1] = 0
                self.assertTrue(s2.index.is_floating())
                self.assertTrue(idxr(s2)[0.1] == 0)

                s2 = s.copy()
                idxr(s2)[0.0] = 0
                tm.assert_index_equal(s2.index, s.index)

                s2 = s.copy()
                idxr(s2)['0'] = 0
                self.assertTrue(s2.index.is_object())

    def test_float_index_to_mixed(self):
        df = DataFrame({0.0: np.random.rand(10), 1.0: np.random.rand(10)})
        df['a'] = 10
        tm.assert_frame_equal(DataFrame({0.0: df[0.0],
                                         1.0: df[1.0],
                                         'a': [10] * 10}),
                              df)

    def test_duplicate_ix_returns_series(self):
        df = DataFrame(np.random.randn(3, 3), index=[0.1, 0.2, 0.2],
                       columns=list('abc'))
        r = df.ix[0.2, 'a']
        e = df.loc[0.2, 'a']
        tm.assert_series_equal(r, e)

    def test_float_index_non_scalar_assignment(self):
        df = DataFrame({'a': [1, 2, 3], 'b': [3, 4, 5]}, index=[1., 2., 3.])
        df.loc[df.index[:2]] = 1
        expected = DataFrame({'a': [1, 1, 3], 'b': [1, 1, 5]}, index=df.index)
        tm.assert_frame_equal(expected, df)

        df = DataFrame({'a': [1, 2, 3], 'b': [3, 4, 5]}, index=[1., 2., 3.])
        df2 = df.copy()
        df.loc[df.index] = df.loc[df.index]
        tm.assert_frame_equal(df, df2)

    def test_float_index_at_iat(self):
        s = pd.Series([1, 2, 3], index=[0.1, 0.2, 0.3])
        for el, item in s.iteritems():
            self.assertEqual(s.at[el], item)
        for i in range(len(s)):
            self.assertEqual(s.iat[i], i + 1)

    def test_rhs_alignment(self):
        # GH8258, tests that both rows & columns are aligned to what is
        # assigned to. covers both uniform data-type & multi-type cases
        def run_tests(df, rhs, right):
            # label, index, slice
            r, i, s = list('bcd'), [1, 2, 3], slice(1, 4)
            c, j, l = ['joe', 'jolie'], [1, 2], slice(1, 3)

            left = df.copy()
            left.loc[r, c] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.iloc[i, j] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.ix[s, l] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.ix[i, j] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.ix[r, c] = rhs
            tm.assert_frame_equal(left, right)

        xs = np.arange(20).reshape(5, 4)
        cols = ['jim', 'joe', 'jolie', 'joline']
        df = pd.DataFrame(xs, columns=cols, index=list('abcde'))

        # right hand side; permute the indices and multiplpy by -2
        rhs = -2 * df.iloc[3:0:-1, 2:0:-1]

        # expected `right` result; just multiply by -2
        right = df.copy()
        right.iloc[1:4, 1:3] *= -2

        # run tests with uniform dtypes
        run_tests(df, rhs, right)

        # make frames multi-type & re-run tests
        for frame in [df, rhs, right]:
            frame['joe'] = frame['joe'].astype('float64')
            frame['jolie'] = frame['jolie'].map('@{0}'.format)

        run_tests(df, rhs, right)

    def test_str_label_slicing_with_negative_step(self):
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            tm.assert_series_equal(s.loc[l_slc], s.iloc[i_slc])

            if not idx.is_integer:
                # For integer indices, ix and plain getitem are position-based.
                tm.assert_series_equal(s[l_slc], s.iloc[i_slc])
                tm.assert_series_equal(s.ix[l_slc], s.iloc[i_slc])

        for idx in [_mklbl('A', 20), np.arange(20) + 100,
                    np.linspace(100, 150, 20)]:
            idx = Index(idx)
            s = Series(np.arange(20), index=idx)
            assert_slices_equivalent(SLC[idx[9]::-1], SLC[9::-1])
            assert_slices_equivalent(SLC[:idx[9]:-1], SLC[:8:-1])
            assert_slices_equivalent(SLC[idx[13]:idx[9]:-1], SLC[13:8:-1])
            assert_slices_equivalent(SLC[idx[9]:idx[13]:-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        s = Series(np.arange(20), index=_mklbl('A', 20))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s.ix[::0])

    def test_indexing_assignment_dict_already_exists(self):
        df = pd.DataFrame({'x': [1, 2, 6],
                           'y': [2, 2, 8],
                           'z': [-5, 0, 5]}).set_index('z')
        expected = df.copy()
        rhs = dict(x=9, y=99)
        df.loc[5] = rhs
        expected.loc[5] = [9, 99]
        tm.assert_frame_equal(df, expected)

    def test_indexing_dtypes_on_empty(self):
        # Check that .iloc and .ix return correct dtypes GH9983
        df = DataFrame({'a': [1, 2, 3], 'b': ['b', 'b2', 'b3']})
        df2 = df.ix[[], :]

        self.assertEqual(df2.loc[:, 'a'].dtype, np.int64)
        tm.assert_series_equal(df2.loc[:, 'a'], df2.iloc[:, 0])
        tm.assert_series_equal(df2.loc[:, 'a'], df2.ix[:, 0])

    def test_range_in_series_indexing(self):
        # range can cause an indexing error
        # GH 11652
        for x in [5, 999999, 1000000]:
            s = pd.Series(index=range(x))
            s.loc[range(1)] = 42
            tm.assert_series_equal(s.loc[range(1)], Series(42.0, index=[0]))

            s.loc[range(2)] = 43
            tm.assert_series_equal(s.loc[range(2)], Series(43.0, index=[0, 1]))

    def test_non_reducing_slice(self):
        df = pd.DataFrame([[0, 1], [2, 3]])

        slices = [
            # pd.IndexSlice[:, :],
            pd.IndexSlice[:, 1],
            pd.IndexSlice[1, :],
            pd.IndexSlice[[1], [1]],
            pd.IndexSlice[1, [1]],
            pd.IndexSlice[[1], 1],
            pd.IndexSlice[1],
            pd.IndexSlice[1, 1],
            slice(None, None, None),
            [0, 1],
            np.array([0, 1]),
            pd.Series([0, 1])
        ]
        for slice_ in slices:
            tslice_ = _non_reducing_slice(slice_)
            self.assertTrue(isinstance(df.loc[tslice_], DataFrame))

    def test_list_slice(self):
        # like dataframe getitem
        slices = [['A'], pd.Series(['A']), np.array(['A'])]
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4]}, index=['A', 'B'])
        expected = pd.IndexSlice[:, ['A']]
        for subset in slices:
            result = _non_reducing_slice(subset)
            tm.assert_frame_equal(df.loc[result], df.loc[expected])

    def test_maybe_numeric_slice(self):
        df = pd.DataFrame({'A': [1, 2], 'B': ['c', 'd'], 'C': [True, False]})
        result = _maybe_numeric_slice(df, slice_=None)
        expected = pd.IndexSlice[:, ['A']]
        self.assertEqual(result, expected)

        result = _maybe_numeric_slice(df, None, include_bool=True)
        expected = pd.IndexSlice[:, ['A', 'C']]
        result = _maybe_numeric_slice(df, [1])
        expected = [1]
        self.assertEqual(result, expected)


class TestSeriesNoneCoercion(tm.TestCase):
    EXPECTED_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),

        # For datetime series, we should coerce to NaT.
        ([datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
         [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)]),

        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[0] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_loc_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[0] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_setitem_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[start_series == start_series[0]] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_loc_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[start_series == start_series[0]] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)


class TestDataframeNoneCoercion(tm.TestCase):
    EXPECTED_SINGLE_ROW_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),

        # For datetime series, we should coerce to NaT.
        ([datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
         [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)]),

        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_loc(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe.loc[0, ['foo']] = None

            expected_dataframe = DataFrame({'foo': expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_coercion_with_setitem_and_dataframe(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe[start_dataframe['foo'] == start_dataframe['foo'][
                0]] = None

            expected_dataframe = DataFrame({'foo': expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_none_coercion_loc_and_dataframe(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe.loc[start_dataframe['foo'] == start_dataframe[
                'foo'][0]] = None

            expected_dataframe = DataFrame({'foo': expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_none_coercion_mixed_dtypes(self):
        start_dataframe = DataFrame({
            'a': [1, 2, 3],
            'b': [1.0, 2.0, 3.0],
            'c': [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1,
                                                                       3)],
            'd': ['a', 'b', 'c']
        })
        start_dataframe.iloc[0] = None

        exp = DataFrame({'a': [np.nan, 2, 3],
                         'b': [np.nan, 2.0, 3.0],
                         'c': [NaT, datetime(2000, 1, 2),
                               datetime(2000, 1, 3)],
                         'd': [None, 'b', 'c']})
        tm.assert_frame_equal(start_dataframe, exp)
