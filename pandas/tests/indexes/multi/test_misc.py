# -*- coding: utf-8 -*-

import re
import warnings

from datetime import timedelta
from itertools import product

import pytest

import numpy as np

import pandas as pd

from pandas import (CategoricalIndex, DataFrame, Index, MultiIndex,
                    compat, date_range, period_range)
from pandas.compat import PY3, long, lrange, lzip, range, u, PYPY
from pandas.errors import PerformanceWarning, UnsortedIndexError
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.indexes.base import InvalidIndexError
from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike
from pandas._libs.tslib import Timestamp

import pandas.util.testing as tm

from pandas.util.testing import assert_almost_equal, assert_copy

from .common import Base


class TestMultiIndex(Base):
    _holder = MultiIndex
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setup_method(self, method):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.indices = dict(index=MultiIndex(levels=[major_axis, minor_axis],
                                             labels=[major_labels, minor_labels
                                                     ], names=self.index_names,
                                             verify_integrity=False))
        self.setup_indices()

    def create_index(self):
        return self.index

    def test_boolean_context_compat2(self):

        # boolean context compat
        # GH7897
        i1 = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        i2 = MultiIndex.from_tuples([('A', 1), ('A', 3)])
        common = i1.intersection(i2)

        def f():
            if common:
                pass

        tm.assert_raises_regex(ValueError, 'The truth value of a', f)

    def test_labels_dtypes(self):

        # GH 8456
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        assert i.labels[0].dtype == 'int8'
        assert i.labels[1].dtype == 'int8'

        i = MultiIndex.from_product([['a'], range(40)])
        assert i.labels[1].dtype == 'int8'
        i = MultiIndex.from_product([['a'], range(400)])
        assert i.labels[1].dtype == 'int16'
        i = MultiIndex.from_product([['a'], range(40000)])
        assert i.labels[1].dtype == 'int32'

        i = pd.MultiIndex.from_product([['a'], range(1000)])
        assert (i.labels[0] >= 0).all()
        assert (i.labels[1] >= 0).all()

    def test_where(self):
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])

        def f():
            i.where(True)

        pytest.raises(NotImplementedError, f)

    def test_where_array_like(self):
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        klasses = [list, tuple, np.array, pd.Series]
        cond = [False, True]

        for klass in klasses:
            f = lambda: i.where(klass(cond))
            pytest.raises(NotImplementedError, f)

    def test_repeat(self):
        reps = 2
        numbers = [1, 2, 3]
        names = np.array(['foo', 'bar'])

        m = MultiIndex.from_product([
            numbers, names], names=names)
        expected = MultiIndex.from_product([
            numbers, names.repeat(reps)], names=names)
        tm.assert_index_equal(m.repeat(reps), expected)

        with tm.assert_produces_warning(FutureWarning):
            result = m.repeat(n=reps)
            tm.assert_index_equal(result, expected)

    def test_numpy_repeat(self):
        reps = 2
        numbers = [1, 2, 3]
        names = np.array(['foo', 'bar'])

        m = MultiIndex.from_product([
            numbers, names], names=names)
        expected = MultiIndex.from_product([
            numbers, names.repeat(reps)], names=names)
        tm.assert_index_equal(np.repeat(m, reps), expected)

        msg = "the 'axis' parameter is not supported"
        tm.assert_raises_regex(
            ValueError, msg, np.repeat, m, reps, axis=1)

    def test_metadata_immutable(self):
        levels, labels = self.index.levels, self.index.labels
        # shouldn't be able to set at either the top level or base level
        mutable_regex = re.compile('does not support mutable operations')
        with tm.assert_raises_regex(TypeError, mutable_regex):
            levels[0] = levels[0]
        with tm.assert_raises_regex(TypeError, mutable_regex):
            levels[0][0] = levels[0][0]
        # ditto for labels
        with tm.assert_raises_regex(TypeError, mutable_regex):
            labels[0] = labels[0]
        with tm.assert_raises_regex(TypeError, mutable_regex):
            labels[0][0] = labels[0][0]
        # and for names
        names = self.index.names
        with tm.assert_raises_regex(TypeError, mutable_regex):
            names[0] = names[0]

    def test_inplace_mutation_resets_values(self):
        levels = [['a', 'b', 'c'], [4]]
        levels2 = [[1, 2, 3], ['a']]
        labels = [[0, 1, 0, 2, 2, 0], [0, 0, 0, 0, 0, 0]]

        mi1 = MultiIndex(levels=levels, labels=labels)
        mi2 = MultiIndex(levels=levels2, labels=labels)
        vals = mi1.values.copy()
        vals2 = mi2.values.copy()

        assert mi1._tuples is not None

        # Make sure level setting works
        new_vals = mi1.set_levels(levels2).values
        tm.assert_almost_equal(vals2, new_vals)

        # Non-inplace doesn't kill _tuples [implementation detail]
        tm.assert_almost_equal(mi1._tuples, vals)

        # ...and values is still same too
        tm.assert_almost_equal(mi1.values, vals)

        # Inplace should kill _tuples
        mi1.set_levels(levels2, inplace=True)
        tm.assert_almost_equal(mi1.values, vals2)

        # Make sure label setting works too
        labels2 = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        exp_values = np.empty((6,), dtype=object)
        exp_values[:] = [(long(1), 'a')] * 6

        # Must be 1d array of tuples
        assert exp_values.shape == (6,)
        new_values = mi2.set_labels(labels2).values

        # Not inplace shouldn't change
        tm.assert_almost_equal(mi2._tuples, vals2)

        # Should have correct values
        tm.assert_almost_equal(exp_values, new_values)

        # ...and again setting inplace should kill _tuples, etc
        mi2.set_labels(labels2, inplace=True)
        tm.assert_almost_equal(mi2.values, new_values)

    def test_astype(self):
        expected = self.index.copy()
        actual = self.index.astype('O')
        assert_copy(actual.levels, expected.levels)
        assert_copy(actual.labels, expected.labels)
        self.check_level_names(actual, expected.names)

        with tm.assert_raises_regex(TypeError, "^Setting.*dtype.*object"):
            self.index.astype(np.dtype(int))

    @pytest.mark.parametrize('ordered', [True, False])
    def test_astype_category(self, ordered):
        # GH 18630
        msg = '> 1 ndim Categorical are not supported at this time'
        with tm.assert_raises_regex(NotImplementedError, msg):
            self.index.astype(CategoricalDtype(ordered=ordered))

        if ordered is False:
            # dtype='category' defaults to ordered=False, so only test once
            with tm.assert_raises_regex(NotImplementedError, msg):
                self.index.astype('category')

    def assert_multiindex_copied(self, copy, original):
        # Levels should be (at least, shallow copied)
        tm.assert_copy(copy.levels, original.levels)
        tm.assert_almost_equal(copy.labels, original.labels)

        # Labels doesn't matter which way copied
        tm.assert_almost_equal(copy.labels, original.labels)
        assert copy.labels is not original.labels

        # Names doesn't matter which way copied
        assert copy.names == original.names
        assert copy.names is not original.names

        # Sort order should be copied
        assert copy.sortorder == original.sortorder

    def test_copy(self):
        i_copy = self.index.copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_shallow_copy(self):
        i_copy = self.index._shallow_copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_view(self):
        i_view = self.index.view()

        self.assert_multiindex_copied(i_view, self.index)

    def test_values_boxed(self):
        tuples = [(1, pd.Timestamp('2000-01-01')), (2, pd.NaT),
                  (3, pd.Timestamp('2000-01-03')),
                  (1, pd.Timestamp('2000-01-04')),
                  (2, pd.Timestamp('2000-01-02')),
                  (3, pd.Timestamp('2000-01-03'))]
        result = pd.MultiIndex.from_tuples(tuples)
        expected = construct_1d_object_array_from_listlike(tuples)
        tm.assert_numpy_array_equal(result.values, expected)
        # Check that code branches for boxed values produce identical results
        tm.assert_numpy_array_equal(result.values[:4], result[:4].values)

    def test_values_multiindex_datetimeindex(self):
        # Test to ensure we hit the boxing / nobox part of MI.values
        ints = np.arange(10 ** 18, 10 ** 18 + 5)
        naive = pd.DatetimeIndex(ints)
        aware = pd.DatetimeIndex(ints, tz='US/Central')

        idx = pd.MultiIndex.from_arrays([naive, aware])
        result = idx.values

        outer = pd.DatetimeIndex([x[0] for x in result])
        tm.assert_index_equal(outer, naive)

        inner = pd.DatetimeIndex([x[1] for x in result])
        tm.assert_index_equal(inner, aware)

        # n_lev > n_lab
        result = idx[:2].values

        outer = pd.DatetimeIndex([x[0] for x in result])
        tm.assert_index_equal(outer, naive[:2])

        inner = pd.DatetimeIndex([x[1] for x in result])
        tm.assert_index_equal(inner, aware[:2])

    def test_values_multiindex_periodindex(self):
        # Test to ensure we hit the boxing / nobox part of MI.values
        ints = np.arange(2007, 2012)
        pidx = pd.PeriodIndex(ints, freq='D')

        idx = pd.MultiIndex.from_arrays([ints, pidx])
        result = idx.values

        outer = pd.Int64Index([x[0] for x in result])
        tm.assert_index_equal(outer, pd.Int64Index(ints))

        inner = pd.PeriodIndex([x[1] for x in result])
        tm.assert_index_equal(inner, pidx)

        # n_lev > n_lab
        result = idx[:2].values

        outer = pd.Int64Index([x[0] for x in result])
        tm.assert_index_equal(outer, pd.Int64Index(ints[:2]))

        inner = pd.PeriodIndex([x[1] for x in result])
        tm.assert_index_equal(inner, pidx[:2])

    def test_append(self):
        result = self.index[:3].append(self.index[3:])
        assert result.equals(self.index)

        foos = [self.index[:1], self.index[1:3], self.index[3:]]
        result = foos[0].append(foos[1:])
        assert result.equals(self.index)

        # empty
        result = self.index.append([])
        assert result.equals(self.index)

    def test_append_mixed_dtypes(self):
        # GH 13660
        dti = date_range('2011-01-01', freq='M', periods=3, )
        dti_tz = date_range('2011-01-01', freq='M', periods=3, tz='US/Eastern')
        pi = period_range('2011-01', freq='M', periods=3)

        mi = MultiIndex.from_arrays([[1, 2, 3],
                                     [1.1, np.nan, 3.3],
                                     ['a', 'b', 'c'],
                                     dti, dti_tz, pi])
        assert mi.nlevels == 6

        res = mi.append(mi)
        exp = MultiIndex.from_arrays([[1, 2, 3, 1, 2, 3],
                                      [1.1, np.nan, 3.3, 1.1, np.nan, 3.3],
                                      ['a', 'b', 'c', 'a', 'b', 'c'],
                                      dti.append(dti),
                                      dti_tz.append(dti_tz),
                                      pi.append(pi)])
        tm.assert_index_equal(res, exp)

        other = MultiIndex.from_arrays([['x', 'y', 'z'], ['x', 'y', 'z'],
                                        ['x', 'y', 'z'], ['x', 'y', 'z'],
                                        ['x', 'y', 'z'], ['x', 'y', 'z']])

        res = mi.append(other)
        exp = MultiIndex.from_arrays([[1, 2, 3, 'x', 'y', 'z'],
                                      [1.1, np.nan, 3.3, 'x', 'y', 'z'],
                                      ['a', 'b', 'c', 'x', 'y', 'z'],
                                      dti.append(pd.Index(['x', 'y', 'z'])),
                                      dti_tz.append(pd.Index(['x', 'y', 'z'])),
                                      pi.append(pd.Index(['x', 'y', 'z']))])
        tm.assert_index_equal(res, exp)

    def test_reorder_levels(self):
        # this blows up
        tm.assert_raises_regex(IndexError, '^Too many levels',
                               self.index.reorder_levels, [2, 1, 0])

    def test_nlevels(self):
        assert self.index.nlevels == 2

    def test_iter(self):
        result = list(self.index)
        expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                    ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
        assert result == expected

    def test_legacy_pickle(self):
        if PY3:
            pytest.skip("testing for legacy pickles not "
                        "support on py3")

        path = tm.get_data_path('multiindex_v1.pickle')
        obj = pd.read_pickle(path)

        obj2 = MultiIndex.from_tuples(obj.values)
        assert obj.equals(obj2)

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj), dtype=np.intp)
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_legacy_v2_unpickle(self):

        # 0.7.3 -> 0.8.0 format manage
        path = tm.get_data_path('mindex_073.pickle')
        obj = pd.read_pickle(path)

        obj2 = MultiIndex.from_tuples(obj.values)
        assert obj.equals(obj2)

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj), dtype=np.intp)
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = MultiIndex.from_product(
            [[1, 2], ['a', 'b'], date_range('20130101', periods=3,
                                            tz='US/Eastern')
             ], names=['one', 'two', 'three'])
        unpickled = tm.round_trip_pickle(index)
        assert index.equal_levels(unpickled)

    def test_from_tuples_index_values(self):
        result = MultiIndex.from_tuples(self.index)
        assert (result.values == self.index.values).all()

    def test_is_all_dates(self):
        assert not self.index.is_all_dates

    def test_is_numeric(self):
        # MultiIndex is never numeric
        assert not self.index.is_numeric()

    def test_consistency(self):
        # need to construct an overflow
        major_axis = lrange(70000)
        minor_axis = lrange(10)

        major_labels = np.arange(70000)
        minor_labels = np.repeat(lrange(10), 7000)

        # the fact that is works means it's consistent
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        # inconsistent
        major_labels = np.array([0, 0, 1, 1, 1, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1])
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        assert not index.is_unique

    def test_truncate(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        result = index.truncate(before=1)
        assert 'foo' not in result.levels[0]
        assert 1 in result.levels[0]

        result = index.truncate(after=1)
        assert 2 not in result.levels[0]
        assert 1 in result.levels[0]

        result = index.truncate(before=1, after=2)
        assert len(result.levels[0]) == 2

        # after < before
        pytest.raises(ValueError, index.truncate, 3, 1)

    def test_hash_collisions(self):
        # non-smoke test that we don't get hash collisions

        index = MultiIndex.from_product([np.arange(1000), np.arange(1000)],
                                        names=['one', 'two'])
        result = index.get_indexer(index.values)
        tm.assert_numpy_array_equal(result, np.arange(
            len(index), dtype='intp'))

        for i in [0, 1, len(index) - 2, len(index) - 1]:
            result = index.get_loc(index[i])
            assert result == i

    def test_to_frame(self):
        tuples = [(1, 'one'), (1, 'two'), (2, 'one'), (2, 'two')]

        index = MultiIndex.from_tuples(tuples)
        result = index.to_frame(index=False)
        expected = DataFrame(tuples)
        tm.assert_frame_equal(result, expected)

        result = index.to_frame()
        expected.index = index
        tm.assert_frame_equal(result, expected)

        tuples = [(1, 'one'), (1, 'two'), (2, 'one'), (2, 'two')]
        index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
        result = index.to_frame(index=False)
        expected = DataFrame(tuples)
        expected.columns = ['first', 'second']
        tm.assert_frame_equal(result, expected)

        result = index.to_frame()
        expected.index = index
        tm.assert_frame_equal(result, expected)

        index = MultiIndex.from_product([range(5),
                                         pd.date_range('20130101', periods=3)])
        result = index.to_frame(index=False)
        expected = DataFrame(
            {0: np.repeat(np.arange(5, dtype='int64'), 3),
             1: np.tile(pd.date_range('20130101', periods=3), 5)})
        tm.assert_frame_equal(result, expected)

        index = MultiIndex.from_product([range(5),
                                         pd.date_range('20130101', periods=3)])
        result = index.to_frame()
        expected.index = index
        tm.assert_frame_equal(result, expected)

    def test_to_hierarchical(self):
        index = MultiIndex.from_tuples([(1, 'one'), (1, 'two'), (2, 'one'), (
            2, 'two')])
        result = index.to_hierarchical(3)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                      [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        assert result.names == index.names

        # K > 1
        result = index.to_hierarchical(3, 2)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                                      [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        assert result.names == index.names

        # non-sorted
        index = MultiIndex.from_tuples([(2, 'c'), (1, 'b'),
                                        (2, 'a'), (2, 'b')],
                                       names=['N1', 'N2'])

        result = index.to_hierarchical(2)
        expected = MultiIndex.from_tuples([(2, 'c'), (2, 'c'), (1, 'b'),
                                           (1, 'b'),
                                           (2, 'a'), (2, 'a'),
                                           (2, 'b'), (2, 'b')],
                                          names=['N1', 'N2'])
        tm.assert_index_equal(result, expected)
        assert result.names == index.names

    def test_bounds(self):
        self.index._bounds

    def test_equals_multi(self):
        assert self.index.equals(self.index)
        assert not self.index.equals(self.index.values)
        assert self.index.equals(Index(self.index.values))

        assert self.index.equal_levels(self.index)
        assert not self.index.equals(self.index[:-1])
        assert not self.index.equals(self.index[-1])

        # different number of levels
        index = MultiIndex(levels=[Index(lrange(4)), Index(lrange(4)), Index(
            lrange(4))], labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]), np.array(
                [0, 1, 0, 0, 0, 1, 0, 1]), np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        index2 = MultiIndex(levels=index.levels[:-1], labels=index.labels[:-1])
        assert not index.equals(index2)
        assert not index.equal_levels(index2)

        # levels are different
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        assert not self.index.equals(index)
        assert not self.index.equal_levels(index)

        # some of the labels are different
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        assert not self.index.equals(index)

    def test_equals_missing_values(self):
        # make sure take is not using -1
        i = pd.MultiIndex.from_tuples([(0, pd.NaT),
                                       (0, pd.Timestamp('20130101'))])
        result = i[0:1].equals(i[0])
        assert not result
        result = i[1:2].equals(i[1])
        assert not result

    def test_identical(self):
        mi = self.index.copy()
        mi2 = self.index.copy()
        assert mi.identical(mi2)

        mi = mi.set_names(['new1', 'new2'])
        assert mi.equals(mi2)
        assert not mi.identical(mi2)

        mi2 = mi2.set_names(['new1', 'new2'])
        assert mi.identical(mi2)

        mi3 = Index(mi.tolist(), names=mi.names)
        mi4 = Index(mi.tolist(), names=mi.names, tupleize_cols=False)
        assert mi.identical(mi3)
        assert not mi.identical(mi4)
        assert mi.equals(mi4)

    def test_is_(self):
        mi = MultiIndex.from_tuples(lzip(range(10), range(10)))
        assert mi.is_(mi)
        assert mi.is_(mi.view())
        assert mi.is_(mi.view().view().view().view())
        mi2 = mi.view()
        # names are metadata, they don't change id
        mi2.names = ["A", "B"]
        assert mi2.is_(mi)
        assert mi.is_(mi2)

        assert mi.is_(mi.set_names(["C", "D"]))
        mi2 = mi.view()
        mi2.set_names(["E", "F"], inplace=True)
        assert mi.is_(mi2)
        # levels are inherent properties, they change identity
        mi3 = mi2.set_levels([lrange(10), lrange(10)])
        assert not mi3.is_(mi2)
        # shouldn't change
        assert mi2.is_(mi)
        mi4 = mi3.view()

        # GH 17464 - Remove duplicate MultiIndex levels
        mi4.set_levels([lrange(10), lrange(10)], inplace=True)
        assert not mi4.is_(mi3)
        mi5 = mi.view()
        mi5.set_levels(mi5.levels, inplace=True)
        assert not mi5.is_(mi)

    def test_union(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_union = piece1 | piece2

        tups = sorted(self.index.values)
        expected = MultiIndex.from_tuples(tups)

        assert the_union.equals(expected)

        # corner case, pass self or empty thing:
        the_union = self.index.union(self.index)
        assert the_union is self.index

        the_union = self.index.union(self.index[:0])
        assert the_union is self.index

        # won't work in python 3
        # tuples = self.index.values
        # result = self.index[:4] | tuples[4:]
        # assert result.equals(tuples)

        # not valid for python 3
        # def test_union_with_regular_index(self):
        #     other = Index(['A', 'B', 'C'])

        #     result = other.union(self.index)
        #     assert ('foo', 'one') in result
        #     assert 'B' in result

        #     result2 = self.index.union(other)
        #     assert result.equals(result2)

    def test_intersection(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_int = piece1 & piece2
        tups = sorted(self.index[3:5].values)
        expected = MultiIndex.from_tuples(tups)
        assert the_int.equals(expected)

        # corner case, pass self
        the_int = self.index.intersection(self.index)
        assert the_int is self.index

        # empty intersection: disjoint
        empty = self.index[:2] & self.index[2:]
        expected = self.index[:0]
        assert empty.equals(expected)

        # can't do in python 3
        # tuples = self.index.values
        # result = self.index & tuples
        # assert result.equals(tuples)

    def test_sub(self):

        first = self.index

        # - now raises (previously was set op difference)
        with pytest.raises(TypeError):
            first - self.index[-3:]
        with pytest.raises(TypeError):
            self.index[-3:] - first
        with pytest.raises(TypeError):
            self.index[-3:] - first.tolist()
        with pytest.raises(TypeError):
            first.tolist() - self.index[-3:]

    def test_difference(self):

        first = self.index
        result = first.difference(self.index[-3:])
        expected = MultiIndex.from_tuples(sorted(self.index[:-3].values),
                                          sortorder=0,
                                          names=self.index.names)

        assert isinstance(result, MultiIndex)
        assert result.equals(expected)
        assert result.names == self.index.names

        # empty difference: reflexive
        result = self.index.difference(self.index)
        expected = self.index[:0]
        assert result.equals(expected)
        assert result.names == self.index.names

        # empty difference: superset
        result = self.index[-3:].difference(self.index)
        expected = self.index[:0]
        assert result.equals(expected)
        assert result.names == self.index.names

        # empty difference: degenerate
        result = self.index[:0].difference(self.index)
        expected = self.index[:0]
        assert result.equals(expected)
        assert result.names == self.index.names

        # names not the same
        chunklet = self.index[-3:]
        chunklet.names = ['foo', 'baz']
        result = first.difference(chunklet)
        assert result.names == (None, None)

        # empty, but non-equal
        result = self.index.difference(self.index.sortlevel(1)[0])
        assert len(result) == 0

        # raise Exception called with non-MultiIndex
        result = first.difference(first.values)
        assert result.equals(first[:0])

        # name from empty array
        result = first.difference([])
        assert first.equals(result)
        assert first.names == result.names

        # name from non-empty array
        result = first.difference([('foo', 'one')])
        expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), (
            'foo', 'two'), ('qux', 'one'), ('qux', 'two')])
        expected.names = first.names
        assert first.names == result.names
        tm.assert_raises_regex(TypeError, "other must be a MultiIndex "
                                          "or a list of tuples",
                               first.difference, [1, 2, 3, 4, 5])

    def test_argsort(self):
        result = self.index.argsort()
        expected = self.index.values.argsort()
        tm.assert_numpy_array_equal(result, expected)

    def test_dims(self):
        pass

    def test_insert(self):
        # key contained in all levels
        new_index = self.index.insert(0, ('bar', 'two'))
        assert new_index.equal_levels(self.index)
        assert new_index[0] == ('bar', 'two')

        # key not contained in all levels
        new_index = self.index.insert(0, ('abc', 'three'))

        exp0 = Index(list(self.index.levels[0]) + ['abc'], name='first')
        tm.assert_index_equal(new_index.levels[0], exp0)

        exp1 = Index(list(self.index.levels[1]) + ['three'], name='second')
        tm.assert_index_equal(new_index.levels[1], exp1)
        assert new_index[0] == ('abc', 'three')

        # key wrong length
        msg = "Item must have length equal to number of levels"
        with tm.assert_raises_regex(ValueError, msg):
            self.index.insert(0, ('foo2',))

        left = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1]],
                            columns=['1st', '2nd', '3rd'])
        left.set_index(['1st', '2nd'], inplace=True)
        ts = left['3rd'].copy(deep=True)

        left.loc[('b', 'x'), '3rd'] = 2
        left.loc[('b', 'a'), '3rd'] = -1
        left.loc[('b', 'b'), '3rd'] = 3
        left.loc[('a', 'x'), '3rd'] = 4
        left.loc[('a', 'w'), '3rd'] = 5
        left.loc[('a', 'a'), '3rd'] = 6

        ts.loc[('b', 'x')] = 2
        ts.loc['b', 'a'] = -1
        ts.loc[('b', 'b')] = 3
        ts.loc['a', 'x'] = 4
        ts.loc[('a', 'w')] = 5
        ts.loc['a', 'a'] = 6

        right = pd.DataFrame([['a', 'b', 0], ['b', 'd', 1], ['b', 'x', 2],
                              ['b', 'a', -1], ['b', 'b', 3], ['a', 'x', 4],
                              ['a', 'w', 5], ['a', 'a', 6]],
                             columns=['1st', '2nd', '3rd'])
        right.set_index(['1st', '2nd'], inplace=True)
        # FIXME data types changes to float because
        # of intermediate nan insertion;
        tm.assert_frame_equal(left, right, check_dtype=False)
        tm.assert_series_equal(ts, right['3rd'])

        # GH9250
        idx = [('test1', i) for i in range(5)] + \
              [('test2', i) for i in range(6)] + \
              [('test', 17), ('test', 18)]

        left = pd.Series(np.linspace(0, 10, 11),
                         pd.MultiIndex.from_tuples(idx[:-2]))

        left.loc[('test', 17)] = 11
        left.loc[('test', 18)] = 12

        right = pd.Series(np.linspace(0, 12, 13),
                          pd.MultiIndex.from_tuples(idx))

        tm.assert_series_equal(left, right)

    def test_take_fill_value(self):
        # GH 12631
        vals = [['A', 'B'],
                [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
        idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])

        result = idx.take(np.array([1, 0, -1]))
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    ('B', pd.Timestamp('2011-01-02'))]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    (np.nan, pd.NaT)]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        exp_vals = [('A', pd.Timestamp('2011-01-02')),
                    ('A', pd.Timestamp('2011-01-01')),
                    ('B', pd.Timestamp('2011-01-02'))]
        expected = pd.MultiIndex.from_tuples(exp_vals, names=['str', 'dt'])
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with pytest.raises(IndexError):
            idx.take(np.array([1, -5]))

    def take_invalid_kwargs(self):
        vals = [['A', 'B'],
                [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02')]]
        idx = pd.MultiIndex.from_product(vals, names=['str', 'dt'])
        indices = [1, 2]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        tm.assert_raises_regex(TypeError, msg, idx.take,
                               indices, foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, idx.take,
                               indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, idx.take,
                               indices, mode='clip')

    def test_tolist(self):
        result = self.index.tolist()
        exp = list(self.index.values)
        assert result == exp

    def test_str(self):
        # tested elsewhere
        pass

    def test_isna_behavior(self):
        # should not segfault GH5123
        # NOTE: if MI representation changes, may make sense to allow
        # isna(MI)
        with pytest.raises(NotImplementedError):
            pd.isna(self.index)

    def test_level_setting_resets_attributes(self):
        ind = pd.MultiIndex.from_arrays([
            ['A', 'A', 'B', 'B', 'B'], [1, 2, 1, 2, 3]
        ])
        assert ind.is_monotonic
        ind.set_levels([['A', 'B'], [1, 3, 2]], inplace=True)
        # if this fails, probably didn't reset the cache correctly.
        assert not ind.is_monotonic

    def test_reconstruct_sort(self):

        # starts off lexsorted & monotonic
        mi = MultiIndex.from_arrays([
            ['A', 'A', 'B', 'B', 'B'], [1, 2, 1, 2, 3]
        ])
        assert mi.is_lexsorted()
        assert mi.is_monotonic

        recons = mi._sort_levels_monotonic()
        assert recons.is_lexsorted()
        assert recons.is_monotonic
        assert mi is recons

        assert mi.equals(recons)
        assert Index(mi.values).equals(Index(recons.values))

        # cannot convert to lexsorted
        mi = pd.MultiIndex.from_tuples([('z', 'a'), ('x', 'a'), ('y', 'b'),
                                        ('x', 'b'), ('y', 'a'), ('z', 'b')],
                                       names=['one', 'two'])
        assert not mi.is_lexsorted()
        assert not mi.is_monotonic

        recons = mi._sort_levels_monotonic()
        assert not recons.is_lexsorted()
        assert not recons.is_monotonic

        assert mi.equals(recons)
        assert Index(mi.values).equals(Index(recons.values))

        # cannot convert to lexsorted
        mi = MultiIndex(levels=[['b', 'd', 'a'], [1, 2, 3]],
                        labels=[[0, 1, 0, 2], [2, 0, 0, 1]],
                        names=['col1', 'col2'])
        assert not mi.is_lexsorted()
        assert not mi.is_monotonic

        recons = mi._sort_levels_monotonic()
        assert not recons.is_lexsorted()
        assert not recons.is_monotonic

        assert mi.equals(recons)
        assert Index(mi.values).equals(Index(recons.values))

    def test_reconstruct_remove_unused(self):
        # xref to GH 2770
        df = DataFrame([['deleteMe', 1, 9],
                        ['keepMe', 2, 9],
                        ['keepMeToo', 3, 9]],
                       columns=['first', 'second', 'third'])
        df2 = df.set_index(['first', 'second'], drop=False)
        df2 = df2[df2['first'] != 'deleteMe']

        # removed levels are there
        expected = MultiIndex(levels=[['deleteMe', 'keepMe', 'keepMeToo'],
                                      [1, 2, 3]],
                              labels=[[1, 2], [1, 2]],
                              names=['first', 'second'])
        result = df2.index
        tm.assert_index_equal(result, expected)

        expected = MultiIndex(levels=[['keepMe', 'keepMeToo'],
                                      [2, 3]],
                              labels=[[0, 1], [0, 1]],
                              names=['first', 'second'])
        result = df2.index.remove_unused_levels()
        tm.assert_index_equal(result, expected)

        # idempotent
        result2 = result.remove_unused_levels()
        tm.assert_index_equal(result2, expected)
        assert result2.is_(result)

    @pytest.mark.parametrize('level0', [['a', 'd', 'b'],
                                        ['a', 'd', 'b', 'unused']])
    @pytest.mark.parametrize('level1', [['w', 'x', 'y', 'z'],
                                        ['w', 'x', 'y', 'z', 'unused']])
    def test_remove_unused_nan(self, level0, level1):
        # GH 18417
        mi = pd.MultiIndex(levels=[level0, level1],
                           labels=[[0, 2, -1, 1, -1], [0, 1, 2, 3, 2]])

        result = mi.remove_unused_levels()
        tm.assert_index_equal(result, mi)
        for level in 0, 1:
            assert('unused' not in result.levels[level])

    @pytest.mark.parametrize('first_type,second_type', [
        ('int64', 'int64'),
        ('datetime64[D]', 'str')])
    def test_remove_unused_levels_large(self, first_type, second_type):
        # GH16556

        # because tests should be deterministic (and this test in particular
        # checks that levels are removed, which is not the case for every
        # random input):
        rng = np.random.RandomState(4)  # seed is arbitrary value that works

        size = 1 << 16
        df = DataFrame(dict(
            first=rng.randint(0, 1 << 13, size).astype(first_type),
            second=rng.randint(0, 1 << 10, size).astype(second_type),
            third=rng.rand(size)))
        df = df.groupby(['first', 'second']).sum()
        df = df[df.third < 0.1]

        result = df.index.remove_unused_levels()
        assert len(result.levels[0]) < len(df.index.levels[0])
        assert len(result.levels[1]) < len(df.index.levels[1])
        assert result.equals(df.index)

        expected = df.reset_index().set_index(['first', 'second']).index
        tm.assert_index_equal(result, expected)

    def test_groupby(self):
        groups = self.index.groupby(np.array([1, 1, 1, 2, 2, 2]))
        labels = self.index.get_values().tolist()
        exp = {1: labels[:3], 2: labels[3:]}
        tm.assert_dict_equal(groups, exp)

        # GH5620
        groups = self.index.groupby(self.index)
        exp = {key: [key] for key in self.index}
        tm.assert_dict_equal(groups, exp)

    def test_equals_operator(self):
        # GH9785
        assert (self.index == self.index).all()

    def test_large_multiindex_error(self):
        # GH12527
        df_below_1000000 = pd.DataFrame(
            1, index=pd.MultiIndex.from_product([[1, 2], range(499999)]),
            columns=['dest'])
        with pytest.raises(KeyError):
            df_below_1000000.loc[(-1, 0), 'dest']
        with pytest.raises(KeyError):
            df_below_1000000.loc[(3, 0), 'dest']
        df_above_1000000 = pd.DataFrame(
            1, index=pd.MultiIndex.from_product([[1, 2], range(500001)]),
            columns=['dest'])
        with pytest.raises(KeyError):
            df_above_1000000.loc[(-1, 0), 'dest']
        with pytest.raises(KeyError):
            df_above_1000000.loc[(3, 0), 'dest']

    def test_partial_string_timestamp_multiindex(self):
        # GH10331
        dr = pd.date_range('2016-01-01', '2016-01-03', freq='12H')
        abc = ['a', 'b', 'c']
        ix = pd.MultiIndex.from_product([dr, abc])
        df = pd.DataFrame({'c1': range(0, 15)}, index=ix)
        idx = pd.IndexSlice

        #                        c1
        # 2016-01-01 00:00:00 a   0
        #                     b   1
        #                     c   2
        # 2016-01-01 12:00:00 a   3
        #                     b   4
        #                     c   5
        # 2016-01-02 00:00:00 a   6
        #                     b   7
        #                     c   8
        # 2016-01-02 12:00:00 a   9
        #                     b  10
        #                     c  11
        # 2016-01-03 00:00:00 a  12
        #                     b  13
        #                     c  14

        # partial string matching on a single index
        for df_swap in (df.swaplevel(),
                        df.swaplevel(0),
                        df.swaplevel(0, 1)):
            df_swap = df_swap.sort_index()
            just_a = df_swap.loc['a']
            result = just_a.loc['2016-01-01']
            expected = df.loc[idx[:, 'a'], :].iloc[0:2]
            expected.index = expected.index.droplevel(1)
            tm.assert_frame_equal(result, expected)

        # indexing with IndexSlice
        result = df.loc[idx['2016-01-01':'2016-02-01', :], :]
        expected = df
        tm.assert_frame_equal(result, expected)

        # match on secondary index
        result = df_swap.loc[idx[:, '2016-01-01':'2016-01-01'], :]
        expected = df_swap.iloc[[0, 1, 5, 6, 10, 11]]
        tm.assert_frame_equal(result, expected)

        # Even though this syntax works on a single index, this is somewhat
        # ambiguous and we don't want to extend this behavior forward to work
        # in multi-indexes. This would amount to selecting a scalar from a
        # column.
        with pytest.raises(KeyError):
            df['2016-01-01']

        # partial string match on year only
        result = df.loc['2016']
        expected = df
        tm.assert_frame_equal(result, expected)

        # partial string match on date
        result = df.loc['2016-01-01']
        expected = df.iloc[0:6]
        tm.assert_frame_equal(result, expected)

        # partial string match on date and hour, from middle
        result = df.loc['2016-01-02 12']
        expected = df.iloc[9:12]
        tm.assert_frame_equal(result, expected)

        # partial string match on secondary index
        result = df_swap.loc[idx[:, '2016-01-02'], :]
        expected = df_swap.iloc[[2, 3, 7, 8, 12, 13]]
        tm.assert_frame_equal(result, expected)

        # tuple selector with partial string match on date
        result = df.loc[('2016-01-01', 'a'), :]
        expected = df.iloc[[0, 3]]
        tm.assert_frame_equal(result, expected)

        # Slicing date on first level should break (of course)
        with pytest.raises(KeyError):
            df_swap.loc['2016-01-01']

        # GH12685 (partial string with daily resolution or below)
        dr = date_range('2013-01-01', periods=100, freq='D')
        ix = MultiIndex.from_product([dr, ['a', 'b']])
        df = DataFrame(np.random.randn(200, 1), columns=['A'], index=ix)

        result = df.loc[idx['2013-03':'2013-03', :], :]
        expected = df.iloc[118:180]
        tm.assert_frame_equal(result, expected)

    def test_rangeindex_fallback_coercion_bug(self):
        # GH 12893
        foo = pd.DataFrame(np.arange(100).reshape((10, 10)))
        bar = pd.DataFrame(np.arange(100).reshape((10, 10)))
        df = pd.concat({'foo': foo.stack(), 'bar': bar.stack()}, axis=1)
        df.index.names = ['fizz', 'buzz']

        str(df)
        expected = pd.DataFrame({'bar': np.arange(100),
                                 'foo': np.arange(100)},
                                index=pd.MultiIndex.from_product(
                                    [range(10), range(10)],
                                    names=['fizz', 'buzz']))
        tm.assert_frame_equal(df, expected, check_like=True)

        result = df.index.get_level_values('fizz')
        expected = pd.Int64Index(np.arange(10), name='fizz').repeat(10)
        tm.assert_index_equal(result, expected)

        result = df.index.get_level_values('buzz')
        expected = pd.Int64Index(np.tile(np.arange(10), 10), name='buzz')
        tm.assert_index_equal(result, expected)

    def test_dropna(self):
        # GH 6194
        idx = pd.MultiIndex.from_arrays([[1, np.nan, 3, np.nan, 5],
                                         [1, 2, np.nan, np.nan, 5],
                                         ['a', 'b', 'c', np.nan, 'e']])

        exp = pd.MultiIndex.from_arrays([[1, 5],
                                         [1, 5],
                                         ['a', 'e']])
        tm.assert_index_equal(idx.dropna(), exp)
        tm.assert_index_equal(idx.dropna(how='any'), exp)

        exp = pd.MultiIndex.from_arrays([[1, np.nan, 3, 5],
                                         [1, 2, np.nan, 5],
                                         ['a', 'b', 'c', 'e']])
        tm.assert_index_equal(idx.dropna(how='all'), exp)

        msg = "invalid how option: xxx"
        with tm.assert_raises_regex(ValueError, msg):
            idx.dropna(how='xxx')

    def test_unsortedindex(self):
        # GH 11897
        mi = pd.MultiIndex.from_tuples([('z', 'a'), ('x', 'a'), ('y', 'b'),
                                        ('x', 'b'), ('y', 'a'), ('z', 'b')],
                                       names=['one', 'two'])
        df = pd.DataFrame([[i, 10 * i] for i in lrange(6)], index=mi,
                          columns=['one', 'two'])

        # GH 16734: not sorted, but no real slicing
        result = df.loc(axis=0)['z', 'a']
        expected = df.iloc[0]
        tm.assert_series_equal(result, expected)

        with pytest.raises(UnsortedIndexError):
            df.loc(axis=0)['z', slice('a')]
        df.sort_index(inplace=True)
        assert len(df.loc(axis=0)['z', :]) == 2

        with pytest.raises(KeyError):
            df.loc(axis=0)['q', :]

    def test_unsortedindex_doc_examples(self):
        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#sorting-a-multiindex  # noqa
        dfm = DataFrame({'jim': [0, 0, 1, 1],
                         'joe': ['x', 'x', 'z', 'y'],
                         'jolie': np.random.rand(4)})

        dfm = dfm.set_index(['jim', 'joe'])
        with tm.assert_produces_warning(PerformanceWarning):
            dfm.loc[(1, 'z')]

        with pytest.raises(UnsortedIndexError):
            dfm.loc[(0, 'y'):(1, 'z')]

        assert not dfm.index.is_lexsorted()
        assert dfm.index.lexsort_depth == 1

        # sort it
        dfm = dfm.sort_index()
        dfm.loc[(1, 'z')]
        dfm.loc[(0, 'y'):(1, 'z')]

        assert dfm.index.is_lexsorted()
        assert dfm.index.lexsort_depth == 2

    def test_nan_stays_float(self):

        # GH 7031
        idx0 = pd.MultiIndex(levels=[["A", "B"], []],
                             labels=[[1, 0], [-1, -1]],
                             names=[0, 1])
        idx1 = pd.MultiIndex(levels=[["C"], ["D"]],
                             labels=[[0], [0]],
                             names=[0, 1])
        idxm = idx0.join(idx1, how='outer')
        assert pd.isna(idx0.get_level_values(1)).all()
        # the following failed in 0.14.1
        assert pd.isna(idxm.get_level_values(1)[:-1]).all()

        df0 = pd.DataFrame([[1, 2]], index=idx0)
        df1 = pd.DataFrame([[3, 4]], index=idx1)
        dfm = df0 - df1
        assert pd.isna(df0.index.get_level_values(1)).all()
        # the following failed in 0.14.1
        assert pd.isna(dfm.index.get_level_values(1)[:-1]).all()

    def test_million_record_attribute_error(self):
        # GH 18165
        r = list(range(1000000))
        df = pd.DataFrame({'a': r, 'b': r},
                          index=pd.MultiIndex.from_tuples([(x, x) for x in r]))

        with tm.assert_raises_regex(AttributeError,
                                    "'Series' object has no attribute 'foo'"):
            df['a'].foo()
