# -*- coding: utf-8 -*-
# pylint: disable=E1101,E1103,W0232

from datetime import datetime, timedelta, time
from pandas.compat import range, lrange, lzip, u, zip, PY3
import operator
import re
import nose
import warnings
import os

import numpy as np

from pandas import (period_range, date_range, Categorical, Series,
                    Index, Float64Index, Int64Index, MultiIndex,
                    CategoricalIndex, DatetimeIndex, TimedeltaIndex, PeriodIndex)
from pandas.core.index import InvalidIndexError, NumericIndex
from pandas.util.testing import (assert_almost_equal, assertRaisesRegexp,
                                 assert_copy)
from pandas import compat
from pandas.compat import long, is_platform_windows

import pandas.util.testing as tm
import pandas.core.config as cf

from pandas.tseries.index import _to_m8
import pandas.tseries.offsets as offsets

import pandas as pd
from pandas.lib import Timestamp
from itertools import product


class Base(object):
    """ base class for index sub-class tests """
    _holder = None
    _compat_props = ['shape', 'ndim', 'size', 'itemsize', 'nbytes']

    def setup_indices(self):
        # setup the test indices in the self.indicies dict
        for name, ind in self.indices.items():
            setattr(self, name, ind)

    def verify_pickle(self,index):
        unpickled = self.round_trip_pickle(index)
        self.assertTrue(index.equals(unpickled))

    def test_pickle_compat_construction(self):
        # this is testing for pickle compat
        if self._holder is None:
            return

        # need an object to create with
        self.assertRaises(TypeError, self._holder)

    def test_shift(self):

        # GH8083 test the base class for shift
        idx = self.create_index()
        self.assertRaises(NotImplementedError, idx.shift, 1)
        self.assertRaises(NotImplementedError, idx.shift, 1, 2)

    def test_numeric_compat(self):

        idx = self.create_index()
        tm.assertRaisesRegexp(TypeError,
                              "cannot perform __mul__",
                              lambda : idx * 1)
        tm.assertRaisesRegexp(TypeError,
                              "cannot perform __mul__",
                              lambda : 1 * idx)

        div_err = "cannot perform __truediv__" if compat.PY3 else "cannot perform __div__"
        tm.assertRaisesRegexp(TypeError,
                              div_err,
                              lambda : idx / 1)
        tm.assertRaisesRegexp(TypeError,
                              div_err,
                              lambda : 1 / idx)
        tm.assertRaisesRegexp(TypeError,
                              "cannot perform __floordiv__",
                              lambda : idx // 1)
        tm.assertRaisesRegexp(TypeError,
                              "cannot perform __floordiv__",
                              lambda : 1 // idx)

    def test_logical_compat(self):
        idx = self.create_index()
        tm.assertRaisesRegexp(TypeError,
                              'cannot perform all',
                              lambda : idx.all())
        tm.assertRaisesRegexp(TypeError,
                              'cannot perform any',
                              lambda : idx.any())

    def test_boolean_context_compat(self):

        # boolean context compat
        idx = self.create_index()
        def f():
            if idx:
                pass
        tm.assertRaisesRegexp(ValueError,'The truth value of a',f)

    def test_reindex_base(self):
        idx = self.create_index()
        expected = np.arange(idx.size)

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with tm.assertRaisesRegexp(ValueError, 'Invalid fill method'):
            idx.get_indexer(idx, method='invalid')

    def test_ndarray_compat_properties(self):

        idx = self.create_index()
        self.assertTrue(idx.T.equals(idx))
        self.assertTrue(idx.transpose().equals(idx))

        values = idx.values
        for prop in self._compat_props:
            self.assertEqual(getattr(idx, prop), getattr(values, prop))

        # test for validity
        idx.nbytes
        idx.values.nbytes

    def test_repr_roundtrip(self):

        idx = self.create_index()
        tm.assert_index_equal(eval(repr(idx)),idx)

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = 'foo'
        self.assertTrue("'foo'" in str(idx))
        self.assertTrue(idx.__class__.__name__ in str(idx))

    def test_dtype_str(self):
        for idx in self.indices.values():
            dtype = idx.dtype_str
            self.assertIsInstance(dtype, compat.string_types)
            if isinstance(idx, PeriodIndex):
                self.assertEqual(dtype, 'period')
            else:
                self.assertEqual(dtype, str(idx.dtype))

    def test_repr_max_seq_item_setting(self):
        # GH10182
        idx = self.create_index()
        idx = idx.repeat(50)
        with pd.option_context("display.max_seq_items", None):
            repr(idx)
            self.assertFalse('...' in str(idx))

    def test_wrong_number_names(self):
        def testit(ind):
            ind.names = ["apple", "banana", "carrot"]

        for ind in self.indices.values():
            assertRaisesRegexp(ValueError, "^Length", testit, ind)

    def test_set_name_methods(self):
        new_name = "This is the new name for this index"
        for ind in self.indices.values():

            # don't tests a MultiIndex here (as its tested separated)
            if isinstance(ind, MultiIndex):
                continue

            original_name = ind.name
            new_ind = ind.set_names([new_name])
            self.assertEqual(new_ind.name, new_name)
            self.assertEqual(ind.name, original_name)
            res = ind.rename(new_name, inplace=True)

            # should return None
            self.assertIsNone(res)
            self.assertEqual(ind.name, new_name)
            self.assertEqual(ind.names, [new_name])
            #with assertRaisesRegexp(TypeError, "list-like"):
            #    # should still fail even if it would be the right length
            #    ind.set_names("a")
            with assertRaisesRegexp(ValueError, "Level must be None"):
                ind.set_names("a", level=0)

            # rename in place just leaves tuples and other containers alone
            name = ('A', 'B')
            ind.rename(name, inplace=True)
            self.assertEqual(ind.name, name)
            self.assertEqual(ind.names, [name])

    def test_hash_error(self):
        for ind in self.indices.values():
            with tm.assertRaisesRegexp(TypeError,
                                       "unhashable type: %r" %
                                       type(ind).__name__):
                hash(ind)

    def test_copy_and_deepcopy(self):
        from copy import copy, deepcopy

        for ind in self.indices.values():

            # don't tests a MultiIndex here (as its tested separated)
            if isinstance(ind, MultiIndex):
                continue

            for func in (copy, deepcopy):
                idx_copy = func(ind)
                self.assertIsNot(idx_copy, ind)
                self.assertTrue(idx_copy.equals(ind))

            new_copy = ind.copy(deep=True, name="banana")
            self.assertEqual(new_copy.name, "banana")

    def test_duplicates(self):
        for ind in self.indices.values():

            if not len(ind):
                continue
            if isinstance(ind, MultiIndex):
                continue
            idx = self._holder([ind[0]]*5)
            self.assertFalse(idx.is_unique)
            self.assertTrue(idx.has_duplicates)

            # GH 10115
            # preserve names
            idx.name = 'foo'
            result = idx.drop_duplicates()
            self.assertEqual(result.name, 'foo')
            self.assert_index_equal(result, Index([ind[0]],name='foo'))

    def test_sort(self):
        for ind in self.indices.values():
            self.assertRaises(TypeError, ind.sort)

    def test_order(self):
        for ind in self.indices.values():
            # 9816 deprecated
            with tm.assert_produces_warning(FutureWarning):
                ind.order()

    def test_mutability(self):
        for ind in self.indices.values():
            if not len(ind):
                continue
            self.assertRaises(TypeError, ind.__setitem__, 0, ind[0])

    def test_view(self):
        for ind in self.indices.values():
            i_view = ind.view()
            self.assertEqual(i_view.name, ind.name)

    def test_compat(self):
        for ind in self.indices.values():
            self.assertEqual(ind.tolist(),list(ind))

    def test_argsort(self):
        for k, ind in self.indices.items():

            # sep teststed
            if k in ['catIndex']:
                continue

            result = ind.argsort()
            expected = np.array(ind).argsort()
            tm.assert_numpy_array_equal(result, expected)

    def test_pickle(self):
        for ind in self.indices.values():
            self.verify_pickle(ind)
            ind.name = 'foo'
            self.verify_pickle(ind)

    def test_take(self):
        indexer = [4, 3, 0, 2]
        for k, ind in self.indices.items():

            # separate
            if k in ['boolIndex','tuples','empty']:
                continue

            result = ind.take(indexer)
            expected = ind[indexer]
            self.assertTrue(result.equals(expected))

            if not isinstance(ind, (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
                # GH 10791
                with tm.assertRaises(AttributeError):
                    ind.freq

    def test_setops_errorcases(self):
        for name, idx in compat.iteritems(self.indices):
            # # non-iterable input
            cases = [0.5, 'xxx']
            methods = [idx.intersection, idx.union, idx.difference, idx.sym_diff]

            for method in methods:
                for case in cases:
                    assertRaisesRegexp(TypeError,
                                       "Input must be Index or array-like",
                                       method, case)

    def test_intersection_base(self):
        for name, idx in compat.iteritems(self.indices):
            first = idx[:5]
            second = idx[:3]
            intersect = first.intersection(second)

            if isinstance(idx, CategoricalIndex):
                pass
            else:
                self.assertTrue(tm.equalContents(intersect, second))

            # GH 10149
            cases = [klass(second.values) for klass in [np.array, Series, list]]
            for case in cases:
                if isinstance(idx, PeriodIndex):
                    msg = "can only call with other PeriodIndex-ed objects"
                    with tm.assertRaisesRegexp(ValueError, msg):
                        result = first.intersection(case)
                elif isinstance(idx, CategoricalIndex):
                    pass
                else:
                    result = first.intersection(case)
                    self.assertTrue(tm.equalContents(result, second))

            if isinstance(idx, MultiIndex):
                msg = "other must be a MultiIndex or a list of tuples"
                with tm.assertRaisesRegexp(TypeError, msg):
                    result = first.intersection([1, 2, 3])

    def test_union_base(self):
        for name, idx in compat.iteritems(self.indices):
            first = idx[3:]
            second = idx[:5]
            everything = idx
            union = first.union(second)
            self.assertTrue(tm.equalContents(union, everything))

            # GH 10149
            cases = [klass(second.values) for klass in [np.array, Series, list]]
            for case in cases:
                if isinstance(idx, PeriodIndex):
                    msg = "can only call with other PeriodIndex-ed objects"
                    with tm.assertRaisesRegexp(ValueError, msg):
                        result = first.union(case)
                elif isinstance(idx, CategoricalIndex):
                    pass
                else:
                    result = first.union(case)
                    self.assertTrue(tm.equalContents(result, everything))

            if isinstance(idx, MultiIndex):
                msg = "other must be a MultiIndex or a list of tuples"
                with tm.assertRaisesRegexp(TypeError, msg):
                    result = first.union([1, 2, 3])

    def test_difference_base(self):
        for name, idx in compat.iteritems(self.indices):
            first = idx[2:]
            second = idx[:4]
            answer = idx[4:]
            result = first.difference(second)

            if isinstance(idx, CategoricalIndex):
                pass
            else:
                self.assertTrue(tm.equalContents(result, answer))

            # GH 10149
            cases = [klass(second.values) for klass in [np.array, Series, list]]
            for case in cases:
                if isinstance(idx, PeriodIndex):
                    msg = "can only call with other PeriodIndex-ed objects"
                    with tm.assertRaisesRegexp(ValueError, msg):
                        result = first.difference(case)
                elif isinstance(idx, CategoricalIndex):
                    pass
                elif isinstance(idx, (DatetimeIndex, TimedeltaIndex)):
                    self.assertEqual(result.__class__, answer.__class__)
                    tm.assert_numpy_array_equal(result.asi8, answer.asi8)
                else:
                    result = first.difference(case)
                    self.assertTrue(tm.equalContents(result, answer))

            if isinstance(idx, MultiIndex):
                msg = "other must be a MultiIndex or a list of tuples"
                with tm.assertRaisesRegexp(TypeError, msg):
                    result = first.difference([1, 2, 3])

    def test_symmetric_diff(self):
        for name, idx in compat.iteritems(self.indices):
            first = idx[1:]
            second = idx[:-1]
            if isinstance(idx, CategoricalIndex):
                pass
            else:
                answer = idx[[0, -1]]
                result = first.sym_diff(second)
                self.assertTrue(tm.equalContents(result, answer))

            # GH 10149
            cases = [klass(second.values) for klass in [np.array, Series, list]]
            for case in cases:
                if isinstance(idx, PeriodIndex):
                    msg = "can only call with other PeriodIndex-ed objects"
                    with tm.assertRaisesRegexp(ValueError, msg):
                        result = first.sym_diff(case)
                elif isinstance(idx, CategoricalIndex):
                    pass
                else:
                    result = first.sym_diff(case)
                    self.assertTrue(tm.equalContents(result, answer))

            if isinstance(idx, MultiIndex):
                msg = "other must be a MultiIndex or a list of tuples"
                with tm.assertRaisesRegexp(TypeError, msg):
                    result = first.sym_diff([1, 2, 3])

    def test_insert_base(self):

        for name, idx in compat.iteritems(self.indices):
            result = idx[1:4]

            if not len(idx):
                continue

            #test 0th element
            self.assertTrue(idx[0:4].equals(
                result.insert(0, idx[0])))

    def test_delete_base(self):

        for name, idx in compat.iteritems(self.indices):

            if not len(idx):
                continue

            expected = idx[1:]
            result = idx.delete(0)
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)

            expected = idx[:-1]
            result = idx.delete(-1)
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)

            with tm.assertRaises((IndexError, ValueError)):
                # either depending on numpy version
                result = idx.delete(len(idx))

    def test_equals_op(self):
        # GH9947, GH10637
        index_a = self.create_index()
        if isinstance(index_a, PeriodIndex):
            return

        n = len(index_a)
        index_b = index_a[0:-1]
        index_c = index_a[0:-1].append(index_a[-2:-1])
        index_d = index_a[0:1]
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == index_b
        expected1 = np.array([True] * n)
        expected2 = np.array([True] * (n - 1) + [False])
        tm.assert_numpy_array_equal(index_a == index_a, expected1)
        tm.assert_numpy_array_equal(index_a == index_c, expected2)

        # test comparisons with numpy arrays
        array_a = np.array(index_a)
        array_b = np.array(index_a[0:-1])
        array_c = np.array(index_a[0:-1].append(index_a[-2:-1]))
        array_d = np.array(index_a[0:1])
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == array_b
        tm.assert_numpy_array_equal(index_a == array_a, expected1)
        tm.assert_numpy_array_equal(index_a == array_c, expected2)

        # test comparisons with Series
        series_a = Series(array_a)
        series_b = Series(array_b)
        series_c = Series(array_c)
        series_d = Series(array_d)
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == series_b
        tm.assert_numpy_array_equal(index_a == series_a, expected1)
        tm.assert_numpy_array_equal(index_a == series_c, expected2)

        # cases where length is 1 for one of them
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == index_d
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == series_d
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            index_a == array_d
        with tm.assertRaisesRegexp(ValueError, "Series lengths must match"):
            series_a == series_d
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            series_a == array_d

        # comparing with a scalar should broadcast; note that we are excluding
        # MultiIndex because in this case each item in the index is a tuple of
        # length 2, and therefore is considered an array of length 2 in the
        # comparison instead of a scalar
        if not isinstance(index_a, MultiIndex):
            expected3 = np.array([False] * (len(index_a) - 2) + [True, False])
            # assuming the 2nd to last item is unique in the data
            item = index_a[-2]
            tm.assert_numpy_array_equal(index_a == item, expected3)
            tm.assert_numpy_array_equal(series_a == item, expected3)

    def test_numpy_ufuncs(self):
        # test ufuncs of numpy 1.9.2. see:
        # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

        # some functions are skipped because it may return different result
        # for unicode input depending on numpy version

        for name, idx in compat.iteritems(self.indices):
            for func in [np.exp, np.exp2, np.expm1, np.log, np.log2, np.log10,
                         np.log1p, np.sqrt, np.sin, np.cos,
                         np.tan, np.arcsin, np.arccos, np.arctan,
                         np.sinh, np.cosh, np.tanh, np.arcsinh, np.arccosh,
                         np.arctanh, np.deg2rad, np.rad2deg]:
                if isinstance(idx, pd.tseries.base.DatetimeIndexOpsMixin):
                    # raise TypeError or ValueError (PeriodIndex)
                    # PeriodIndex behavior should be changed in future version
                    with tm.assertRaises(Exception):
                        func(idx)
                elif isinstance(idx, (Float64Index, Int64Index)):
                    # coerces to float (e.g. np.sin)
                    result = func(idx)
                    exp = Index(func(idx.values), name=idx.name)
                    self.assert_index_equal(result, exp)
                    self.assertIsInstance(result, pd.Float64Index)
                else:
                    # raise AttributeError or TypeError
                    if len(idx) == 0:
                        continue
                    else:
                        with tm.assertRaises(Exception):
                            func(idx)

            for func in [np.isfinite, np.isinf, np.isnan, np.signbit]:
                if isinstance(idx, pd.tseries.base.DatetimeIndexOpsMixin):
                    # raise TypeError or ValueError (PeriodIndex)
                    with tm.assertRaises(Exception):
                        func(idx)
                elif isinstance(idx, (Float64Index, Int64Index)):
                    # results in bool array
                    result = func(idx)
                    exp = func(idx.values)
                    self.assertIsInstance(result, np.ndarray)
                    tm.assertNotIsInstance(result, Index)
                else:
                    if len(idx) == 0:
                        continue
                    else:
                        with tm.assertRaises(Exception):
                            func(idx)

    def test_hasnans_isnans(self):
        # GH 11343, added tests for hasnans / isnans
        for name, index in self.indices.items():
            if isinstance(index, MultiIndex):
                pass
            else:
                idx = index.copy()

                # cases in indices doesn't include NaN
                expected = np.array([False] * len(idx), dtype=bool)
                self.assert_numpy_array_equal(idx._isnan, expected)
                self.assertFalse(idx.hasnans)

                idx = index.copy()
                values = idx.values

                if len(index) == 0:
                    continue
                elif isinstance(index, pd.tseries.base.DatetimeIndexOpsMixin):
                    values[1] = pd.tslib.iNaT
                elif isinstance(index, Int64Index):
                    continue
                else:
                    values[1] = np.nan

                if isinstance(index, PeriodIndex):
                    idx = index.__class__(values, freq=index.freq)
                else:
                    idx = index.__class__(values)

                expected = np.array([False] * len(idx), dtype=bool)
                expected[1] = True
                self.assert_numpy_array_equal(idx._isnan, expected)
                self.assertTrue(idx.hasnans)

    def test_fillna(self):
        # GH 11343
        for name, index in self.indices.items():
            if len(index) == 0:
                pass
            elif isinstance(index, MultiIndex):
                idx = index.copy()
                msg = "isnull is not defined for MultiIndex"
                with self.assertRaisesRegexp(NotImplementedError, msg):
                    idx.fillna(idx[0])
            else:
                idx = index.copy()
                result = idx.fillna(idx[0])
                self.assert_index_equal(result, idx)
                self.assertFalse(result is idx)

                msg = "'value' must be a scalar, passed: "
                with self.assertRaisesRegexp(TypeError, msg):
                    idx.fillna([idx[0]])

                idx = index.copy()
                values = idx.values

                if isinstance(index, pd.tseries.base.DatetimeIndexOpsMixin):
                    values[1] = pd.tslib.iNaT
                elif isinstance(index, Int64Index):
                    continue
                else:
                    values[1] = np.nan

                if isinstance(index, PeriodIndex):
                    idx = index.__class__(values, freq=index.freq)
                else:
                    idx = index.__class__(values)

                expected = np.array([False] * len(idx), dtype=bool)
                expected[1] = True
                self.assert_numpy_array_equal(idx._isnan, expected)
                self.assertTrue(idx.hasnans)


class TestIndex(Base, tm.TestCase):
    _holder = Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(
            unicodeIndex = tm.makeUnicodeIndex(100),
            strIndex = tm.makeStringIndex(100),
            dateIndex = tm.makeDateIndex(100),
            periodIndex = tm.makePeriodIndex(100),
            tdIndex = tm.makeTimedeltaIndex(100),
            intIndex = tm.makeIntIndex(100),
            floatIndex = tm.makeFloatIndex(100),
            boolIndex = Index([True,False]),
            catIndex = tm.makeCategoricalIndex(100),
            empty = Index([]),
            tuples = MultiIndex.from_tuples(lzip(['foo', 'bar', 'baz'],
                                                 [1, 2, 3]))
        )
        self.setup_indices()

    def create_index(self):
        return Index(list('abcde'))

    def test_new_axis(self):
        new_index = self.dateIndex[None, :]
        self.assertEqual(new_index.ndim, 2)
        tm.assertIsInstance(new_index, np.ndarray)

    def test_copy_and_deepcopy(self):
        super(TestIndex, self).test_copy_and_deepcopy()

        new_copy2 = self.intIndex.copy(dtype=int)
        self.assertEqual(new_copy2.dtype.kind, 'i')

    def test_constructor(self):
        # regular instance creation
        tm.assert_contains_all(self.strIndex, self.strIndex)
        tm.assert_contains_all(self.dateIndex, self.dateIndex)

        # casting
        arr = np.array(self.strIndex)
        index = Index(arr)
        tm.assert_contains_all(arr, index)
        tm.assert_numpy_array_equal(self.strIndex, index)

        # copy
        arr = np.array(self.strIndex)
        index = Index(arr, copy=True, name='name')
        tm.assertIsInstance(index, Index)
        self.assertEqual(index.name, 'name')
        tm.assert_numpy_array_equal(arr, index)
        arr[0] = "SOMEBIGLONGSTRING"
        self.assertNotEqual(index[0], "SOMEBIGLONGSTRING")

        # what to do here?
        # arr = np.array(5.)
        # self.assertRaises(Exception, arr.view, Index)

    def test_constructor_corner(self):
        # corner case
        self.assertRaises(TypeError, Index, 0)

    def test_construction_list_mixed_tuples(self):
        # 10697
        # if we are constructing from a mixed list of tuples, make sure that we
        # are independent of the sorting order
        idx1 = Index([('A',1),'B'])
        self.assertIsInstance(idx1, Index) and self.assertNotInstance(idx1, MultiIndex)
        idx2 = Index(['B',('A',1)])
        self.assertIsInstance(idx2, Index) and self.assertNotInstance(idx2, MultiIndex)

    def test_constructor_from_series(self):

        expected = DatetimeIndex([Timestamp('20110101'),Timestamp('20120101'),Timestamp('20130101')])
        s = Series([Timestamp('20110101'),Timestamp('20120101'),Timestamp('20130101')])
        result = Index(s)
        self.assertTrue(result.equals(expected))
        result = DatetimeIndex(s)
        self.assertTrue(result.equals(expected))

        # GH 6273
        # create from a series, passing a freq
        s = Series(pd.to_datetime(['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990']))
        result = DatetimeIndex(s, freq='MS')
        expected = DatetimeIndex(['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990'],freq='MS')
        self.assertTrue(result.equals(expected))

        df = pd.DataFrame(np.random.rand(5,3))
        df['date'] = ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990']
        result = DatetimeIndex(df['date'], freq='MS')
        self.assertTrue(result.equals(expected))
        self.assertEqual(df['date'].dtype, object)

        exp = pd.Series(['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990'], name='date')
        self.assert_series_equal(df['date'], exp)

        # GH 6274
        # infer freq of same
        result = pd.infer_freq(df['date'])
        self.assertEqual(result,'MS')

    def test_constructor_ndarray_like(self):
        # GH 5460#issuecomment-44474502
        # it should be possible to convert any object that satisfies the numpy
        # ndarray interface directly into an Index
        class ArrayLike(object):
            def __init__(self, array):
                self.array = array
            def __array__(self, dtype=None):
                return self.array

        for array in [np.arange(5),
                      np.array(['a', 'b', 'c']),
                      date_range('2000-01-01', periods=3).values]:
            expected = pd.Index(array)
            result = pd.Index(ArrayLike(array))
            self.assertTrue(result.equals(expected))

    def test_index_ctor_infer_periodindex(self):
        xp = period_range('2012-1-1', freq='M', periods=3)
        rs = Index(xp)
        tm.assert_numpy_array_equal(rs, xp)
        tm.assertIsInstance(rs, PeriodIndex)

    def test_constructor_simple_new(self):
        idx = Index([1, 2, 3, 4, 5], name='int')
        result = idx._simple_new(idx, 'int')
        self.assertTrue(result.equals(idx))

        idx = Index([1.1, np.nan, 2.2, 3.0], name='float')
        result = idx._simple_new(idx, 'float')
        self.assertTrue(result.equals(idx))

        idx = Index(['A', 'B', 'C', np.nan], name='obj')
        result = idx._simple_new(idx, 'obj')
        self.assertTrue(result.equals(idx))

    def test_constructor_dtypes(self):

        for idx in [Index(np.array([1, 2, 3], dtype=int)),
                    Index(np.array([1, 2, 3], dtype=int), dtype=int),
                    Index(np.array([1., 2., 3.], dtype=float), dtype=int),
                    Index([1, 2, 3], dtype=int),
                    Index([1., 2., 3.], dtype=int)]:
            self.assertIsInstance(idx, Int64Index)

        for idx in [Index(np.array([1., 2., 3.], dtype=float)),
                    Index(np.array([1, 2, 3], dtype=int), dtype=float),
                    Index(np.array([1., 2., 3.], dtype=float), dtype=float),
                    Index([1, 2, 3], dtype=float),
                    Index([1., 2., 3.], dtype=float)]:
            self.assertIsInstance(idx, Float64Index)

        for idx in [Index(np.array([True, False, True], dtype=bool)),
                    Index([True, False, True]),
                    Index(np.array([True, False, True], dtype=bool), dtype=bool),
                    Index([True, False, True], dtype=bool)]:
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

        for idx in [Index(np.array([1, 2, 3], dtype=int), dtype='category'),
                    Index([1, 2, 3], dtype='category'),
                    Index(np.array([np.datetime64('2011-01-01'), np.datetime64('2011-01-02')]), dtype='category'),
                    Index([datetime(2011, 1, 1), datetime(2011, 1, 2)], dtype='category')]:
            self.assertIsInstance(idx, CategoricalIndex)

        for idx in [Index(np.array([np.datetime64('2011-01-01'), np.datetime64('2011-01-02')])),
                    Index([datetime(2011, 1, 1), datetime(2011, 1, 2)])]:
            self.assertIsInstance(idx, DatetimeIndex)

        for idx in [Index(np.array([np.datetime64('2011-01-01'), np.datetime64('2011-01-02')]), dtype=object),
                    Index([datetime(2011, 1, 1), datetime(2011, 1, 2)], dtype=object)]:
            self.assertNotIsInstance(idx, DatetimeIndex)
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

        for idx in [Index(np.array([np.timedelta64(1, 'D'), np.timedelta64(1, 'D')])),
                    Index([timedelta(1), timedelta(1)])]:
            self.assertIsInstance(idx, TimedeltaIndex)

        for idx in [Index(np.array([np.timedelta64(1, 'D'), np.timedelta64(1, 'D')]), dtype=object),
                    Index([timedelta(1), timedelta(1)], dtype=object)]:
            self.assertNotIsInstance(idx, TimedeltaIndex)
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

    def test_view_with_args(self):

        restricted = ['unicodeIndex','strIndex','catIndex','boolIndex','empty']

        for i in restricted:
            ind = self.indices[i]

            # with arguments
            self.assertRaises(TypeError, lambda : ind.view('i8'))

        # these are ok
        for i in list(set(self.indices.keys())-set(restricted)):
            ind = self.indices[i]

            # with arguments
            ind.view('i8')

    def test_legacy_pickle_identity(self):

        # GH 8431
        pth = tm.get_data_path()
        s1 = pd.read_pickle(os.path.join(pth,'s1-0.12.0.pickle'))
        s2 = pd.read_pickle(os.path.join(pth,'s2-0.12.0.pickle'))
        self.assertFalse(s1.index.identical(s2.index))
        self.assertFalse(s1.index.equals(s2.index))

    def test_astype(self):
        casted = self.intIndex.astype('i8')

        # it works!
        casted.get_loc(5)

        # pass on name
        self.intIndex.name = 'foobar'
        casted = self.intIndex.astype('i8')
        self.assertEqual(casted.name, 'foobar')

    def test_equals(self):
        # same
        self.assertTrue(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'c'])))

        # different length
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b'])))

        # same length, different values
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'd'])))

        # Must also be an Index
        self.assertFalse(Index(['a', 'b', 'c']).equals(['a', 'b', 'c']))

    def test_insert(self):

        # GH 7256
        # validate neg/pos inserts
        result = Index(['b', 'c', 'd'])

        #test 0th element
        self.assertTrue(Index(['a', 'b', 'c', 'd']).equals(
            result.insert(0, 'a')))

        #test Nth element that follows Python list behavior
        self.assertTrue(Index(['b', 'c', 'e', 'd']).equals(
            result.insert(-1, 'e')))

        #test loc +/- neq (0, -1)
        self.assertTrue(result.insert(1, 'z').equals(
            result.insert(-2, 'z')))

        #test empty
        null_index = Index([])
        self.assertTrue(Index(['a']).equals(
            null_index.insert(0, 'a')))

    def test_delete(self):
        idx = Index(['a', 'b', 'c', 'd'], name='idx')

        expected = Index(['b', 'c', 'd'], name='idx')
        result = idx.delete(0)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)

        expected = Index(['a', 'b', 'c'], name='idx')
        result = idx.delete(-1)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

    def test_identical(self):

        # index
        i1 = Index(['a', 'b', 'c'])
        i2 = Index(['a', 'b', 'c'])

        self.assertTrue(i1.identical(i2))

        i1 = i1.rename('foo')
        self.assertTrue(i1.equals(i2))
        self.assertFalse(i1.identical(i2))

        i2 = i2.rename('foo')
        self.assertTrue(i1.identical(i2))

        i3 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')])
        i4 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')], tupleize_cols=False)
        self.assertFalse(i3.identical(i4))

    def test_is_(self):
        ind = Index(range(10))
        self.assertTrue(ind.is_(ind))
        self.assertTrue(ind.is_(ind.view().view().view().view()))
        self.assertFalse(ind.is_(Index(range(10))))
        self.assertFalse(ind.is_(ind.copy()))
        self.assertFalse(ind.is_(ind.copy(deep=False)))
        self.assertFalse(ind.is_(ind[:]))
        self.assertFalse(ind.is_(ind.view(np.ndarray).view(Index)))
        self.assertFalse(ind.is_(np.array(range(10))))

        # quasi-implementation dependent
        self.assertTrue(ind.is_(ind.view()))
        ind2 = ind.view()
        ind2.name = 'bob'
        self.assertTrue(ind.is_(ind2))
        self.assertTrue(ind2.is_(ind))
        # doesn't matter if Indices are *actually* views of underlying data,
        self.assertFalse(ind.is_(Index(ind.values)))
        arr = np.array(range(1, 11))
        ind1 = Index(arr, copy=False)
        ind2 = Index(arr, copy=False)
        self.assertFalse(ind1.is_(ind2))

    def test_asof(self):
        d = self.dateIndex[0]
        self.assertEqual(self.dateIndex.asof(d), d)
        self.assertTrue(np.isnan(self.dateIndex.asof(d - timedelta(1))))

        d = self.dateIndex[-1]
        self.assertEqual(self.dateIndex.asof(d + timedelta(1)), d)

        d = self.dateIndex[0].to_datetime()
        tm.assertIsInstance(self.dateIndex.asof(d), Timestamp)

    def test_asof_datetime_partial(self):
        idx = pd.date_range('2010-01-01', periods=2, freq='m')
        expected = Timestamp('2010-02-28')
        result = idx.asof('2010-02')
        self.assertEqual(result, expected)
        self.assertFalse(isinstance(result, Index))

    def test_nanosecond_index_access(self):
        s = Series([Timestamp('20130101')]).values.view('i8')[0]
        r = DatetimeIndex([s + 50 + i for i in range(100)])
        x = Series(np.random.randn(100), index=r)

        first_value = x.asof(x.index[0])

        # this does not yet work, as parsing strings is done via dateutil
        #self.assertEqual(first_value, x['2013-01-01 00:00:00.000000050+0000'])

        self.assertEqual(first_value, x[Timestamp(np.datetime64('2013-01-01 00:00:00.000000050+0000', 'ns'))])

    def test_comparators(self):
        index = self.dateIndex
        element = index[len(index) // 2]
        element = _to_m8(element)

        arr = np.array(index)

        def _check(op):
            arr_result = op(arr, element)
            index_result = op(index, element)

            self.assertIsInstance(index_result, np.ndarray)
            tm.assert_numpy_array_equal(arr_result, index_result)

        _check(operator.eq)
        _check(operator.ne)
        _check(operator.gt)
        _check(operator.lt)
        _check(operator.ge)
        _check(operator.le)

    def test_booleanindex(self):
        boolIdx = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIdx[5:30:2] = False

        subIndex = self.strIndex[boolIdx]

        for i, val in enumerate(subIndex):
            self.assertEqual(subIndex.get_loc(val), i)

        subIndex = self.strIndex[list(boolIdx)]
        for i, val in enumerate(subIndex):
            self.assertEqual(subIndex.get_loc(val), i)

    def test_fancy(self):
        sl = self.strIndex[[1, 2, 3]]
        for i in sl:
            self.assertEqual(i, sl[sl.get_loc(i)])

    def test_empty_fancy(self):
        empty_farr = np.array([], dtype=np.float_)
        empty_iarr = np.array([], dtype=np.int_)
        empty_barr = np.array([], dtype=np.bool_)

        # pd.DatetimeIndex is excluded, because it overrides getitem and should
        # be tested separately.
        for idx in [self.strIndex, self.intIndex, self.floatIndex]:
            empty_idx = idx.__class__([])
            values = idx.values

            self.assertTrue(idx[[]].identical(empty_idx))
            self.assertTrue(idx[empty_iarr].identical(empty_idx))
            self.assertTrue(idx[empty_barr].identical(empty_idx))

            # np.ndarray only accepts ndarray of int & bool dtypes, so should
            # Index.
            self.assertRaises(IndexError, idx.__getitem__, empty_farr)

    def test_getitem(self):
        arr = np.array(self.dateIndex)
        exp = self.dateIndex[5]
        exp = _to_m8(exp)

        self.assertEqual(exp, arr[5])

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)
        self.assertTrue(tm.equalContents(intersect, second))

        # Corner cases
        inter = first.intersection(first)
        self.assertIs(inter, first)

        idx1 = Index([1, 2, 3, 4, 5], name='idx')
        # if target has the same name, it is preserved
        idx2 = Index([3, 4, 5, 6, 7], name='idx')
        expected2 = Index([3, 4, 5], name='idx')
        result2 = idx1.intersection(idx2)
        self.assertTrue(result2.equals(expected2))
        self.assertEqual(result2.name, expected2.name)

        # if target name is different, it will be reset
        idx3 = Index([3, 4, 5, 6, 7], name='other')
        expected3 = Index([3, 4, 5], name=None)
        result3 = idx1.intersection(idx3)
        self.assertTrue(result3.equals(expected3))
        self.assertEqual(result3.name, expected3.name)

        # non monotonic
        idx1 = Index([5, 3, 2, 4, 1], name='idx')
        idx2 = Index([4, 7, 6, 5, 3], name='idx')
        result2 = idx1.intersection(idx2)
        self.assertTrue(tm.equalContents(result2, expected2))
        self.assertEqual(result2.name, expected2.name)

        idx3 = Index([4, 7, 6, 5, 3], name='other')
        result3 = idx1.intersection(idx3)
        self.assertTrue(tm.equalContents(result3, expected3))
        self.assertEqual(result3.name, expected3.name)

        # non-monotonic non-unique
        idx1 = Index(['A','B','A','C'])
        idx2 = Index(['B','D'])
        expected = Index(['B'], dtype='object')
        result = idx1.intersection(idx2)
        self.assertTrue(result.equals(expected))

    def test_union(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        everything = self.strIndex[:20]
        union = first.union(second)
        self.assertTrue(tm.equalContents(union, everything))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.union(case)
            self.assertTrue(tm.equalContents(result, everything))

        # Corner cases
        union = first.union(first)
        self.assertIs(union, first)

        union = first.union([])
        self.assertIs(union, first)

        union = Index([]).union(first)
        self.assertIs(union, first)

        # preserve names
        first.name = 'A'
        second.name = 'A'
        union = first.union(second)
        self.assertEqual(union.name, 'A')

        second.name = 'B'
        union = first.union(second)
        self.assertIsNone(union.name)

    def test_add(self):

        # - API change GH 8226
        with tm.assert_produces_warning():
            self.strIndex + self.strIndex
        with tm.assert_produces_warning():
            self.strIndex + self.strIndex.tolist()
        with tm.assert_produces_warning():
            self.strIndex.tolist() + self.strIndex

        with tm.assert_produces_warning(RuntimeWarning):
            firstCat = self.strIndex.union(self.dateIndex)
        secondCat = self.strIndex.union(self.strIndex)

        if self.dateIndex.dtype == np.object_:
            appended = np.append(self.strIndex, self.dateIndex)
        else:
            appended = np.append(self.strIndex, self.dateIndex.astype('O'))

        self.assertTrue(tm.equalContents(firstCat, appended))
        self.assertTrue(tm.equalContents(secondCat, self.strIndex))
        tm.assert_contains_all(self.strIndex, firstCat)
        tm.assert_contains_all(self.strIndex, secondCat)
        tm.assert_contains_all(self.dateIndex, firstCat)

        # test add and radd
        idx = Index(list('abc'))
        expected = Index(['a1', 'b1', 'c1'])
        self.assert_index_equal(idx + '1', expected)
        expected = Index(['1a', '1b', '1c'])
        self.assert_index_equal('1' + idx, expected)

    def test_append_multiple(self):
        index = Index(['a', 'b', 'c', 'd', 'e', 'f'])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        self.assertTrue(result.equals(index))

        # empty
        result = index.append([])
        self.assertTrue(result.equals(index))

    def test_append_empty_preserve_name(self):
        left = Index([], name='foo')
        right = Index([1, 2, 3], name='foo')

        result = left.append(right)
        self.assertEqual(result.name, 'foo')

        left = Index([], name='foo')
        right = Index([1, 2, 3], name='bar')

        result = left.append(right)
        self.assertIsNone(result.name)

    def test_add_string(self):
        # from bug report
        index = Index(['a', 'b', 'c'])
        index2 = index + 'foo'

        self.assertNotIn('a', index2)
        self.assertIn('afoo', index2)

    def test_iadd_string(self):
        index = pd.Index(['a', 'b', 'c'])
        # doesn't fail test unless there is a check before `+=`
        self.assertIn('a', index)

        index += '_x'
        self.assertIn('a_x', index)

    def test_difference(self):

        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]
        first.name = 'name'
        # different names
        result = first.difference(second)

        self.assertTrue(tm.equalContents(result, answer))
        self.assertEqual(result.name, None)

        # same names
        second.name = 'name'
        result = first.difference(second)
        self.assertEqual(result.name, 'name')

        # with empty
        result = first.difference([])
        self.assertTrue(tm.equalContents(result, first))
        self.assertEqual(result.name, first.name)

        # with everythin
        result = first.difference(first)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.name, first.name)

    def test_symmetric_diff(self):
        # smoke
        idx1 = Index([1, 2, 3, 4], name='idx1')
        idx2 = Index([2, 3, 4, 5])
        result = idx1.sym_diff(idx2)
        expected = Index([1, 5])
        self.assertTrue(tm.equalContents(result, expected))
        self.assertIsNone(result.name)

        # __xor__ syntax
        expected = idx1 ^ idx2
        self.assertTrue(tm.equalContents(result, expected))
        self.assertIsNone(result.name)

        # multiIndex
        idx1 = MultiIndex.from_tuples(self.tuples)
        idx2 = MultiIndex.from_tuples([('foo', 1), ('bar', 3)])
        result = idx1.sym_diff(idx2)
        expected = MultiIndex.from_tuples([('bar', 2), ('baz', 3), ('bar', 3)])
        self.assertTrue(tm.equalContents(result, expected))

        # nans:
        # GH #6444, sorting of nans. Make sure the number of nans is right
        # and the correct non-nan values are there. punt on sorting.
        idx1 = Index([1, 2, 3, np.nan])
        idx2 = Index([0, 1, np.nan])
        result = idx1.sym_diff(idx2)
        # expected = Index([0.0, np.nan, 2.0, 3.0, np.nan])

        nans = pd.isnull(result)
        self.assertEqual(nans.sum(), 1)
        self.assertEqual((~nans).sum(), 3)
        [self.assertIn(x, result) for x in [0.0, 2.0, 3.0]]

        # other not an Index:
        idx1 = Index([1, 2, 3, 4], name='idx1')
        idx2 = np.array([2, 3, 4, 5])
        expected = Index([1, 5])
        result = idx1.sym_diff(idx2)
        self.assertTrue(tm.equalContents(result, expected))
        self.assertEqual(result.name, 'idx1')

        result = idx1.sym_diff(idx2, result_name='new_name')
        self.assertTrue(tm.equalContents(result, expected))
        self.assertEqual(result.name, 'new_name')

    def test_is_numeric(self):
        self.assertFalse(self.dateIndex.is_numeric())
        self.assertFalse(self.strIndex.is_numeric())
        self.assertTrue(self.intIndex.is_numeric())
        self.assertTrue(self.floatIndex.is_numeric())
        self.assertFalse(self.catIndex.is_numeric())

    def test_is_object(self):
        self.assertTrue(self.strIndex.is_object())
        self.assertTrue(self.boolIndex.is_object())
        self.assertFalse(self.catIndex.is_object())
        self.assertFalse(self.intIndex.is_object())
        self.assertFalse(self.dateIndex.is_object())
        self.assertFalse(self.floatIndex.is_object())

    def test_is_all_dates(self):
        self.assertTrue(self.dateIndex.is_all_dates)
        self.assertFalse(self.strIndex.is_all_dates)
        self.assertFalse(self.intIndex.is_all_dates)

    def test_summary(self):
        self._check_method_works(Index.summary)
        # GH3869
        ind = Index(['{other}%s', "~:{range}:0"], name='A')
        result = ind.summary()
        # shouldn't be formatted accidentally.
        self.assertIn('~:{range}:0', result)
        self.assertIn('{other}%s', result)

    def test_format(self):
        self._check_method_works(Index.format)

        index = Index([datetime.now()])


        # windows has different precision on datetime.datetime.now (it doesn't include us
        # since the default for Timestamp shows these but Index formating does not
        # we are skipping
        if not is_platform_windows():
            formatted = index.format()
            expected = [str(index[0])]
            self.assertEqual(formatted, expected)

        # 2845
        index = Index([1, 2.0+3.0j, np.nan])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEqual(formatted, expected)

        # is this really allowed?
        index = Index([1, 2.0+3.0j, None])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEqual(formatted, expected)

        self.strIndex[:0].format()

    def test_format_with_name_time_info(self):
        # bug I fixed 12/20/2011
        inc = timedelta(hours=4)
        dates = Index([dt + inc for dt in self.dateIndex], name='something')

        formatted = dates.format(name=True)
        self.assertEqual(formatted[0], 'something')

    def test_format_datetime_with_time(self):
        t = Index([datetime(2012, 2, 7), datetime(2012, 2, 7, 23)])

        result = t.format()
        expected = ['2012-02-07 00:00:00', '2012-02-07 23:00:00']
        self.assertEqual(len(result), 2)
        self.assertEqual(result, expected)

    def test_format_none(self):
        values = ['a', 'b', 'c', None]

        idx = Index(values)
        idx.format()
        self.assertIsNone(idx[3])

    def test_logical_compat(self):
        idx = self.create_index()
        self.assertEqual(idx.all(), idx.values.all())
        self.assertEqual(idx.any(), idx.values.any())

    def _check_method_works(self, method):
        method(self.empty)
        method(self.dateIndex)
        method(self.unicodeIndex)
        method(self.strIndex)
        method(self.intIndex)
        method(self.tuples)
        method(self.catIndex)

    def test_get_indexer(self):
        idx1 = Index([1, 2, 3, 4, 5])
        idx2 = Index([2, 4, 6])

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, [1, 3, -1])

        r1 = idx2.get_indexer(idx1, method='pad')
        e1 = [-1, 0, 0, 1, 1]
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='pad')
        assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        e1 = [0, 0, 1, 1, 2]
        assert_almost_equal(r1, e1)

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        r2 = idx2.get_indexer(idx1[::-1], method='backfill')
        assert_almost_equal(r2, e1[::-1])

    def test_get_indexer_invalid(self):
        # GH10411
        idx = Index(np.arange(10))

        with tm.assertRaisesRegexp(ValueError, 'tolerance argument'):
            idx.get_indexer([1, 0], tolerance=1)

        with tm.assertRaisesRegexp(ValueError, 'limit argument'):
            idx.get_indexer([1, 0], limit=1)

    def test_get_indexer_nearest(self):
        idx = Index(np.arange(10))

        all_methods = ['pad', 'backfill', 'nearest']
        for method in all_methods:
            actual = idx.get_indexer([0, 5, 9], method=method)
            tm.assert_numpy_array_equal(actual, [0, 5, 9])

            actual = idx.get_indexer([0, 5, 9], method=method, tolerance=0)
            tm.assert_numpy_array_equal(actual, [0, 5, 9])

        for method, expected in zip(all_methods, [[0, 1, 8], [1, 2, 9], [0, 2, 9]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method)
            tm.assert_numpy_array_equal(actual, expected)

            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method, tolerance=1)
            tm.assert_numpy_array_equal(actual, expected)

        for method, expected in zip(all_methods, [[0, -1, -1], [-1, 2, -1], [0, 2, -1]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method, tolerance=0.2)
            tm.assert_numpy_array_equal(actual, expected)

        with tm.assertRaisesRegexp(ValueError, 'limit argument'):
            idx.get_indexer([1, 0], method='nearest', limit=1)

    def test_get_indexer_nearest_decreasing(self):
        idx = Index(np.arange(10))[::-1]

        all_methods = ['pad', 'backfill', 'nearest']
        for method in all_methods:
            actual = idx.get_indexer([0, 5, 9], method=method)
            tm.assert_numpy_array_equal(actual, [9, 4, 0])

        for method, expected in zip(all_methods, [[8, 7, 0], [9, 8, 1], [9, 7, 0]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method)
            tm.assert_numpy_array_equal(actual, expected)

    def test_get_indexer_strings(self):
        idx = pd.Index(['b', 'c'])

        actual = idx.get_indexer(['a', 'b', 'c', 'd'], method='pad')
        expected = [-1, 0, 1, 1]
        tm.assert_numpy_array_equal(actual, expected)

        actual = idx.get_indexer(['a', 'b', 'c', 'd'], method='backfill')
        expected = [0, 0, 1, -1]
        tm.assert_numpy_array_equal(actual, expected)

        with tm.assertRaises(TypeError):
            idx.get_indexer(['a', 'b', 'c', 'd'], method='nearest')

        with tm.assertRaises(TypeError):
            idx.get_indexer(['a', 'b', 'c', 'd'], method='pad', tolerance=2)

    def test_get_loc(self):
        idx = pd.Index([0, 1, 2])
        all_methods = [None, 'pad', 'backfill', 'nearest']
        for method in all_methods:
            self.assertEqual(idx.get_loc(1, method=method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(1, method=method, tolerance=0), 1)
            with tm.assertRaises(TypeError):
                idx.get_loc([1, 2], method=method)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method), loc)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method, tolerance=1), loc)

        for method in ['pad', 'backfill', 'nearest']:
            with tm.assertRaises(KeyError):
                idx.get_loc(1.1, method, tolerance=0.05)

        with tm.assertRaisesRegexp(ValueError, 'must be numeric'):
            idx.get_loc(1.1, 'nearest', tolerance='invalid')
        with tm.assertRaisesRegexp(ValueError, 'tolerance .* valid if'):
            idx.get_loc(1.1, tolerance=1)

        idx = pd.Index(['a', 'c'])
        with tm.assertRaises(TypeError):
            idx.get_loc('a', method='nearest')
        with tm.assertRaises(TypeError):
            idx.get_loc('a', method='pad', tolerance='invalid')

    def test_slice_locs(self):
        for dtype in [int, float]:
            idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=dtype))
            n = len(idx)

            self.assertEqual(idx.slice_locs(start=2), (2, n))
            self.assertEqual(idx.slice_locs(start=3), (3, n))
            self.assertEqual(idx.slice_locs(3, 8), (3, 6))
            self.assertEqual(idx.slice_locs(5, 10), (3, n))
            self.assertEqual(idx.slice_locs(end=8), (0, 6))
            self.assertEqual(idx.slice_locs(end=9), (0, 7))

            # reversed
            idx2 = idx[::-1]
            self.assertEqual(idx2.slice_locs(8, 2), (2, 6))
            self.assertEqual(idx2.slice_locs(7, 3), (2, 5))

        # float slicing
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=float))
        n = len(idx)
        self.assertEqual(idx.slice_locs(5.0, 10.0), (3, n))
        self.assertEqual(idx.slice_locs(4.5, 10.5), (3, 8))
        idx2 = idx[::-1]
        self.assertEqual(idx2.slice_locs(8.5, 1.5), (2, 6))
        self.assertEqual(idx2.slice_locs(10.5, -1), (0, n))

        # int slicing with floats
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=int))
        self.assertEqual(idx.slice_locs(5.0, 10.0), (3, n))
        self.assertEqual(idx.slice_locs(4.5, 10.5), (3, 8))
        idx2 = idx[::-1]
        self.assertEqual(idx2.slice_locs(8.5, 1.5), (2, 6))
        self.assertEqual(idx2.slice_locs(10.5, -1), (0, n))

    def test_slice_locs_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        self.assertEqual(idx.slice_locs('a', 'd'), (0, 6))
        self.assertEqual(idx.slice_locs(end='d'), (0, 6))
        self.assertEqual(idx.slice_locs('a', 'c'), (0, 4))
        self.assertEqual(idx.slice_locs('b', 'd'), (2, 6))

        idx2 = idx[::-1]
        self.assertEqual(idx2.slice_locs('d', 'a'), (0, 6))
        self.assertEqual(idx2.slice_locs(end='a'), (0, 6))
        self.assertEqual(idx2.slice_locs('d', 'b'), (0, 4))
        self.assertEqual(idx2.slice_locs('c', 'a'), (2, 6))

        for dtype in [int, float]:
            idx = Index(np.array([10, 12, 12, 14], dtype=dtype))
            self.assertEqual(idx.slice_locs(12, 12), (1, 3))
            self.assertEqual(idx.slice_locs(11, 13), (1, 3))

            idx2 = idx[::-1]
            self.assertEqual(idx2.slice_locs(12, 12), (1, 3))
            self.assertEqual(idx2.slice_locs(13, 11), (1, 3))

    def test_slice_locs_na(self):
        idx = Index([np.nan, 1, 2])
        self.assertRaises(KeyError, idx.slice_locs, start=1.5)
        self.assertRaises(KeyError, idx.slice_locs, end=1.5)
        self.assertEqual(idx.slice_locs(1), (1, 3))
        self.assertEqual(idx.slice_locs(np.nan), (0, 3))

        idx = Index([0, np.nan, np.nan, 1, 2])
        self.assertEqual(idx.slice_locs(np.nan), (1, 5))

    def test_slice_locs_negative_step(self):
        idx = Index(list('bcdxy'))

        SLC = pd.IndexSlice

        def check_slice(in_slice, expected):
            s_start, s_stop = idx.slice_locs(in_slice.start, in_slice.stop,
                                             in_slice.step)
            result = idx[s_start:s_stop:in_slice.step]
            expected = pd.Index(list(expected))
            self.assertTrue(result.equals(expected))

        for in_slice, expected in [
                (SLC[::-1], 'yxdcb'), (SLC['b':'y':-1], ''),
                (SLC['b'::-1], 'b'), (SLC[:'b':-1], 'yxdcb'),
                (SLC[:'y':-1], 'y'), (SLC['y'::-1], 'yxdcb'),
                (SLC['y'::-4], 'yb'),
                # absent labels
                (SLC[:'a':-1], 'yxdcb'), (SLC[:'a':-2], 'ydb'),
                (SLC['z'::-1], 'yxdcb'), (SLC['z'::-3], 'yc'),
                (SLC['m'::-1], 'dcb'), (SLC[:'m':-1], 'yx'),
                (SLC['a':'a':-1], ''), (SLC['z':'z':-1], ''),
                (SLC['m':'m':-1], '')
        ]:
            check_slice(in_slice, expected)

    def test_drop(self):
        n = len(self.strIndex)

        drop = self.strIndex[lrange(5, 10)]
        dropped = self.strIndex.drop(drop)
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        self.assertTrue(dropped.equals(expected))

        self.assertRaises(ValueError, self.strIndex.drop, ['foo', 'bar'])
        self.assertRaises(ValueError, self.strIndex.drop, ['1', 'bar'])

        # errors='ignore'
        mixed = drop.tolist() + ['foo']
        dropped = self.strIndex.drop(mixed, errors='ignore')
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        self.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(['foo', 'bar'], errors='ignore')
        expected = self.strIndex[lrange(n)]
        self.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(self.strIndex[0])
        expected = self.strIndex[1:]
        self.assert_index_equal(dropped, expected)

        ser = Index([1, 2, 3])
        dropped = ser.drop(1)
        expected = Index([2, 3])
        self.assert_index_equal(dropped, expected)

        # errors='ignore'
        self.assertRaises(ValueError, ser.drop, [3, 4])

        dropped = ser.drop(4, errors='ignore')
        expected = Index([1, 2, 3])
        self.assert_index_equal(dropped, expected)

        dropped = ser.drop([3, 4, 5], errors='ignore')
        expected = Index([1, 2])
        self.assert_index_equal(dropped, expected)

    def test_tuple_union_bug(self):
        import pandas
        import numpy as np

        aidx1 = np.array([(1, 'A'), (2, 'A'), (1, 'B'), (2, 'B')],
                         dtype=[('num', int), ('let', 'a1')])
        aidx2 = np.array([(1, 'A'), (2, 'A'), (1, 'B'), (2, 'B'), (1, 'C'), (2,
                         'C')], dtype=[('num', int), ('let', 'a1')])

        idx1 = pandas.Index(aidx1)
        idx2 = pandas.Index(aidx2)

        # intersection broken?
        int_idx = idx1.intersection(idx2)
        # needs to be 1d like idx1 and idx2
        expected = idx1[:4]  # pandas.Index(sorted(set(idx1) & set(idx2)))
        self.assertEqual(int_idx.ndim, 1)
        self.assertTrue(int_idx.equals(expected))

        # union broken
        union_idx = idx1.union(idx2)
        expected = idx2
        self.assertEqual(union_idx.ndim, 1)
        self.assertTrue(union_idx.equals(expected))

    def test_is_monotonic_incomparable(self):
        index = Index([5, datetime.now(), 7])
        self.assertFalse(index.is_monotonic)
        self.assertFalse(index.is_monotonic_decreasing)

    def test_get_set_value(self):
        values = np.random.randn(100)
        date = self.dateIndex[67]

        assert_almost_equal(self.dateIndex.get_value(values, date),
                            values[67])

        self.dateIndex.set_value(values, date, 10)
        self.assertEqual(values[67], 10)

    def test_isin(self):
        values = ['foo', 'bar', 'quux']

        idx = Index(['qux', 'baz', 'foo', 'bar'])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = Index([])
        result = idx.isin(values)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.dtype, np.bool_)

    def test_isin_nan(self):
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([np.nan]), [False, True])
        tm.assert_numpy_array_equal(
            Index(['a', pd.NaT]).isin([pd.NaT]), [False, True])
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([float('nan')]), [False, False])
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([pd.NaT]), [False, False])
        # Float64Index overrides isin, so must be checked separately
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([np.nan]), [False, True])
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([float('nan')]), [False, True])
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([pd.NaT]), [False, True])

    def test_isin_level_kwarg(self):
        def check_idx(idx):
            values = idx.tolist()[-2:] + ['nonexisting']

            expected = np.array([False, False, True, True])
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=0))
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=-1))

            self.assertRaises(IndexError, idx.isin, values, level=1)
            self.assertRaises(IndexError, idx.isin, values, level=10)
            self.assertRaises(IndexError, idx.isin, values, level=-2)

            self.assertRaises(KeyError, idx.isin, values, level=1.0)
            self.assertRaises(KeyError, idx.isin, values, level='foobar')

            idx.name = 'foobar'
            tm.assert_numpy_array_equal(expected,
                                          idx.isin(values, level='foobar'))

            self.assertRaises(KeyError, idx.isin, values, level='xyzzy')
            self.assertRaises(KeyError, idx.isin, values, level=np.nan)

        check_idx(Index(['qux', 'baz', 'foo', 'bar']))
        # Float64Index overrides isin, so must be checked separately
        check_idx(Float64Index([1.0, 2.0, 3.0, 4.0]))

    def test_boolean_cmp(self):
        values = [1, 2, 3, 4]

        idx = Index(values)
        res = (idx == values)

        tm.assert_numpy_array_equal(res,np.array([True,True,True,True],dtype=bool))

    def test_get_level_values(self):
        result = self.strIndex.get_level_values(0)
        self.assertTrue(result.equals(self.strIndex))

    def test_slice_keep_name(self):
        idx = Index(['a', 'b'], name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

    def test_join_self(self):
        # instance attributes of the form self.<name>Index
        indices = 'unicode', 'str', 'date', 'int', 'float'
        kinds = 'outer', 'inner', 'left', 'right'
        for index_kind in indices:
            res = getattr(self, '{0}Index'.format(index_kind))

            for kind in kinds:
                joined = res.join(res, how=kind)
                self.assertIs(res, joined)

    def test_str_attribute(self):
        # GH9068
        methods = ['strip', 'rstrip', 'lstrip']
        idx = Index([' jack', 'jill ', ' jesse ', 'frank'])
        for method in methods:
            expected = Index([getattr(str, method)(x) for x in idx.values])
            tm.assert_index_equal(getattr(Index.str, method)(idx.str), expected)

        # create a few instances that are not able to use .str accessor
        indices = [Index(range(5)),
                   tm.makeDateIndex(10),
                   MultiIndex.from_tuples([('foo', '1'), ('bar', '3')]),
                   PeriodIndex(start='2000', end='2010', freq='A')]
        for idx in indices:
            with self.assertRaisesRegexp(AttributeError, 'only use .str accessor'):
                idx.str.repeat(2)

        idx = Index(['a b c', 'd e', 'f'])
        expected = Index([['a', 'b', 'c'], ['d', 'e'], ['f']])
        tm.assert_index_equal(idx.str.split(), expected)
        tm.assert_index_equal(idx.str.split(expand=False), expected)

        expected = MultiIndex.from_tuples([('a', 'b', 'c'),
                                           ('d', 'e', np.nan),
                                           ('f', np.nan, np.nan)])
        tm.assert_index_equal(idx.str.split(expand=True), expected)

        # test boolean case, should return np.array instead of boolean Index
        idx = Index(['a1', 'a2', 'b1', 'b2'])
        expected = np.array([True, True, False, False])
        tm.assert_numpy_array_equal(idx.str.startswith('a'), expected)
        self.assertIsInstance(idx.str.startswith('a'), np.ndarray)
        s = Series(range(4), index=idx)
        expected = Series(range(2), index=['a1', 'a2'])
        tm.assert_series_equal(s[s.index.str.startswith('a')], expected)

    def test_tab_completion(self):
        # GH 9910
        idx = Index(list('abcd'))
        self.assertTrue('str' in dir(idx))

        idx = Index(range(4))
        self.assertTrue('str' not in dir(idx))

    def test_indexing_doesnt_change_class(self):
        idx = Index([1, 2, 3, 'a', 'b', 'c'])

        self.assertTrue(idx[1:3].identical(
            pd.Index([2, 3], dtype=np.object_)))
        self.assertTrue(idx[[0,1]].identical(
            pd.Index([1, 2], dtype=np.object_)))

    def test_outer_join_sort(self):
        left_idx = Index(np.random.permutation(15))
        right_idx = tm.makeDateIndex(10)

        with tm.assert_produces_warning(RuntimeWarning):
            joined = left_idx.join(right_idx, how='outer')

        # right_idx in this case because DatetimeIndex has join precedence over
        # Int64Index
        with tm.assert_produces_warning(RuntimeWarning):
            expected = right_idx.astype(object).union(left_idx.astype(object))
        tm.assert_index_equal(joined, expected)

    def test_nan_first_take_datetime(self):
        idx = Index([pd.NaT, Timestamp('20130101'), Timestamp('20130102')])
        res = idx.take([-1, 0, 1])
        exp = Index([idx[-1], idx[0], idx[1]])
        tm.assert_index_equal(res, exp)

    def test_reindex_preserves_name_if_target_is_list_or_ndarray(self):
        # GH6552
        idx = pd.Index([0, 1, 2])

        dt_idx = pd.date_range('20130101', periods=3)

        idx.name = None
        self.assertEqual(idx.reindex([])[0].name, None)
        self.assertEqual(idx.reindex(np.array([]))[0].name, None)
        self.assertEqual(idx.reindex(idx.tolist())[0].name, None)
        self.assertEqual(idx.reindex(idx.tolist()[:-1])[0].name, None)
        self.assertEqual(idx.reindex(idx.values)[0].name, None)
        self.assertEqual(idx.reindex(idx.values[:-1])[0].name, None)

        # Must preserve name even if dtype changes.
        self.assertEqual(idx.reindex(dt_idx.values)[0].name, None)
        self.assertEqual(idx.reindex(dt_idx.tolist())[0].name, None)

        idx.name = 'foobar'
        self.assertEqual(idx.reindex([])[0].name, 'foobar')
        self.assertEqual(idx.reindex(np.array([]))[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.tolist())[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.tolist()[:-1])[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.values)[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.values[:-1])[0].name, 'foobar')

        # Must preserve name even if dtype changes.
        self.assertEqual(idx.reindex(dt_idx.values)[0].name, 'foobar')
        self.assertEqual(idx.reindex(dt_idx.tolist())[0].name, 'foobar')

    def test_reindex_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.Index(list('abc'))
        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        self.assertEqual(get_reindex_type([]), np.object_)
        self.assertEqual(get_reindex_type(np.array([])), np.object_)
        self.assertEqual(get_reindex_type(np.array([], dtype=np.int64)),
                         np.object_)

    def test_reindex_doesnt_preserve_type_if_target_is_empty_index(self):
        # GH7774
        idx = pd.Index(list('abc'))
        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        self.assertEqual(get_reindex_type(pd.Int64Index([])), np.int64)
        self.assertEqual(get_reindex_type(pd.Float64Index([])), np.float64)
        self.assertEqual(get_reindex_type(pd.DatetimeIndex([])), np.datetime64)

        reindexed = idx.reindex(pd.MultiIndex([pd.Int64Index([]),
                                               pd.Float64Index([])],
                                              [[], []]))[0]
        self.assertEqual(reindexed.levels[0].dtype.type, np.int64)
        self.assertEqual(reindexed.levels[1].dtype.type, np.float64)

    def test_groupby(self):
        idx = Index(range(5))
        groups = idx.groupby(np.array([1,1,2,2,2]))
        exp = {1: [0, 1], 2: [2, 3, 4]}
        tm.assert_dict_equal(groups, exp)

    def test_equals_op_multiindex(self):
        # GH9785
        # test comparisons of multiindex
        from pandas.compat import StringIO
        df = pd.read_csv(StringIO('a,b,c\n1,2,3\n4,5,6'), index_col=[0, 1])
        tm.assert_numpy_array_equal(df.index == df.index, np.array([True, True]))

        mi1 = MultiIndex.from_tuples([(1, 2), (4, 5)])
        tm.assert_numpy_array_equal(df.index == mi1, np.array([True, True]))
        mi2 = MultiIndex.from_tuples([(1, 2), (4, 6)])
        tm.assert_numpy_array_equal(df.index == mi2, np.array([True, False]))
        mi3 = MultiIndex.from_tuples([(1, 2), (4, 5), (8, 9)])
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            df.index == mi3

        index_a = Index(['foo', 'bar', 'baz'])
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            df.index == index_a
        tm.assert_numpy_array_equal(index_a == mi3, np.array([False, False, False]))

    def test_conversion_preserves_name(self):
        #GH 10875
        i = pd.Index(['01:02:03', '01:02:04'], name='label')
        self.assertEqual(i.name, pd.to_datetime(i).name)
        self.assertEqual(i.name, pd.to_timedelta(i).name)

    def test_string_index_repr(self):
        # py3/py2 repr can differ because of "u" prefix
        # which also affects to displayed element size

        # short
        idx = pd.Index(['a', 'bb', 'ccc'])
        if PY3:
            expected = u"""Index(['a', 'bb', 'ccc'], dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'a', u'bb', u'ccc'], dtype='object')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.Index(['a', 'bb', 'ccc'] * 10)
        if PY3:
            expected = u"""Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object')"""
            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.Index(['a', 'bb', 'ccc'] * 100)
        if PY3:
            expected = u"""Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
       ...
       'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object', length=300)"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       ...
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object', length=300)"""
            self.assertEqual(unicode(idx), expected)

        # short
        idx = pd.Index([u'あ', u'いい', u'ううう'])
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう'], dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう'], dtype='object')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
      dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object')"""
            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
       ...
       'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
      dtype='object', length=300)"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       ...
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object', length=300)"""
            self.assertEqual(unicode(idx), expected)

        # Emable Unicode option -----------------------------------------
        with cf.option_context('display.unicode.east_asian_width', True):

            # short
            idx = pd.Index([u'あ', u'いい', u'ううう'])
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう'], dtype='object')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう'], dtype='object')"""
                self.assertEqual(unicode(idx), expected)

            # multiple lines
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう'],
      dtype='object')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object')"""
                self.assertEqual(unicode(idx), expected)

            # truncated
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ',
       ...
       'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
       'ううう'],
      dtype='object', length=300)"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ',
       ...
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう'],
      dtype='object', length=300)"""
                self.assertEqual(unicode(idx), expected)


class TestCategoricalIndex(Base, tm.TestCase):
    _holder = CategoricalIndex

    def setUp(self):
        self.indices = dict(catIndex = tm.makeCategoricalIndex(100))
        self.setup_indices()

    def create_index(self, categories=None, ordered=False):
        if categories is None:
            categories = list('cab')
        return CategoricalIndex(list('aabbca'), categories=categories, ordered=ordered)

    def test_construction(self):

        ci = self.create_index(categories=list('abcd'))
        categories = ci.categories

        result = Index(ci)
        tm.assert_index_equal(result,ci,exact=True)
        self.assertFalse(result.ordered)

        result = Index(ci.values)
        tm.assert_index_equal(result,ci,exact=True)
        self.assertFalse(result.ordered)

        # empty
        result = CategoricalIndex(categories=categories)
        self.assertTrue(result.categories.equals(Index(categories)))
        tm.assert_numpy_array_equal(result.codes, np.array([],dtype='int8'))
        self.assertFalse(result.ordered)

        # passing categories
        result = CategoricalIndex(list('aabbca'),categories=categories)
        self.assertTrue(result.categories.equals(Index(categories)))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,2,0],dtype='int8'))

        c = pd.Categorical(list('aabbca'))
        result = CategoricalIndex(c)
        self.assertTrue(result.categories.equals(Index(list('abc'))))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,2,0],dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(c,categories=categories)
        self.assertTrue(result.categories.equals(Index(categories)))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,2,0],dtype='int8'))
        self.assertFalse(result.ordered)

        ci = CategoricalIndex(c,categories=list('abcd'))
        result = CategoricalIndex(ci)
        self.assertTrue(result.categories.equals(Index(categories)))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,2,0],dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(ci, categories=list('ab'))
        self.assertTrue(result.categories.equals(Index(list('ab'))))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,-1,0],dtype='int8'))
        self.assertFalse(result.ordered)

        result = CategoricalIndex(ci, categories=list('ab'), ordered=True)
        self.assertTrue(result.categories.equals(Index(list('ab'))))
        tm.assert_numpy_array_equal(result.codes,np.array([0,0,1,1,-1,0],dtype='int8'))
        self.assertTrue(result.ordered)

        # turn me to an Index
        result = Index(np.array(ci))
        self.assertIsInstance(result, Index)
        self.assertNotIsInstance(result, CategoricalIndex)

    def test_construction_with_dtype(self):

        # specify dtype
        ci = self.create_index(categories=list('abc'))

        result = Index(np.array(ci), dtype='category')
        tm.assert_index_equal(result,ci,exact=True)

        result = Index(np.array(ci).tolist(), dtype='category')
        tm.assert_index_equal(result,ci,exact=True)

        # these are generally only equal when the categories are reordered
        ci = self.create_index()

        result = Index(np.array(ci), dtype='category').reorder_categories(ci.categories)
        tm.assert_index_equal(result,ci,exact=True)

        # make sure indexes are handled
        expected = CategoricalIndex([0,1,2], categories=[0,1,2], ordered=True)
        idx = Index(range(3))
        result = CategoricalIndex(idx, categories=idx, ordered=True)
        tm.assert_index_equal(result, expected, exact=True)

    def test_disallow_set_ops(self):

        # GH 10039
        # set ops (+/-) raise TypeError
        idx = pd.Index(pd.Categorical(['a', 'b']))

        self.assertRaises(TypeError, lambda : idx - idx)
        self.assertRaises(TypeError, lambda : idx + idx)
        self.assertRaises(TypeError, lambda : idx - ['a','b'])
        self.assertRaises(TypeError, lambda : idx + ['a','b'])
        self.assertRaises(TypeError, lambda : ['a','b'] - idx)
        self.assertRaises(TypeError, lambda : ['a','b'] + idx)

    def test_method_delegation(self):

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.set_categories(list('cab'))
        tm.assert_index_equal(result, CategoricalIndex(list('aabbca'), categories=list('cab')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.rename_categories(list('efg'))
        tm.assert_index_equal(result, CategoricalIndex(list('ffggef'), categories=list('efg')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.add_categories(['d'])
        tm.assert_index_equal(result, CategoricalIndex(list('aabbca'), categories=list('cabd')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cab'))
        result = ci.remove_categories(['c'])
        tm.assert_index_equal(result, CategoricalIndex(list('aabb') + [np.nan] + ['a'], categories=list('ab')))

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.as_unordered()
        tm.assert_index_equal(result, ci)

        ci = CategoricalIndex(list('aabbca'), categories=list('cabdef'))
        result = ci.as_ordered()
        tm.assert_index_equal(result, CategoricalIndex(list('aabbca'), categories=list('cabdef'), ordered=True))

        # invalid
        self.assertRaises(ValueError, lambda : ci.set_categories(list('cab'), inplace=True))

    def test_contains(self):

        ci = self.create_index(categories=list('cabdef'))

        self.assertTrue('a' in ci)
        self.assertTrue('z' not in ci)
        self.assertTrue('e' not in ci)
        self.assertTrue(np.nan not in ci)

        # assert codes NOT in index
        self.assertFalse(0 in ci)
        self.assertFalse(1 in ci)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            ci = CategoricalIndex(list('aabbca'), categories=list('cabdef') + [np.nan])
        self.assertFalse(np.nan in ci)

        ci = CategoricalIndex(list('aabbca') + [np.nan], categories=list('cabdef'))
        self.assertTrue(np.nan in ci)

    def test_min_max(self):

        ci = self.create_index(ordered=False)
        self.assertRaises(TypeError, lambda : ci.min())
        self.assertRaises(TypeError, lambda : ci.max())

        ci = self.create_index(ordered=True)

        self.assertEqual(ci.min(),'c')
        self.assertEqual(ci.max(),'b')

    def test_append(self):

        ci = self.create_index()
        categories = ci.categories

        # append cats with the same categories
        result = ci[:3].append(ci[3:])
        tm.assert_index_equal(result,ci,exact=True)

        foos = [ci[:1], ci[1:3], ci[3:]]
        result = foos[0].append(foos[1:])
        tm.assert_index_equal(result,ci,exact=True)

        # empty
        result = ci.append([])
        tm.assert_index_equal(result,ci,exact=True)

        # appending with different categories or reoreded is not ok
        self.assertRaises(TypeError, lambda : ci.append(ci.values.set_categories(list('abcd'))))
        self.assertRaises(TypeError, lambda : ci.append(ci.values.reorder_categories(list('abc'))))

        # with objects
        result = ci.append(['c','a'])
        expected = CategoricalIndex(list('aabbcaca'), categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        # invalid objects
        self.assertRaises(TypeError, lambda : ci.append(['a','d']))

    def test_insert(self):

        ci = self.create_index()
        categories = ci.categories

        #test 0th element
        result = ci.insert(0, 'a')
        expected = CategoricalIndex(list('aaabbca'),categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        #test Nth element that follows Python list behavior
        result = ci.insert(-1, 'a')
        expected = CategoricalIndex(list('aabbcaa'),categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        #test empty
        result = CategoricalIndex(categories=categories).insert(0, 'a')
        expected = CategoricalIndex(['a'],categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        # invalid
        self.assertRaises(TypeError, lambda : ci.insert(0,'d'))

    def test_delete(self):

        ci = self.create_index()
        categories = ci.categories

        result = ci.delete(0)
        expected = CategoricalIndex(list('abbca'),categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        result = ci.delete(-1)
        expected = CategoricalIndex(list('aabbc'),categories=categories)
        tm.assert_index_equal(result,expected,exact=True)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = ci.delete(10)

    def test_astype(self):

        ci = self.create_index()
        result = ci.astype('category')
        tm.assert_index_equal(result,ci,exact=True)

        result = ci.astype(object)
        self.assertTrue(result.equals(Index(np.array(ci))))

        # this IS equal, but not the same class
        self.assertTrue(result.equals(ci))
        self.assertIsInstance(result, Index)
        self.assertNotIsInstance(result, CategoricalIndex)

    def test_reindex_base(self):

        # determined by cat ordering
        idx = self.create_index()
        expected = np.array([4,0,1,5,2,3])

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with tm.assertRaisesRegexp(ValueError, 'Invalid fill method'):
            idx.get_indexer(idx, method='invalid')

    def test_reindexing(self):

        ci = self.create_index()
        oidx = Index(np.array(ci))

        for n in [1,2,5,len(ci)]:
            finder = oidx[np.random.randint(0,len(ci),size=n)]
            expected = oidx.get_indexer_non_unique(finder)[0]

            actual = ci.get_indexer(finder)
            tm.assert_numpy_array_equal(expected, actual)

    def test_duplicates(self):

        idx = CategoricalIndex([0, 0, 0], name='foo')
        self.assertFalse(idx.is_unique)
        self.assertTrue(idx.has_duplicates)

        expected = CategoricalIndex([0], name='foo')
        self.assert_index_equal(idx.drop_duplicates(), expected)

    def test_get_indexer(self):

        idx1 = CategoricalIndex(list('aabcde'),categories=list('edabc'))
        idx2 = CategoricalIndex(list('abf'))

        for indexer in [idx2, list('abf'), Index(list('abf'))]:
            r1 = idx1.get_indexer(idx2)
            assert_almost_equal(r1, [0, 1, 2, -1])

        self.assertRaises(NotImplementedError, lambda : idx2.get_indexer(idx1, method='pad'))
        self.assertRaises(NotImplementedError, lambda : idx2.get_indexer(idx1, method='backfill'))
        self.assertRaises(NotImplementedError, lambda : idx2.get_indexer(idx1, method='nearest'))

    def test_repr_roundtrip(self):

        ci = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        str(ci)
        tm.assert_index_equal(eval(repr(ci)),ci,exact=True)

        # formatting
        if compat.PY3:
            str(ci)
        else:
            compat.text_type(ci)

        # long format
        # this is not reprable
        ci = CategoricalIndex(np.random.randint(0,5,size=100))
        if compat.PY3:
            str(ci)
        else:
            compat.text_type(ci)

    def test_isin(self):

        ci = CategoricalIndex(list('aabca') + [np.nan],categories=['c','a','b'])
        tm.assert_numpy_array_equal(ci.isin(['c']),np.array([False,False,False,True,False,False]))
        tm.assert_numpy_array_equal(ci.isin(['c','a','b']),np.array([True]*5 + [False]))
        tm.assert_numpy_array_equal(ci.isin(['c','a','b',np.nan]),np.array([True]*6))

        # mismatched categorical -> coerced to ndarray so doesn't matter
        tm.assert_numpy_array_equal(ci.isin(ci.set_categories(list('abcdefghi'))),np.array([True]*6))
        tm.assert_numpy_array_equal(ci.isin(ci.set_categories(list('defghi'))),np.array([False]*5 + [True]))

    def test_identical(self):

        ci1 = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        ci2 = CategoricalIndex(['a', 'b'], categories=['a', 'b', 'c'], ordered=True)
        self.assertTrue(ci1.identical(ci1))
        self.assertTrue(ci1.identical(ci1.copy()))
        self.assertFalse(ci1.identical(ci2))

    def test_equals(self):

        ci1 = CategoricalIndex(['a', 'b'], categories=['a', 'b'], ordered=True)
        ci2 = CategoricalIndex(['a', 'b'], categories=['a', 'b', 'c'], ordered=True)

        self.assertTrue(ci1.equals(ci1))
        self.assertFalse(ci1.equals(ci2))
        self.assertTrue(ci1.equals(ci1.astype(object)))
        self.assertTrue(ci1.astype(object).equals(ci1))

        self.assertTrue((ci1 == ci1).all())
        self.assertFalse((ci1 != ci1).all())
        self.assertFalse((ci1 > ci1).all())
        self.assertFalse((ci1 < ci1).all())
        self.assertTrue((ci1 <= ci1).all())
        self.assertTrue((ci1 >= ci1).all())

        self.assertFalse((ci1 == 1).all())
        self.assertTrue((ci1 == Index(['a','b'])).all())
        self.assertTrue((ci1 == ci1.values).all())

        # invalid comparisons
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            ci1 == Index(['a','b','c'])
        self.assertRaises(TypeError, lambda : ci1 == ci2)
        self.assertRaises(TypeError, lambda : ci1 == Categorical(ci1.values, ordered=False))
        self.assertRaises(TypeError, lambda : ci1 == Categorical(ci1.values, categories=list('abc')))

        # tests
        # make sure that we are testing for category inclusion properly
        self.assertTrue(CategoricalIndex(list('aabca'),categories=['c','a','b']).equals(list('aabca')))
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertTrue(CategoricalIndex(list('aabca'),categories=['c','a','b',np.nan]).equals(list('aabca')))

        self.assertFalse(CategoricalIndex(list('aabca') + [np.nan],categories=['c','a','b']).equals(list('aabca')))
        self.assertTrue(CategoricalIndex(list('aabca') + [np.nan],categories=['c','a','b']).equals(list('aabca') + [np.nan]))

    def test_string_categorical_index_repr(self):
        # short
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'])
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc'], categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc'], categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'] * 10)
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb',
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
                  u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
                  u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc',
                  u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
                 categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.CategoricalIndex(['a', 'bb', 'ccc'] * 100)
        if PY3:
            expected = u"""CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  ...
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category', length=300)"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
                  u'ccc', u'a',
                  ...
                  u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
                  u'bb', u'ccc'],
                 categories=[u'a', u'bb', u'ccc'], ordered=False, dtype='category', length=300)"""
            self.assertEqual(unicode(idx), expected)

        # larger categories
        idx = pd.CategoricalIndex(list('abcdefghijklmmo'))
        if PY3:
            expected = u"""CategoricalIndex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                  'm', 'm', 'o'],
                 categories=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', ...], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'a', u'b', u'c', u'd', u'e', u'f', u'g', u'h', u'i', u'j',
                  u'k', u'l', u'm', u'm', u'o'],
                 categories=[u'a', u'b', u'c', u'd', u'e', u'f', u'g', u'h', ...], ordered=False, dtype='category')"""

            self.assertEqual(unicode(idx), expected)

        # short
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'])
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう'], categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # multiple lines
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 10)
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう',
                  u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # truncated
        idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 100)
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ',
                  ...
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category', length=300)"""
            self.assertEqual(unicode(idx), expected)

        # larger categories
        idx = pd.CategoricalIndex(list(u'あいうえおかきくけこさしすせそ'))
        if PY3:
            expected = u"""CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ', 'さ', 'し',
                  'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""CategoricalIndex([u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', u'け', u'こ',
                  u'さ', u'し', u'す', u'せ', u'そ'],
                 categories=[u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', ...], ordered=False, dtype='category')"""
            self.assertEqual(unicode(idx), expected)

        # Emable Unicode option -----------------------------------------
        with cf.option_context('display.unicode.east_asian_width', True):

            # short
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'])
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう'], categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(unicode(idx), expected)

            # multiple lines
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 10)
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category')"""
                self.assertEqual(unicode(idx), expected)

            # truncated
            idx = pd.CategoricalIndex([u'あ', u'いい', u'ううう'] * 100)
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
                  u'いい', u'ううう', u'あ',
                  ...
                  u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
                  u'ううう', u'あ', u'いい', u'ううう'],
                 categories=[u'あ', u'いい', u'ううう'], ordered=False, dtype='category', length=300)"""
                self.assertEqual(unicode(idx), expected)

            # larger categories
            idx = pd.CategoricalIndex(list(u'あいうえおかきくけこさしすせそ'))
            if PY3:
                expected = u"""CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ',
                  'さ', 'し', 'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""CategoricalIndex([u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く',
                  u'け', u'こ', u'さ', u'し', u'す', u'せ', u'そ'],
                 categories=[u'あ', u'い', u'う', u'え', u'お', u'か', u'き', u'く', ...], ordered=False, dtype='category')"""
                self.assertEqual(unicode(idx), expected)

    def test_fillna_categorical(self):
        # GH 11343
        idx = CategoricalIndex([1.0, np.nan, 3.0, 1.0], name='x')
        # fill by value in categories
        exp = CategoricalIndex([1.0, 1.0, 3.0, 1.0], name='x')
        self.assert_index_equal(idx.fillna(1.0), exp)

        # fill by value not in categories raises ValueError
        with tm.assertRaisesRegexp(ValueError, 'fill value must be in categories'):
            idx.fillna(2.0)


class Numeric(Base):

    def test_numeric_compat(self):

        idx = self._holder(np.arange(5,dtype='int64'))
        didx = self._holder(np.arange(5,dtype='int64')**2
                            )
        result = idx * 1
        tm.assert_index_equal(result, idx)

        result = 1 * idx
        tm.assert_index_equal(result, idx)

        result = idx * idx
        tm.assert_index_equal(result, didx)

        result = idx / 1
        tm.assert_index_equal(result, idx)

        result = idx // 1
        tm.assert_index_equal(result, idx)

        result = idx * np.array(5,dtype='int64')
        tm.assert_index_equal(result, self._holder(np.arange(5,dtype='int64')*5))

        result = idx * np.arange(5,dtype='int64')
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5,dtype='int64'))
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5,dtype='float64')+0.1)
        tm.assert_index_equal(result,
                              Float64Index(np.arange(5,dtype='float64')*(np.arange(5,dtype='float64')+0.1)))

        # invalid
        self.assertRaises(TypeError, lambda : idx * date_range('20130101',periods=5))
        self.assertRaises(ValueError, lambda : idx * self._holder(np.arange(3)))
        self.assertRaises(ValueError, lambda : idx * np.array([1,2]))


    def test_explicit_conversions(self):

        # GH 8608
        # add/sub are overriden explicity for Float/Int Index
        idx = self._holder(np.arange(5,dtype='int64'))

        # float conversions
        arr = np.arange(5,dtype='int64')*3.2
        expected = Float64Index(arr)
        fidx = idx * 3.2
        tm.assert_index_equal(fidx,expected)
        fidx = 3.2 * idx
        tm.assert_index_equal(fidx,expected)

        # interops with numpy arrays
        expected = Float64Index(arr)
        a = np.zeros(5,dtype='float64')
        result = fidx - a
        tm.assert_index_equal(result,expected)

        expected = Float64Index(-arr)
        a = np.zeros(5,dtype='float64')
        result = a - fidx
        tm.assert_index_equal(result,expected)

    def test_ufunc_compat(self):
        idx = self._holder(np.arange(5,dtype='int64'))
        result = np.sin(idx)
        expected = Float64Index(np.sin(np.arange(5,dtype='int64')))
        tm.assert_index_equal(result, expected)

    def test_index_groupby(self):
        int_idx = Index(range(6))
        float_idx = Index(np.arange(0, 0.6, 0.1))
        obj_idx = Index('A B C D E F'.split())
        dt_idx = pd.date_range('2013-01-01', freq='M', periods=6)

        for idx in [int_idx, float_idx, obj_idx, dt_idx]:
            to_groupby = np.array([1, 2, np.nan, np.nan, 2, 1])
            self.assertEqual(idx.groupby(to_groupby),
                             {1.0: [idx[0], idx[5]], 2.0: [idx[1], idx[4]]})

            to_groupby = Index([datetime(2011, 11, 1), datetime(2011, 12, 1),
                                pd.NaT, pd.NaT,
                                datetime(2011, 12, 1), datetime(2011, 11, 1)], tz='UTC').values

            ex_keys = pd.tslib.datetime_to_datetime64(np.array([Timestamp('2011-11-01'), Timestamp('2011-12-01')]))
            expected = {ex_keys[0][0]: [idx[0], idx[5]], ex_keys[0][1]: [idx[1], idx[4]]}
            self.assertEqual(idx.groupby(to_groupby), expected)


class TestFloat64Index(Numeric, tm.TestCase):
    _holder = Float64Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(mixed = Float64Index([1.5, 2, 3, 4, 5]),
                            float = Float64Index(np.arange(5) * 2.5))
        self.setup_indices()

    def create_index(self):
        return Float64Index(np.arange(5, dtype='float64'))

    def test_repr_roundtrip(self):
        for ind in (self.mixed, self.float):
            tm.assert_index_equal(eval(repr(ind)), ind)

    def check_is_index(self, i):
        self.assertIsInstance(i, Index)
        self.assertNotIsInstance(i, Float64Index)

    def check_coerce(self, a, b, is_float_index=True):
        self.assertTrue(a.equals(b))
        if is_float_index:
            self.assertIsInstance(b, Float64Index)
        else:
            self.check_is_index(b)

    def test_constructor(self):

        # explicit construction
        index = Float64Index([1,2,3,4,5])
        self.assertIsInstance(index, Float64Index)
        self.assertTrue((index.values == np.array([1,2,3,4,5],dtype='float64')).all())
        index = Float64Index(np.array([1,2,3,4,5]))
        self.assertIsInstance(index, Float64Index)
        index = Float64Index([1.,2,3,4,5])
        self.assertIsInstance(index, Float64Index)
        index = Float64Index(np.array([1.,2,3,4,5]))
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, float)

        index = Float64Index(np.array([1.,2,3,4,5]),dtype=np.float32)
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, np.float64)

        index = Float64Index(np.array([1,2,3,4,5]),dtype=np.float32)
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, np.float64)

        # nan handling
        result = Float64Index([np.nan, np.nan])
        self.assertTrue(pd.isnull(result.values).all())
        result = Float64Index(np.array([np.nan]))
        self.assertTrue(pd.isnull(result.values).all())
        result = Index(np.array([np.nan]))
        self.assertTrue(pd.isnull(result.values).all())

    def test_constructor_invalid(self):

        # invalid
        self.assertRaises(TypeError, Float64Index, 0.)
        self.assertRaises(TypeError, Float64Index, ['a','b',0.])
        self.assertRaises(TypeError, Float64Index, [Timestamp('20130101')])

    def test_constructor_coerce(self):

        self.check_coerce(self.mixed,Index([1.5, 2, 3, 4, 5]))
        self.check_coerce(self.float,Index(np.arange(5) * 2.5))
        self.check_coerce(self.float,Index(np.array(np.arange(5) * 2.5, dtype=object)))

    def test_constructor_explicit(self):

        # these don't auto convert
        self.check_coerce(self.float,Index((np.arange(5) * 2.5), dtype=object),
                          is_float_index=False)
        self.check_coerce(self.mixed,Index([1.5, 2, 3, 4, 5],dtype=object),
                          is_float_index=False)

    def test_astype(self):

        result = self.float.astype(object)
        self.assertTrue(result.equals(self.float))
        self.assertTrue(self.float.equals(result))
        self.check_is_index(result)

        i = self.mixed.copy()
        i.name = 'foo'
        result = i.astype(object)
        self.assertTrue(result.equals(i))
        self.assertTrue(i.equals(result))
        self.check_is_index(result)

    def test_equals(self):

        i = Float64Index([1.0,2.0])
        self.assertTrue(i.equals(i))
        self.assertTrue(i.identical(i))

        i2 = Float64Index([1.0,2.0])
        self.assertTrue(i.equals(i2))

        i = Float64Index([1.0,np.nan])
        self.assertTrue(i.equals(i))
        self.assertTrue(i.identical(i))

        i2 = Float64Index([1.0,np.nan])
        self.assertTrue(i.equals(i2))

    def test_get_indexer(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = [-0.1, 0.5, 1.1]
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'), [0, 1, 1])

    def test_get_loc(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(1, method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(1, method, tolerance=0), 1)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method), loc)
            self.assertEqual(idx.get_loc(1.1, method, tolerance=0.9), loc)

        self.assertRaises(KeyError, idx.get_loc, 'foo')
        self.assertRaises(KeyError, idx.get_loc, 1.5)
        self.assertRaises(KeyError, idx.get_loc, 1.5,
                          method='pad', tolerance=0.1)

        with tm.assertRaisesRegexp(ValueError, 'must be numeric'):
            idx.get_loc(1.4, method='nearest', tolerance='foo')

    def test_get_loc_na(self):
        idx = Float64Index([np.nan, 1, 2])
        self.assertEqual(idx.get_loc(1), 1)
        self.assertEqual(idx.get_loc(np.nan), 0)

        idx = Float64Index([np.nan, 1, np.nan])
        self.assertEqual(idx.get_loc(1), 1)

        # representable by slice [0:2:2]
        # self.assertRaises(KeyError, idx.slice_locs, np.nan)
        sliced = idx.slice_locs(np.nan)
        self.assertTrue(isinstance(sliced, tuple))
        self.assertEqual(sliced, (0, 3))

        # not representable by slice
        idx = Float64Index([np.nan, 1, np.nan, np.nan])
        self.assertEqual(idx.get_loc(1), 1)
        self.assertRaises(KeyError, idx.slice_locs, np.nan)

    def test_contains_nans(self):
        i = Float64Index([1.0, 2.0, np.nan])
        self.assertTrue(np.nan in i)

    def test_contains_not_nans(self):
        i = Float64Index([1.0, 2.0, np.nan])
        self.assertTrue(1.0 in i)

    def test_doesnt_contain_all_the_things(self):
        i = Float64Index([np.nan])
        self.assertFalse(i.isin([0]).item())
        self.assertFalse(i.isin([1]).item())
        self.assertTrue(i.isin([np.nan]).item())

    def test_nan_multiple_containment(self):
        i = Float64Index([1.0, np.nan])
        tm.assert_numpy_array_equal(i.isin([1.0]), np.array([True, False]))
        tm.assert_numpy_array_equal(i.isin([2.0, np.pi]),
                                    np.array([False, False]))
        tm.assert_numpy_array_equal(i.isin([np.nan]),
                                    np.array([False, True]))
        tm.assert_numpy_array_equal(i.isin([1.0, np.nan]),
                                    np.array([True, True]))
        i = Float64Index([1.0, 2.0])
        tm.assert_numpy_array_equal(i.isin([np.nan]),
                                    np.array([False, False]))

    def test_astype_from_object(self):
        index = Index([1.0, np.nan, 0.2], dtype='object')
        result = index.astype(float)
        expected = Float64Index([1.0, np.nan, 0.2])
        tm.assert_equal(result.dtype, expected.dtype)
        tm.assert_index_equal(result, expected)

    def test_fillna_float64(self):
        # GH 11343
        idx = Index([1.0, np.nan, 3.0], dtype=float, name='x')
        # can't downcast
        exp = Index([1.0, 0.1, 3.0], name='x')
        self.assert_index_equal(idx.fillna(0.1), exp)

        # downcast
        exp = Int64Index([1, 2, 3], name='x')
        self.assert_index_equal(idx.fillna(2), exp)

        # object
        exp = Index([1, 'obj', 3], name='x')
        self.assert_index_equal(idx.fillna('obj'), exp)


class TestInt64Index(Numeric, tm.TestCase):
    _holder = Int64Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index = Int64Index(np.arange(0, 20, 2)))
        self.setup_indices()

    def create_index(self):
        return Int64Index(np.arange(5, dtype='int64'))

    def test_too_many_names(self):
        def testit():
            self.index.names = ["roger", "harold"]
        assertRaisesRegexp(ValueError, "^Length", testit)

    def test_constructor(self):
        # pass list, coerce fine
        index = Int64Index([-5, 0, 1, 2])
        expected = np.array([-5, 0, 1, 2], dtype=np.int64)
        tm.assert_numpy_array_equal(index, expected)

        # from iterable
        index = Int64Index(iter([-5, 0, 1, 2]))
        tm.assert_numpy_array_equal(index, expected)

        # scalar raise Exception
        self.assertRaises(TypeError, Int64Index, 5)

        # copy
        arr = self.index.values
        new_index = Int64Index(arr, copy=True)
        tm.assert_numpy_array_equal(new_index, self.index)
        val = arr[0] + 3000
        # this should not change index
        arr[0] = val
        self.assertNotEqual(new_index[0], val)

    def test_constructor_corner(self):
        arr = np.array([1, 2, 3, 4], dtype=object)
        index = Int64Index(arr)
        self.assertEqual(index.values.dtype, np.int64)
        self.assertTrue(index.equals(arr))

        # preventing casting
        arr = np.array([1, '2', 3, '4'], dtype=object)
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr)

        arr_with_floats = [0, 2, 3, 4, 5, 1.25, 3, -1]
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr_with_floats)

    def test_copy(self):
        i = Int64Index([], name='Foo')
        i_copy = i.copy()
        self.assertEqual(i_copy.name, 'Foo')

    def test_view(self):
        super(TestInt64Index, self).test_view()

        i = Int64Index([], name='Foo')
        i_view = i.view()
        self.assertEqual(i_view.name, 'Foo')

        i_view = i.view('i8')
        tm.assert_index_equal(i, Int64Index(i_view, name='Foo'))

        i_view = i.view(Int64Index)
        tm.assert_index_equal(i, Int64Index(i_view, name='Foo'))

    def test_coerce_list(self):
        # coerce things
        arr = Index([1, 2, 3, 4])
        tm.assertIsInstance(arr, Int64Index)

        # but not if explicit dtype passed
        arr = Index([1, 2, 3, 4], dtype=object)
        tm.assertIsInstance(arr, Index)

    def test_dtype(self):
        self.assertEqual(self.index.dtype, np.int64)

    def test_is_monotonic(self):
        self.assertTrue(self.index.is_monotonic)
        self.assertTrue(self.index.is_monotonic_increasing)
        self.assertFalse(self.index.is_monotonic_decreasing)

        index = Int64Index([4, 3, 2, 1])
        self.assertFalse(index.is_monotonic)
        self.assertTrue(index.is_monotonic_decreasing)

        index = Int64Index([1])
        self.assertTrue(index.is_monotonic)
        self.assertTrue(index.is_monotonic_increasing)
        self.assertTrue(index.is_monotonic_decreasing)

    def test_is_monotonic_na(self):
        examples = [Index([np.nan]),
                    Index([np.nan, 1]),
                    Index([1, 2, np.nan]),
                    Index(['a', 'b', np.nan]),
                    pd.to_datetime(['NaT']),
                    pd.to_datetime(['NaT', '2000-01-01']),
                    pd.to_datetime(['2000-01-01', 'NaT', '2000-01-02']),
                    pd.to_timedelta(['1 day', 'NaT']),
                   ]
        for index in examples:
            self.assertFalse(index.is_monotonic_increasing)
            self.assertFalse(index.is_monotonic_decreasing)

    def test_equals(self):
        same_values = Index(self.index, dtype=object)
        self.assertTrue(self.index.equals(same_values))
        self.assertTrue(same_values.equals(self.index))

    def test_logical_compat(self):
        idx = self.create_index()
        self.assertEqual(idx.all(), idx.values.all())
        self.assertEqual(idx.any(), idx.values.any())

    def test_identical(self):
        i = Index(self.index.copy())
        self.assertTrue(i.identical(self.index))

        same_values_different_type = Index(i, dtype=object)
        self.assertFalse(i.identical(same_values_different_type))

        i = self.index.copy(dtype=object)
        i = i.rename('foo')
        same_values = Index(i, dtype=object)
        self.assertTrue(same_values.identical(self.index.copy(dtype=object)))

        self.assertFalse(i.identical(self.index))
        self.assertTrue(Index(same_values, name='foo', dtype=object
                              ).identical(i))

        self.assertFalse(
            self.index.copy(dtype=object)
            .identical(self.index.copy(dtype='int64')))

    def test_get_indexer(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target)
        expected = np.array([0, -1, 1, -1, 2, -1, 3, -1, 4, -1])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_pad(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='pad')
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_backfill(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='backfill')
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_join_outer(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        # guarantee of sortedness
        res, lidx, ridx = self.index.join(other, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other, how='outer')
        self.assertTrue(res.equals(noidx_res))

        eres = Int64Index([0, 1, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 25])
        elidx = np.array([0, -1, 1, 2, -1, 3, -1, 4, 5, 6, 7, 8, 9, -1],
                         dtype=np.int64)
        eridx = np.array([-1, 3, 4, -1, 5, -1, 0, -1, -1, 1, -1, -1, -1, 2],
                         dtype=np.int64)

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other_mono, how='outer')
        self.assertTrue(res.equals(noidx_res))

        eridx = np.array([-1, 0, 1, -1, 2, -1, 3, -1, -1, 4, -1, -1, -1, 5],
                         dtype=np.int64)
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_inner(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='inner',
                                          return_indexers=True)

        # no guarantee of sortedness, so sort for comparison purposes
        ind = res.argsort()
        res = res.take(ind)
        lidx = lidx.take(ind)
        ridx = ridx.take(ind)

        eres = Int64Index([2, 12])
        elidx = np.array([1, 6])
        eridx = np.array([4, 1])

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='inner',
                                          return_indexers=True)

        res2 = self.index.intersection(other_mono)
        self.assertTrue(res.equals(res2))

        eridx = np.array([1, 4])
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_left(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='left',
                                          return_indexers=True)
        eres = self.index
        eridx = np.array([-1, 4, -1, -1, -1, -1, 1, -1, -1, -1],
                         dtype=np.int64)

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        self.assertIsNone(lidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='left',
                                          return_indexers=True)
        eridx = np.array([-1, 1, -1, -1, -1, -1, 4, -1, -1, -1],
                         dtype=np.int64)
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        self.assertIsNone(lidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # non-unique
        idx = Index([1, 1, 2, 5])
        idx2 = Index([1, 2, 5, 7, 9])
        res, lidx, ridx = idx2.join(idx, how='left', return_indexers=True)
        eres = Index([1, 1, 2, 5, 7, 9])  # 1 is in idx2, so it should be x2
        eridx = np.array([0, 1, 2, 3, -1, -1])
        elidx = np.array([0, 0, 1, 2, 3, 4])
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_right(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='right',
                                          return_indexers=True)
        eres = other
        elidx = np.array([-1, 6, -1, -1, 1, -1],
                         dtype=np.int64)

        tm.assertIsInstance(other, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        self.assertIsNone(ridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='right',
                                          return_indexers=True)
        eres = other_mono
        elidx = np.array([-1, 1, -1, -1, 6, -1],
                         dtype=np.int64)
        tm.assertIsInstance(other, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        self.assertIsNone(ridx)

        # non-unique
        idx = Index([1, 1, 2, 5])
        idx2 = Index([1, 2, 5, 7, 9])
        res, lidx, ridx = idx.join(idx2, how='right', return_indexers=True)
        eres = Index([1, 1, 2, 5, 7, 9])  # 1 is in idx2, so it should be x2
        elidx = np.array([0, 1, 2, 3, -1, -1])
        eridx = np.array([0, 0, 1, 2, 3, 4])
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_non_int_index(self):
        other = Index([3, 6, 7, 8, 10], dtype=object)

        outer = self.index.join(other, how='outer')
        outer2 = other.join(self.index, how='outer')
        expected = Index([0, 2, 3, 4, 6, 7, 8, 10, 12, 14,
                          16, 18], dtype=object)
        self.assertTrue(outer.equals(outer2))
        self.assertTrue(outer.equals(expected))

        inner = self.index.join(other, how='inner')
        inner2 = other.join(self.index, how='inner')
        expected = Index([6, 8, 10], dtype=object)
        self.assertTrue(inner.equals(inner2))
        self.assertTrue(inner.equals(expected))

        left = self.index.join(other, how='left')
        self.assertTrue(left.equals(self.index))

        left2 = other.join(self.index, how='left')
        self.assertTrue(left2.equals(other))

        right = self.index.join(other, how='right')
        self.assertTrue(right.equals(other))

        right2 = other.join(self.index, how='right')
        self.assertTrue(right2.equals(self.index))

    def test_join_non_unique(self):
        left = Index([4, 4, 3, 3])

        joined, lidx, ridx = left.join(left, return_indexers=True)

        exp_joined = Index([3, 3, 3, 3, 4, 4, 4, 4])
        self.assertTrue(joined.equals(exp_joined))

        exp_lidx = np.array([2, 2, 3, 3, 0, 0, 1, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(lidx, exp_lidx)

        exp_ridx = np.array([2, 3, 2, 3, 0, 1, 0, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(ridx, exp_ridx)

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = self.index.join(self.index, how=kind)
            self.assertIs(self.index, joined)

    def test_intersection(self):
        other = Index([1, 2, 3, 4, 5])
        result = self.index.intersection(other)
        expected = np.sort(np.intersect1d(self.index.values, other.values))
        tm.assert_numpy_array_equal(result, expected)

        result = other.intersection(self.index)
        expected = np.sort(np.asarray(np.intersect1d(self.index.values,
                                                     other.values)))
        tm.assert_numpy_array_equal(result, expected)

    def test_intersect_str_dates(self):
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        i1 = Index(dt_dates, dtype=object)
        i2 = Index(['aa'], dtype=object)
        res = i2.intersection(i1)

        self.assertEqual(len(res), 0)

    def test_union_noncomparable(self):
        from datetime import datetime, timedelta
        # corner case, non-Int64Index
        now = datetime.now()
        other = Index([now + timedelta(i) for i in range(4)], dtype=object)
        result = self.index.union(other)
        expected = np.concatenate((self.index, other))
        tm.assert_numpy_array_equal(result, expected)

        result = other.union(self.index)
        expected = np.concatenate((other, self.index))
        tm.assert_numpy_array_equal(result, expected)

    def test_cant_or_shouldnt_cast(self):
        # can't
        data = ['foo', 'bar', 'baz']
        self.assertRaises(TypeError, Int64Index, data)

        # shouldn't
        data = ['0', '1', '2']
        self.assertRaises(TypeError, Int64Index, data)

    def test_view_Index(self):
        self.index.view(Index)

    def test_prevent_casting(self):
        result = self.index.astype('O')
        self.assertEqual(result.dtype, np.object_)

    def test_take_preserve_name(self):
        index = Int64Index([1, 2, 3, 4], name='foo')
        taken = index.take([3, 0, 1])
        self.assertEqual(index.name, taken.name)

    def test_int_name_format(self):
        from pandas import Series, DataFrame
        index = Index(['a', 'b', 'c'], name=0)
        s = Series(lrange(3), index)
        df = DataFrame(lrange(3), index=index)
        repr(s)
        repr(df)

    def test_print_unicode_columns(self):
        df = pd.DataFrame(
            {u("\u05d0"): [1, 2, 3], "\u05d1": [4, 5, 6], "c": [7, 8, 9]})
        repr(df.columns)  # should not raise UnicodeDecodeError

    def test_repr_summary(self):
        with cf.option_context('display.max_seq_items', 10):
            r = repr(pd.Index(np.arange(1000)))
            self.assertTrue(len(r) < 200)
            self.assertTrue("..." in r)

    def test_repr_roundtrip(self):
        tm.assert_index_equal(eval(repr(self.index)), self.index)

    def test_unicode_string_with_unicode(self):
        idx = Index(lrange(1000))

        if compat.PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        idx = Index(lrange(1000))
        if compat.PY3:
            bytes(idx)
        else:
            str(idx)

    def test_slice_keep_name(self):
        idx = Int64Index([1, 2], name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

    def test_ufunc_coercions(self):
        idx = pd.Int64Index([1, 2, 3, 4, 5], name='x')

        result = np.sqrt(idx)
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index(np.sqrt(np.array([1, 2, 3, 4, 5])), name='x')
        tm.assert_index_equal(result, exp)

        result = np.divide(idx, 2.)
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        tm.assert_index_equal(result, exp)

        # _evaluate_numeric_binop
        result = idx + 2.
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index([3., 4., 5., 6., 7.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx - 2.
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index([-1., 0., 1., 2., 3.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx * 1.
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index([1., 2., 3., 4., 5.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx / 2.
        tm.assertIsInstance(result, Float64Index)
        exp = pd.Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        tm.assert_index_equal(result, exp)


class DatetimeLike(Base):

    def test_shift_identity(self):

        idx = self.create_index()
        self.assert_index_equal(idx, idx.shift(0))

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = 'foo'
        self.assertFalse("length=%s" % len(idx) in str(idx))
        self.assertTrue("'foo'" in str(idx))
        self.assertTrue(idx.__class__.__name__ in str(idx))

        if hasattr(idx,'tz'):
            if idx.tz is not None:
                self.assertTrue(idx.tz in str(idx))
        if hasattr(idx,'freq'):
            self.assertTrue("freq='%s'" % idx.freqstr in str(idx))

    def test_view(self):
        super(DatetimeLike, self).test_view()

        i = self.create_index()

        i_view = i.view('i8')
        result = self._holder(i)
        tm.assert_index_equal(result, i)

        i_view = i.view(self._holder)
        result = self._holder(i)
        tm.assert_index_equal(result, i)

class TestDatetimeIndex(DatetimeLike, tm.TestCase):
    _holder = DatetimeIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index = tm.makeDateIndex(10))
        self.setup_indices()

    def create_index(self):
        return date_range('20130101', periods=5)

    def test_shift(self):

        # test shift for datetimeIndex and non datetimeIndex
        # GH8083

        drange = self.create_index()
        result = drange.shift(1)
        expected = DatetimeIndex(['2013-01-02', '2013-01-03', '2013-01-04', '2013-01-05',
               '2013-01-06'], freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(-1)
        expected = DatetimeIndex(['2012-12-31','2013-01-01', '2013-01-02', '2013-01-03', '2013-01-04'],
                                 freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D')
        expected = DatetimeIndex(['2013-01-07', '2013-01-08', '2013-01-09', '2013-01-10',
               '2013-01-11'],freq='D')
        self.assert_index_equal(result, expected)


    def test_construction_with_alt(self):

        i = pd.date_range('20130101',periods=5,freq='H',tz='US/Eastern')
        i2 = DatetimeIndex(i, dtype=i.dtype)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(i.tz_localize(None).asi8, tz=i.dtype.tz)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(i.tz_localize(None).asi8, dtype=i.dtype)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(i.tz_localize(None).asi8, dtype=i.dtype, tz=i.dtype.tz)
        self.assert_index_equal(i, i2)

        # localize into the provided tz
        i2 = DatetimeIndex(i.tz_localize(None).asi8, tz='UTC')
        expected = i.tz_localize(None).tz_localize('UTC')
        self.assert_index_equal(i2, expected)

        i2 = DatetimeIndex(i, tz='UTC')
        expected = i.tz_convert('UTC')
        self.assert_index_equal(i2, expected)

        # incompat tz/dtype
        self.assertRaises(ValueError, lambda : DatetimeIndex(i.tz_localize(None).asi8, dtype=i.dtype, tz='US/Pacific'))

    def test_pickle_compat_construction(self):
        pass

    def test_get_loc(self):
        idx = pd.date_range('2000-01-01', periods=3)

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_pydatetime(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(idx[1], method,
                                             tolerance=pd.Timedelta('0 days')),
                                 1)

        self.assertEqual(idx.get_loc('2000-01-01', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest'), 1)

        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance='1 day'), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=pd.Timedelta('1D')), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=np.timedelta64(1, 'D')), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=timedelta(1)), 1)
        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc('2000-01-01T12', method='nearest', tolerance='foo')
        with tm.assertRaises(KeyError):
            idx.get_loc('2000-01-01T03', method='nearest',
                        tolerance='2 hours')

        self.assertEqual(idx.get_loc('2000', method='nearest'), slice(0, 3))
        self.assertEqual(idx.get_loc('2000-01', method='nearest'), slice(0, 3))

        self.assertEqual(idx.get_loc('1999', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2001', method='nearest'), 2)

        with tm.assertRaises(KeyError):
            idx.get_loc('1999', method='pad')
        with tm.assertRaises(KeyError):
            idx.get_loc('2001', method='backfill')

        with tm.assertRaises(KeyError):
            idx.get_loc('foobar')
        with tm.assertRaises(TypeError):
            idx.get_loc(slice(2))

        idx = pd.to_datetime(['2000-01-01', '2000-01-04'])
        self.assertEqual(idx.get_loc('2000-01-02', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2000-01-03', method='nearest'), 1)
        self.assertEqual(idx.get_loc('2000-01', method='nearest'), slice(0, 2))

        # time indexing
        idx = pd.date_range('2000-01-01', periods=24, freq='H')
        tm.assert_numpy_array_equal(idx.get_loc(time(12)), [12])
        tm.assert_numpy_array_equal(idx.get_loc(time(12, 30)), [])
        with tm.assertRaises(NotImplementedError):
            idx.get_loc(time(12, 30), method='pad')

    def test_get_indexer(self):
        idx = pd.date_range('2000-01-01', periods=3)
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = idx[0] + pd.to_timedelta(['-1 hour', '12 hours', '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest', tolerance=pd.Timedelta('1 hour')),
            [0, -1, 1])
        with tm.assertRaises(ValueError):
            idx.get_indexer(idx[[0]], method='nearest', tolerance='foo')

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index=date_range('20130101',periods=3,tz='US/Eastern',name='foo')
        unpickled = self.round_trip_pickle(index)
        self.assertTrue(index.equals(unpickled))

    def test_reindex_preserves_tz_if_target_is_empty_list_or_array(self):
        # GH7774
        index = date_range('20130101', periods=3, tz='US/Eastern')
        self.assertEqual(str(index.reindex([])[0].tz), 'US/Eastern')
        self.assertEqual(str(index.reindex(np.array([]))[0].tz), 'US/Eastern')

    def test_time_loc(self):  # GH8667
        from datetime import time
        from pandas.index import _SIZE_CUTOFF

        ns = _SIZE_CUTOFF + np.array([-100, 100],dtype=np.int64)
        key = time(15, 11, 30)
        start = key.hour * 3600 + key.minute * 60 + key.second
        step = 24 * 3600

        for n in ns:
            idx = pd.date_range('2014-11-26', periods=n, freq='S')
            ts = pd.Series(np.random.randn(n), index=idx)
            i = np.arange(start, n, step)

            tm.assert_numpy_array_equal(ts.index.get_loc(key), i)
            tm.assert_series_equal(ts[key], ts.iloc[i])

            left, right = ts.copy(), ts.copy()
            left[key] *= -10
            right.iloc[i] *= -10
            tm.assert_series_equal(left, right)

    def test_time_overflow_for_32bit_machines(self):
        # GH8943.  On some machines NumPy defaults to np.int32 (for example,
        # 32-bit Linux machines).  In the function _generate_regular_range
        # found in tseries/index.py, `periods` gets multiplied by `strides`
        # (which has value 1e9) and since the max value for np.int32 is ~2e9,
        # and since those machines won't promote np.int32 to np.int64, we get
        # overflow.
        periods = np.int_(1000)

        idx1 = pd.date_range(start='2000', periods=periods, freq='S')
        self.assertEqual(len(idx1), periods)

        idx2 = pd.date_range(end='2000', periods=periods, freq='S')
        self.assertEqual(len(idx2), periods)

    def test_intersection(self):
        first = self.index
        second = self.index[5:]
        intersect = first.intersection(second)
        self.assertTrue(tm.equalContents(intersect, second))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            self.assertTrue(tm.equalContents(result, second))

        third = Index(['a', 'b', 'c'])
        result = first.intersection(third)
        expected = pd.Index([], dtype=object)
        self.assert_index_equal(result, expected)

    def test_union(self):
        first = self.index[:5]
        second = self.index[5:]
        everything = self.index
        union = first.union(second)
        self.assertTrue(tm.equalContents(union, everything))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.union(case)
            self.assertTrue(tm.equalContents(result, everything))

    def test_nat(self):
        self.assertIs(DatetimeIndex([np.nan])[0], pd.NaT)


    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '2D')

        for result in [idx - delta, np.subtract(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '2D')

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '3D')

        for result in [idx - delta, np.subtract(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, 'D')

    def test_fillna_datetime64(self):
        # GH 11343
        for tz in ['US/Eastern', 'Asia/Tokyo']:
            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT, '2011-01-01 11:00'])

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'])
            self.assert_index_equal(idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # tz mismatch
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'), pd.Timestamp('2011-01-01 10:00', tz=tz),
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'), 'x',
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)


            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT, '2011-01-01 11:00'], tz=tz)

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'], tz=tz)
            self.assert_index_equal(idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz), pd.Timestamp('2011-01-01 10:00'),
                            pd.Timestamp('2011-01-01 11:00', tz=tz)], dtype=object)
            self.assert_index_equal(idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz), 'x',
                            pd.Timestamp('2011-01-01 11:00', tz=tz)], dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)


class TestPeriodIndex(DatetimeLike, tm.TestCase):
    _holder = PeriodIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index = tm.makePeriodIndex(10))
        self.setup_indices()

    def create_index(self):
        return period_range('20130101', periods=5, freq='D')

    def test_shift(self):

        # test shift for PeriodIndex
        # GH8083
        drange = self.create_index()
        result = drange.shift(1)
        expected = PeriodIndex(['2013-01-02', '2013-01-03', '2013-01-04', '2013-01-05',
             '2013-01-06'], freq='D')
        self.assert_index_equal(result, expected)

    def test_pickle_compat_construction(self):
        pass

    def test_get_loc(self):
        idx = pd.period_range('2000-01-01', periods=3)

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(idx.get_loc(idx[1].asfreq('H', how='start'), method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_timestamp(), method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_timestamp().to_pydatetime(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)

        idx = pd.period_range('2000-01-01', periods=5)[::2]
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance='1 day'), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=pd.Timedelta('1D')), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=np.timedelta64(1, 'D')), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=timedelta(1)), 1)
        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc('2000-01-10', method='nearest', tolerance='foo')

        msg = 'Input has different freq from PeriodIndex\\(freq=D\\)'
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 hour')
        with tm.assertRaises(KeyError):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 day')

    def test_get_indexer(self):
        idx = pd.period_range('2000-01-01', periods=3).asfreq('H', how='start')
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = pd.PeriodIndex(['1999-12-31T23', '2000-01-01T12',
                                 '2000-01-02T01'], freq='H')
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest', tolerance='1 hour'),
            [0, -1, 1])

        msg = 'Input has different freq from PeriodIndex\\(freq=H\\)'
        with self.assertRaisesRegexp(ValueError, msg):
            idx.get_indexer(target, 'nearest', tolerance='1 minute')

        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest', tolerance='1 day'), [0, 1, 1])

    def test_repeat(self):
        # GH10183
        idx = pd.period_range('2000-01-01', periods=3, freq='D')
        res = idx.repeat(3)
        exp = PeriodIndex(idx.values.repeat(3), freq='D')
        self.assert_index_equal(res, exp)
        self.assertEqual(res.freqstr, 'D')

    def test_period_index_indexer(self):

        #GH4125
        idx = pd.period_range('2002-01','2003-12', freq='M')
        df = pd.DataFrame(pd.np.random.randn(24,10), index=idx)
        self.assert_frame_equal(df, df.ix[idx])
        self.assert_frame_equal(df, df.ix[list(idx)])
        self.assert_frame_equal(df, df.loc[list(idx)])
        self.assert_frame_equal(df.iloc[0:5], df.loc[idx[0:5]])
        self.assert_frame_equal(df, df.loc[list(idx)])

    def test_fillna_period(self):
        # GH 11343
        idx = pd.PeriodIndex(['2011-01-01 09:00', pd.NaT, '2011-01-01 11:00'], freq='H')

        exp = pd.PeriodIndex(['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'], freq='H')
        self.assert_index_equal(idx.fillna(pd.Period('2011-01-01 10:00', freq='H')), exp)

        exp = pd.Index([pd.Period('2011-01-01 09:00', freq='H'), 'x',
                        pd.Period('2011-01-01 11:00', freq='H')], dtype=object)
        self.assert_index_equal(idx.fillna('x'), exp)

        with tm.assertRaisesRegexp(ValueError, 'Input has different freq=D from PeriodIndex\\(freq=H\\)'):
            idx.fillna(pd.Period('2011-01-01', freq='D'))


class TestTimedeltaIndex(DatetimeLike, tm.TestCase):
    _holder = TimedeltaIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index = tm.makeTimedeltaIndex(10))
        self.setup_indices()

    def create_index(self):
        return pd.to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)

    def test_shift(self):
        # test shift for TimedeltaIndex
        # err8083

        drange = self.create_index()
        result = drange.shift(1)
        expected = TimedeltaIndex(['1 days 01:00:00', '2 days 01:00:00', '3 days 01:00:00',
                '4 days 01:00:00', '5 days 01:00:00'],freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D 1s')
        expected = TimedeltaIndex(['6 days 01:00:03', '7 days 01:00:03', '8 days 01:00:03',
                '9 days 01:00:03', '10 days 01:00:03'],freq='D')
        self.assert_index_equal(result, expected)

    def test_get_loc(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_pytimedelta(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)

        self.assertEqual(idx.get_loc(idx[1], 'pad', tolerance=pd.Timedelta(0)), 1)
        self.assertEqual(idx.get_loc(idx[1], 'pad', tolerance=np.timedelta64(0, 's')), 1)
        self.assertEqual(idx.get_loc(idx[1], 'pad', tolerance=timedelta(0)), 1)

        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc(idx[1], method='nearest', tolerance='foo')

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc('1 day 1 hour', method), loc)

    def test_get_indexer(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = pd.to_timedelta(['-1 hour', '12 hours', '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest',
                            tolerance=pd.Timedelta('1 hour')),
            [0, -1, 1])

    def test_numeric_compat(self):

        idx = self._holder(np.arange(5,dtype='int64'))
        didx = self._holder(np.arange(5,dtype='int64')**2
                            )
        result = idx * 1
        tm.assert_index_equal(result, idx)

        result = 1 * idx
        tm.assert_index_equal(result, idx)

        result = idx / 1
        tm.assert_index_equal(result, idx)

        result = idx // 1
        tm.assert_index_equal(result, idx)

        result = idx * np.array(5,dtype='int64')
        tm.assert_index_equal(result, self._holder(np.arange(5,dtype='int64')*5))

        result = idx * np.arange(5,dtype='int64')
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5,dtype='int64'))
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5,dtype='float64')+0.1)
        tm.assert_index_equal(result,
                              Float64Index(np.arange(5,dtype='float64')*(np.arange(5,dtype='float64')+0.1)))


        # invalid
        self.assertRaises(TypeError, lambda : idx * idx)
        self.assertRaises(ValueError, lambda : idx * self._holder(np.arange(3)))
        self.assertRaises(ValueError, lambda : idx * np.array([1,2]))

    def test_pickle_compat_construction(self):
        pass

    def test_ufunc_coercions(self):
        # normal ops are also tested in tseries/test_timedeltas.py
        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                              freq='2H', name='x')

        for result in [idx * 2, np.multiply(idx, 2)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['4H', '8H', '12H', '16H', '20H'],
                                 freq='4H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '4H')

        for result in [idx / 2, np.divide(idx, 2)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['1H', '2H', '3H', '4H', '5H'],
                                 freq='H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, 'H')

        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                             freq='2H', name='x')
        for result in [ - idx, np.negative(idx)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['-2H', '-4H', '-6H', '-8H', '-10H'],
                                 freq='-2H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '-2H')

        idx = TimedeltaIndex(['-2H', '-1H', '0H', '1H', '2H'],
                             freq='H', name='x')
        for result in [ abs(idx), np.absolute(idx)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['2H', '1H', '0H', '1H', '2H'],
                                 freq=None, name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, None)

    def test_fillna_timedelta(self):
        # GH 11343
        idx = pd.TimedeltaIndex(['1 day', pd.NaT, '3 day'])

        exp = pd.TimedeltaIndex(['1 day', '2 day', '3 day'])
        self.assert_index_equal(idx.fillna(pd.Timedelta('2 day')), exp)

        exp = pd.TimedeltaIndex(['1 day', '3 hour', '3 day'])
        idx.fillna(pd.Timedelta('3 hour'))

        exp = pd.Index([pd.Timedelta('1 day'), 'x', pd.Timedelta('3 day')], dtype=object)
        self.assert_index_equal(idx.fillna('x'), exp)


class TestMultiIndex(Base, tm.TestCase):
    _holder = MultiIndex
    _multiprocess_can_split_ = True
    _compat_props = ['shape', 'ndim', 'size', 'itemsize']

    def setUp(self):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.indices = dict(index = MultiIndex(levels=[major_axis, minor_axis],
                                               labels=[major_labels, minor_labels],
                                               names=self.index_names, verify_integrity=False))
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
        tm.assertRaisesRegexp(ValueError,'The truth value of a',f)

    def test_labels_dtypes(self):

        # GH 8456
        i = MultiIndex.from_tuples([('A', 1), ('A', 2)])
        self.assertTrue(i.labels[0].dtype == 'int8')
        self.assertTrue(i.labels[1].dtype == 'int8')

        i = MultiIndex.from_product([['a'],range(40)])
        self.assertTrue(i.labels[1].dtype == 'int8')
        i = MultiIndex.from_product([['a'],range(400)])
        self.assertTrue(i.labels[1].dtype == 'int16')
        i = MultiIndex.from_product([['a'],range(40000)])
        self.assertTrue(i.labels[1].dtype == 'int32')

        i = pd.MultiIndex.from_product([['a'],range(1000)])
        self.assertTrue((i.labels[0]>=0).all())
        self.assertTrue((i.labels[1]>=0).all())

    def test_set_name_methods(self):
        # so long as these are synonyms, we don't need to test set_names
        self.assertEqual(self.index.rename, self.index.set_names)
        new_names = [name + "SUFFIX" for name in self.index_names]
        ind = self.index.set_names(new_names)
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, new_names)
        with assertRaisesRegexp(ValueError, "^Length"):
            ind.set_names(new_names + new_names)
        new_names2 = [name + "SUFFIX2" for name in new_names]
        res = ind.set_names(new_names2, inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, new_names2)

        # set names for specific level (# GH7792)
        ind = self.index.set_names(new_names[0], level=0)
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, [new_names[0], self.index_names[1]])

        res = ind.set_names(new_names2[0], level=0, inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, [new_names2[0], self.index_names[1]])

        # set names for multiple levels
        ind = self.index.set_names(new_names, level=[0, 1])
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, new_names)

        res = ind.set_names(new_names2, level=[0, 1], inplace=True)
        self.assertIsNone(res)
        self.assertEqual(ind.names, new_names2)


    def test_set_levels(self):

        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        levels, labels = self.index.levels, self.index.labels
        new_levels = [[lev + 'a' for lev in level] for level in levels]

        def assert_matching(actual, expected):
            # avoid specifying internal representation
            # as much as possible
            self.assertEqual(len(actual), len(expected))
            for act, exp in zip(actual, expected):
                act = np.asarray(act)
                exp = np.asarray(exp)
                assert_almost_equal(act, exp)

        # level changing [w/o mutation]
        ind2 = self.index.set_levels(new_levels)
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # level changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, new_levels)

        # level changing specific level [w/o mutation]
        ind2 = self.index.set_levels(new_levels[0], level=0)
        assert_matching(ind2.levels, [new_levels[0], levels[1]])
        assert_matching(self.index.levels, levels)

        ind2 = self.index.set_levels(new_levels[1], level=1)
        assert_matching(ind2.levels, [levels[0], new_levels[1]])
        assert_matching(self.index.levels, levels)

        # level changing multiple levels [w/o mutation]
        ind2 = self.index.set_levels(new_levels, level=[0, 1])
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # level changing specific level [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels[0], level=0, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, [new_levels[0], levels[1]])
        assert_matching(self.index.levels, levels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels[1], level=1, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, [levels[0], new_levels[1]])
        assert_matching(self.index.levels, levels)

        # level changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, level=[0, 1], inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

    def test_set_labels(self):
        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        levels, labels = self.index.levels, self.index.labels
        major_labels, minor_labels = labels
        major_labels = [(x + 1) % 3 for x in major_labels]
        minor_labels = [(x + 1) % 1 for x in minor_labels]
        new_labels = [major_labels, minor_labels]

        def assert_matching(actual, expected):
            # avoid specifying internal representation
            # as much as possible
            self.assertEqual(len(actual), len(expected))
            for act, exp in zip(actual, expected):
                act = np.asarray(act)
                exp = np.asarray(exp)
                assert_almost_equal(act, exp)

        # label changing [w/o mutation]
        ind2 = self.index.set_labels(new_labels)
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, new_labels)

        # label changing specific level [w/o mutation]
        ind2 = self.index.set_labels(new_labels[0], level=0)
        assert_matching(ind2.labels, [new_labels[0], labels[1]])
        assert_matching(self.index.labels, labels)

        ind2 = self.index.set_labels(new_labels[1], level=1)
        assert_matching(ind2.labels, [labels[0], new_labels[1]])
        assert_matching(self.index.labels, labels)

        # label changing multiple levels [w/o mutation]
        ind2 = self.index.set_labels(new_labels, level=[0, 1])
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing specific level [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels[0], level=0, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, [new_labels[0], labels[1]])
        assert_matching(self.index.labels, labels)

        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels[1], level=1, inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, [labels[0], new_labels[1]])
        assert_matching(self.index.labels, labels)

        # label changing multiple levels [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, level=[0, 1], inplace=True)
        self.assertIsNone(inplace_return)
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

    def test_set_levels_labels_names_bad_input(self):
        levels, labels = self.index.levels, self.index.labels
        names = self.index.names

        with tm.assertRaisesRegexp(ValueError, 'Length of levels'):
            self.index.set_levels([levels[0]])

        with tm.assertRaisesRegexp(ValueError, 'Length of labels'):
            self.index.set_labels([labels[0]])

        with tm.assertRaisesRegexp(ValueError, 'Length of names'):
            self.index.set_names([names[0]])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0])

        # shouldn't scalar data error, instead should demand list-like
        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_names(names[0])

        # should have equal lengths
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_levels(levels[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_levels(levels, level=0)

        # should have equal lengths
        with tm.assertRaisesRegexp(TypeError, 'list of lists-like'):
            self.index.set_labels(labels[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'list-like'):
            self.index.set_labels(labels, level=0)

        # should have equal lengths
        with tm.assertRaisesRegexp(ValueError, 'Length of names'):
            self.index.set_names(names[0], level=[0, 1])

        with tm.assertRaisesRegexp(TypeError, 'string'):
            self.index.set_names(names, level=0)

    def test_metadata_immutable(self):
        levels, labels = self.index.levels, self.index.labels
        # shouldn't be able to set at either the top level or base level
        mutable_regex = re.compile('does not support mutable operations')
        with assertRaisesRegexp(TypeError, mutable_regex):
            levels[0] = levels[0]
        with assertRaisesRegexp(TypeError, mutable_regex):
            levels[0][0] = levels[0][0]
        # ditto for labels
        with assertRaisesRegexp(TypeError, mutable_regex):
            labels[0] = labels[0]
        with assertRaisesRegexp(TypeError, mutable_regex):
            labels[0][0] = labels[0][0]
        # and for names
        names = self.index.names
        with assertRaisesRegexp(TypeError, mutable_regex):
            names[0] = names[0]

    def test_inplace_mutation_resets_values(self):
        levels = [['a', 'b', 'c'], [4]]
        levels2 = [[1, 2, 3], ['a']]
        labels = [[0, 1, 0, 2, 2, 0], [0, 0, 0, 0, 0, 0]]
        mi1 = MultiIndex(levels=levels, labels=labels)
        mi2 = MultiIndex(levels=levels2, labels=labels)
        vals = mi1.values.copy()
        vals2 = mi2.values.copy()
        self.assertIsNotNone(mi1._tuples)

        # make sure level setting works
        new_vals = mi1.set_levels(levels2).values
        assert_almost_equal(vals2, new_vals)
        # non-inplace doesn't kill _tuples [implementation detail]
        assert_almost_equal(mi1._tuples, vals)
        # and values is still same too
        assert_almost_equal(mi1.values, vals)

        # inplace should kill _tuples
        mi1.set_levels(levels2, inplace=True)
        assert_almost_equal(mi1.values, vals2)

        # make sure label setting works too
        labels2 = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        exp_values = np.empty((6, ), dtype=object)
        exp_values[:] = [(long(1), 'a')] * 6
        # must be 1d array of tuples
        self.assertEqual(exp_values.shape, (6, ))
        new_values = mi2.set_labels(labels2).values
        # not inplace shouldn't change
        assert_almost_equal(mi2._tuples, vals2)
        # should have correct values
        assert_almost_equal(exp_values, new_values)

        # and again setting inplace should kill _tuples, etc
        mi2.set_labels(labels2, inplace=True)
        assert_almost_equal(mi2.values, new_values)

    def test_copy_in_constructor(self):
        levels = np.array(["a", "b", "c"])
        labels = np.array([1, 1, 2, 0, 0, 1, 1])
        val = labels[0]
        mi = MultiIndex(levels=[levels, levels], labels=[labels, labels],
                        copy=True)
        self.assertEqual(mi.labels[0][0], val)
        labels[0] = 15
        self.assertEqual(mi.labels[0][0], val)
        val = levels[0]
        levels[0] = "PANDA"
        self.assertEqual(mi.levels[0][0], val)

    def test_set_value_keeps_names(self):
        # motivating example from #3742
        lev1 = ['hans', 'hans', 'hans', 'grethe', 'grethe', 'grethe']
        lev2 = ['1', '2', '3'] * 2
        idx = pd.MultiIndex.from_arrays(
            [lev1, lev2],
            names=['Name', 'Number'])
        df = pd.DataFrame(
            np.random.randn(6, 4),
            columns=['one', 'two', 'three', 'four'],
            index=idx)
        df = df.sortlevel()
        self.assertIsNone(df.is_copy)
        self.assertEqual(df.index.names, ('Name', 'Number'))
        df = df.set_value(('grethe', '4'), 'one', 99.34)
        self.assertIsNone(df.is_copy)
        self.assertEqual(df.index.names, ('Name', 'Number'))

    def test_names(self):

        # names are assigned in __init__
        names = self.index_names
        level_names = [level.name for level in self.index.levels]
        self.assertEqual(names, level_names)

        # setting bad names on existing
        index = self.index
        assertRaisesRegexp(ValueError, "^Length of names", setattr, index,
                           "names", list(index.names) + ["third"])
        assertRaisesRegexp(ValueError, "^Length of names", setattr, index,
                           "names", [])

        # initializing with bad names (should always be equivalent)
        major_axis, minor_axis = self.index.levels
        major_labels, minor_labels = self.index.labels
        assertRaisesRegexp(ValueError, "^Length of names", MultiIndex,
                           levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels],
                           names=['first'])
        assertRaisesRegexp(ValueError, "^Length of names", MultiIndex,
                           levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels],
                           names=['first', 'second', 'third'])

        # names are assigned
        index.names = ["a", "b"]
        ind_names = list(index.names)
        level_names = [level.name for level in index.levels]
        self.assertEqual(ind_names, level_names)

    def test_reference_duplicate_name(self):
        idx = MultiIndex.from_tuples([('a', 'b'), ('c', 'd')], names=['x', 'x'])
        self.assertTrue(idx._reference_duplicate_name('x'))

        idx = MultiIndex.from_tuples([('a', 'b'), ('c', 'd')], names=['x', 'y'])
        self.assertFalse(idx._reference_duplicate_name('x'))

    def test_astype(self):
        expected = self.index.copy()
        actual = self.index.astype('O')
        assert_copy(actual.levels, expected.levels)
        assert_copy(actual.labels, expected.labels)
        self.check_level_names(actual, expected.names)

        with assertRaisesRegexp(TypeError, "^Setting.*dtype.*object"):
            self.index.astype(np.dtype(int))

    def test_constructor_single_level(self):
        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]],
                                  names=['first'])
        tm.assertIsInstance(single_level, Index)
        self.assertNotIsInstance(single_level, MultiIndex)
        self.assertEqual(single_level.name, 'first')

        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]])
        self.assertIsNone(single_level.name)

    def test_constructor_no_levels(self):
        assertRaisesRegexp(ValueError, "non-zero number of levels/labels",
                           MultiIndex, levels=[], labels=[])
        both_re = re.compile('Must pass both levels and labels')
        with tm.assertRaisesRegexp(TypeError, both_re):
            MultiIndex(levels=[])
        with tm.assertRaisesRegexp(TypeError, both_re):
            MultiIndex(labels=[])

    def test_constructor_mismatched_label_levels(self):
        labels = [np.array([1]), np.array([2]), np.array([3])]
        levels = ["a"]
        assertRaisesRegexp(ValueError, "Length of levels and labels must be"
                           " the same", MultiIndex, levels=levels,
                           labels=labels)
        length_error = re.compile('>= length of level')
        label_error = re.compile(r'Unequal label lengths: \[4, 2\]')

        # important to check that it's looking at the right thing.
        with tm.assertRaisesRegexp(ValueError, length_error):
            MultiIndex(levels=[['a'], ['b']], labels=[[0, 1, 2, 3], [0, 3, 4, 1]])

        with tm.assertRaisesRegexp(ValueError, label_error):
            MultiIndex(levels=[['a'], ['b']], labels=[[0, 0, 0, 0], [0, 0]])

        # external API
        with tm.assertRaisesRegexp(ValueError, length_error):
            self.index.copy().set_levels([['a'], ['b']])

        with tm.assertRaisesRegexp(ValueError, label_error):
            self.index.copy().set_labels([[0, 0, 0, 0], [0, 0]])

        # deprecated properties
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            with tm.assertRaisesRegexp(ValueError, length_error):
                self.index.copy().levels = [['a'], ['b']]

            with tm.assertRaisesRegexp(ValueError, label_error):
                self.index.copy().labels = [[0, 0, 0, 0], [0, 0]]


    def assert_multiindex_copied(self, copy, original):
        # levels shoudl be (at least, shallow copied)
        assert_copy(copy.levels, original.levels)

        assert_almost_equal(copy.labels, original.labels)

        # labels doesn't matter which way copied
        assert_almost_equal(copy.labels, original.labels)
        self.assertIsNot(copy.labels, original.labels)

        # names doesn't matter which way copied
        self.assertEqual(copy.names, original.names)
        self.assertIsNot(copy.names, original.names)

        # sort order should be copied
        self.assertEqual(copy.sortorder, original.sortorder)

    def test_copy(self):
        i_copy = self.index.copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_shallow_copy(self):
        i_copy = self.index._shallow_copy()

        self.assert_multiindex_copied(i_copy, self.index)

    def test_view(self):
        i_view = self.index.view()

        self.assert_multiindex_copied(i_view, self.index)

    def check_level_names(self, index, names):
        self.assertEqual([level.name for level in index.levels], list(names))

    def test_changing_names(self):

        # names should be applied to levels
        level_names = [level.name for level in self.index.levels]
        self.check_level_names(self.index, self.index.names)

        view = self.index.view()
        copy = self.index.copy()
        shallow_copy = self.index._shallow_copy()

        # changing names should change level names on object
        new_names = [name + "a" for name in self.index.names]
        self.index.names = new_names
        self.check_level_names(self.index, new_names)

        # but not on copies
        self.check_level_names(view, level_names)
        self.check_level_names(copy, level_names)
        self.check_level_names(shallow_copy, level_names)

        # and copies shouldn't change original
        shallow_copy.names = [name + "c" for name in shallow_copy.names]
        self.check_level_names(self.index, new_names)

    def test_duplicate_names(self):
        self.index.names = ['foo', 'foo']
        assertRaisesRegexp(KeyError, 'Level foo not found',
                           self.index._get_level_number, 'foo')

    def test_get_level_number_integer(self):
        self.index.names = [1, 0]
        self.assertEqual(self.index._get_level_number(1), 0)
        self.assertEqual(self.index._get_level_number(0), 1)
        self.assertRaises(IndexError, self.index._get_level_number, 2)
        assertRaisesRegexp(KeyError, 'Level fourth not found',
                           self.index._get_level_number, 'fourth')

    def test_from_arrays(self):
        arrays = []
        for lev, lab in zip(self.index.levels, self.index.labels):
            arrays.append(np.asarray(lev).take(lab))

        result = MultiIndex.from_arrays(arrays)
        self.assertEqual(list(result), list(self.index))

        # infer correctly
        result = MultiIndex.from_arrays([[pd.NaT, Timestamp('20130101')], ['a', 'b']])
        self.assertTrue(result.levels[0].equals(Index([Timestamp('20130101')])))
        self.assertTrue(result.levels[1].equals(Index(['a','b'])))

    def test_from_product(self):

        first = ['foo', 'bar', 'buz']
        second = ['a', 'b', 'c']
        names = ['first', 'second']
        result = MultiIndex.from_product([first, second], names=names)

        tuples = [('foo', 'a'), ('foo', 'b'), ('foo', 'c'),
                  ('bar', 'a'), ('bar', 'b'), ('bar', 'c'),
                  ('buz', 'a'), ('buz', 'b'), ('buz', 'c')]
        expected = MultiIndex.from_tuples(tuples, names=names)

        tm.assert_numpy_array_equal(result, expected)
        self.assertEqual(result.names, names)

    def test_from_product_datetimeindex(self):
        dt_index = date_range('2000-01-01', periods=2)
        mi = pd.MultiIndex.from_product([[1, 2], dt_index])
        etalon = pd.lib.list_to_object_array([(1, pd.Timestamp('2000-01-01')),
                                              (1, pd.Timestamp('2000-01-02')),
                                              (2, pd.Timestamp('2000-01-01')),
                                              (2, pd.Timestamp('2000-01-02'))])
        tm.assert_numpy_array_equal(mi.values, etalon)

    def test_values_boxed(self):
        tuples = [(1, pd.Timestamp('2000-01-01')),
                  (2, pd.NaT),
                  (3, pd.Timestamp('2000-01-03')),
                  (1, pd.Timestamp('2000-01-04')),
                  (2, pd.Timestamp('2000-01-02')),
                  (3, pd.Timestamp('2000-01-03'))]
        mi = pd.MultiIndex.from_tuples(tuples)
        tm.assert_numpy_array_equal(mi.values, pd.lib.list_to_object_array(tuples))
        # Check that code branches for boxed values produce identical results
        tm.assert_numpy_array_equal(mi.values[:4], mi[:4].values)

    def test_append(self):
        result = self.index[:3].append(self.index[3:])
        self.assertTrue(result.equals(self.index))

        foos = [self.index[:1], self.index[1:3], self.index[3:]]
        result = foos[0].append(foos[1:])
        self.assertTrue(result.equals(self.index))

        # empty
        result = self.index.append([])
        self.assertTrue(result.equals(self.index))

    def test_get_level_values(self):
        result = self.index.get_level_values(0)
        expected = ['foo', 'foo', 'bar', 'baz', 'qux', 'qux']
        tm.assert_numpy_array_equal(result, expected)

        self.assertEqual(result.name, 'first')

        result = self.index.get_level_values('first')
        expected = self.index.get_level_values(0)
        tm.assert_numpy_array_equal(result, expected)

        # GH 10460
        index = MultiIndex(levels=[CategoricalIndex(['A', 'B']),
                                   CategoricalIndex([1, 2, 3])],
                           labels=[np.array([0, 0, 0, 1, 1, 1]),
                                   np.array([0, 1, 2, 0, 1, 2])])
        exp = CategoricalIndex(['A', 'A', 'A', 'B', 'B', 'B'])
        self.assert_index_equal(index.get_level_values(0), exp)
        exp = CategoricalIndex([1, 2 ,3, 1, 2, 3])
        self.assert_index_equal(index.get_level_values(1), exp)

    def test_get_level_values_na(self):
        arrays = [['a', 'b', 'b'], [1, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = [1, np.nan, 2]
        tm.assert_numpy_array_equal(values.values.astype(float), expected)

        arrays = [['a', 'b', 'b'], [np.nan, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = [np.nan, np.nan, 2]
        tm.assert_numpy_array_equal(values.values.astype(float), expected)

        arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        expected = [np.nan, np.nan, np.nan]
        tm.assert_numpy_array_equal(values.values.astype(float), expected)
        values = index.get_level_values(1)
        expected = np.array(['a', np.nan, 1],dtype=object)
        tm.assert_numpy_array_equal(values.values, expected)

        arrays = [['a', 'b', 'b'], pd.DatetimeIndex([0, 1, pd.NaT])]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = pd.DatetimeIndex([0, 1, pd.NaT])
        tm.assert_numpy_array_equal(values.values, expected.values)

        arrays = [[], []]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        self.assertEqual(values.shape, (0,))

    def test_reorder_levels(self):
        # this blows up
        assertRaisesRegexp(IndexError, '^Too many levels',
                           self.index.reorder_levels, [2, 1, 0])

    def test_nlevels(self):
        self.assertEqual(self.index.nlevels, 2)

    def test_iter(self):
        result = list(self.index)
        expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                    ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
        self.assertEqual(result, expected)

    def test_legacy_pickle(self):
        if compat.PY3:
            raise nose.SkipTest("testing for legacy pickles not support on py3")

        path = tm.get_data_path('multiindex_v1.pickle')
        obj = pd.read_pickle(path)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assertTrue(obj.equals(obj2))

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj))
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
        self.assertTrue(obj.equals(obj2))

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj))
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index=MultiIndex.from_product([[1,2],['a','b'],date_range('20130101',periods=3,tz='US/Eastern')],names=['one','two','three'])
        unpickled = self.round_trip_pickle(index)
        self.assertTrue(index.equal_levels(unpickled))

    def test_from_tuples_index_values(self):
        result = MultiIndex.from_tuples(self.index)
        self.assertTrue((result.values == self.index.values).all())

    def test_contains(self):
        self.assertIn(('foo', 'two'), self.index)
        self.assertNotIn(('bar', 'two'), self.index)
        self.assertNotIn(None, self.index)

    def test_is_all_dates(self):
        self.assertFalse(self.index.is_all_dates)

    def test_is_numeric(self):
        # MultiIndex is never numeric
        self.assertFalse(self.index.is_numeric())

    def test_getitem(self):
        # scalar
        self.assertEqual(self.index[2], ('bar', 'one'))

        # slice
        result = self.index[2:5]
        expected = self.index[[2, 3, 4]]
        self.assertTrue(result.equals(expected))

        # boolean
        result = self.index[[True, False, True, False, True, True]]
        result2 = self.index[np.array([True, False, True, False, True, True])]
        expected = self.index[[0, 2, 4, 5]]
        self.assertTrue(result.equals(expected))
        self.assertTrue(result2.equals(expected))

    def test_getitem_group_select(self):
        sorted_idx, _ = self.index.sortlevel(0)
        self.assertEqual(sorted_idx.get_loc('baz'), slice(3, 4))
        self.assertEqual(sorted_idx.get_loc('foo'), slice(0, 2))

    def test_get_loc(self):
        self.assertEqual(self.index.get_loc(('foo', 'two')), 1)
        self.assertEqual(self.index.get_loc(('baz', 'two')), 3)
        self.assertRaises(KeyError, self.index.get_loc, ('bar', 'two'))
        self.assertRaises(KeyError, self.index.get_loc, 'quux')

        self.assertRaises(NotImplementedError, self.index.get_loc, 'foo',
                          method='nearest')

        # 3 levels
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])
        self.assertRaises(KeyError, index.get_loc, (1, 1))
        self.assertEqual(index.get_loc((2, 0)), slice(3, 5))

    def test_get_loc_duplicates(self):
        index = Index([2, 2, 2, 2])
        result = index.get_loc(2)
        expected = slice(0, 4)
        self.assertEqual(result, expected)
        # self.assertRaises(Exception, index.get_loc, 2)

        index = Index(['c', 'a', 'a', 'b', 'b'])
        rs = index.get_loc('c')
        xp = 0
        assert(rs == xp)

    def test_get_loc_level(self):
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        loc, new_index = index.get_loc_level((0, 1))
        expected = slice(1, 2)
        exp_index = index[expected].droplevel(0).droplevel(0)
        self.assertEqual(loc, expected)
        self.assertTrue(new_index.equals(exp_index))

        loc, new_index = index.get_loc_level((0, 1, 0))
        expected = 1
        self.assertEqual(loc, expected)
        self.assertIsNone(new_index)

        self.assertRaises(KeyError, index.get_loc_level, (2, 2))

        index = MultiIndex(levels=[[2000], lrange(4)],
                           labels=[np.array([0, 0, 0, 0]),
                                   np.array([0, 1, 2, 3])])
        result, new_index = index.get_loc_level((2000, slice(None, None)))
        expected = slice(None, None)
        self.assertEqual(result, expected)
        self.assertTrue(new_index.equals(index.droplevel(0)))

    def test_slice_locs(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index

        slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
        sliced = stacked[slob]
        expected = df[5:16].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

        slob = slice(*idx.slice_locs(df.index[5] + timedelta(seconds=30),
                                     df.index[15] - timedelta(seconds=30)))
        sliced = stacked[slob]
        expected = df[6:15].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

    def test_slice_locs_with_type_mismatch(self):
        df = tm.makeTimeDataFrame()
        stacked = df.stack()
        idx = stacked.index
        assertRaisesRegexp(TypeError, '^Level type mismatch', idx.slice_locs,
                           (1, 3))
        assertRaisesRegexp(TypeError, '^Level type mismatch', idx.slice_locs,
                           df.index[5] + timedelta(seconds=30), (5, 2))
        df = tm.makeCustomDataframe(5, 5)
        stacked = df.stack()
        idx = stacked.index
        with assertRaisesRegexp(TypeError, '^Level type mismatch'):
            idx.slice_locs(timedelta(seconds=30))
        # TODO: Try creating a UnicodeDecodeError in exception message
        with assertRaisesRegexp(TypeError, '^Level type mismatch'):
            idx.slice_locs(df.index[1], (16, "a"))

    def test_slice_locs_not_sorted(self):
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        assertRaisesRegexp(KeyError, "[Kk]ey length.*greater than MultiIndex"
                           " lexsort depth", index.slice_locs, (1, 0, 1),
                           (2, 1, 0))

        # works
        sorted_index, _ = index.sortlevel(0)
        # should there be a test case here???
        sorted_index.slice_locs((1, 0, 1), (2, 1, 0))

    def test_slice_locs_partial(self):
        sorted_idx, _ = self.index.sortlevel(0)

        result = sorted_idx.slice_locs(('foo', 'two'), ('qux', 'one'))
        self.assertEqual(result, (1, 5))

        result = sorted_idx.slice_locs(None, ('qux', 'one'))
        self.assertEqual(result, (0, 5))

        result = sorted_idx.slice_locs(('foo', 'two'), None)
        self.assertEqual(result, (1, len(sorted_idx)))

        result = sorted_idx.slice_locs('bar', 'baz')
        self.assertEqual(result, (2, 4))

    def test_slice_locs_not_contained(self):
        # some searchsorted action

        index = MultiIndex(levels=[[0, 2, 4, 6], [0, 2, 4]],
                           labels=[[0, 0, 0, 1, 1, 2, 3, 3, 3],
                                   [0, 1, 2, 1, 2, 2, 0, 1, 2]],
                           sortorder=0)

        result = index.slice_locs((1, 0), (5, 2))
        self.assertEqual(result, (3, 6))

        result = index.slice_locs(1, 5)
        self.assertEqual(result, (3, 6))

        result = index.slice_locs((2, 2), (5, 2))
        self.assertEqual(result, (3, 6))

        result = index.slice_locs(2, 5)
        self.assertEqual(result, (3, 6))

        result = index.slice_locs((1, 0), (6, 3))
        self.assertEqual(result, (3, 8))

        result = index.slice_locs(-1, 10)
        self.assertEqual(result, (0, len(index)))

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

        self.assertFalse(index.is_unique)

    def test_truncate(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        result = index.truncate(before=1)
        self.assertNotIn('foo', result.levels[0])
        self.assertIn(1, result.levels[0])

        result = index.truncate(after=1)
        self.assertNotIn(2, result.levels[0])
        self.assertIn(1, result.levels[0])

        result = index.truncate(before=1, after=2)
        self.assertEqual(len(result.levels[0]), 2)

        # after < before
        self.assertRaises(ValueError, index.truncate, 3, 1)

    def test_get_indexer(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        idx1 = index[:5]
        idx2 = index[[1, 3, 5]]

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, [1, 3, -1])

        r1 = idx2.get_indexer(idx1, method='pad')
        e1 = [-1, 0, 0, 1, 1]
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='pad')
        assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        e1 = [0, 0, 1, 1, 2]
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='backfill')
        assert_almost_equal(r2, e1[::-1])

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        # pass non-MultiIndex
        r1 = idx1.get_indexer(idx2._tuple_index)
        rexp1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1, 2, 3])
        self.assertTrue((r1 == [-1, -1, -1]).all())

        # create index with duplicates
        idx1 = Index(lrange(10) + lrange(10))
        idx2 = Index(lrange(20))
        assertRaisesRegexp(InvalidIndexError, "Reindexing only valid with"
                           " uniquely valued Index objects",
                           idx1.get_indexer, idx2)

    def test_get_indexer_nearest(self):
        midx = MultiIndex.from_tuples([('a', 1), ('b', 2)])
        with tm.assertRaises(NotImplementedError):
            midx.get_indexer(['a'], method='nearest')
        with tm.assertRaises(NotImplementedError):
            midx.get_indexer(['a'], method='pad', tolerance=2)

    def test_format(self):
        self.index.format()
        self.index[:0].format()

    def test_format_integer_names(self):
        index = MultiIndex(levels=[[0, 1], [0, 1]],
                           labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                           names=[0, 1])
        index.format(names=True)

    def test_format_sparse_display(self):
        index = MultiIndex(levels=[[0, 1], [0, 1], [0, 1], [0]],
                           labels=[[0, 0, 0, 1, 1, 1],
                                   [0, 0, 1, 0, 0, 1],
                                   [0, 1, 0, 0, 1, 0],
                                   [0, 0, 0, 0, 0, 0]])

        result = index.format()
        self.assertEqual(result[3], '1  0  0  0')

    def test_format_sparse_config(self):
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")
        # GH1538
        pd.set_option('display.multi_sparse', False)

        result = self.index.format()
        self.assertEqual(result[1], 'foo  two')

        self.reset_display_options()

        warnings.filters = warn_filters

    def test_to_hierarchical(self):
        index = MultiIndex.from_tuples([(1, 'one'), (1, 'two'),
                                        (2, 'one'), (2, 'two')])
        result = index.to_hierarchical(3)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                      [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

        # K > 1
        result = index.to_hierarchical(3, 2)
        expected = MultiIndex(levels=[[1, 2], ['one', 'two']],
                              labels=[[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                                      [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

        # non-sorted
        index = MultiIndex.from_tuples([(2, 'c'), (1, 'b'),
                                        (2, 'a'), (2, 'b')],
                                       names=['N1', 'N2'])

        result = index.to_hierarchical(2)
        expected = MultiIndex.from_tuples([(2, 'c'), (2, 'c'), (1, 'b'), (1, 'b'),
                                           (2, 'a'), (2, 'a'), (2, 'b'), (2, 'b')],
                                          names=['N1', 'N2'])
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.names, index.names)

    def test_bounds(self):
        self.index._bounds

    def test_equals(self):
        self.assertTrue(self.index.equals(self.index))
        self.assertTrue(self.index.equal_levels(self.index))

        self.assertFalse(self.index.equals(self.index[:-1]))

        self.assertTrue(self.index.equals(self.index._tuple_index))

        # different number of levels
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        index2 = MultiIndex(levels=index.levels[:-1],
                            labels=index.labels[:-1])
        self.assertFalse(index.equals(index2))
        self.assertFalse(index.equal_levels(index2))

        # levels are different
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assertFalse(self.index.equals(index))
        self.assertFalse(self.index.equal_levels(index))

        # some of the labels are different
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assertFalse(self.index.equals(index))

    def test_identical(self):
        mi = self.index.copy()
        mi2 = self.index.copy()
        self.assertTrue(mi.identical(mi2))

        mi = mi.set_names(['new1', 'new2'])
        self.assertTrue(mi.equals(mi2))
        self.assertFalse(mi.identical(mi2))

        mi2 = mi2.set_names(['new1', 'new2'])
        self.assertTrue(mi.identical(mi2))

        mi3 = Index(mi.tolist(), names=mi.names)
        mi4 = Index(mi.tolist(), names=mi.names, tupleize_cols=False)
        self.assertTrue(mi.identical(mi3))
        self.assertFalse(mi.identical(mi4))
        self.assertTrue(mi.equals(mi4))

    def test_is_(self):

        mi = MultiIndex.from_tuples(lzip(range(10), range(10)))
        self.assertTrue(mi.is_(mi))
        self.assertTrue(mi.is_(mi.view()))
        self.assertTrue(mi.is_(mi.view().view().view().view()))
        mi2 = mi.view()
        # names are metadata, they don't change id
        mi2.names = ["A", "B"]
        self.assertTrue(mi2.is_(mi))
        self.assertTrue(mi.is_(mi2))

        self.assertTrue(mi.is_(mi.set_names(["C", "D"])))
        mi2 = mi.view()
        mi2.set_names(["E", "F"], inplace=True)
        self.assertTrue(mi.is_(mi2))
        # levels are inherent properties, they change identity
        mi3 = mi2.set_levels([lrange(10), lrange(10)])
        self.assertFalse(mi3.is_(mi2))
        # shouldn't change
        self.assertTrue(mi2.is_(mi))
        mi4 = mi3.view()
        mi4.set_levels([[1 for _ in range(10)], lrange(10)], inplace=True)
        self.assertFalse(mi4.is_(mi3))
        mi5 = mi.view()
        mi5.set_levels(mi5.levels, inplace=True)
        self.assertFalse(mi5.is_(mi))

    def test_union(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_union = piece1 | piece2

        tups = sorted(self.index._tuple_index)
        expected = MultiIndex.from_tuples(tups)

        self.assertTrue(the_union.equals(expected))

        # corner case, pass self or empty thing:
        the_union = self.index.union(self.index)
        self.assertIs(the_union, self.index)

        the_union = self.index.union(self.index[:0])
        self.assertIs(the_union, self.index)

        # won't work in python 3
        # tuples = self.index._tuple_index
        # result = self.index[:4] | tuples[4:]
        # self.assertTrue(result.equals(tuples))

    # not valid for python 3
    # def test_union_with_regular_index(self):
    #     other = Index(['A', 'B', 'C'])

    #     result = other.union(self.index)
    #     self.assertIn(('foo', 'one'), result)
    #     self.assertIn('B', result)

    #     result2 = self.index.union(other)
    #     self.assertTrue(result.equals(result2))

    def test_intersection(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_int = piece1 & piece2
        tups = sorted(self.index[3:5]._tuple_index)
        expected = MultiIndex.from_tuples(tups)
        self.assertTrue(the_int.equals(expected))

        # corner case, pass self
        the_int = self.index.intersection(self.index)
        self.assertIs(the_int, self.index)

        # empty intersection: disjoint
        empty = self.index[:2] & self.index[2:]
        expected = self.index[:0]
        self.assertTrue(empty.equals(expected))

        # can't do in python 3
        # tuples = self.index._tuple_index
        # result = self.index & tuples
        # self.assertTrue(result.equals(tuples))

    def test_difference(self):

        first = self.index
        result = first.difference(self.index[-3:])

        # - API change GH 8226
        with tm.assert_produces_warning():
            first - self.index[-3:]
        with tm.assert_produces_warning():
            self.index[-3:] - first
        with tm.assert_produces_warning():
            self.index[-3:] - first.tolist()

        self.assertRaises(TypeError, lambda : first.tolist() - self.index[-3:])

        expected = MultiIndex.from_tuples(sorted(self.index[:-3].values),
                                          sortorder=0,
                                          names=self.index.names)

        tm.assertIsInstance(result, MultiIndex)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: reflexive
        result = self.index.difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: superset
        result = self.index[-3:].difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: degenerate
        result = self.index[:0].difference(self.index)
        expected = self.index[:0]
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # names not the same
        chunklet = self.index[-3:]
        chunklet.names = ['foo', 'baz']
        result = first.difference(chunklet)
        self.assertEqual(result.names, (None, None))

        # empty, but non-equal
        result = self.index.difference(self.index.sortlevel(1)[0])
        self.assertEqual(len(result), 0)

        # raise Exception called with non-MultiIndex
        result = first.difference(first._tuple_index)
        self.assertTrue(result.equals(first[:0]))

        # name from empty array
        result = first.difference([])
        self.assertTrue(first.equals(result))
        self.assertEqual(first.names, result.names)

        # name from non-empty array
        result = first.difference([('foo', 'one')])
        expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'),
                                            ('foo', 'two'), ('qux', 'one'),
                                            ('qux', 'two')])
        expected.names = first.names
        self.assertEqual(first.names, result.names)
        assertRaisesRegexp(TypeError, "other must be a MultiIndex or a list"
                           " of tuples", first.difference, [1, 2, 3, 4, 5])

    def test_from_tuples(self):
        assertRaisesRegexp(TypeError, 'Cannot infer number of levels from'
                           ' empty list', MultiIndex.from_tuples, [])

        idx = MultiIndex.from_tuples(((1, 2), (3, 4)), names=['a', 'b'])
        self.assertEqual(len(idx), 2)

    def test_argsort(self):
        result = self.index.argsort()
        expected = self.index._tuple_index.argsort()
        tm.assert_numpy_array_equal(result, expected)

    def test_sortlevel(self):
        import random

        tuples = list(self.index)
        random.shuffle(tuples)

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

    def test_sortlevel_not_sort_remaining(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        sorted_idx, _ = mi.sortlevel('A', sort_remaining=False)
        self.assertTrue(sorted_idx.equals(mi))

    def test_sortlevel_deterministic(self):
        tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
                  ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assertTrue(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assertTrue(sorted_idx.equals(expected[::-1]))

    def test_dims(self):
        pass

    def test_drop(self):
        dropped = self.index.drop([('foo', 'two'), ('qux', 'one')])

        index = MultiIndex.from_tuples([('foo', 'two'), ('qux', 'one')])
        dropped2 = self.index.drop(index)

        expected = self.index[[0, 2, 3, 5]]
        self.assert_index_equal(dropped, expected)
        self.assert_index_equal(dropped2, expected)

        dropped = self.index.drop(['bar'])
        expected = self.index[[0, 1, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop('foo')
        expected = self.index[[2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        index = MultiIndex.from_tuples([('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, [('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, index)
        self.assertRaises(KeyError, self.index.drop, ['foo', 'two'])

        # partially correct argument
        mixed_index = MultiIndex.from_tuples([('qux', 'one'), ('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, mixed_index)

        # error='ignore'
        dropped = self.index.drop(index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[0, 1, 2, 3, 5]]
        self.assert_index_equal(dropped, expected)

        dropped = self.index.drop(['foo', 'two'], errors='ignore')
        expected = self.index[[2, 3, 4, 5]]
        self.assert_index_equal(dropped, expected)

        # mixed partial / full drop
        dropped = self.index.drop(['foo', ('qux', 'one')])
        expected = self.index[[2, 3, 5]]
        self.assert_index_equal(dropped, expected)

        # mixed partial / full drop / error='ignore'
        mixed_index = ['foo', ('qux', 'one'), 'two']
        self.assertRaises(KeyError, self.index.drop, mixed_index)
        dropped = self.index.drop(mixed_index, errors='ignore')
        expected = self.index[[2, 3, 5]]
        self.assert_index_equal(dropped, expected)

    def test_droplevel_with_names(self):
        index = self.index[self.index.get_loc('foo')]
        dropped = index.droplevel(0)
        self.assertEqual(dropped.name, 'second')

        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])],
                           names=['one', 'two', 'three'])
        dropped = index.droplevel(0)
        self.assertEqual(dropped.names, ('two', 'three'))

        dropped = index.droplevel('two')
        expected = index.droplevel(1)
        self.assertTrue(dropped.equals(expected))

    def test_droplevel_multiple(self):
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])],
                           names=['one', 'two', 'three'])

        dropped = index[:2].droplevel(['three', 'one'])
        expected = index[:2].droplevel(2).droplevel(0)
        self.assertTrue(dropped.equals(expected))

    def test_insert(self):
        # key contained in all levels
        new_index = self.index.insert(0, ('bar', 'two'))
        self.assertTrue(new_index.equal_levels(self.index))
        self.assertEqual(new_index[0], ('bar', 'two'))

        # key not contained in all levels
        new_index = self.index.insert(0, ('abc', 'three'))
        tm.assert_numpy_array_equal(new_index.levels[0],
                                      list(self.index.levels[0]) + ['abc'])
        tm.assert_numpy_array_equal(new_index.levels[1],
                                      list(self.index.levels[1]) + ['three'])
        self.assertEqual(new_index[0], ('abc', 'three'))

        # key wrong length
        assertRaisesRegexp(ValueError, "Item must have length equal to number"
                           " of levels", self.index.insert, 0, ('foo2',))

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

        right = pd.DataFrame([['a', 'b',  0],
                              ['b', 'd',  1],
                              ['b', 'x',  2],
                              ['b', 'a', -1],
                              ['b', 'b',  3],
                              ['a', 'x',  4],
                              ['a', 'w',  5],
                              ['a', 'a',  6]],
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
        left.ix[('test', 18)] = 12

        right = pd.Series(np.linspace(0, 12, 13),
                          pd.MultiIndex.from_tuples(idx))

        tm.assert_series_equal(left, right)

    def test_take_preserve_name(self):
        taken = self.index.take([3, 0, 1])
        self.assertEqual(taken.names, self.index.names)

    def test_join_level(self):
        def _check_how(other, how):
            join_index, lidx, ridx = other.join(self.index, how=how,
                                                level='second',
                                                return_indexers=True)

            exp_level = other.join(self.index.levels[1], how=how)
            self.assertTrue(join_index.levels[0].equals(self.index.levels[0]))
            self.assertTrue(join_index.levels[1].equals(exp_level))

            # pare down levels
            mask = np.array(
                [x[1] in exp_level for x in self.index], dtype=bool)
            exp_values = self.index.values[mask]
            tm.assert_numpy_array_equal(join_index.values, exp_values)

            if how in ('outer', 'inner'):
                join_index2, ridx2, lidx2 = \
                    self.index.join(other, how=how, level='second',
                                    return_indexers=True)

                self.assertTrue(join_index.equals(join_index2))
                tm.assert_numpy_array_equal(lidx, lidx2)
                tm.assert_numpy_array_equal(ridx, ridx2)
                tm.assert_numpy_array_equal(join_index2.values, exp_values)

        def _check_all(other):
            _check_how(other, 'outer')
            _check_how(other, 'inner')
            _check_how(other, 'left')
            _check_how(other, 'right')

        _check_all(Index(['three', 'one', 'two']))
        _check_all(Index(['one']))
        _check_all(Index(['one', 'three']))

        # some corner cases
        idx = Index(['three', 'one', 'two'])
        result = idx.join(self.index, level='second')
        tm.assertIsInstance(result, MultiIndex)

        assertRaisesRegexp(TypeError, "Join.*MultiIndex.*ambiguous",
                           self.index.join, self.index, level=1)

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            res = self.index
            joined = res.join(res, how=kind)
            self.assertIs(res, joined)

    def test_join_multi(self):
        # GH 10665
        midx = pd.MultiIndex.from_product([np.arange(4), np.arange(4)], names=['a', 'b'])
        idx = pd.Index([1, 2, 5], name='b')

        # inner
        jidx, lidx, ridx = midx.join(idx, how='inner', return_indexers=True)
        exp_idx = pd.MultiIndex.from_product([np.arange(4), [1, 2]], names=['a', 'b'])
        exp_lidx = np.array([1, 2, 5, 6, 9, 10, 13, 14])
        exp_ridx = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        self.assert_index_equal(jidx, exp_idx)
        self.assert_numpy_array_equal(lidx, exp_lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='inner', return_indexers=True)
        self.assert_index_equal(jidx, exp_idx)
        self.assert_numpy_array_equal(lidx, exp_lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)

        # keep MultiIndex
        jidx, lidx, ridx = midx.join(idx, how='left', return_indexers=True)
        exp_ridx = np.array([-1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1])
        self.assert_index_equal(jidx, midx)
        self.assertIsNone(lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)
        # flip
        jidx, ridx, lidx = idx.join(midx, how='right', return_indexers=True)
        self.assert_index_equal(jidx, midx)
        self.assertIsNone(lidx)
        self.assert_numpy_array_equal(ridx, exp_ridx)

    def test_reindex(self):
        result, indexer = self.index.reindex(list(self.index[:4]))
        tm.assertIsInstance(result, MultiIndex)
        self.check_level_names(result, self.index[:4].names)

        result, indexer = self.index.reindex(list(self.index))
        tm.assertIsInstance(result, MultiIndex)
        self.assertIsNone(indexer)
        self.check_level_names(result, self.index.names)

    def test_reindex_level(self):
        idx = Index(['one'])

        target, indexer = self.index.reindex(idx, level='second')
        target2, indexer2 = idx.reindex(self.index, level='second')

        exp_index = self.index.join(idx, level='second', how='right')
        exp_index2 = self.index.join(idx, level='second', how='left')

        self.assertTrue(target.equals(exp_index))
        exp_indexer = np.array([0, 2, 4])
        tm.assert_numpy_array_equal(indexer, exp_indexer)

        self.assertTrue(target2.equals(exp_index2))
        exp_indexer2 = np.array([0, -1, 0, -1, 0, -1])
        tm.assert_numpy_array_equal(indexer2, exp_indexer2)

        assertRaisesRegexp(TypeError, "Fill method not supported",
                           self.index.reindex, self.index, method='pad',
                           level='second')

        assertRaisesRegexp(TypeError, "Fill method not supported",
                           idx.reindex, idx, method='bfill', level='first')

    def test_duplicates(self):
        self.assertFalse(self.index.has_duplicates)
        self.assertTrue(self.index.append(self.index).has_duplicates)

        index = MultiIndex(levels=[[0, 1], [0, 1, 2]],
                           labels=[[0, 0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 0, 1, 2]])
        self.assertTrue(index.has_duplicates)

        # GH 9075
        t = [(u('x'), u('out'), u('z'), 5, u('y'), u('in'), u('z'), 169),
             (u('x'), u('out'), u('z'), 7, u('y'), u('in'), u('z'), 119),
             (u('x'), u('out'), u('z'), 9, u('y'), u('in'), u('z'), 135),
             (u('x'), u('out'), u('z'), 13, u('y'), u('in'), u('z'), 145),
             (u('x'), u('out'), u('z'), 14, u('y'), u('in'), u('z'), 158),
             (u('x'), u('out'), u('z'), 16, u('y'), u('in'), u('z'), 122),
             (u('x'), u('out'), u('z'), 17, u('y'), u('in'), u('z'), 160),
             (u('x'), u('out'), u('z'), 18, u('y'), u('in'), u('z'), 180),
             (u('x'), u('out'), u('z'), 20, u('y'), u('in'), u('z'), 143),
             (u('x'), u('out'), u('z'), 21, u('y'), u('in'), u('z'), 128),
             (u('x'), u('out'), u('z'), 22, u('y'), u('in'), u('z'), 129),
             (u('x'), u('out'), u('z'), 25, u('y'), u('in'), u('z'), 111),
             (u('x'), u('out'), u('z'), 28, u('y'), u('in'), u('z'), 114),
             (u('x'), u('out'), u('z'), 29, u('y'), u('in'), u('z'), 121),
             (u('x'), u('out'), u('z'), 31, u('y'), u('in'), u('z'), 126),
             (u('x'), u('out'), u('z'), 32, u('y'), u('in'), u('z'), 155),
             (u('x'), u('out'), u('z'), 33, u('y'), u('in'), u('z'), 123),
             (u('x'), u('out'), u('z'), 12, u('y'), u('in'), u('z'), 144)]

        index = pd.MultiIndex.from_tuples(t)
        self.assertFalse(index.has_duplicates)

        # handle int64 overflow if possible
        def check(nlevels, with_nulls):
            labels = np.tile(np.arange(500), 2)
            level = np.arange(500)

            if with_nulls:  # inject some null values
                labels[500] = -1  # common nan value
                labels = list(labels.copy() for i in range(nlevels))
                for i in range(nlevels):
                    labels[i][500 + i - nlevels // 2 ] = -1

                labels += [np.array([-1, 1]).repeat(500)]
            else:
                labels = [labels] * nlevels + [np.arange(2).repeat(500)]

            levels = [level] * nlevels + [[0, 1]]

            # no dups
            index = MultiIndex(levels=levels, labels=labels)
            self.assertFalse(index.has_duplicates)

            # with a dup
            if with_nulls:
                f = lambda a: np.insert(a, 1000, a[0])
                labels = list(map(f, labels))
                index = MultiIndex(levels=levels, labels=labels)
            else:
                values = index.values.tolist()
                index = MultiIndex.from_tuples(values + [values[0]])

            self.assertTrue(index.has_duplicates)

        # no overflow
        check(4, False)
        check(4, True)

        # overflow possible
        check(8, False)
        check(8, True)

        # GH 9125
        n, k = 200, 5000
        levels = [np.arange(n), tm.makeStringIndex(n), 1000 + np.arange(n)]
        labels = [np.random.choice(n, k * n) for lev in levels]
        mi = MultiIndex(levels=levels, labels=labels)

        for keep in ['first', 'last', False]:
            left = mi.duplicated(keep=keep)
            right = pd.lib.duplicated(mi.values, keep=keep)
            tm.assert_numpy_array_equal(left, right)

        # GH5873
        for a in [101, 102]:
            mi = MultiIndex.from_arrays([[101, a], [3.5, np.nan]])
            self.assertFalse(mi.has_duplicates)
            self.assertEqual(mi.get_duplicates(), [])
            tm.assert_numpy_array_equal(mi.duplicated(), np.zeros(2, dtype='bool'))

        for n in range(1, 6):  # 1st level shape
            for m in range(1, 5):  # 2nd level shape
                # all possible unique combinations, including nan
                lab = product(range(-1, n), range(-1, m))
                mi = MultiIndex(levels=[list('abcde')[:n], list('WXYZ')[:m]],
                                labels=np.random.permutation(list(lab)).T)
                self.assertEqual(len(mi), (n + 1) * (m + 1))
                self.assertFalse(mi.has_duplicates)
                self.assertEqual(mi.get_duplicates(), [])
                tm.assert_numpy_array_equal(mi.duplicated(),
                                            np.zeros(len(mi), dtype='bool'))

    def test_duplicate_meta_data(self):
        # GH 10115
        index = MultiIndex(levels=[[0, 1], [0, 1, 2]],
                           labels=[[0, 0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 0, 1, 2]])
        for idx in [index,
                    index.set_names([None, None]),
                    index.set_names([None, 'Num']),
                    index.set_names(['Upper','Num']),
                   ]:
            self.assertTrue(idx.has_duplicates)
            self.assertEqual(idx.drop_duplicates().names, idx.names)

    def test_tolist(self):
        result = self.index.tolist()
        exp = list(self.index.values)
        self.assertEqual(result, exp)

    def test_repr_with_unicode_data(self):
        with pd.core.config.option_context("display.encoding",'UTF-8'):
            d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
            index = pd.DataFrame(d).set_index(["a", "b"]).index
            self.assertFalse("\\u" in repr(index))  # we don't want unicode-escaped

    def test_repr_roundtrip(self):

        mi = MultiIndex.from_product([list('ab'),range(3)],names=['first','second'])
        str(mi)

        if compat.PY3:
            tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
        else:
            result = eval(repr(mi))
            # string coerces to unicode
            tm.assert_index_equal(result, mi, exact=False)
            self.assertEqual(mi.get_level_values('first').inferred_type, 'string')
            self.assertEqual(result.get_level_values('first').inferred_type, 'unicode')

        mi_u = MultiIndex.from_product([list(u'ab'),range(3)],names=['first','second'])
        result = eval(repr(mi_u))
        tm.assert_index_equal(result, mi_u, exact=True)

        # formatting
        if compat.PY3:
            str(mi)
        else:
            compat.text_type(mi)

        # long format
        mi = MultiIndex.from_product([list('abcdefg'),range(10)],names=['first','second'])
        result = str(mi)

        if compat.PY3:
            tm.assert_index_equal(eval(repr(mi)), mi, exact=True)
        else:
            result = eval(repr(mi))
            # string coerces to unicode
            tm.assert_index_equal(result, mi, exact=False)
            self.assertEqual(mi.get_level_values('first').inferred_type, 'string')
            self.assertEqual(result.get_level_values('first').inferred_type, 'unicode')

        mi = MultiIndex.from_product([list(u'abcdefg'),range(10)],names=['first','second'])
        result = eval(repr(mi_u))
        tm.assert_index_equal(result, mi_u, exact=True)

    def test_str(self):
        # tested elsewhere
        pass

    def test_unicode_string_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if compat.PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        d = {"a": [u("\u05d0"), 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}
        idx = pd.DataFrame(d).set_index(["a", "b"]).index

        if compat.PY3:
            bytes(idx)
        else:
            str(idx)

    def test_slice_keep_name(self):
        x = MultiIndex.from_tuples([('a', 'b'), (1, 2), ('c', 'd')],
                                   names=['x', 'y'])
        self.assertEqual(x[1:].names, x.names)

    def test_isnull_behavior(self):
        # should not segfault GH5123
        # NOTE: if MI representation changes, may make sense to allow
        # isnull(MI)
        with tm.assertRaises(NotImplementedError):
            pd.isnull(self.index)

    def test_level_setting_resets_attributes(self):
        ind = MultiIndex.from_arrays([
            ['A', 'A', 'B', 'B', 'B'],
            [1, 2, 1, 2, 3]])
        assert ind.is_monotonic
        ind.set_levels([['A', 'B', 'A', 'A', 'B'], [2, 1, 3, -2, 5]],
                       inplace=True)
        # if this fails, probably didn't reset the cache correctly.
        assert not ind.is_monotonic

    def test_isin(self):
        values = [('foo', 2), ('bar', 3), ('quux', 4)]

        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'],
                                      np.arange(4)])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = MultiIndex.from_arrays([[], []])
        result = idx.isin(values)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.dtype, np.bool_)

    def test_isin_nan(self):
        idx = MultiIndex.from_arrays([['foo', 'bar'], [1.0, np.nan]])
        tm.assert_numpy_array_equal(idx.isin([('bar', np.nan)]),
                                      [False, False])
        tm.assert_numpy_array_equal(idx.isin([('bar', float('nan'))]),
                                      [False, False])

    def test_isin_level_kwarg(self):
        idx = MultiIndex.from_arrays([['qux', 'baz', 'foo', 'bar'],
                                      np.arange(4)])

        vals_0 = ['foo', 'bar', 'quux']
        vals_1 = [2, 3, 10]

        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=0))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level=-2))

        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=1))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level=-1))

        self.assertRaises(IndexError, idx.isin, vals_0, level=5)
        self.assertRaises(IndexError, idx.isin, vals_0, level=-5)

        self.assertRaises(KeyError, idx.isin, vals_0, level=1.0)
        self.assertRaises(KeyError, idx.isin, vals_1, level=-1.0)
        self.assertRaises(KeyError, idx.isin, vals_1, level='A')

        idx.names = ['A', 'B']
        tm.assert_numpy_array_equal(expected, idx.isin(vals_0, level='A'))
        tm.assert_numpy_array_equal(expected, idx.isin(vals_1, level='B'))

        self.assertRaises(KeyError, idx.isin, vals_1, level='C')

    def test_reindex_preserves_names_when_target_is_list_or_ndarray(self):
        # GH6552
        idx = self.index.copy()
        target = idx.copy()
        idx.names = target.names = [None, None]

        other_dtype = pd.MultiIndex.from_product([[1, 2], [3, 4]])

        # list & ndarray cases
        self.assertEqual(idx.reindex([])[0].names, [None, None])
        self.assertEqual(idx.reindex(np.array([]))[0].names, [None, None])
        self.assertEqual(idx.reindex(target.tolist())[0].names, [None, None])
        self.assertEqual(idx.reindex(target.values)[0].names, [None, None])
        self.assertEqual(idx.reindex(other_dtype.tolist())[0].names, [None, None])
        self.assertEqual(idx.reindex(other_dtype.values)[0].names, [None, None])

        idx.names = ['foo', 'bar']
        self.assertEqual(idx.reindex([])[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(np.array([]))[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(target.tolist())[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(target.values)[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(other_dtype.tolist())[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex(other_dtype.values)[0].names, ['foo', 'bar'])

    def test_reindex_lvl_preserves_names_when_target_is_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']],
                                         names=['foo', 'bar'])
        self.assertEqual(idx.reindex([], level=0)[0].names, ['foo', 'bar'])
        self.assertEqual(idx.reindex([], level=1)[0].names, ['foo', 'bar'])

    def test_reindex_lvl_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
        self.assertEqual(idx.reindex([], level=0)[0].levels[0].dtype.type,
                         np.int64)
        self.assertEqual(idx.reindex([], level=1)[0].levels[1].dtype.type,
                         np.object_)

    def test_groupby(self):
        groups = self.index.groupby(np.array([1, 1, 1, 2, 2, 2]))
        labels = self.index.get_values().tolist()
        exp = {1: labels[:3], 2: labels[3:]}
        tm.assert_dict_equal(groups, exp)

        # GH5620
        groups = self.index.groupby(self.index)
        exp = dict((key, [key]) for key in self.index)
        tm.assert_dict_equal(groups, exp)

    def test_index_name_retained(self):
        # GH9857
        result = pd.DataFrame({'x': [1, 2, 6],
                               'y': [2, 2, 8],
                               'z': [-5, 0, 5]})
        result = result.set_index('z')
        result.loc[10] = [9, 10]
        df_expected = pd.DataFrame({'x': [1, 2, 6, 9],
                                    'y': [2, 2, 8, 10],
                                    'z': [-5, 0, 5, 10]})
        df_expected = df_expected.set_index('z')
        tm.assert_frame_equal(result, df_expected)

    def test_equals_operator(self):
        # GH9785
        self.assertTrue((self.index == self.index).all())


def test_get_combined_index():
    from pandas.core.index import _get_combined_index
    result = _get_combined_index([])
    assert(result.equals(Index([])))



if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
