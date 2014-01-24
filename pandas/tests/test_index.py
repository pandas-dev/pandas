# pylint: disable=E1101,E1103,W0232

from datetime import datetime, timedelta
from pandas.compat import range, lrange, lzip, u, zip
import operator
import pickle
import re
import nose
import warnings
import os

import numpy as np
from numpy.testing import assert_array_equal

from pandas.core.index import (Index, Float64Index, Int64Index, MultiIndex,
                               InvalidIndexError)
from pandas.tseries.index import DatetimeIndex
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.util.testing import (assert_almost_equal, assertRaisesRegexp,
                                 assert_copy)
from pandas import compat
from pandas.compat import long

import pandas.util.testing as tm
import pandas.core.config as cf

from pandas.tseries.index import _to_m8
import pandas.tseries.offsets as offsets

import pandas as pd
from pandas.lib import Timestamp

from pandas import _np_version_under1p7

def _skip_if_need_numpy_1_7():
    if _np_version_under1p7:
        raise nose.SkipTest('numpy >= 1.7 required')

class TestIndex(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(
            unicodeIndex = tm.makeUnicodeIndex(100),
            strIndex = tm.makeStringIndex(100),
            dateIndex = tm.makeDateIndex(100),
            intIndex = tm.makeIntIndex(100),
            floatIndex = tm.makeFloatIndex(100),
            empty = Index([]),
            tuples = Index(lzip(['foo', 'bar', 'baz'], [1, 2, 3])),
        )
        for name, ind in self.indices.items():
            setattr(self, name, ind)

    def test_wrong_number_names(self):
        def testit(ind):
            ind.names = ["apple", "banana", "carrot"]

        for ind in self.indices.values():
            assertRaisesRegexp(ValueError, "^Length", testit, ind)

    def test_set_name_methods(self):
        new_name = "This is the new name for this index"
        indices = (self.dateIndex, self.intIndex, self.unicodeIndex,
                   self.empty)
        for ind in indices:
            original_name = ind.name
            new_ind = ind.set_names([new_name])
            self.assertEqual(new_ind.name, new_name)
            self.assertEqual(ind.name, original_name)
            res = ind.rename(new_name, inplace=True)
            # should return None
            self.assert_(res is None)
            self.assertEqual(ind.name, new_name)
            self.assertEqual(ind.names, [new_name])
            with assertRaisesRegexp(TypeError, "list-like"):
                # should still fail even if it would be the right length
                ind.set_names("a")
        # rename in place just leaves tuples and other containers alone
        name = ('A', 'B')
        ind = self.intIndex
        ind.rename(name, inplace=True)
        self.assertEqual(ind.name, name)
        self.assertEqual(ind.names, [name])

    def test_hash_error(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(self.strIndex).__name__):
            hash(self.strIndex)

    def test_new_axis(self):
        new_index = self.dateIndex[None, :]
        self.assert_(new_index.ndim == 2)
        tm.assert_isinstance(new_index, np.ndarray)

    def test_copy_and_deepcopy(self):
        from copy import copy, deepcopy

        for func in (copy, deepcopy):
            idx_copy = func(self.strIndex)
            self.assert_(idx_copy is not self.strIndex)
            self.assert_(idx_copy.equals(self.strIndex))

        new_copy = self.strIndex.copy(deep=True, name="banana")
        self.assertEqual(new_copy.name, "banana")
        new_copy2 = self.intIndex.copy(dtype=int)
        self.assertEqual(new_copy2.dtype.kind, 'i')

    def test_duplicates(self):
        idx = Index([0, 0, 0])
        self.assert_(not idx.is_unique)

    def test_sort(self):
        self.assertRaises(TypeError, self.strIndex.sort)

    def test_mutability(self):
        self.assertRaises(TypeError, self.strIndex.__setitem__, 0, 'foo')

    def test_constructor(self):
        # regular instance creation
        tm.assert_contains_all(self.strIndex, self.strIndex)
        tm.assert_contains_all(self.dateIndex, self.dateIndex)

        # casting
        arr = np.array(self.strIndex)
        index = arr.view(Index)
        tm.assert_contains_all(arr, index)
        self.assert_(np.array_equal(self.strIndex, index))

        # copy
        arr = np.array(self.strIndex)
        index = Index(arr, copy=True, name='name')
        tm.assert_isinstance(index, Index)
        self.assert_(index.name == 'name')
        assert_array_equal(arr, index)
        arr[0] = "SOMEBIGLONGSTRING"
        self.assertNotEqual(index[0], "SOMEBIGLONGSTRING")

        # what to do here?
        # arr = np.array(5.)
        # self.assertRaises(Exception, arr.view, Index)

    def test_constructor_corner(self):
        # corner case
        self.assertRaises(TypeError, Index, 0)

    def test_index_ctor_infer_periodindex(self):
        from pandas import period_range, PeriodIndex
        xp = period_range('2012-1-1', freq='M', periods=3)
        rs = Index(xp)
        assert_array_equal(rs, xp)
        tm.assert_isinstance(rs, PeriodIndex)

    def test_copy(self):
        i = Index([], name='Foo')
        i_copy = i.copy()
        self.assert_(i_copy.name == 'Foo')

    def test_view(self):
        i = Index([], name='Foo')
        i_view = i.view()
        self.assert_(i_view.name == 'Foo')

    def test_astype(self):
        casted = self.intIndex.astype('i8')

        # it works!
        casted.get_loc(5)

        # pass on name
        self.intIndex.name = 'foobar'
        casted = self.intIndex.astype('i8')
        self.assertEqual(casted.name, 'foobar')

    def test_compat(self):
        self.strIndex.tolist()

    def test_equals(self):
        # same
        self.assert_(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'c'])))

        # different length
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b'])))

        # same length, different values
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'd'])))

        # Must also be an Index
        self.assertFalse(Index(['a', 'b', 'c']).equals(['a', 'b', 'c']))

    def test_identical(self):

        # index
        i1 = Index(['a', 'b', 'c'])
        i2 = Index(['a', 'b', 'c'])

        self.assert_(i1.identical(i2))

        i1 = i1.rename('foo')
        self.assert_(i1.equals(i2))
        self.assert_(not i1.identical(i2))

        i2 = i2.rename('foo')
        self.assert_(i1.identical(i2))

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
        self.assertTrue(ind.is_(ind.view().base))
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
        self.assert_(self.dateIndex.asof(d) is d)
        self.assert_(np.isnan(self.dateIndex.asof(d - timedelta(1))))

        d = self.dateIndex[-1]
        self.assert_(self.dateIndex.asof(d + timedelta(1)) == d)

        d = self.dateIndex[0].to_datetime()
        tm.assert_isinstance(self.dateIndex.asof(d), Timestamp)

    def test_nanosecond_index_access(self):
        _skip_if_need_numpy_1_7()

        s = Series([Timestamp('20130101')]).values.view('i8')[0]
        r = DatetimeIndex([s + 50 + i for i in range(100)])
        x = Series(np.random.randn(100), index=r)

        first_value = x.asof(x.index[0])

        # this does not yet work, as parsing strings is done via dateutil
        #self.assertEqual(first_value, x['2013-01-01 00:00:00.000000050+0000'])

        self.assertEqual(first_value, x[Timestamp(np.datetime64('2013-01-01 00:00:00.000000050+0000', 'ns'))])

    def test_argsort(self):
        result = self.strIndex.argsort()
        expected = np.array(self.strIndex).argsort()
        self.assert_(np.array_equal(result, expected))

    def test_comparators(self):
        index = self.dateIndex
        element = index[len(index) // 2]
        element = _to_m8(element)

        arr = np.array(index)

        def _check(op):
            arr_result = op(arr, element)
            index_result = op(index, element)

            tm.assert_isinstance(index_result, np.ndarray)
            self.assert_(not isinstance(index_result, Index))
            self.assert_(np.array_equal(arr_result, index_result))

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

    def test_getitem(self):
        arr = np.array(self.dateIndex)
        exp = self.dateIndex[5]
        exp = _to_m8(exp)

        self.assertEquals(exp, arr[5])

    def test_shift(self):
        shifted = self.dateIndex.shift(0, timedelta(1))
        self.assert_(shifted is self.dateIndex)

        shifted = self.dateIndex.shift(5, timedelta(1))
        self.assert_(np.array_equal(shifted, self.dateIndex + timedelta(5)))

        shifted = self.dateIndex.shift(1, 'B')
        self.assert_(np.array_equal(shifted, self.dateIndex + offsets.BDay()))

        shifted.name = 'shifted'
        self.assertEqual(shifted.name, shifted.shift(1, 'D').name)

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)

        self.assert_(tm.equalContents(intersect, second))

        # Corner cases
        inter = first.intersection(first)
        self.assert_(inter is first)

        # non-iterable input
        assertRaisesRegexp(TypeError, "iterable", first.intersection, 0.5)

    def test_union(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        everything = self.strIndex[:20]
        union = first.union(second)
        self.assert_(tm.equalContents(union, everything))

        # Corner cases
        union = first.union(first)
        self.assert_(union is first)

        union = first.union([])
        self.assert_(union is first)

        union = Index([]).union(first)
        self.assert_(union is first)

        # non-iterable input
        assertRaisesRegexp(TypeError, "iterable", first.union, 0.5)

        # preserve names
        first.name = 'A'
        second.name = 'A'
        union = first.union(second)
        self.assert_(union.name == 'A')

        second.name = 'B'
        union = first.union(second)
        self.assert_(union.name is None)

    def test_add(self):
        firstCat = self.strIndex + self.dateIndex
        secondCat = self.strIndex + self.strIndex

        if self.dateIndex.dtype == np.object_:
            appended = np.append(self.strIndex, self.dateIndex)
        else:
            appended = np.append(self.strIndex, self.dateIndex.astype('O'))

        self.assert_(tm.equalContents(firstCat, appended))
        self.assert_(tm.equalContents(secondCat, self.strIndex))
        tm.assert_contains_all(self.strIndex, firstCat)
        tm.assert_contains_all(self.strIndex, secondCat)
        tm.assert_contains_all(self.dateIndex, firstCat)

    def test_append_multiple(self):
        index = Index(['a', 'b', 'c', 'd', 'e', 'f'])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        self.assert_(result.equals(index))

        # empty
        result = index.append([])
        self.assert_(result.equals(index))

    def test_append_empty_preserve_name(self):
        left = Index([], name='foo')
        right = Index([1, 2, 3], name='foo')

        result = left.append(right)
        self.assert_(result.name == 'foo')

        left = Index([], name='foo')
        right = Index([1, 2, 3], name='bar')

        result = left.append(right)
        self.assert_(result.name is None)

    def test_add_string(self):
        # from bug report
        index = Index(['a', 'b', 'c'])
        index2 = index + 'foo'

        self.assert_('a' not in index2)
        self.assert_('afoo' in index2)

    def test_iadd_string(self):
        index = pd.Index(['a', 'b', 'c'])
        # doesn't fail test unless there is a check before `+=`
        self.assert_('a' in index)

        index += '_x'
        self.assert_('a_x' in index)

    def test_diff(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]
        first.name = 'name'
        # different names
        result = first - second

        self.assert_(tm.equalContents(result, answer))
        self.assertEqual(result.name, None)

        # same names
        second.name = 'name'
        result = first - second
        self.assertEqual(result.name, 'name')

        # with empty
        result = first.diff([])
        self.assert_(tm.equalContents(result, first))
        self.assertEqual(result.name, first.name)

        # with everythin
        result = first.diff(first)
        self.assert_(len(result) == 0)
        self.assertEqual(result.name, first.name)

        # non-iterable input
        assertRaisesRegexp(TypeError, "iterable", first.diff, 0.5)

    def test_pickle(self):
        def testit(index):
            pickled = pickle.dumps(index)
            unpickled = pickle.loads(pickled)

            tm.assert_isinstance(unpickled, Index)
            self.assert_(np.array_equal(unpickled, index))
            self.assertEquals(unpickled.name, index.name)

            # tm.assert_dict_equal(unpickled.indexMap, index.indexMap)

        testit(self.strIndex)
        self.strIndex.name = 'foo'
        testit(self.strIndex)

        testit(self.dateIndex)

    def test_is_numeric(self):
        self.assert_(not self.dateIndex.is_numeric())
        self.assert_(not self.strIndex.is_numeric())
        self.assert_(self.intIndex.is_numeric())
        self.assert_(self.floatIndex.is_numeric())

    def test_is_all_dates(self):
        self.assert_(self.dateIndex.is_all_dates)
        self.assert_(not self.strIndex.is_all_dates)
        self.assert_(not self.intIndex.is_all_dates)

    def test_summary(self):
        self._check_method_works(Index.summary)
        # GH3869
        ind = Index(['{other}%s', "~:{range}:0"], name='A')
        result = ind.summary()
        # shouldn't be formatted accidentally.
        self.assert_('~:{range}:0' in result)
        self.assert_('{other}%s' in result)

    def test_format(self):
        self._check_method_works(Index.format)

        index = Index([datetime.now()])
        formatted = index.format()
        expected = [str(index[0])]
        self.assertEquals(formatted, expected)

        # 2845
        index = Index([1, 2.0+3.0j, np.nan])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEquals(formatted, expected)

        # is this really allowed?
        index = Index([1, 2.0+3.0j, None])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEquals(formatted, expected)

        self.strIndex[:0].format()

    def test_format_with_name_time_info(self):
        # bug I fixed 12/20/2011
        inc = timedelta(hours=4)
        dates = Index([dt + inc for dt in self.dateIndex], name='something')

        formatted = dates.format(name=True)
        self.assert_(formatted[0] == 'something')

    def test_format_datetime_with_time(self):
        t = Index([datetime(2012, 2, 7), datetime(2012, 2, 7, 23)])

        result = t.format()
        expected = ['2012-02-07 00:00:00', '2012-02-07 23:00:00']
        self.assert_(len(result) == 2)
        self.assertEquals(result, expected)

    def test_format_none(self):
        values = ['a', 'b', 'c', None]

        idx = Index(values)
        idx.format()
        self.assert_(idx[3] is None)

    def test_take(self):
        indexer = [4, 3, 0, 2]
        result = self.dateIndex.take(indexer)
        expected = self.dateIndex[indexer]
        self.assert_(result.equals(expected))

    def _check_method_works(self, method):
        method(self.empty)
        method(self.dateIndex)
        method(self.unicodeIndex)
        method(self.strIndex)
        method(self.intIndex)
        method(self.tuples)

    def test_get_indexer(self):
        idx1 = Index([1, 2, 3, 4, 5])
        idx2 = Index([2, 4, 6])

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, [1, 3, -1])

        r1 = idx2.get_indexer(idx1, method='pad')
        assert_almost_equal(r1, [-1, 0, 0, 1, 1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        assert_almost_equal(r1, [0, 0, 1, 1, 2])

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

    def test_slice_locs(self):
        idx = Index([0, 1, 2, 5, 6, 7, 9, 10])
        n = len(idx)

        self.assertEquals(idx.slice_locs(start=2), (2, n))
        self.assertEquals(idx.slice_locs(start=3), (3, n))
        self.assertEquals(idx.slice_locs(3, 8), (3, 6))
        self.assertEquals(idx.slice_locs(5, 10), (3, n))
        self.assertEquals(idx.slice_locs(end=8), (0, 6))
        self.assertEquals(idx.slice_locs(end=9), (0, 7))

        idx2 = idx[::-1]
        self.assertRaises(KeyError, idx2.slice_locs, 8, 2)
        self.assertRaises(KeyError, idx2.slice_locs, 7, 3)

    def test_slice_locs_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        rs = idx.slice_locs('a', 'd')
        self.assert_(rs == (0, 6))

        rs = idx.slice_locs(end='d')
        self.assert_(rs == (0, 6))

        rs = idx.slice_locs('a', 'c')
        self.assert_(rs == (0, 4))

        rs = idx.slice_locs('b', 'd')
        self.assert_(rs == (2, 6))

    def test_drop(self):
        n = len(self.strIndex)

        dropped = self.strIndex.drop(self.strIndex[lrange(5, 10)])
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        self.assert_(dropped.equals(expected))

        self.assertRaises(ValueError, self.strIndex.drop, ['foo', 'bar'])

        dropped = self.strIndex.drop(self.strIndex[0])
        expected = self.strIndex[1:]
        self.assert_(dropped.equals(expected))

        ser = Index([1, 2, 3])
        dropped = ser.drop(1)
        expected = Index([2, 3])
        self.assert_(dropped.equals(expected))

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
        self.assert_(int_idx.ndim == 1)
        self.assert_(int_idx.equals(expected))

        # union broken
        union_idx = idx1.union(idx2)
        expected = pandas.Index(sorted(set(idx1) | set(idx2)))
        self.assert_(union_idx.ndim == 1)
        self.assert_(union_idx.equals(expected))

    def test_is_monotonic_incomparable(self):
        index = Index([5, datetime.now(), 7])
        self.assert_(not index.is_monotonic)

    def test_get_set_value(self):
        values = np.random.randn(100)
        date = self.dateIndex[67]

        assert_almost_equal(self.dateIndex.get_value(values, date),
                            values[67])

        self.dateIndex.set_value(values, date, 10)
        self.assertEquals(values[67], 10)

    def test_isin(self):
        values = ['foo', 'bar']

        idx = Index(['qux', 'baz', 'foo', 'bar'])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        self.assert_(np.array_equal(result, expected))

        # empty, return dtype bool
        idx = Index([])
        result = idx.isin(values)
        self.assert_(len(result) == 0)
        self.assert_(result.dtype == np.bool_)

    def test_boolean_cmp(self):
        values = [1, 2, 3, 4]

        idx = Index(values)
        res = (idx == values)

        self.assert_(res.all())
        self.assert_(res.dtype == 'bool')
        self.assert_(not isinstance(res, Index))

    def test_get_level_values(self):
        result = self.strIndex.get_level_values(0)
        self.assert_(result.equals(self.strIndex))

    def test_slice_keep_name(self):
        idx = Index(['a', 'b'], name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

    def test_join_self(self):
        indices = 'unicode', 'str', 'date', 'int', 'float'
        kinds = 'outer', 'inner', 'left', 'right'
        for index_kind in indices:
            for kind in kinds:
                res = getattr(self, '{0}Index'.format(index_kind))
                joined = res.join(res, how=kind)
                self.assert_(res is joined)


class TestFloat64Index(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.mixed = Float64Index([1.5, 2, 3, 4, 5])
        self.float = Float64Index(np.arange(5) * 2.5)

    def test_hash_error(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(self.float).__name__):
            hash(self.float)

    def test_repr_roundtrip(self):
        for ind in (self.mixed, self.float):
            tm.assert_index_equal(eval(repr(ind)), ind)

    def check_is_index(self, i):
        self.assert_(isinstance(i, Index) and not isinstance(i, Float64Index))

    def check_coerce(self, a, b, is_float_index=True):
        self.assert_(a.equals(b))
        if is_float_index:
            self.assert_(isinstance(b, Float64Index))
        else:
            self.check_is_index(b)

    def test_constructor(self):

        # explicit construction
        index = Float64Index([1,2,3,4,5])
        self.assert_(isinstance(index, Float64Index))
        self.assert_((index.values == np.array([1,2,3,4,5],dtype='float64')).all())
        index = Float64Index(np.array([1,2,3,4,5]))
        self.assert_(isinstance(index, Float64Index))
        index = Float64Index([1.,2,3,4,5])
        self.assert_(isinstance(index, Float64Index))
        index = Float64Index(np.array([1.,2,3,4,5]))
        self.assert_(isinstance(index, Float64Index))
        self.assert_(index.dtype == object)

        index = Float64Index(np.array([1.,2,3,4,5]),dtype=np.float32)
        self.assert_(isinstance(index, Float64Index))
        self.assert_(index.dtype == object)

        index = Float64Index(np.array([1,2,3,4,5]),dtype=np.float32)
        self.assert_(isinstance(index, Float64Index))
        self.assert_(index.dtype == object)

        # nan handling
        result = Float64Index([np.nan, np.nan])
        self.assert_(pd.isnull(result.values).all())
        result = Float64Index(np.array([np.nan]))
        self.assert_(pd.isnull(result.values).all())
        result = Index(np.array([np.nan]))
        self.assert_(pd.isnull(result.values).all())

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
        self.assert_(result.equals(self.float))
        self.assert_(self.float.equals(result))
        self.check_is_index(result)

        i = self.mixed.copy()
        i.name = 'foo'
        result = i.astype(object)
        self.assert_(result.equals(i))
        self.assert_(i.equals(result))
        self.check_is_index(result)


class TestInt64Index(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.index = Int64Index(np.arange(0, 20, 2))

    def test_too_many_names(self):
        def testit():
            self.index.names = ["roger", "harold"]
        assertRaisesRegexp(ValueError, "^Length", testit)

    def test_constructor(self):
        # pass list, coerce fine
        index = Int64Index([-5, 0, 1, 2])
        expected = np.array([-5, 0, 1, 2], dtype=np.int64)
        self.assert_(np.array_equal(index, expected))

        # from iterable
        index = Int64Index(iter([-5, 0, 1, 2]))
        self.assert_(np.array_equal(index, expected))

        # scalar raise Exception
        self.assertRaises(TypeError, Int64Index, 5)

        # copy
        arr = self.index.values
        new_index = Int64Index(arr, copy=True)
        self.assert_(np.array_equal(new_index, self.index))
        val = arr[0] + 3000
        # this should not change index
        arr[0] = val
        self.assertNotEqual(new_index[0], val)

    def test_constructor_corner(self):
        arr = np.array([1, 2, 3, 4], dtype=object)
        index = Int64Index(arr)
        self.assert_(index.values.dtype == np.int64)
        self.assert_(index.equals(arr))

        # preventing casting
        arr = np.array([1, '2', 3, '4'], dtype=object)
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr)

        arr_with_floats = [0, 2, 3, 4, 5, 1.25, 3, -1]
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr_with_floats)

    def test_hash_error(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(self.index).__name__):
            hash(self.index)

    def test_copy(self):
        i = Int64Index([], name='Foo')
        i_copy = i.copy()
        self.assert_(i_copy.name == 'Foo')

    def test_view(self):
        i = Int64Index([], name='Foo')
        i_view = i.view()
        self.assert_(i_view.name == 'Foo')

    def test_coerce_list(self):
        # coerce things
        arr = Index([1, 2, 3, 4])
        tm.assert_isinstance(arr, Int64Index)

        # but not if explicit dtype passed
        arr = Index([1, 2, 3, 4], dtype=object)
        tm.assert_isinstance(arr, Index)

    def test_dtype(self):
        self.assert_(self.index.dtype == np.int64)

    def test_is_monotonic(self):
        self.assert_(self.index.is_monotonic)

        index = Int64Index([4, 3, 2, 1])
        self.assert_(not index.is_monotonic)

    def test_equals(self):
        same_values = Index(self.index, dtype=object)
        self.assert_(self.index.equals(same_values))
        self.assert_(same_values.equals(self.index))

    def test_identical(self):

        i = self.index.copy()
        same_values = Index(i, dtype=object)
        self.assert_(i.identical(same_values))

        i = self.index.copy()
        i = i.rename('foo')
        same_values = Index(i, dtype=object)
        self.assert_(same_values.identical(self.index))

        self.assertFalse(i.identical(self.index))
        self.assert_(Index(same_values, name='foo').identical(i))

    def test_get_indexer(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target)
        expected = np.array([0, -1, 1, -1, 2, -1, 3, -1, 4, -1])
        self.assert_(np.array_equal(indexer, expected))

    def test_get_indexer_pad(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='pad')
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        self.assert_(np.array_equal(indexer, expected))

    def test_get_indexer_backfill(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='backfill')
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5])
        self.assert_(np.array_equal(indexer, expected))

    def test_join_outer(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        # guarantee of sortedness
        res, lidx, ridx = self.index.join(other, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other, how='outer')
        self.assert_(res.equals(noidx_res))

        eres = Int64Index([0, 1, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 25])
        elidx = np.array([0, -1, 1, 2, -1, 3, -1, 4, 5, 6, 7, 8, 9, -1],
                         dtype=np.int64)
        eridx = np.array([-1, 3, 4, -1, 5, -1, 0, -1, -1, 1, -1, -1, -1, 2],
                         dtype=np.int64)

        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other_mono, how='outer')
        self.assert_(res.equals(noidx_res))

        eridx = np.array([-1, 0, 1, -1, 2, -1, 3, -1, -1, 4, -1, -1, -1, 5],
                         dtype=np.int64)
        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

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

        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='inner',
                                          return_indexers=True)

        res2 = self.index.intersection(other_mono)
        self.assert_(res.equals(res2))

        eridx = np.array([1, 4])
        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

    def test_join_left(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='left',
                                          return_indexers=True)
        eres = self.index
        eridx = np.array([-1, 4, -1, -1, -1, -1, 1, -1, -1, -1],
                         dtype=np.int64)

        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(lidx is None)
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='left',
                                          return_indexers=True)
        eridx = np.array([-1, 1, -1, -1, -1, -1, 4, -1, -1, -1],
                         dtype=np.int64)
        tm.assert_isinstance(res, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(lidx is None)
        self.assert_(np.array_equal(ridx, eridx))

        # non-unique
        """
        idx = Index([1,1,2,5])
        idx2 = Index([1,2,5,7,9])
        res, lidx, ridx = idx2.join(idx, how='left', return_indexers=True)
        eres = idx2
        eridx = np.array([0, 2, 3, -1, -1])
        elidx = np.array([0, 1, 2, 3, 4])
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))
        """

    def test_join_right(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='right',
                                          return_indexers=True)
        eres = other
        elidx = np.array([-1, 6, -1, -1, 1, -1],
                         dtype=np.int64)

        tm.assert_isinstance(other, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(ridx is None)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='right',
                                          return_indexers=True)
        eres = other_mono
        elidx = np.array([-1, 1, -1, -1, 6, -1],
                         dtype=np.int64)
        tm.assert_isinstance(other, Int64Index)
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(ridx is None)

        # non-unique
        """
        idx = Index([1,1,2,5])
        idx2 = Index([1,2,5,7,9])
        res, lidx, ridx = idx.join(idx2, how='right', return_indexers=True)
        eres = idx2
        elidx = np.array([0, 2, 3, -1, -1])
        eridx = np.array([0, 1, 2, 3, 4])
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

        idx = Index([1,1,2,5])
        idx2 = Index([1,2,5,9,7])
        res = idx.join(idx2, how='right', return_indexers=False)
        eres = idx2
        self.assert(res.equals(eres))
        """

    def test_join_non_int_index(self):
        other = Index([3, 6, 7, 8, 10], dtype=object)

        outer = self.index.join(other, how='outer')
        outer2 = other.join(self.index, how='outer')
        expected = Index([0, 2, 3, 4, 6, 7, 8, 10, 12, 14,
                          16, 18], dtype=object)
        self.assert_(outer.equals(outer2))
        self.assert_(outer.equals(expected))

        inner = self.index.join(other, how='inner')
        inner2 = other.join(self.index, how='inner')
        expected = Index([6, 8, 10], dtype=object)
        self.assert_(inner.equals(inner2))
        self.assert_(inner.equals(expected))

        left = self.index.join(other, how='left')
        self.assert_(left.equals(self.index))

        left2 = other.join(self.index, how='left')
        self.assert_(left2.equals(other))

        right = self.index.join(other, how='right')
        self.assert_(right.equals(other))

        right2 = other.join(self.index, how='right')
        self.assert_(right2.equals(self.index))

    def test_join_non_unique(self):
        left = Index([4, 4, 3, 3])

        joined, lidx, ridx = left.join(left, return_indexers=True)

        exp_joined = Index([3, 3, 3, 3, 4, 4, 4, 4])
        self.assert_(joined.equals(exp_joined))

        exp_lidx = np.array([2, 2, 3, 3, 0, 0, 1, 1], dtype=np.int64)
        self.assert_(np.array_equal(lidx, exp_lidx))

        exp_ridx = np.array([2, 3, 2, 3, 0, 1, 0, 1], dtype=np.int64)
        self.assert_(np.array_equal(ridx, exp_ridx))

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = self.index.join(self.index, how=kind)
            self.assert_(self.index is joined)

    def test_intersection(self):
        other = Index([1, 2, 3, 4, 5])
        result = self.index.intersection(other)
        expected = np.sort(np.intersect1d(self.index.values, other.values))
        self.assert_(np.array_equal(result, expected))

        result = other.intersection(self.index)
        expected = np.sort(np.asarray(np.intersect1d(self.index.values,
                                                     other.values)))
        self.assert_(np.array_equal(result, expected))

    def test_intersect_str_dates(self):
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        i1 = Index(dt_dates, dtype=object)
        i2 = Index(['aa'], dtype=object)
        res = i2.intersection(i1)

        self.assert_(len(res) == 0)

    def test_union_noncomparable(self):
        from datetime import datetime, timedelta
        # corner case, non-Int64Index
        now = datetime.now()
        other = Index([now + timedelta(i) for i in range(4)], dtype=object)
        result = self.index.union(other)
        expected = np.concatenate((self.index, other))
        self.assert_(np.array_equal(result, expected))

        result = other.union(self.index)
        expected = np.concatenate((other, self.index))
        self.assert_(np.array_equal(result, expected))

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
        self.assert_(result.dtype == np.object_)

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
            self.assertTrue(len(r) < 100)
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


class TestMultiIndex(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])
        self.index_names = ['first', 'second']
        self.index = MultiIndex(levels=[major_axis, minor_axis],
                                labels=[major_labels, minor_labels],
                                names=self.index_names, verify_integrity=False)

    def test_hash_error(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(self.index).__name__):
            hash(self.index)

    def test_set_names_and_rename(self):
        # so long as these are synonyms, we don't need to test set_names
        self.assert_(self.index.rename == self.index.set_names)
        new_names = [name + "SUFFIX" for name in self.index_names]
        ind = self.index.set_names(new_names)
        self.assertEqual(self.index.names, self.index_names)
        self.assertEqual(ind.names, new_names)
        with assertRaisesRegexp(ValueError, "^Length"):
            ind.set_names(new_names + new_names)
        new_names2 = [name + "SUFFIX2" for name in new_names]
        res = ind.set_names(new_names2, inplace=True)
        self.assert_(res is None)
        self.assertEqual(ind.names, new_names2)

    def test_set_levels_and_set_labels(self):
        # side note - you probably wouldn't want to use levels and labels
        # directly like this - but it is possible.
        levels, labels = self.index.levels, self.index.labels
        new_levels = [[lev + 'a' for lev in level] for level in levels]
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

        # level changing [w/o mutation]
        ind2 = self.index.set_levels(new_levels)
        assert_matching(ind2.levels, new_levels)
        assert_matching(self.index.levels, levels)

        # level changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_levels(new_levels, inplace=True)
        self.assert_(inplace_return is None)
        assert_matching(ind2.levels, new_levels)

        # label changing [w/o mutation]
        ind2 = self.index.set_labels(new_labels)
        assert_matching(ind2.labels, new_labels)
        assert_matching(self.index.labels, labels)

        # label changing [w/ mutation]
        ind2 = self.index.copy()
        inplace_return = ind2.set_labels(new_labels, inplace=True)
        self.assert_(inplace_return is None)
        assert_matching(ind2.labels, new_labels)

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
        self.assert_(mi1._tuples is not None)

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
        exp_values = np.array([(long(1), 'a')] * 6, dtype=object)
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
        self.assert_(df.is_copy is None)
        self.assertEqual(df.index.names, ('Name', 'Number'))
        df = df.set_value(('grethe', '4'), 'one', 99.34)
        self.assert_(df.is_copy is None)
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
        tm.assert_isinstance(single_level, Index)
        self.assert_(not isinstance(single_level, MultiIndex))
        self.assert_(single_level.name == 'first')

        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]])
        self.assert_(single_level.name is None)

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
        self.assert_(copy.labels is not original.labels)

        # names doesn't matter which way copied
        self.assertEqual(copy.names, original.names)
        self.assert_(copy.names is not original.names)

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
        self.assertEquals(list(result), list(self.index))

    def test_from_product(self):
        first = ['foo', 'bar', 'buz']
        second = ['a', 'b', 'c']
        names = ['first', 'second']
        result = MultiIndex.from_product([first, second], names=names)

        tuples = [('foo', 'a'), ('foo', 'b'), ('foo', 'c'),
                  ('bar', 'a'), ('bar', 'b'), ('bar', 'c'),
                  ('buz', 'a'), ('buz', 'b'), ('buz', 'c')]
        expected = MultiIndex.from_tuples(tuples, names=names)

        assert_array_equal(result, expected)
        self.assertEquals(result.names, names)

    def test_append(self):
        result = self.index[:3].append(self.index[3:])
        self.assert_(result.equals(self.index))

        foos = [self.index[:1], self.index[1:3], self.index[3:]]
        result = foos[0].append(foos[1:])
        self.assert_(result.equals(self.index))

        # empty
        result = self.index.append([])
        self.assert_(result.equals(self.index))

    def test_get_level_values(self):
        result = self.index.get_level_values(0)
        expected = ['foo', 'foo', 'bar', 'baz', 'qux', 'qux']
        self.assert_(np.array_equal(result, expected))

        self.assertEquals(result.name, 'first')

        result = self.index.get_level_values('first')
        expected = self.index.get_level_values(0)
        self.assert_(np.array_equal(result, expected))

    def test_get_level_values_na(self):
        arrays = [['a', 'b', 'b'], [1, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = [1, np.nan, 2]
        assert_array_equal(values.values.astype(float), expected)

        arrays = [['a', 'b', 'b'], [np.nan, np.nan, 2]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(1)
        expected = [np.nan, np.nan, 2]
        assert_array_equal(values.values.astype(float), expected)

        arrays = [[np.nan, np.nan, np.nan], ['a', np.nan, 1]]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        expected = [np.nan, np.nan, np.nan]
        assert_array_equal(values.values.astype(float), expected)
        values = index.get_level_values(1)
        expected = ['a', np.nan, 1]
        assert_array_equal(values.values, expected)

        if not _np_version_under1p7:
            arrays = [['a', 'b', 'b'], pd.DatetimeIndex([0, 1, pd.NaT])]
            index = pd.MultiIndex.from_arrays(arrays)
            values = index.get_level_values(1)
            expected = pd.DatetimeIndex([0, 1, pd.NaT])
            assert_array_equal(values.values, expected.values)

        arrays = [[], []]
        index = pd.MultiIndex.from_arrays(arrays)
        values = index.get_level_values(0)
        self.assertEqual(values.shape, (0,))

    def test_reorder_levels(self):
        # this blows up
        assertRaisesRegexp(IndexError, '^Too many levels',
                           self.index.reorder_levels, [2, 1, 0])

    def test_nlevels(self):
        self.assertEquals(self.index.nlevels, 2)

    def test_iter(self):
        result = list(self.index)
        expected = [('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                    ('baz', 'two'), ('qux', 'one'), ('qux', 'two')]
        self.assert_(result == expected)

    def test_pickle(self):
        pickled = pickle.dumps(self.index)
        unpickled = pickle.loads(pickled)
        self.assert_(self.index.equals(unpickled))

    def test_legacy_pickle(self):
        if compat.PY3:
            raise nose.SkipTest("doesn't work on Python 3")

        def curpath():
            pth, _ = os.path.split(os.path.abspath(__file__))
            return pth

        ppath = os.path.join(curpath(), 'data/multiindex_v1.pickle')
        obj = pickle.load(open(ppath, 'r'))

        self.assert_(obj._is_v1)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assert_(obj.equals(obj2))

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
        pth, _ = os.path.split(os.path.abspath(__file__))
        filepath = os.path.join(pth, 'data', 'mindex_073.pickle')

        obj = pd.read_pickle(filepath)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assert_(obj.equals(obj2))

        res = obj.get_indexer(obj)
        exp = np.arange(len(obj))
        assert_almost_equal(res, exp)

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_from_tuples_index_values(self):
        result = MultiIndex.from_tuples(self.index)
        self.assert_((result.values == self.index.values).all())

    def test_contains(self):
        self.assert_(('foo', 'two') in self.index)
        self.assert_(('bar', 'two') not in self.index)
        self.assert_(None not in self.index)

    def test_is_all_dates(self):
        self.assert_(not self.index.is_all_dates)

    def test_is_numeric(self):
        # MultiIndex is never numeric
        self.assert_(not self.index.is_numeric())

    def test_getitem(self):
        # scalar
        self.assertEquals(self.index[2], ('bar', 'one'))

        # slice
        result = self.index[2:5]
        expected = self.index[[2, 3, 4]]
        self.assert_(result.equals(expected))

        # boolean
        result = self.index[[True, False, True, False, True, True]]
        result2 = self.index[np.array([True, False, True, False, True, True])]
        expected = self.index[[0, 2, 4, 5]]
        self.assert_(result.equals(expected))
        self.assert_(result2.equals(expected))

    def test_getitem_group_select(self):
        sorted_idx, _ = self.index.sortlevel(0)
        self.assertEquals(sorted_idx.get_loc('baz'), slice(3, 4))
        self.assertEquals(sorted_idx.get_loc('foo'), slice(0, 2))

    def test_get_loc(self):
        self.assert_(self.index.get_loc(('foo', 'two')) == 1)
        self.assert_(self.index.get_loc(('baz', 'two')) == 3)
        self.assertRaises(KeyError, self.index.get_loc, ('bar', 'two'))
        self.assertRaises(KeyError, self.index.get_loc, 'quux')

        # 3 levels
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])
        self.assertRaises(KeyError, index.get_loc, (1, 1))
        self.assert_(index.get_loc((2, 0)) == slice(3, 5))

    def test_get_loc_duplicates(self):
        index = Index([2, 2, 2, 2])
        result = index.get_loc(2)
        expected = slice(0, 4)
        assert(result == expected)
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
        self.assert_(new_index.equals(exp_index))

        loc, new_index = index.get_loc_level((0, 1, 0))
        expected = 1
        self.assertEqual(loc, expected)
        self.assert_(new_index is None)

        self.assertRaises(KeyError, index.get_loc_level, (2, 2))

        index = MultiIndex(levels=[[2000], lrange(4)],
                           labels=[np.array([0, 0, 0, 0]),
                                   np.array([0, 1, 2, 3])])
        result, new_index = index.get_loc_level((2000, slice(None, None)))
        expected = slice(None, None)
        self.assertEqual(result, expected)
        self.assert_(new_index.equals(index.droplevel(0)))

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
        self.assertEquals(result, (1, 5))

        result = sorted_idx.slice_locs(None, ('qux', 'one'))
        self.assertEquals(result, (0, 5))

        result = sorted_idx.slice_locs(('foo', 'two'), None)
        self.assertEquals(result, (1, len(sorted_idx)))

        result = sorted_idx.slice_locs('bar', 'baz')
        self.assertEquals(result, (2, 4))

    def test_slice_locs_not_contained(self):
        # some searchsorted action

        index = MultiIndex(levels=[[0, 2, 4, 6], [0, 2, 4]],
                           labels=[[0, 0, 0, 1, 1, 2, 3, 3, 3],
                                   [0, 1, 2, 1, 2, 2, 0, 1, 2]],
                           sortorder=0)

        result = index.slice_locs((1, 0), (5, 2))
        self.assertEquals(result, (3, 6))

        result = index.slice_locs(1, 5)
        self.assertEquals(result, (3, 6))

        result = index.slice_locs((2, 2), (5, 2))
        self.assertEquals(result, (3, 6))

        result = index.slice_locs(2, 5)
        self.assertEquals(result, (3, 6))

        result = index.slice_locs((1, 0), (6, 3))
        self.assertEquals(result, (3, 8))

        result = index.slice_locs(-1, 10)
        self.assertEquals(result, (0, len(index)))

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

        self.assert_(not index.is_unique)

    def test_truncate(self):
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        result = index.truncate(before=1)
        self.assert_('foo' not in result.levels[0])
        self.assert_(1 in result.levels[0])

        result = index.truncate(after=1)
        self.assert_(2 not in result.levels[0])
        self.assert_(1 in result.levels[0])

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
        assert_almost_equal(r1, [-1, 0, 0, 1, 1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        assert_almost_equal(r1, [0, 0, 1, 1, 2])

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        # pass non-MultiIndex
        r1 = idx1.get_indexer(idx2._tuple_index)
        rexp1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1, 2, 3])
        self.assert_((r1 == [-1, -1, -1]).all())

        # create index with duplicates
        idx1 = Index(lrange(10) + lrange(10))
        idx2 = Index(lrange(20))
        assertRaisesRegexp(InvalidIndexError, "Reindexing only valid with"
                           " uniquely valued Index objects",
                           idx1.get_indexer, idx2)

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
        import warnings
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")
        # GH1538
        pd.set_option('display.multi_sparse', False)

        result = self.index.format()
        self.assertEqual(result[1], 'foo  two')

        pd.reset_option("^display\.")

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
        self.assert_(self.index.equals(self.index))
        self.assert_(self.index.equal_levels(self.index))

        self.assert_(not self.index.equals(self.index[:-1]))

        self.assert_(self.index.equals(self.index._tuple_index))

        # different number of levels
        index = MultiIndex(levels=[Index(lrange(4)),
                                   Index(lrange(4)),
                                   Index(lrange(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        index2 = MultiIndex(levels=index.levels[:-1],
                            labels=index.labels[:-1])
        self.assert_(not index.equals(index2))
        self.assert_(not index.equal_levels(index2))

        # levels are different
        major_axis = Index(lrange(4))
        minor_axis = Index(lrange(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assert_(not self.index.equals(index))
        self.assert_(not self.index.equal_levels(index))

        # some of the labels are different
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        self.assert_(not self.index.equals(index))

    def test_identical(self):
        mi = self.index.copy()
        mi2 = self.index.copy()
        self.assert_(mi.identical(mi2))

        mi = mi.set_names(['new1', 'new2'])
        self.assert_(mi.equals(mi2))
        self.assert_(not mi.identical(mi2))

        mi2 = mi2.set_names(['new1', 'new2'])
        self.assert_(mi.identical(mi2))

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

        self.assert_(the_union.equals(expected))

        # corner case, pass self or empty thing:
        the_union = self.index.union(self.index)
        self.assert_(the_union is self.index)

        the_union = self.index.union(self.index[:0])
        self.assert_(the_union is self.index)

        # won't work in python 3
        # tuples = self.index._tuple_index
        # result = self.index[:4] | tuples[4:]
        # self.assert_(result.equals(tuples))

    # not valid for python 3
    # def test_union_with_regular_index(self):
    #     other = Index(['A', 'B', 'C'])

    #     result = other.union(self.index)
    #     self.assert_(('foo', 'one') in result)
    #     self.assert_('B' in result)

    #     result2 = self.index.union(other)
    #     self.assert_(result.equals(result2))

    def test_intersection(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_int = piece1 & piece2
        tups = sorted(self.index[3:5]._tuple_index)
        expected = MultiIndex.from_tuples(tups)
        self.assert_(the_int.equals(expected))

        # corner case, pass self
        the_int = self.index.intersection(self.index)
        self.assert_(the_int is self.index)

        # empty intersection: disjoint
        empty = self.index[:2] & self.index[2:]
        expected = self.index[:0]
        self.assert_(empty.equals(expected))

        # can't do in python 3
        # tuples = self.index._tuple_index
        # result = self.index & tuples
        # self.assert_(result.equals(tuples))

    def test_diff(self):
        first = self.index
        result = first - self.index[-3:]
        expected = MultiIndex.from_tuples(sorted(self.index[:-3].values),
                                          sortorder=0,
                                          names=self.index.names)

        tm.assert_isinstance(result, MultiIndex)
        self.assert_(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: reflexive
        result = self.index - self.index
        expected = self.index[:0]
        self.assert_(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: superset
        result = self.index[-3:] - self.index
        expected = self.index[:0]
        self.assert_(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # empty difference: degenerate
        result = self.index[:0] - self.index
        expected = self.index[:0]
        self.assert_(result.equals(expected))
        self.assertEqual(result.names, self.index.names)

        # names not the same
        chunklet = self.index[-3:]
        chunklet.names = ['foo', 'baz']
        result = first - chunklet
        self.assertEqual(result.names, (None, None))

        # empty, but non-equal
        result = self.index - self.index.sortlevel(1)[0]
        self.assert_(len(result) == 0)

        # raise Exception called with non-MultiIndex
        result = first.diff(first._tuple_index)
        self.assertTrue(result.equals(first[:0]))

        # name from empty array
        result = first.diff([])
        self.assert_(first.equals(result))
        self.assertEqual(first.names, result.names)

        # name from non-empty array
        result = first.diff([('foo', 'one')])
        expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'),
                                            ('foo', 'two'), ('qux', 'one'),
                                            ('qux', 'two')])
        expected.names = first.names
        self.assertEqual(first.names, result.names)
        assertRaisesRegexp(TypeError, "other must be a MultiIndex or a list"
                           " of tuples", first.diff, [1, 2, 3, 4, 5])

    def test_from_tuples(self):
        assertRaisesRegexp(TypeError, 'Cannot infer number of levels from'
                           ' empty list', MultiIndex.from_tuples, [])

        idx = MultiIndex.from_tuples(((1, 2), (3, 4)), names=['a', 'b'])
        self.assertEquals(len(idx), 2)

    def test_argsort(self):
        result = self.index.argsort()
        expected = self.index._tuple_index.argsort()
        self.assert_(np.array_equal(result, expected))

    def test_sortlevel(self):
        import random

        tuples = list(self.index)
        random.shuffle(tuples)

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assert_(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assert_(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assert_(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assert_(sorted_idx.equals(expected[::-1]))

    def test_sortlevel_deterministic(self):
        tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
                  ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

        index = MultiIndex.from_tuples(tuples)

        sorted_idx, _ = index.sortlevel(0)
        expected = MultiIndex.from_tuples(sorted(tuples))
        self.assert_(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(0, ascending=False)
        self.assert_(sorted_idx.equals(expected[::-1]))

        sorted_idx, _ = index.sortlevel(1)
        by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
        expected = MultiIndex.from_tuples(by1)
        self.assert_(sorted_idx.equals(expected))

        sorted_idx, _ = index.sortlevel(1, ascending=False)
        self.assert_(sorted_idx.equals(expected[::-1]))

    def test_dims(self):
        pass

    def test_drop(self):
        dropped = self.index.drop([('foo', 'two'), ('qux', 'one')])

        index = MultiIndex.from_tuples([('foo', 'two'), ('qux', 'one')])
        dropped2 = self.index.drop(index)

        expected = self.index[[0, 2, 3, 5]]
        self.assert_(dropped.equals(expected))
        self.assert_(dropped2.equals(expected))

        dropped = self.index.drop(['bar'])
        expected = self.index[[0, 1, 3, 4, 5]]
        self.assert_(dropped.equals(expected))

        index = MultiIndex.from_tuples([('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, [('bar', 'two')])
        self.assertRaises(KeyError, self.index.drop, index)

        # mixed partial / full drop
        dropped = self.index.drop(['foo', ('qux', 'one')])
        expected = self.index[[2, 3, 5]]
        self.assert_(dropped.equals(expected))

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
        self.assert_(dropped.equals(expected))

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
        self.assert_(dropped.equals(expected))

    def test_insert(self):
        # key contained in all levels
        new_index = self.index.insert(0, ('bar', 'two'))
        self.assert_(new_index.equal_levels(self.index))
        self.assert_(new_index[0] == ('bar', 'two'))

        # key not contained in all levels
        new_index = self.index.insert(0, ('abc', 'three'))
        self.assert_(np.array_equal(new_index.levels[0],
                                    list(self.index.levels[0]) + ['abc']))
        self.assert_(np.array_equal(new_index.levels[1],
                                    list(self.index.levels[1]) + ['three']))
        self.assert_(new_index[0] == ('abc', 'three'))

        # key wrong length
        assertRaisesRegexp(ValueError, "Item must have length equal to number"
                           " of levels", self.index.insert, 0, ('foo2',))

    def test_take_preserve_name(self):
        taken = self.index.take([3, 0, 1])
        self.assertEqual(taken.names, self.index.names)

    def test_join_level(self):
        def _check_how(other, how):
            join_index, lidx, ridx = other.join(self.index, how=how,
                                                level='second',
                                                return_indexers=True)

            exp_level = other.join(self.index.levels[1], how=how)
            self.assert_(join_index.levels[0].equals(self.index.levels[0]))
            self.assert_(join_index.levels[1].equals(exp_level))

            # pare down levels
            mask = np.array(
                [x[1] in exp_level for x in self.index], dtype=bool)
            exp_values = self.index.values[mask]
            self.assert_(np.array_equal(join_index.values, exp_values))

            if how in ('outer', 'inner'):
                join_index2, ridx2, lidx2 = \
                    self.index.join(other, how=how, level='second',
                                    return_indexers=True)

                self.assert_(join_index.equals(join_index2))
                self.assert_(np.array_equal(lidx, lidx2))
                self.assert_(np.array_equal(ridx, ridx2))
                self.assert_(np.array_equal(join_index2.values, exp_values))

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
        tm.assert_isinstance(result, MultiIndex)

        assertRaisesRegexp(TypeError, "Join.*MultiIndex.*ambiguous",
                           self.index.join, self.index, level=1)

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            res = self.index
            joined = res.join(res, how=kind)
            self.assert_(res is joined)

    def test_reindex(self):
        result, indexer = self.index.reindex(list(self.index[:4]))
        tm.assert_isinstance(result, MultiIndex)
        self.check_level_names(result, self.index[:4].names)

        result, indexer = self.index.reindex(list(self.index))
        tm.assert_isinstance(result, MultiIndex)
        self.assert_(indexer is None)
        self.check_level_names(result, self.index.names)

    def test_reindex_level(self):
        idx = Index(['one'])

        target, indexer = self.index.reindex(idx, level='second')
        target2, indexer2 = idx.reindex(self.index, level='second')

        exp_index = self.index.join(idx, level='second', how='right')
        exp_index2 = self.index.join(idx, level='second', how='left')

        self.assert_(target.equals(exp_index))
        exp_indexer = np.array([0, 2, 4])
        self.assert_(np.array_equal(indexer, exp_indexer))

        self.assert_(target2.equals(exp_index2))
        exp_indexer2 = np.array([0, -1, 0, -1, 0, -1])
        self.assert_(np.array_equal(indexer2, exp_indexer2))

        assertRaisesRegexp(TypeError, "Fill method not supported",
                           self.index.reindex, self.index, method='pad',
                           level='second')

        assertRaisesRegexp(TypeError, "Fill method not supported",
                           idx.reindex, idx, method='bfill', level='first')

    def test_has_duplicates(self):
        self.assert_(not self.index.has_duplicates)
        self.assert_(self.index.append(self.index).has_duplicates)

        index = MultiIndex(levels=[[0, 1], [0, 1, 2]],
                           labels=[[0, 0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 0, 1, 2]])
        self.assert_(index.has_duplicates)

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
        tm.assert_index_equal(eval(repr(self.index)), self.index)

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


def test_get_combined_index():
    from pandas.core.index import _get_combined_index
    result = _get_combined_index([])
    assert(result.equals(Index([])))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
