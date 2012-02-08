# pylint: disable=E1101,E1103,W0232

from datetime import datetime, timedelta
import operator
import pickle
import unittest
import nose

import numpy as np

from pandas.core.index import Index, Int64Index, Factor, MultiIndex
from pandas.util.testing import assert_almost_equal
from pandas.util import py3compat

import pandas.util.testing as tm
import pandas._tseries as tseries

class TestIndex(unittest.TestCase):

    def setUp(self):
        self.strIndex = tm.makeStringIndex(100)
        self.dateIndex = tm.makeDateIndex(100)
        self.intIndex = tm.makeIntIndex(100)
        self.empty = Index([])
        self.tuples = Index(zip(['foo', 'bar', 'baz'], [1, 2, 3]))

    def test_hash_error(self):
        self.assertRaises(TypeError, hash, self.strIndex)

    def test_new_axis(self):
        new_index = self.dateIndex[None, :]
        self.assert_(new_index.ndim == 2)
        self.assert_(type(new_index) == np.ndarray)

    def test_deepcopy(self):
        from copy import deepcopy

        copy = deepcopy(self.strIndex)
        self.assert_(copy is self.strIndex)

    def test_duplicates(self):
        idx = Index([0, 0, 0])
        self.assert_(not idx._verify_integrity())
        self.assertRaises(Exception, getattr, idx, 'indexMap')

    def test_sort(self):
        self.assertRaises(Exception, self.strIndex.sort)

    def test_mutability(self):
        self.assertRaises(Exception, self.strIndex.__setitem__, 5, 0)
        self.assertRaises(Exception, self.strIndex.__setitem__, slice(1,5), 0)

    def test_constructor(self):
        # regular instance creation
        tm.assert_contains_all(self.strIndex, self.strIndex)
        tm.assert_contains_all(self.dateIndex, self.dateIndex)

        # casting
        arr = np.array(self.strIndex)
        index = arr.view(Index)
        tm.assert_contains_all(arr, index)
        self.assert_(np.array_equal(self.strIndex, index))

        # what to do here?
        # arr = np.array(5.)
        # self.assertRaises(Exception, arr.view, Index)

    def test_constructor_corner(self):
        # corner case
        self.assertRaises(Exception, Index, 0)

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

    def test_asof(self):
        d = self.dateIndex[0]
        self.assert_(self.dateIndex.asof(d) is d)
        self.assert_(np.isnan(self.dateIndex.asof(d - timedelta(1))))

        d = self.dateIndex[-1]
        self.assert_(self.dateIndex.asof(d + timedelta(1)) is d)

    def test_argsort(self):
        result = self.strIndex.argsort()
        expected = np.array(self.strIndex).argsort()
        self.assert_(np.array_equal(result, expected))

    def test_comparators(self):
        index = self.dateIndex
        element = index[len(index) // 2]
        arr = np.array(index)

        def _check(op):
            arr_result = op(arr, element)
            index_result = op(index, element)

            self.assert_(isinstance(index_result, np.ndarray))
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
        tm.assert_dict_equal(tseries.map_indices_object(subIndex),
                             subIndex.indexMap)

        subIndex = self.strIndex[list(boolIdx)]
        tm.assert_dict_equal(tseries.map_indices_object(subIndex),
                             subIndex.indexMap)

    def test_fancy(self):
        sl = self.strIndex[[1,2,3]]
        for i in sl:
            self.assertEqual(i, sl[sl.indexMap[i]])

    def test_getitem(self):
        arr = np.array(self.dateIndex)
        self.assertEquals(self.dateIndex[5], arr[5])

    def test_shift(self):
        shifted = self.dateIndex.shift(0, timedelta(1))
        self.assert_(shifted is self.dateIndex)

        shifted = self.dateIndex.shift(5, timedelta(1))
        self.assert_(np.array_equal(shifted, self.dateIndex + timedelta(5)))

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)

        self.assert_(tm.equalContents(intersect, second))

        # Corner cases
        inter = first.intersection(first)
        self.assert_(inter is first)

        # non-iterable input
        self.assertRaises(Exception, first.intersection, 0.5)

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

        # non-iterable input
        self.assertRaises(Exception, first.union, 0.5)

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

        self.assert_(tm.equalContents(np.append(self.strIndex,
                                                self.dateIndex), firstCat))
        self.assert_(tm.equalContents(secondCat, self.strIndex))
        tm.assert_contains_all(self.strIndex, firstCat.indexMap)
        tm.assert_contains_all(self.strIndex, secondCat.indexMap)
        tm.assert_contains_all(self.dateIndex, firstCat.indexMap)

        # this is valid too
        shifted = self.dateIndex + timedelta(1)

    def test_append_multiple(self):
        index = Index(['a', 'b', 'c', 'd', 'e', 'f'])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        self.assert_(result.equals(index))

        # empty
        result = index.append([])
        self.assert_(result.equals(index))

    def test_add_string(self):
        # from bug report
        index = Index(['a', 'b', 'c'])
        index2 = index + 'foo'

        self.assert_('a' not in index2.indexMap)
        self.assert_('afoo' in index2.indexMap)

    def test_diff(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]
        result = first - second

        self.assert_(tm.equalContents(result, answer))

        diff = first.diff(first)
        self.assert_(len(diff) == 0)

        # non-iterable input
        self.assertRaises(Exception, first.diff, 0.5)

    def test_pickle(self):
        def testit(index):
            pickled = pickle.dumps(index)
            unpickled = pickle.loads(pickled)

            self.assert_(isinstance(unpickled, Index))
            self.assert_(np.array_equal(unpickled, index))
            self.assertEquals(unpickled.name, index.name)

            tm.assert_dict_equal(unpickled.indexMap, index.indexMap)

        testit(self.strIndex)
        self.strIndex.name = 'foo'
        testit(self.strIndex)

        testit(self.dateIndex)

    # def test_always_get_null_index(self):
    #     empty = Index([])
    #     self.assert_(empty is NULL_INDEX)
    #     self.assert_(self.dateIndex[15:15] is NULL_INDEX)

    def test_is_all_dates(self):
        self.assert_(self.dateIndex.is_all_dates)
        self.assert_(not self.strIndex.is_all_dates)
        self.assert_(not self.intIndex.is_all_dates)

    def test_summary(self):
        self._check_method_works(Index.summary)

    def test_format(self):
        self._check_method_works(Index.format)

        index = Index([datetime.now()])
        formatted = index.format()
        expected = [str(index[0])]
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

    def test_drop(self):
        n = len(self.strIndex)

        dropped = self.strIndex.drop(self.strIndex[range(5, 10)])
        expected = self.strIndex[range(5) + range(10, n)]
        self.assert_(dropped.equals(expected))

        self.assertRaises(ValueError, self.strIndex.drop, ['foo', 'bar'])

    def test_tuple_union_bug(self):
        import pandas
        import numpy as np

        aidx1 = np.array([(1, 'A'),(2, 'A'),(1, 'B'),(2, 'B')], dtype=[('num',
        int),('let', 'a1')])
        aidx2 = np.array([(1, 'A'),(2, 'A'),(1, 'B'),(2, 'B'),(1,'C'),(2,
        'C')], dtype=[('num', int),('let', 'a1')])

        idx1 = pandas.Index(aidx1)
        idx2 = pandas.Index(aidx2)

        # intersection broken?
        int_idx = idx1.intersection(idx2)
        # needs to be 1d like idx1 and idx2
        expected = idx1[:4] # pandas.Index(sorted(set(idx1) & set(idx2)))
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

class TestInt64Index(unittest.TestCase):

    def setUp(self):
        self.index = Int64Index(np.arange(0, 20, 2))

    def test_constructor(self):
        # pass list, coerce fine
        index = Int64Index([-5, 0, 1, 2])
        expected = np.array([-5, 0, 1, 2], dtype=np.int64)
        self.assert_(np.array_equal(index, expected))

        # from iterable
        index = Int64Index(iter([-5, 0, 1, 2]))
        self.assert_(np.array_equal(index, expected))

        # scalar raise Exception
        self.assertRaises(ValueError, Int64Index, 5)

    def test_constructor_corner(self):
        arr = np.array([1, 2, 3, 4], dtype=object)
        index = Int64Index(arr)
        self.assert_(index.values.dtype == np.int64)
        self.assert_(index.equals(arr))

        # preventing casting
        arr = np.array([1, '2', 3, '4'], dtype=object)
        self.assertRaises(TypeError, Int64Index, arr)

    def test_coerce_list(self):
        # coerce things
        arr = Index([1, 2, 3, 4])
        self.assert_(type(arr) == Int64Index)

        # but not if explicit dtype passed
        arr = Index([1, 2, 3, 4], dtype=object)
        self.assert_(type(arr) == Index)

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
                         dtype='i4')
        eridx = np.array([-1, 3, 4, -1, 5, -1, 0, -1, -1, 1, -1, -1, -1, 2],
                         dtype='i4')

        self.assert_(isinstance(res, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other_mono, how='outer')
        self.assert_(res.equals(noidx_res))

        eridx = np.array([-1, 0, 1, -1, 2, -1, 3, -1, -1, 4, -1, -1, -1, 5],
                         dtype='i4')
        self.assert_(isinstance(res, Int64Index))
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

        self.assert_(isinstance(res, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='inner',
                                          return_indexers=True)

        res2 = self.index.intersection(other_mono)
        self.assert_(res.equals(res2))

        eridx = np.array([1, 4])
        self.assert_(isinstance(res, Int64Index))
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
                         dtype='i4')

        self.assert_(isinstance(res, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(lidx is None)
        self.assert_(np.array_equal(ridx, eridx))

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='left',
                                          return_indexers=True)
        eridx = np.array([-1, 1, -1, -1, -1, -1, 4, -1, -1, -1],
                         dtype='i4')
        self.assert_(isinstance(res, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(lidx is None)
        self.assert_(np.array_equal(ridx, eridx))

    def test_join_right(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='right',
                                          return_indexers=True)
        eres = other
        elidx = np.array([-1, 6, -1, -1, 1, -1],
                         dtype='i4')

        self.assert_(isinstance(other, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(ridx is None)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='right',
                                          return_indexers=True)
        eres = other_mono
        elidx = np.array([-1, 1, -1, -1, 6, -1],
                         dtype='i4')
        self.assert_(isinstance(other, Int64Index))
        self.assert_(res.equals(eres))
        self.assert_(np.array_equal(lidx, elidx))
        self.assert_(ridx is None)

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

    def test_intersection(self):
        other = Index([1, 2, 3, 4, 5])
        result = self.index.intersection(other)
        expected = np.sort(np.intersect1d(self.index.values, other.values))
        self.assert_(np.array_equal(result, expected))

        result = other.intersection(self.index)
        expected = np.sort(np.asarray(np.intersect1d(self.index.values,
                                                     other.values)))
        self.assert_(np.array_equal(result, expected))

    def test_union_noncomparable(self):
        from datetime import datetime, timedelta
        # corner case, non-Int64Index
        now = datetime.now()
        other = Index([now + timedelta(i) for i in xrange(4)])
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
        index = Int64Index([1,2,3,4], name='foo')
        taken = index.take([3,0,1])
        self.assertEqual(index.name, taken.name)

    def test_int_name_format(self):
        from pandas import Series, DataFrame
        index = Index(['a', 'b', 'c'], name=0)
        s = Series(range(3), index)
        df = DataFrame(range(3), index=index)
        repr(s)
        repr(df)

class TestMultiIndex(unittest.TestCase):

    def setUp(self):
        major_axis = Index(['foo', 'bar', 'baz', 'qux'])
        minor_axis = Index(['one', 'two'])

        major_labels = np.array([0, 0, 1, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 0, 1])

        self.index = MultiIndex(levels=[major_axis, minor_axis],
                                labels=[major_labels, minor_labels],
                                names=['first', 'second'])

    def test_constructor_single_level(self):
        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]],
                                  names=['first'])
        self.assert_(isinstance(single_level, Index))
        self.assert_(not isinstance(single_level, MultiIndex))
        self.assert_(single_level.name == 'first')

        single_level = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux']],
                                  labels=[[0, 1, 2, 3]])
        self.assert_(single_level.name is None)

    def test_constructor_no_levels(self):
        self.assertRaises(Exception, MultiIndex, levels=[], labels=[])

    def test_duplicate_names(self):
        self.index.names = ['foo', 'foo']
        self.assertRaises(Exception, self.index._get_level_number, 'foo')

    def test_get_level_number_integer(self):
        self.index.names = [1, 0]
        self.assertEqual(self.index._get_level_number(1), 0)
        self.assertEqual(self.index._get_level_number(0), 1)
        self.assertRaises(Exception, self.index._get_level_number, 2)

        self.assertRaises(Exception, self.index._get_level_number, 'fourth')

    def test_from_arrays(self):
        arrays = []
        for lev, lab in zip(self.index.levels, self.index.labels):
            arrays.append(np.asarray(lev).take(lab))

        result = MultiIndex.from_arrays(arrays)
        self.assertEquals(list(result), list(self.index))

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

        result = self.index.get_level_values('first')
        expected = self.index.get_level_values(0)
        self.assert_(np.array_equal(result, expected))

    def test_reorder_levels(self):
        # this blows up
        self.assertRaises(Exception, self.index.reorder_levels,
                          [2, 1, 0])

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
        if py3compat.PY3:
            raise nose.SkipTest

        import os
        def curpath():
            pth, _ = os.path.split(os.path.abspath(__file__))
            return pth

        ppath = os.path.join(curpath(), 'data/multiindex_v1.pickle')
        obj = pickle.load(open(ppath, 'r'))

        self.assert_(obj._is_legacy_format)

        obj2 = MultiIndex.from_tuples(obj.values)
        self.assert_(obj.equals(obj2))

        res = obj.get_indexer(obj2[::-1])
        exp = obj.get_indexer(obj[::-1])
        exp2 = obj2.get_indexer(obj2[::-1])
        assert_almost_equal(res, exp)
        assert_almost_equal(exp, exp2)

    def test_contains(self):
        self.assert_(('foo', 'two') in self.index)
        self.assert_(('bar', 'two') not in self.index)
        self.assert_(None not in self.index)

    def test_is_all_dates(self):
        self.assert_(not self.index.is_all_dates)

    def test_getitem(self):
        # scalar
        self.assertEquals(self.index[2], ('bar', 'one'))

        # slice
        result = self.index[2:5]
        expected = self.index[[2,3,4]]
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
        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])
        self.assertRaises(KeyError, index.get_loc, (1, 1))
        self.assert_(index.get_loc((2, 0)) == slice(3, 5))

    def test_get_loc_duplicates(self):
        index = Index([2, 2, 2, 2])
        self.assertRaises(Exception, index.get_loc, 2)

    def test_get_loc_level(self):
        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
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

        index = MultiIndex(levels=[[2000], range(4)],
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

    def test_slice_locs_not_sorted(self):
        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        self.assertRaises(Exception, index.slice_locs, (1, 0, 1),
                          (2, 1, 0))

        # works
        sorted_index, _ = index.sortlevel(0)
        result = sorted_index.slice_locs((1, 0, 1), (2, 1, 0))

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
        major_axis = range(70000)
        minor_axis = range(10)

        major_labels = np.arange(70000)
        minor_labels = np.repeat(range(10), 7000)

        # the fact that is works means it's consistent
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        # inconsistent
        major_labels = np.array([0, 0, 1, 1, 1, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1])
        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])

        self.assertRaises(Exception, getattr, index, 'indexMap')

    def test_truncate(self):
        major_axis = Index(range(4))
        minor_axis = Index(range(2))

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
        major_axis = Index(range(4))
        minor_axis = Index(range(2))

        major_labels = np.array([0, 0, 1, 2, 2, 3, 3])
        minor_labels = np.array([0, 1, 0, 0, 1, 0, 1])

        index = MultiIndex(levels=[major_axis, minor_axis],
                           labels=[major_labels, minor_labels])
        idx1 = index[:5]
        idx2 = index[[1,3,5]]

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
        r1 = idx1.get_indexer(idx2.get_tuple_index())
        rexp1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1,2,3])
        self.assert_( (r1 == [-1, -1, -1]).all() )

        # self.assertRaises(Exception, idx1.get_indexer,
        #                   list(list(zip(*idx2.get_tuple_index()))[0]))

    def test_format(self):
        self.index.format()
        self.index[:0].format()

    def test_format_integer_names(self):
        index = MultiIndex(levels=[[0, 1], [0, 1]],
                           labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                           names=[0, 1])
        index.format(names=True)

    def test_bounds(self):
        self.index._bounds

    def test_equals(self):
        self.assert_(self.index.equals(self.index))
        self.assert_(self.index.equal_levels(self.index))

        self.assert_(not self.index.equals(self.index[:-1]))

        self.assert_(self.index.equals(self.index.get_tuple_index()))

        # different number of levels
        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])])

        index2 = MultiIndex(levels=index.levels[:-1],
                            labels=index.labels[:-1])
        self.assert_(not index.equals(index2))
        self.assert_(not index.equal_levels(index2))

        # levels are different
        major_axis = Index(range(4))
        minor_axis = Index(range(2))

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

    def test_union(self):
        piece1 = self.index[:5][::-1]
        piece2 = self.index[3:]

        the_union = piece1 | piece2

        tups = sorted(self.index.get_tuple_index())
        expected = MultiIndex.from_tuples(tups)

        self.assert_(the_union.equals(expected))

        # corner case, pass self or empty thing:
        the_union = self.index.union(self.index)
        self.assert_(the_union is self.index)

        the_union = self.index.union(self.index[:0])
        self.assert_(the_union is self.index)

        # won't work in python 3
        # tuples = self.index.get_tuple_index()
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
        tups = sorted(self.index[3:5].get_tuple_index())
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
        # tuples = self.index.get_tuple_index()
        # result = self.index & tuples
        # self.assert_(result.equals(tuples))

    def test_diff(self):
        first = self.index
        result = first - self.index[-3:]
        expected = MultiIndex.from_tuples(sorted(self.index[:-3].values),
                                          sortorder=0,
                                          names=self.index.names)

        self.assert_(isinstance(result, MultiIndex))
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
        self.assertEqual(result.names, [None, None])

        # empty, but non-equal
        result = self.index - self.index.sortlevel(1)[0]
        self.assert_(len(result) == 0)

        # raise Exception called with non-MultiIndex
        self.assertRaises(Exception, first.diff, first.get_tuple_index())

    def test_from_tuples(self):
        self.assertRaises(Exception, MultiIndex.from_tuples, [])

    def test_argsort(self):
        result = self.index.argsort()
        expected = self.index.get_tuple_index().argsort()
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
        self.assertRaises(Exception, self.index.drop, [('bar', 'two')])
        self.assertRaises(Exception, self.index.drop, index)

        # mixed partial / full drop
        dropped = self.index.drop(['foo', ('qux', 'one')])
        expected = self.index[[2, 3, 5]]
        self.assert_(dropped.equals(expected))

    def test_droplevel_with_names(self):
        index = self.index[self.index.get_loc('foo')]
        dropped = index.droplevel(0)
        self.assertEqual(dropped.name, 'second')

        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
                           labels=[np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                                   np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                                   np.array([1, 0, 1, 1, 0, 0, 1, 0])],
                           names=['one', 'two', 'three'])
        dropped = index.droplevel(0)
        self.assertEqual(dropped.names, ['two', 'three'])

        dropped = index.droplevel('two')
        expected = index.droplevel(1)
        self.assert_(dropped.equals(expected))

    def test_droplevel_multiple(self):
        index = MultiIndex(levels=[Index(range(4)),
                                   Index(range(4)),
                                   Index(range(4))],
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
        self.assertRaises(Exception, self.index.insert, 0, ('foo2',))

    def test_take_preserve_name(self):
        taken = self.index.take([3,0,1])
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
            mask = np.array([x[1] in exp_level for x in self.index], dtype=bool)
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
        self.assert_(isinstance(result, MultiIndex))

        self.assertRaises(Exception, self.index.join, self.index, level=1)

    def test_reindex(self):
        result, indexer = self.index.reindex(list(self.index[:4]))
        self.assert_(isinstance(result, MultiIndex))

        result, indexer = self.index.reindex(list(self.index))
        self.assert_(isinstance(result, MultiIndex))
        self.assert_(indexer is None)

    def test_reindex_level(self):
        idx = Index(['one'])

        target, indexer = self.index.reindex(idx, level='second')
        target2, indexer2 = idx.reindex(self.index, idx, level='second')

        exp_index = self.index.join(idx, level='second', how='left')
        self.assert_(target.equals(exp_index))
        self.assert_(target2.equals(exp_index))

    def test_has_duplicates(self):
        self.assert_(not self.index.has_duplicates)
        self.assert_(self.index.append(self.index).has_duplicates)

        index = MultiIndex(levels=[[0, 1], [0, 1, 2]],
                           labels=[[0, 0, 0, 0, 1, 1, 1],
                                   [0, 1, 2, 0, 0, 1, 2]])
        self.assert_(index.has_duplicates)

class TestFactor(unittest.TestCase):

    def setUp(self):
        self.factor = Factor(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])

    def test_getitem(self):
        self.assertEqual(self.factor[0], 'a')
        self.assertEqual(self.factor[-1], 'c')

        subf = self.factor[[0, 1, 2]]
        tm.assert_almost_equal(subf.labels, [0, 1, 1])

        subf = self.factor[self.factor.asarray() == 'c']
        tm.assert_almost_equal(subf.labels, [2, 2, 2])

    def test_constructor_unsortable(self):
        arr = np.array([1, 2, 3, datetime.now()], dtype='O')

        # it works!
        factor = Factor(arr)

    def test_factor_agg(self):
        import pandas.core.frame as frame

        arr = np.arange(len(self.factor))

        f = np.sum
        agged = frame.factor_agg(self.factor, arr, f)
        labels = self.factor.labels
        for i, idx in enumerate(self.factor.levels):
            self.assertEqual(f(arr[labels == i]), agged[i])


def test_get_combined_index():
    from pandas.core.index import _get_combined_index, NULL_INDEX

    result = _get_combined_index([])
    assert(result is NULL_INDEX)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                         # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
