from datetime import timedelta
from pandas.core.index import Index
import pandas.util.testing as common
import pandas.lib.tseries as tseries
import numpy as np
import os
import pickle
import unittest

class TestIndex(unittest.TestCase):

    def setUp(self):
        self.strIndex = common.makeStringIndex(100)
        self.dateIndex = common.makeDateIndex(100)
        self.intIndex = common.makeIntIndex(100)

    def test_deepcopy(self):
        from copy import deepcopy

        copy = deepcopy(self.strIndex)
        self.assert_(copy is self.strIndex)

    def test_duplicates(self):
        self.assertRaises(Exception, Index, [0, 0, 0])

    def test_sort(self):
        self.assertRaises(Exception, self.strIndex.sort)

    def test_mutability(self):
        self.assertRaises(Exception, self.strIndex.__setitem__, 5, 0)
        self.assertRaises(Exception, self.strIndex.__setitem__, slice(1,5), 0)

    def test_constructor(self):
        # regular instance creation
        common.assert_contains_all(self.strIndex, self.strIndex)
        common.assert_contains_all(self.dateIndex, self.dateIndex)

        # casting
        arr = np.array(self.strIndex)
        index = arr.view(Index)
        common.assert_contains_all(arr, index)
        self.assert_(np.array_equal(self.strIndex, index))

        # corner case
        self.assertRaises(Exception, Index, 0)

        # arr = np.array(5.)
        # self.assertRaises(Exception, arr.view, Index)

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

    def test_asOfDate(self):
        d = self.dateIndex[0]
        self.assert_(self.dateIndex.asOfDate(d) is d)
        self.assert_(self.dateIndex.asOfDate(d - timedelta(1)) is None)

        d = self.dateIndex[-1]
        self.assert_(self.dateIndex.asOfDate(d + timedelta(1)) is d)

    def test_argsort(self):
        result = self.strIndex.argsort()
        expected = np.array(self.strIndex).argsort()
        self.assert_(np.array_equal(result, expected))

    def test_comparators(self):
        index = self.dateIndex
        element = index[len(index) // 2]
        arr = np.array(index)

        self.assert_(np.array_equal(arr == element, index == element))
        self.assert_(np.array_equal(arr > element, index > element))
        self.assert_(np.array_equal(arr < element, index < element))
        self.assert_(np.array_equal(arr >= element, index >= element))
        self.assert_(np.array_equal(arr <= element, index <= element))

    def test_booleanindex(self):
        boolIdx = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIdx[5:30:2] = False

        subIndex = self.strIndex[boolIdx]
        common.assert_dict_equal(tseries.map_indices(subIndex),
                                 subIndex.indexMap)

    def test_fancy(self):
        sl = self.strIndex[[1,2,3]]
        for i in sl:
            self.assertEqual(i, sl[sl.indexMap[i]])

    def test_getitem(self):
        arr = np.array(self.dateIndex)
        self.assertEquals(self.dateIndex[5], arr[5])

    def test_add(self):
        firstCat = self.strIndex + self.dateIndex
        secondCat = self.strIndex + self.strIndex

        self.assert_(common.equalContents(np.append(self.strIndex,
                                                    self.dateIndex), firstCat))
        self.assert_(common.equalContents(secondCat, self.strIndex))
        common.assert_contains_all(self.strIndex, firstCat.indexMap)
        common.assert_contains_all(self.strIndex, secondCat.indexMap)
        common.assert_contains_all(self.dateIndex, firstCat.indexMap)

        # this is valid too
        shifted = self.dateIndex + timedelta(1)

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)

        self.assert_(common.equalContents(intersect, second))

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
        self.assert_(common.equalContents(union, everything))

        # Corner cases
        union = first.union(first)
        self.assert_(union is first)

        union = first.union([])
        self.assert_(union is first)

        # non-iterable input
        self.assertRaises(Exception, first.union, 0.5)

    def test_diff(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]
        result = first - second

        self.assert_(common.equalContents(result, answer))

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

            common.assert_dict_equal(unpickled.indexMap, index.indexMap)

        testit(self.strIndex)
        testit(self.dateIndex)

