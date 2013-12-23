import unittest

import numpy as np
import pandas.util.testing as tm
from pandas.core.index import Index
from pandas.core.range import RangeIndex

lrange = lambda *args: list(range(*args))

def assert_almost_equal(a, b):
    try:
        tm.assert_almost_equal(a, b)
    except:
        print(a, b)
        raise


def knownfail(f):
    def wrapper():
        try:
            f()
        except Exception as e:
            print("%s: KNOWN FAILURE: %r" % (f.__name__, e))
        else:
            raise AssertionError("Known failure passed! %s" % f.__name__)
    return wrapper


class self(object):

    """Fake for tests!"""

    @staticmethod
    def assertEquals(a, b):
        assert a == b, "%r != %r" % (a, b)

    assertEqual = assertEquals

    @staticmethod
    def assertRaises(exc, f, *args, **kwargs):
        try:
            f(*args, **kwargs)
        except exc:
            return True
        else:
            raise AssertionError(
                "Expected exception of type %s to be raised!" % exc)


class TestRangeIndex(unittest.TestCase):

    def test_basic(self):
        self.assert_(not RangeIndex(1, 0).ascending)
        self.assert_(RangeIndex(0, 100).ascending)
        # make sure conditions work correctly
        # descending
        r = RangeIndex(10, 5)
        self.assertEqual(r.start, 5)
        self.assertEqual(r.stop, 10)
        self.assertEqual(r.left, 10)
        self.assertEqual(r.right, 5)
        self.assertEqual(r.step, -1)

        # ascending
        r2 = RangeIndex(5, 10)
        self.assertEqual(r2.start, 5)
        self.assertEqual(r2.stop, 10)
        self.assertEqual(r2.left, 5)
        self.assertEqual(r2.right, 10)
        self.assertEqual(r2.step, 1)

        # negative values
        r3 = RangeIndex(-10, -9)
        self.assertEqual(r3.start, -10)
        self.assertEqual(r3.stop, -9)
        self.assert_(r3.ascending)
        self.assertEqual(r3.step, 1)

        r4 = RangeIndex(-8, -15)
        self.assertEqual(r4.start, -15)
        self.assertEqual(r4.stop, -8)
        self.assert_(not r4.ascending)
        self.assertEqual(r4.step, -1)

    def test_contains(self):

        r = RangeIndex(10, 5)
        r2 = RangeIndex(5, 10)
        for i in range(5, 10):
            self.assert_(i in r)
            self.assert_(i in r2)

    def test_empty(self):
        assert np.array_equal(RangeIndex(5, 5), Index([], dtype='int64'))

    def test_asarray(self):
        # __array__
        self.assert_(np.array_equal(np.asarray(RangeIndex(1, 0)),
                                    np.array([1])))
        self.assert_(np.array_equal(np.asarray(RangeIndex(0, 100)),
                                    np.arange(0, 100)))
        self.assert_(np.array_equal(RangeIndex(1, 0).values, np.array([1])))
        self.assert_(np.array_equal(RangeIndex(0, 100).values,
                                    np.arange(0, 100)))

    def test_set_ops(self):
        r1 = RangeIndex(1, 10)
        r2 = RangeIndex(5, 10)
        self.assert_(r1._overlaps(r2))
        self.assert_(r2._overlaps(r1))
        # union and intersection - underlying methods)
        self.assert_(r1.intersection(r2).equals(RangeIndex(5, 10)))
        self.assert_(r2.intersection(r1).equals(r1.intersection(r2)))
        self.assert_((r1 & r2).equals(RangeIndex(5, 10)))
        self.assert_((r2 & r1).equals((r1 & r2)))
        self.assert_(r1.union(r2).equals(r2.union(r1)))
        self.assert_(r1.union(r2).equals(r1))
        # union and intersection - with infix operators)
        self.assert_((r1 + r2).equals((r2 + r1)))
        self.assert_((r1 + r2).equals(r1))
        self.assert_((r1 | r2).equals((r2 | r1)))
        self.assert_((r1 | r2).equals(r1))

        # difference - underlying method)
        self.assert_(r1.difference(r2).equals(RangeIndex(1, 5)))
        self.assert_(r2.difference(r1).equals(Index([], dtype='int64')))
        self.assert_(r1.difference(r1).equals(Index([], dtype='int64')))
        self.assert_(r2.difference(r2).equals(Index([], dtype='int64')))
        # difference - with infix operator)
        self.assert_((r1 - r2).equals(RangeIndex(1, 5)))
        self.assert_((r2 - r1).equals(Index([], dtype='int64')))
        self.assert_((r1 - r1).equals(Index([], dtype='int64')))
        self.assert_((r2 - r2).equals(Index([], dtype='int64')))

    def test_getitem_and_iter(self):
        # basic container ops
        pairs = [(RangeIndex(-10, -5), lrange(-10, -5)),
                 (RangeIndex(8, 5), lrange(8, 5, -1)),
                 (RangeIndex(0, 10), lrange(0, 10)),
                 (RangeIndex(-5, 5), lrange(-5, 5)),
                 (RangeIndex(3, -15), lrange(3, -15, -1))]

        for ind, rng in pairs:
            try:
                self.assertEqual(len(ind), len(rng))
                for i in range(len(rng)):
                    self.assertEqual(ind[i], rng[i])
                    self.assertEqual(ind[-i], rng[-i])
            except:
                print(i, ind, ind[i])
                print(i, rng, rng[i])
                raise
            # basic __iter__ test
            assert_almost_equal(list(ind), rng)
            assert np.array_equal(ind.values, np.array(list(rng)))

        cases = 10
        for ind in zip(*pairs)[0]:
            length = len(ind)
            # edges
            self.assertRaises(IndexError, lambda: ind[length])
            self.assertRaises(IndexError, lambda: ind[-length - 1])
            for _ in range(cases):
                i = np.random.randint(1, 100)
                self.assertRaises(IndexError, lambda: ind[length + i])
                self.assertRaises(IndexError, lambda: ind[-length - 1 - i])

    def test_slicing(self):
        pairs = [(RangeIndex(-10, -5), lrange(-10, -5)),  # can remove later
                 (RangeIndex(8, 5), np.arange(8, 5, -1)),
                 (RangeIndex(0, 10), lrange(0, 10)),  # can remove later
                 (RangeIndex(-3, 3), np.arange(-3, 3)),
                 (RangeIndex(3, -2), np.arange(3, -2, -1))]
        # TODO: This is incredibly slow - pick something smaller to work with
        for ind, rng in pairs:
            assert_almost_equal(ind[:], rng[:])
            for i, j in [(i, j) for i in range(len(rng))
                         for j in range(len(rng)) if i >= j]:
                assert_almost_equal(ind[i:], rng[i:])
                assert_almost_equal(ind[:i], rng[:i])
                assert_almost_equal(ind[-i:], rng[-i:])
                assert_almost_equal(ind[:-i], rng[:-i])

                assert_almost_equal(ind[i:j], rng[i:j])
                assert_almost_equal(ind[i:-j], rng[i:-j])
                assert_almost_equal(ind[-i:-j], rng[-i:-j])
                assert_almost_equal(ind[-i:j], rng[-i:j])

                assert_almost_equal(ind[j:i], rng[j:i])
                assert_almost_equal(ind[j:-i], rng[j:-i])
                assert_almost_equal(ind[-j:-i], rng[-j:-i])
                assert_almost_equal(ind[-j:i], rng[-j:i])
        # in range
        # - forward
        # - reversed
        # totally out of range
        # - forward/reversed below
        # - forward/reversed above
        # partial in range
        # - forward/reversed with low
        # - forward/reversed with high
        # [:] yields (shallow copy of) self
        # Empty slice yields Index([], dtype='int64')
        pass

    def test_slicing_with_step_of_1(self):
        # [::-1] yields self but reversed
        r1 = RangeIndex(-5, 5)
        r2 = RangeIndex(20, 10)
        self.assert_(r1[::-1].equals(RangeIndex(5, -5)))
        self.assert_(r2[::-1].equals(RangeIndex(10, 20)))
        self.assert_(r1[::1].equals(r1))
        self.assert_(r2[::1].equals(r2))

    def test_slicing_with_other_steps(self):
        pass

    def test_immutable(self):
        # setitem
        # setslice
        pass

#
# PandasObject properties
#

    def test_copy_and_view(self):
        # shallow / deep copy should be same
        pass

    def test_is__continuity():
        # is should work on views/copies
        # is should not work with two separately constructed indices
        # is should be False when reversed or sliced
        pass

    def test_equals():
        # should work on views/copies
        # should be equal when separately constructed
        # should not be equal when reversed/reduced/etc
        pass

    def test_error_on_bool(self):
        self.assertRaises(ValueError, bool, RangeIndex(1, 5))
        self.assertRaises(ValueError, bool, RangeIndex(-10, -9))
        self.assertRaises(ValueError, bool, RangeIndex(1, 2))


#
# indexing ops
#
@knownfail
def test_get_indexer():
    idx1 = RangeIndex(1, 5)
    # TODO: Consider supporting steps
    idx2 = Index([2, 4, 6])
    idx3 = Index([1, 6, 7, 1, 2])

    r1 = idx1.get_indexer(idx2)
    assert_almost_equal(r1, [1, 3, -1])

    r1 = idx1.get_indexer(idx3)
    assert_almost_equal(r1, np.array([0, -1, -1,  0,  1]))

    r1 = idx1.get_indexer(idx3, method='pad')
    assert_almost_equal(r1, np.array([0,  3,  3, -1, -1]))

    rffill1 = idx1.get_indexer(idx3, method='ffill')
    assert_almost_equal(r1, rffill1)

    r1 = idx1.get_indexer(idx3, method='backfill')
    assert_almost_equal(r1, np.array([0, -1, -1,  0,  1]))

    rbfill1 = idx1.get_indexer(idx3, method='bfill')
    assert_almost_equal(r1, rbfill1)

    # r1 = idx3.get_indexer(idx1, method='pad')
    # assert_almost_equal(r1, [0, 0, 0, 0, 0])

    # rffill1 = idx3.get_indexer(idx1, method='ffill')

    # r1 = idx3.get_indexer(idx1, method='backfill')
    # assert_almost_equal(r1, [0, -1, -1, -1, -1])

    # rbfill1 = idx3.get_indexer(idx1, method='bfill')
    # assert_almost_equal(r1, rbfill1)


@knownfail
def test_range_index_from_range():
    def assert_fails(inpt):
        res = RangeIndex.possibly_convert_array(inpt)
        assert res is None, "Expected %r to return None" % inpt

    def assert_converts(inpt, expected):
        res = RangeIndex.possibly_convert_array(inpt)
        assert expected.equals(res), "With input %r, %r != %r" % (inpt, res,
                                                                  expected)
    assert_converts(range(5), RangeIndex(0, 5))
    assert_fails([1, 3, 7, 5])
    assert_fails([4, 10, 11, 13])
    assert_converts(np.arange(50, 40, -1), RangeIndex(50, 40))
    assert_converts([0], RangeIndex(0, 1))
    assert_fails([])

    # dupe values
    assert_fails([10, 9, 8, 7, 10])
    assert_fails([1, 2, 3, 4, 5, 7])

    # should not try to convert dtype (caller responsibility)
    arr = np.arange(5, 15)
    assert_converts(arr, RangeIndex(5, 15))
    assert_fails(arr.astype(float))

    # works with resort
    assert_fails([-10, -5, -6, -7, -2, -3, -4, -8, -9])
    assert_fails([9, 8, 5, 7, 6])

    # possibilities that *won't* work now but could in the future
    # (i.e., nested ranges, steps)
    assert_fails([15, 13, 11, 9, 7, 5])
    assert_fails([1, 2, 3, 8, 9, 10])
    assert_fails([2, 4, 6, 8, 10, 12])


def test_nonzero():
    r1 = RangeIndex(0, 5)
    a1 = np.arange(0, 5)
    assert_almost_equal(r1.nonzero(), a1.nonzero())
    r2 = RangeIndex(5, 0)
    a2 = np.arange(5, 0, -1)
    assert_almost_equal(r2.nonzero(), a2.nonzero())


def test_get_loc():
    pass


def test_groupby():
    pass


@knownfail
def test_slice_locs():
    idx = RangeIndex(0, 11)
    n = len(idx)

    self.assertEquals(idx.slice_locs(start=2), (2, n))
    self.assertEquals(idx.slice_locs(start=3), (3, n))
    self.assertEquals(idx.slice_locs(3, 8), (3, 8))
    self.assertEquals(idx.slice_locs(5, 10), (3, n))
    self.assertEquals(idx.slice_locs(end=8), (0, 8))
    self.assertEquals(idx.slice_locs(end=9), (0, 9))
    # monotonic *increasing* indexes allow slice_locs that aren't in the Index
    self.assertEquals(idx.slice_locs(-5, 50), (0, 11))

    idx2 = RangeIndex(5, 1)
    self.assertRaises(KeyError, idx2.slice_locs, 8, 2)
    self.assertRaises(KeyError, idx2.slice_locs, 7, 3)

#
# Index inference
#


@knownfail
def test_sorted_index_yields_range():
    ind = Index(range(10))
    assert isinstance(ind, RangeIndex)
    assert ind.equals(RangeIndex(0, 10))
    ind = Index(range(15, -1, -1)),
    assert isinstance(ind, RangeIndex)
    assert ind.equals(RangeIndex(15, -1))
    ind = Index([1, 3, 5, 7])
    assert not isinstance(ind, RangeIndex)
    ind = Index(range(5) + [6])
    assert not isinstance(ind, RangeIndex)
    ind = Index([1, 3, 2, 4, 5])
    assert not isinstance(ind, RangeIndex)
    ind = Index(np.arange(0, 10).astype(float))
    assert not isinstance(ind, RangeIndex)
