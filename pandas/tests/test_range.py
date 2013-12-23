import unittest

import numpy as np
import pandas.util.testing as tm
from pandas.core.index import Index, Int64Index
from pandas.core.range import RangeIndex
import pandas.compat as compat

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
        # self.assertEqual(r.left, 10)
        # self.assertEqual(r.right, 5)
        self.assertEqual(r.step, -1)

        # ascending
        r2 = RangeIndex(5, 10)
        self.assertEqual(r2.start, 5)
        self.assertEqual(r2.stop, 10)
        # self.assertEqual(r2.left, 5)
        # self.assertEqual(r2.right, 10)
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

    def test_bad_input(self):
        with tm.assertRaisesRegexp(TypeError, 'Must be integer'):
            RangeIndex(0, 1.25)

        with tm.assertRaisesRegexp(TypeError, 'invalid literal'):
            RangeIndex(0, 'a')

        with tm.assertRaisesRegexp(TypeError, 'Must be integer'):
            RangeIndex('0', '5')


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
            assert_almost_equal(ind[0:0], Index([], dtype='int64'))
            assert_almost_equal(ind[8:8], Index([], dtype='int64'))

    def test_slicing_with_step_of_1(self):
        # [::-1] yields self but reversed
        rng1 = RangeIndex(-5, 5)
        rev1 = rng1[::-1]
        self.assertEqual(list(rev1), list(range(4, -6, -1)))
        self.assert_(rev1.equals(RangeIndex(4, -6)))
        self.assert_(rev1.equals(Index(np.arange(4, -6, -1))))

        rng2 = RangeIndex(20, 10)
        rev2 = rng2[::-1]
        self.assertEqual(list(rev2), list(range(11, 21, 1)))
        self.assert_(rev2.equals(RangeIndex(11, 21)))
        self.assert_(rev2.equals(Index(np.arange(11, 21, 1))))

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

    def test_is__continuity(self):
        # is should work on views/copies
        # is should not work with two separately constructed indices
        # is should be False when reversed or sliced
        pass

    def test_equals(self):
        # should work on views/copies
        # should be equal when separately constructed
        # should not be equal when reversed/reduced/etc
        pass

    def test_error_on_bool(self):
        self.assertRaises(ValueError, bool, RangeIndex(1, 5))
        self.assertRaises(ValueError, bool, RangeIndex(-10, -9))
        self.assertRaises(ValueError, bool, RangeIndex(1, 2))

    def test_all_and_any(self):
        zero_only = [RangeIndex(0, 1), RangeIndex(0, -1)]
        assert not any(x.any() for x in zero_only)
        assert not any(x.all() for x in zero_only)
        assert RangeIndex(5, 10).any()
        assert RangeIndex(5, 10).all()
        assert not RangeIndex(-5, 5).all()
        assert RangeIndex(-5, 5).any()
        assert RangeIndex(-3, -1).any()
        assert not RangeIndex(-3, 1).all()
        assert RangeIndex(-3, 0).all()

    #
    # indexing ops
    #
    def test_get_indexer(self):
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

    def test_range_index_from_range(self):
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

    def test_nonzero(self):
        r1 = RangeIndex(0, 5)
        a1 = np.arange(0, 5)
        assert_almost_equal(r1.nonzero(), a1.nonzero())
        r2 = RangeIndex(5, 0)
        a2 = np.arange(5, 0, -1)
        assert_almost_equal(r2.nonzero(), a2.nonzero())
        assert_almost_equal(RangeIndex(-10, -5).nonzero(),
                            np.arange(-10, -5).nonzero())

    def test_get_loc(self):
        pass

    def test_groupby(self):
        pass

    def test_slice_locs(self):
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
        self.assertRaises(KeyError, lambda : idx[::-1].slice_locs(-5, 50))

        idx2 = RangeIndex(5, 1)
        self.assertRaises(KeyError, idx2.slice_locs, 8, 2)
        self.assertRaises(KeyError, idx2.slice_locs, 7, 3)

    #
    # Index inference
    #

    @knownfail
    def test_sorted_index_yields_range(self):
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


class TestRangeIndexInt64Compat(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.index = RangeIndex(10, 20)

    def test_too_many_names(self):
        with tm.assertRaisesRegexp(ValueError, "^Length"):
            self.index.names = ["roger", "harold"]

    def test_constructor(self):
        # TODO: Fill this in
        raise AssertionError("Decide what to do here!")
        # scalar raise Exception
        self.assertRaises(TypeError, RangeIndex, 5)

    def test_basic_properties(self):
        self.assertTrue(self.index.is_unique)
        self.assertTrue(self.index.is_monotonic)
        self.assertEqual(self.index.dtype, np.int64)

    def test_basic_functions(self):
        self.assertTrue(self.index.is_numeric())
        self.assertTrue(self.index.is_integer())
        self.assertTrue(self.index.holds_integer())
        self.assertFalse(self.index.is_mixed())
        self.assertFalse(self.index.is_floating())

        self.assertEqual(self.index.nlevels, 1)
        self.assertEqual(self.index.inferred_type, 'integer')
        self.assertEqual(self.get_duplicates(), [])

    def test_hash_error(self):
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(self.index).__name__):
            hash(self.index)

    def test_copy(self):
        i = RangeIndex(0, 1, name='Foo')
        i_copy = i.copy()
        self.assert_(i_copy.name == 'Foo')

    def test_view(self):
        i = RangeIndex(0, 1, name='Foo')
        i_view = i.view()
        self.assert_(i_view.name == 'Foo')

    def test_dtype(self):
        self.assert_(self.index.dtype == np.int64)

    def test_is_monotonic(self):
        # monotonic is monotonically *increasing*
        self.assertTrue(RangeIndex(0, 5).is_monotonic)
        self.assertFalse(RangeIndex(5, 0).is_monotonic)
        self.assertFalse(RangeIndex(-5, 5)[::-1].is_monotonic)
        # TODO: If you have empty Index, need to match regular Index
        # self.assertTrue(RangeIndex(0, 0).is_monotonic)

    def test_equals(self):
        same_values = Index(self.index, dtype=object)
        self.assert_(self.index.equals(same_values))
        self.assert_(same_values.equals(self.index))

    def test_identical(self):
        i = self.index.copy()
        same_values = RangeIndex(i.left, i.right, i.step, name=i.name)
        self.assert_(i.identical(same_values))
        int64_values = Int64Index(list(i), name=i.name)
        self.assertFalse(i.identical(int64_values))
        self.assertFalse(int64_values.identical(i))

        i = self.index.copy()
        i = i.rename('foo')
        same_values = RangeIndex(i.left, i.right, i.step)
        # no name passed through constructor
        self.assert_(same_values.identical(self.index))
        self.assertFalse(i.identical(same_values))

    def test_get_indexer(self):
        def test_indexer(target, expected):
            indexer = self.index.get_indexer(target)
            self.assert_(np.array_equal(indexer, expected))

        test_indexer(RangeIndex(-5, 5),
                     np.array([-1] * 10))

        test_indexer(RangeIndex(5, 15),
                     np.array([-1, -1, -1, -1, -1, 0, 1, 2, 3, 4]))

        test_indexer(Index(list('abcd') + [11]),
                     np.array([-1, -1, -1, -1, 1]))

        test_indexer(Index([0.5, 0.25, 1, 18.0, 17]),
                     np.array([-1, -1, -1, 8, 7]))

    def test_get_indexer_fails_non_monotonic(self):
        ind = RangeIndex(10, 5, -1)

        with tm.assertRaisesRegexp(ValueError, 'monotonic for backward fill'):
            ind.get_indexer(Index([0]), method='bfill')
        with tm.assertRaisesRegexp(ValueError, 'monotonic for backward fill'):
            ind.get_indexer(Index([0]), method='backfill')

        with tm.assertRaisesRegexp(ValueError, 'monotonic for forward fill'):
            ind.get_indexer(Index([0]), method='ffill')
        with tm.assertRaisesRegexp(ValueError, 'monotonic for forward fill'):
            ind.get_indexer(Index([0]), method='pad')

    def test_get_indexer_pad_int64_with_range(self):
        # TODO: Move this over to Int64Index tests instead
        target = RangeIndex(0, 10)
        idx = Index(range(0, 20, 2))
        indexer = idx.get_indexer(target, method='pad')
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        self.assert_(np.array_equal(indexer, expected))

    def test_get_indexer_pad(self):
        idx = RangeIndex(-3, 0)
        target = Index([-2, -1, 0, 1, 2])
        indexer = idx.get_indexer(target, method='pad')
        expected = np.array([1, 2, 2, 2, 2])
        self.assert_(np.array_equal(indexer, expected))

        target2 = Index([-4, -2, 1])
        indexer = idx.get_indexer(target2, method='pad')
        expected = np.array([-1, 1, 2])
        self.assert_(np.array_equal(indexer, expected))

    def test_get_indexer_backfill_int64_with_range(self):
        # TODO: Move this over to Int64Index tests instead
        target = RangeIndex(0, 10)
        idx = Index(range(0, 20, 2))
        indexer = idx.get_indexer(target, method='backfill')
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5])
        self.assert_(np.array_equal(indexer, expected))

    def test_get_indexer_backfill(self):
        idx = RangeIndex(3, 5)
        target = Index([-4, -2, 3, 4, 5, 7])
        indexer = idx.get_indexer(target, method='backfill')
        expected = np.array([0, 0, 1, 2, -1, -1])
        self.assert_(np.array_equal(indexer, expected))

        # # TODO: Decide on ffill, bfill, pad for NON-integer with RangeIndex...

    def test_join_outer(self):
        # TODO: Convert this to RangeIndex formatted
        # 1. Write tests for take
        # 2. Make sure this works with return indexers (which are just args to
        # take).
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
        # TODO: Convert this to RangeIndex formatted
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
        # TODO: Convert this to RangeIndex formatted
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
        # UPDATED
        idx = RangeIndex(0, 15)
        other = Index([3, 6, 7, 8, 10], dtype=object)

        outer = idx.join(other, how='outer')
        outer2 = other.join(idx, how='outer')
        expected = Index(range(0, 15), dtype=object)
        self.assert_(outer.equals(outer2))
        self.assert_(outer.equals(expected))

        inner = idx.join(other, how='inner')
        inner2 = other.join(idx, how='inner')
        expected = other.copy() # avoid is_ stuff
        self.assert_(inner.equals(inner2))
        self.assert_(inner.equals(expected))

        idx2 = RangeIndex(0, 4)
        inner3 = idx2.join(other, how='inner')
        inner4 = other.join(idx2, how='inner')
        expected = Index([3], dtype=object)
        self.assert_(inner3.equals(inner4))
        self.assert_(inner3.equals(expected))


        left = idx.join(other, how='left')
        self.assert_(left.equals(idx))

        left2 = other.join(idx, how='left')
        self.assert_(left2.equals(other))

        right = idx.join(other, how='right')
        self.assert_(right.equals(other))

        right2 = other.join(idx, how='right')
        self.assert_(right2.equals(idx))

    def test_join_self(self):
        # UPDATED
        idx = RangeIndex(-10, -4)
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = idx.join(idx, how=kind)
            self.assert_(self.index is joined)

    def test_intersection(self):
        # TODO: Convert this to RangeIndex formatted
        other = Index([1, 2, 3, 4, 5])
        result = self.index.intersection(other)
        expected = np.sort(np.intersect1d(self.index.values, other.values))
        self.assert_(np.array_equal(result, expected))

        result = other.intersection(self.index)
        expected = np.sort(np.asarray(np.intersect1d(self.index.values,
                                                     other.values)))
        self.assert_(np.array_equal(result, expected))

    def test_union_noncomparable(self):
        # TODO: Convert this to RangeIndex formatted
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

    # def test_view_Index(self):
    #     self.index.view(Index)

    def test_prevent_casting(self):
        # TODO: Convert this to RangeIndex formatted
        result = self.index.astype('O')
        self.assert_(result.dtype == np.object_)

    def test_take_preserve_name(self):
        # TODO: Convert this to RangeIndex formatted
        index = RangeIndex(1, 4, name='foo')
        taken = index.take([3, 0, 1])
        self.assertEqual(index.name, taken.name)

    def test_int_name_format(self):
        from pandas import Series, DataFrame
        index = RangeIndex(3, 0, -1, name=0)
        s = Series(lrange(3), index)
        df = DataFrame(lrange(3), index=index)
        repr(s)
        repr(df)

    def test_repr_roundtrip(self):
        tm.assert_index_equal(eval(repr(self.index)), self.index)

    def test_unicode_string_with_unicode(self):
        idx = RangeIndex(0, 1000)

        if compat.PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        idx = RangeIndex(0, 1000)
        if compat.PY3:
            bytes(idx)
        else:
            str(idx)

    def test_slice_keep_name(self):
        idx = RangeIndex(1, 3, name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs'],
                   exit=False)
