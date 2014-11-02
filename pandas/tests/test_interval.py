import numpy as np

from pandas.core.interval import Interval, IntervalIndex
from pandas.core.index import Index
from pandas.lib import IntervalTree

import pandas.util.testing as tm
import pandas as pd


class TestInterval(tm.TestCase):
    def setUp(self):
        self.interval = Interval(0, 1)

    def test_properties(self):
        self.assertEqual(self.interval.closed, 'right')
        self.assertEqual(self.interval.left, 0)
        self.assertEqual(self.interval.right, 1)
        self.assertEqual(self.interval.mid, 0.5)

    def test_repr(self):
        self.assertEqual(repr(self.interval),
                         "Interval(0, 1, closed='right')")
        self.assertEqual(str(self.interval), "(0, 1]")

        interval_left = Interval(0, 1, closed='left')
        self.assertEqual(repr(interval_left),
                         "Interval(0, 1, closed='left')")
        self.assertEqual(str(interval_left), "[0, 1)")

    def test_contains(self):
        self.assertIn(0.5, self.interval)
        self.assertIn(1, self.interval)
        self.assertNotIn(0, self.interval)
        self.assertRaises(TypeError, lambda: self.interval in self.interval)

        interval = Interval(0, 1, closed='both')
        self.assertIn(0, interval)
        self.assertIn(1, interval)

        interval = Interval(0, 1, closed='neither')
        self.assertNotIn(0, interval)
        self.assertIn(0.5, interval)
        self.assertNotIn(1, interval)

    def test_equal(self):
        self.assertEqual(Interval(0, 1), Interval(0, 1, closed='right'))
        self.assertNotEqual(Interval(0, 1), Interval(0, 1, closed='left'))
        self.assertNotEqual(Interval(0, 1), 0)

    def test_comparison(self):
        with self.assertRaisesRegexp(TypeError, 'unorderable types'):
            Interval(0, 1) < 2

        self.assertTrue(Interval(0, 1) < Interval(1, 2))
        self.assertTrue(Interval(0, 1) < Interval(0, 2))
        self.assertTrue(Interval(0, 1) < Interval(0.5, 1.5))
        self.assertTrue(Interval(0, 1) <= Interval(0, 1))
        self.assertTrue(Interval(0, 1) > Interval(-1, 2))
        self.assertTrue(Interval(0, 1) >= Interval(0, 1))

    def test_hash(self):
        # should not raise
        hash(self.interval)

    # def test_math(self):
    #     expected = Interval(1, 2)
    #     actual = self.interval + 1
    #     self.assertEqual(expected, actual)


class TestIntervalTree(tm.TestCase):
    def setUp(self):
        self.tree = IntervalTree(np.arange(5), np.arange(5) + 2)

    def test_get_loc(self):
        self.assert_numpy_array_equal(self.tree.get_loc(1), [0])
        self.assert_numpy_array_equal(np.sort(self.tree.get_loc(2)), [0, 1])
        with self.assertRaises(KeyError):
            self.tree.get_loc(-1)

    def test_get_indexer(self):
        self.assert_numpy_array_equal(
            self.tree.get_indexer(np.array([1.0, 5.5, 6.5])), [0, 4, -1])
        with self.assertRaises(KeyError):
            self.tree.get_indexer(np.array([3.0]))

    def test_get_indexer_non_unique(self):
        indexer, missing = self.tree.get_indexer_non_unique(
            np.array([1.0, 2.0, 6.5]))
        self.assert_numpy_array_equal(indexer[:1], [0])
        self.assert_numpy_array_equal(np.sort(indexer[1:3]), [0, 1])
        self.assert_numpy_array_equal(np.sort(indexer[3:]), [-1])
        self.assert_numpy_array_equal(missing, [2])

    def test_duplicates(self):
        tree = IntervalTree([0, 0, 0], [1, 1, 1])
        self.assert_numpy_array_equal(np.sort(tree.get_loc(0.5)), [0, 1, 2])

        with self.assertRaises(KeyError):
            tree.get_indexer(np.array([0.5]))

        indexer, missing = tree.get_indexer_non_unique(np.array([0.5]))
        self.assert_numpy_array_equal(np.sort(indexer), [0, 1, 2])
        self.assert_numpy_array_equal(missing, [])

    def test_get_loc_closed(self):
        for closed in ['left', 'right', 'both', 'neither']:
            tree = IntervalTree([0], [1], closed=closed)
            for p, errors in [(0, tree.open_left),
                              (1, tree.open_right)]:
                if errors:
                    with self.assertRaises(KeyError):
                        tree.get_loc(p)
                else:
                    self.assert_numpy_array_equal(tree.get_loc(p),
                                                  np.array([0]))

    def test_get_indexer_closed(self):
        x = np.arange(1000)
        found = x
        not_found = -np.ones(1000)
        for closed in ['left', 'right', 'both', 'neither']:
            tree = IntervalTree(x, x + 0.5, closed=closed)
            self.assert_numpy_array_equal(found, tree.get_indexer(x + 0.25))

            expected = found if tree.closed_left else not_found
            self.assert_numpy_array_equal(expected, tree.get_indexer(x + 0.0))

            expected = found if tree.closed_right else not_found
            self.assert_numpy_array_equal(expected, tree.get_indexer(x + 0.5))


class TestIntervalIndex(tm.TestCase):
    def setUp(self):
        self.index = IntervalIndex([0, 1], [1, 2])

    def test_constructors(self):
        expected = self.index
        actual = IntervalIndex.from_breaks(np.arange(3), closed='right')
        self.assertTrue(expected.equals(actual))

        alternate = IntervalIndex.from_breaks(np.arange(3), closed='left')
        self.assertFalse(expected.equals(alternate))

        actual = IntervalIndex.from_intervals([Interval(0, 1), Interval(1, 2)])
        self.assertTrue(expected.equals(actual))

        self.assertRaises(ValueError, IntervalIndex, [0], [1], closed='invalid')

        # TODO: fix all these commented out tests (here and below)

        intervals = [Interval(0, 1), Interval(1, 2, closed='left')]
        with self.assertRaises(ValueError):
            IntervalIndex.from_intervals(intervals)

        with self.assertRaises(ValueError):
            IntervalIndex([0, 10], [3, 5])

        actual = Index([Interval(0, 1), Interval(1, 2)])
        self.assertIsInstance(actual, IntervalIndex)
        self.assertTrue(expected.equals(actual))

        actual = Index(expected)
        self.assertIsInstance(actual, IntervalIndex)
        self.assertTrue(expected.equals(actual))

        # no point in nesting periods in an IntervalIndex
        # self.assertRaises(ValueError, IntervalIndex.from_breaks,
        #                   pd.period_range('2000-01-01', periods=3))

    def test_properties(self):
        self.assertEqual(len(self.index), 2)
        self.assertEqual(self.index.size, 2)

        self.assert_numpy_array_equal(self.index.left, [0, 1])
        self.assertIsInstance(self.index.left, Index)

        self.assert_numpy_array_equal(self.index.right, [1, 2])
        self.assertIsInstance(self.index.right, Index)

        self.assert_numpy_array_equal(self.index.mid, [0.5, 1.5])
        self.assertIsInstance(self.index.mid, Index)

        self.assertEqual(self.index.closed, 'right')

        expected = np.array([Interval(0, 1), Interval(1, 2)], dtype=object)
        self.assert_numpy_array_equal(np.asarray(self.index), expected)
        self.assert_numpy_array_equal(self.index.values, expected)

    def test_copy(self):
        actual = self.index.copy()
        self.assertTrue(actual.equals(self.index))

        actual = self.index.copy(deep=True)
        self.assertTrue(actual.equals(self.index))
        self.assertIsNot(actual.left, self.index.left)

    def test_delete(self):
        expected = IntervalIndex.from_breaks([1, 2])
        actual = self.index.delete(0)
        self.assertTrue(expected.equals(actual))

    def test_insert(self):
        expected = IntervalIndex.from_breaks(range(4))
        actual = self.index.insert(2, Interval(2, 3))
        self.assertTrue(expected.equals(actual))

        self.assertRaises(ValueError, self.index.insert, 0, 1)
        self.assertRaises(ValueError, self.index.insert, 0,
                          Interval(2, 3, closed='left'))

    def test_take(self):
        actual = self.index.take([0, 1])
        self.assertTrue(self.index.equals(actual))

        expected = IntervalIndex([0, 0, 1], [1, 1, 2])
        actual = self.index.take([0, 0, 1])
        self.assertTrue(expected.equals(actual))

    def test_monotonic_and_unique(self):
        self.assertTrue(self.index.is_monotonic)
        self.assertTrue(self.index.is_unique)

        idx = IntervalIndex.from_tuples([(0, 1), (0.5, 1.5)])
        self.assertTrue(idx.is_monotonic)
        self.assertTrue(idx.is_unique)

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 2)])
        self.assertFalse(idx.is_monotonic)
        self.assertTrue(idx.is_unique)

        idx = IntervalIndex.from_tuples([(0, 2), (0, 2)])
        self.assertFalse(idx.is_unique)
        self.assertTrue(idx.is_monotonic)

    def test_repr(self):
        expected = ("IntervalIndex(left=[0, 1],\n              right=[1, 2],"
                    "\n              closed='right')")
        IntervalIndex((0, 1), (1, 2), closed='right')
        self.assertEqual(repr(self.index), expected)

    def test_get_loc_value(self):
        self.assertRaises(KeyError, self.index.get_loc, 0)
        self.assertEqual(self.index.get_loc(0.5), 0)
        self.assertEqual(self.index.get_loc(1), 0)
        self.assertEqual(self.index.get_loc(1.5), 1)
        self.assertEqual(self.index.get_loc(2), 1)
        self.assertRaises(KeyError, self.index.get_loc, -1)
        self.assertRaises(KeyError, self.index.get_loc, 3)

        idx = IntervalIndex.from_tuples([(0, 2), (1, 3)])
        self.assertEqual(idx.get_loc(0.5), 0)
        self.assertEqual(idx.get_loc(1), 0)
        self.assert_numpy_array_equal(idx.get_loc(1.5), [0, 1])
        self.assert_numpy_array_equal(np.sort(idx.get_loc(2)), [0, 1])
        self.assertEqual(idx.get_loc(3), 1)
        self.assertRaises(KeyError, idx.get_loc, 3.5)

        idx = IntervalIndex([0, 2], [1, 3])
        self.assertRaises(KeyError, idx.get_loc, 1.5)

    def slice_locs_cases(self, breaks):
        # TODO: same tests for more index types
        index = IntervalIndex.from_breaks([0, 1, 2], closed='right')
        self.assertEqual(index.slice_locs(), (0, 2))
        self.assertEqual(index.slice_locs(0, 1), (0, 1))
        self.assertEqual(index.slice_locs(1, 1), (0, 1))
        self.assertEqual(index.slice_locs(0, 2), (0, 2))
        self.assertEqual(index.slice_locs(0.5, 1.5), (0, 2))
        self.assertEqual(index.slice_locs(0, 0.5), (0, 1))
        self.assertEqual(index.slice_locs(start=1), (0, 2))
        self.assertEqual(index.slice_locs(start=1.2), (1, 2))
        self.assertEqual(index.slice_locs(end=1), (0, 1))
        self.assertEqual(index.slice_locs(end=1.1), (0, 2))
        self.assertEqual(index.slice_locs(end=1.0), (0, 1))
        self.assertEqual(*index.slice_locs(-1, -1))

        index = IntervalIndex.from_breaks([0, 1, 2], closed='neither')
        self.assertEqual(index.slice_locs(0, 1), (0, 1))
        self.assertEqual(index.slice_locs(0, 2), (0, 2))
        self.assertEqual(index.slice_locs(0.5, 1.5), (0, 2))
        self.assertEqual(index.slice_locs(1, 1), (1, 1))
        self.assertEqual(index.slice_locs(1, 2), (1, 2))

        index = IntervalIndex.from_breaks([0, 1, 2], closed='both')
        self.assertEqual(index.slice_locs(1, 1), (0, 2))
        self.assertEqual(index.slice_locs(1, 2), (0, 2))

    def test_slice_locs_int64(self):
        self.slice_locs_cases([0, 1, 2])

    def test_slice_locs_float64(self):
        self.slice_locs_cases([0.0, 1.0, 2.0])

    def slice_locs_decreasing_cases(self, tuples):
        index = IntervalIndex.from_tuples(tuples)
        self.assertEqual(index.slice_locs(1.5, 0.5), (1, 3))
        self.assertEqual(index.slice_locs(2, 0), (1, 3))
        self.assertEqual(index.slice_locs(2, 1), (1, 3))
        self.assertEqual(index.slice_locs(3, 1.1), (0, 3))
        self.assertEqual(index.slice_locs(3, 3), (0, 2))
        self.assertEqual(index.slice_locs(3.5, 3.3), (0, 1))
        self.assertEqual(index.slice_locs(1, -3), (2, 3))
        self.assertEqual(*index.slice_locs(-1, -1))

    def test_slice_locs_decreasing_int64(self):
        self.slice_locs_cases([(2, 4), (1, 3), (0, 2)])

    def test_slice_locs_decreasing_float64(self):
        self.slice_locs_cases([(2., 4.), (1., 3.), (0., 2.)])

    def test_slice_locs_fails(self):
        index = IntervalIndex.from_tuples([(1, 2), (0, 1), (2, 3)])
        with self.assertRaises(KeyError):
            index.slice_locs(1, 2)

    def test_get_loc_interval(self):
        self.assertEqual(self.index.get_loc(Interval(0, 1)), 0)
        self.assertEqual(self.index.get_loc(Interval(0, 0.5)), 0)
        self.assertEqual(self.index.get_loc(Interval(0, 1, 'left')), 0)
        self.assertRaises(KeyError, self.index.get_loc, Interval(2, 3))
        self.assertRaises(KeyError, self.index.get_loc, Interval(-1, 0, 'left'))

    def test_get_indexer(self):
        actual = self.index.get_indexer([-1, 0, 0.5, 1, 1.5, 2, 3])
        expected = [-1, -1, 0, 0, 1, 1, -1]
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(self.index)
        expected = [0, 1]
        self.assert_numpy_array_equal(actual, expected)

        index = IntervalIndex.from_breaks([0, 1, 2], closed='left')
        actual = index.get_indexer([-1, 0, 0.5, 1, 1.5, 2, 3])
        expected = [-1, 0, 0, 1, 1, -1, -1]
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(index[:1])
        expected = [0]
        self.assert_numpy_array_equal(actual, expected)

        self.assertRaises(ValueError, self.index.get_indexer, index)

    def test_get_indexer_subintervals(self):
        # return indexers for wholly contained subintervals
        target = IntervalIndex.from_breaks(np.linspace(0, 2, 5))
        actual = self.index.get_indexer(target)
        expected = [0, 0, 1, 1]
        self.assert_numpy_array_equal(actual, expected)

        target = IntervalIndex.from_breaks([0, 0.67, 1.33, 2])
        self.assertRaises(ValueError, self.index.get_indexer, target)

        actual = self.index.get_indexer(target[[0, -1]])
        expected = [0, 1]
        self.assert_numpy_array_equal(actual, expected)

        target = IntervalIndex.from_breaks([0, 0.33, 0.67, 1], closed='left')
        actual = self.index.get_indexer(target)
        expected = [0, 0, 0]
        self.assert_numpy_array_equal(actual, expected)

    def test_contains(self):
        self.assertNotIn(0, self.index)
        self.assertIn(0.5, self.index)
        self.assertIn(2, self.index)

        self.assertIn(Interval(0, 1), self.index)
        self.assertIn(Interval(0, 2), self.index)
        self.assertIn(Interval(0, 0.5), self.index)
        self.assertNotIn(Interval(3, 5), self.index)
        self.assertNotIn(Interval(-1, 0, closed='left'), self.index)

    def test_non_contiguous(self):
        index = IntervalIndex.from_tuples([(0, 1), (2, 3)])
        target = [0.5, 1.5, 2.5]
        actual = index.get_indexer(target)
        expected = [0, -1, 1]
        self.assert_numpy_array_equal(actual, expected)

        self.assertNotIn(1.5, index)

    def test_union(self):
        other = IntervalIndex([2], [3])
        expected = IntervalIndex(range(3), range(1, 4))
        actual = self.index.union(other)
        self.assertTrue(expected.equals(actual))

        actual = other.union(self.index)
        self.assertTrue(expected.equals(actual))

        self.assert_numpy_array_equal(self.index.union(self.index), self.index)
        self.assert_numpy_array_equal(self.index.union(self.index[:1]),
                                      self.index)

    def test_intersection(self):
        other = IntervalIndex.from_breaks([1, 2, 3])
        expected = IntervalIndex.from_breaks([1, 2])
        actual = self.index.intersection(other)
        self.assertTrue(expected.equals(actual))

        self.assert_numpy_array_equal(self.index.intersection(self.index),
                                      self.index)

    def test_difference(self):
        self.assert_numpy_array_equal(self.index.difference(self.index[:1]),
                                      self.index[1:])

    def test_sym_diff(self):
        self.assert_numpy_array_equal(self.index[:1].sym_diff(self.index[1:]),
                                      self.index)

    def test_set_operation_errors(self):
        self.assertRaises(ValueError, self.index.union, self.index.left)

        other = IntervalIndex.from_breaks([0, 1, 2], closed='neither')
        self.assertRaises(ValueError, self.index.union, other)

    def test_isin(self):
        actual = self.index.isin(self.index)
        self.assert_numpy_array_equal([True, True], actual)

        actual = self.index.isin(self.index[:1])
        self.assert_numpy_array_equal([True, False], actual)

    def test_comparison(self):
        actual = Interval(0, 1) < self.index
        expected = [False, True]
        self.assert_numpy_array_equal(actual, expected)

        actual = Interval(0.5, 1.5) < self.index
        expected = [False, True]
        self.assert_numpy_array_equal(actual, expected)
        actual = self.index > Interval(0.5, 1.5)
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index == self.index
        expected = [True, True]
        self.assert_numpy_array_equal(actual, expected)
        actual = self.index <= self.index
        self.assert_numpy_array_equal(actual, expected)
        actual = self.index >= self.index
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index < self.index
        expected = [False, False]
        self.assert_numpy_array_equal(actual, expected)
        actual = self.index > self.index
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index == IntervalIndex.from_breaks([0, 1, 2], 'left')
        self.assert_numpy_array_equal(actual, expected)

        actual = self.index == self.index.values
        self.assert_numpy_array_equal(actual, [True, True])
        actual = self.index.values == self.index
        self.assert_numpy_array_equal(actual, [True, True])
        actual = self.index <= self.index.values
        self.assert_numpy_array_equal(actual, [True, True])
        actual = self.index != self.index.values
        self.assert_numpy_array_equal(actual, [False, False])
        actual = self.index > self.index.values
        self.assert_numpy_array_equal(actual, [False, False])
        actual = self.index.values > self.index
        self.assert_numpy_array_equal(actual, [False, False])

        # invalid comparisons
        actual = self.index == 0
        self.assert_numpy_array_equal(actual, [False, False])
        actual = self.index == self.index.left
        self.assert_numpy_array_equal(actual, [False, False])

        with self.assertRaisesRegexp(TypeError, 'unorderable types'):
            self.index > 0
        with self.assertRaisesRegexp(TypeError, 'unorderable types'):
            self.index <= 0
        with self.assertRaises(TypeError):
            self.index > np.arange(2)
        with self.assertRaises(ValueError):
            self.index > np.arange(3)

    def test_missing_values(self):
        idx = pd.Index([np.nan, pd.Interval(0, 1), pd.Interval(1, 2)])
        idx2 = pd.IntervalIndex([np.nan, 0, 1], [np.nan, 1, 2])
        assert idx.equals(idx2)

        with tm.assertRaisesRegexp(ValueError, 'both left and right sides'):
            pd.IntervalIndex([np.nan, 0, 1], [0, 1, 2])

        self.assert_numpy_array_equal(pd.isnull(idx), [True, False, False])

    def test_order(self):
        expected = IntervalIndex.from_breaks([1, 2, 3, 4])
        actual = IntervalIndex.from_tuples([(3, 4), (1, 2), (2, 3)]).order()
        self.assert_numpy_array_equal(expected, actual)

    def test_datetime(self):
        dates = pd.date_range('2000', periods=3)
        idx = IntervalIndex.from_breaks(dates)

        self.assert_numpy_array_equal(idx.left, dates[:2])
        self.assert_numpy_array_equal(idx.right, dates[-2:])

        expected = pd.date_range('2000-01-01T12:00', periods=2)
        self.assert_numpy_array_equal(idx.mid, expected)

        self.assertIn('2000-01-01T12', idx)

        target = pd.date_range('1999-12-31T12:00', periods=7, freq='12H')
        actual = idx.get_indexer(target)
        expected = [-1, -1, 0, 0, 1, 1, -1]
        self.assert_numpy_array_equal(actual, expected)

    # def test_math(self):
    #     # add, subtract, multiply, divide with scalers should be OK
    #     actual = 2 * self.index + 1
    #     expected = IntervalIndex.from_breaks((2 * np.arange(3) + 1))
    #     self.assertTrue(expected.equals(actual))

    #     actual = self.index / 2.0 - 1
    #     expected = IntervalIndex.from_breaks((np.arange(3) / 2.0 - 1))
    #     self.assertTrue(expected.equals(actual))

    #     with self.assertRaises(TypeError):
    #         # doesn't make sense to add two IntervalIndex objects
    #         self.index + self.index

    # def test_datetime_math(self):

    #     expected = IntervalIndex(pd.date_range('2000-01-02', periods=3))
    #     actual = idx + pd.to_timedelta(1, unit='D')
    #     self.assertTrue(expected.equals(actual))

    # TODO: other set operations (left join, right join, intersection),
    # set operations with conflicting IntervalIndex objects or other dtypes,
    # groupby, cut, reset_index...
