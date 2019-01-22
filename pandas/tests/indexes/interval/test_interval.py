from __future__ import division

from itertools import permutations
import re

import numpy as np
import pytest

from pandas.compat import lzip

import pandas as pd
from pandas import (
    Index, Interval, IntervalIndex, Timedelta, Timestamp, date_range,
    interval_range, isna, notna, timedelta_range)
import pandas.core.common as com
from pandas.tests.indexes.common import Base
import pandas.util.testing as tm


@pytest.fixture(scope='class', params=[None, 'foo'])
def name(request):
    return request.param


class TestIntervalIndex(Base):
    _holder = IntervalIndex

    def setup_method(self, method):
        self.index = IntervalIndex.from_arrays([0, 1], [1, 2])
        self.index_with_nan = IntervalIndex.from_tuples(
            [(0, 1), np.nan, (1, 2)])
        self.indices = dict(intervalIndex=tm.makeIntervalIndex(10))

    def create_index(self, closed='right'):
        return IntervalIndex.from_breaks(range(11), closed=closed)

    def create_index_with_nan(self, closed='right'):
        mask = [True, False] + [True] * 8
        return IntervalIndex.from_arrays(
            np.where(mask, np.arange(10), np.nan),
            np.where(mask, np.arange(1, 11), np.nan), closed=closed)

    def test_properties(self, closed):
        index = self.create_index(closed=closed)
        assert len(index) == 10
        assert index.size == 10
        assert index.shape == (10, )

        tm.assert_index_equal(index.left, Index(np.arange(10)))
        tm.assert_index_equal(index.right, Index(np.arange(1, 11)))
        tm.assert_index_equal(index.mid, Index(np.arange(0.5, 10.5)))

        assert index.closed == closed

        ivs = [Interval(l, r, closed) for l, r in zip(range(10), range(1, 11))]
        expected = np.array(ivs, dtype=object)
        tm.assert_numpy_array_equal(np.asarray(index), expected)

        # with nans
        index = self.create_index_with_nan(closed=closed)
        assert len(index) == 10
        assert index.size == 10
        assert index.shape == (10, )

        expected_left = Index([0, np.nan, 2, 3, 4, 5, 6, 7, 8, 9])
        expected_right = expected_left + 1
        expected_mid = expected_left + 0.5
        tm.assert_index_equal(index.left, expected_left)
        tm.assert_index_equal(index.right, expected_right)
        tm.assert_index_equal(index.mid, expected_mid)

        assert index.closed == closed

        ivs = [Interval(l, r, closed) if notna(l) else np.nan
               for l, r in zip(expected_left, expected_right)]
        expected = np.array(ivs, dtype=object)
        tm.assert_numpy_array_equal(np.asarray(index), expected)

    @pytest.mark.parametrize('breaks', [
        [1, 1, 2, 5, 15, 53, 217, 1014, 5335, 31240, 201608],
        [-np.inf, -100, -10, 0.5, 1, 1.5, 3.8, 101, 202, np.inf],
        pd.to_datetime(['20170101', '20170202', '20170303', '20170404']),
        pd.to_timedelta(['1ns', '2ms', '3s', '4M', '5H', '6D'])])
    def test_length(self, closed, breaks):
        # GH 18789
        index = IntervalIndex.from_breaks(breaks, closed=closed)
        result = index.length
        expected = Index(iv.length for iv in index)
        tm.assert_index_equal(result, expected)

        # with NA
        index = index.insert(1, np.nan)
        result = index.length
        expected = Index(iv.length if notna(iv) else iv for iv in index)
        tm.assert_index_equal(result, expected)

    def test_with_nans(self, closed):
        index = self.create_index(closed=closed)
        assert index.hasnans is False

        result = index.isna()
        expected = np.repeat(False, len(index))
        tm.assert_numpy_array_equal(result, expected)

        result = index.notna()
        expected = np.repeat(True, len(index))
        tm.assert_numpy_array_equal(result, expected)

        index = self.create_index_with_nan(closed=closed)
        assert index.hasnans is True

        result = index.isna()
        expected = np.array([False, True] + [False] * (len(index) - 2))
        tm.assert_numpy_array_equal(result, expected)

        result = index.notna()
        expected = np.array([True, False] + [True] * (len(index) - 2))
        tm.assert_numpy_array_equal(result, expected)

    def test_copy(self, closed):
        expected = self.create_index(closed=closed)

        result = expected.copy()
        assert result.equals(expected)

        result = expected.copy(deep=True)
        assert result.equals(expected)
        assert result.left is not expected.left

    def test_ensure_copied_data(self, closed):
        # exercise the copy flag in the constructor

        # not copying
        index = self.create_index(closed=closed)
        result = IntervalIndex(index, copy=False)
        tm.assert_numpy_array_equal(index.left.values, result.left.values,
                                    check_same='same')
        tm.assert_numpy_array_equal(index.right.values, result.right.values,
                                    check_same='same')

        # by-definition make a copy
        result = IntervalIndex(index._ndarray_values, copy=False)
        tm.assert_numpy_array_equal(index.left.values, result.left.values,
                                    check_same='copy')
        tm.assert_numpy_array_equal(index.right.values, result.right.values,
                                    check_same='copy')

    def test_equals(self, closed):
        expected = IntervalIndex.from_breaks(np.arange(5), closed=closed)
        assert expected.equals(expected)
        assert expected.equals(expected.copy())

        assert not expected.equals(expected.astype(object))
        assert not expected.equals(np.array(expected))
        assert not expected.equals(list(expected))

        assert not expected.equals([1, 2])
        assert not expected.equals(np.array([1, 2]))
        assert not expected.equals(pd.date_range('20130101', periods=2))

        expected_name1 = IntervalIndex.from_breaks(
            np.arange(5), closed=closed, name='foo')
        expected_name2 = IntervalIndex.from_breaks(
            np.arange(5), closed=closed, name='bar')
        assert expected.equals(expected_name1)
        assert expected_name1.equals(expected_name2)

        for other_closed in {'left', 'right', 'both', 'neither'} - {closed}:
            expected_other_closed = IntervalIndex.from_breaks(
                np.arange(5), closed=other_closed)
            assert not expected.equals(expected_other_closed)

    @pytest.mark.parametrize('klass', [list, tuple, np.array, pd.Series])
    def test_where(self, closed, klass):
        idx = self.create_index(closed=closed)
        cond = [True] * len(idx)
        expected = idx
        result = expected.where(klass(cond))
        tm.assert_index_equal(result, expected)

        cond = [False] + [True] * len(idx[1:])
        expected = IntervalIndex([np.nan] + idx[1:].tolist())
        result = idx.where(klass(cond))
        tm.assert_index_equal(result, expected)

    def test_delete(self, closed):
        expected = IntervalIndex.from_breaks(np.arange(1, 11), closed=closed)
        result = self.create_index(closed=closed).delete(0)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('data', [
        interval_range(0, periods=10, closed='neither'),
        interval_range(1.7, periods=8, freq=2.5, closed='both'),
        interval_range(Timestamp('20170101'), periods=12, closed='left'),
        interval_range(Timedelta('1 day'), periods=6, closed='right')])
    def test_insert(self, data):
        item = data[0]
        idx_item = IntervalIndex([item])

        # start
        expected = idx_item.append(data)
        result = data.insert(0, item)
        tm.assert_index_equal(result, expected)

        # end
        expected = data.append(idx_item)
        result = data.insert(len(data), item)
        tm.assert_index_equal(result, expected)

        # mid
        expected = data[:3].append(idx_item).append(data[3:])
        result = data.insert(3, item)
        tm.assert_index_equal(result, expected)

        # invalid type
        msg = 'can only insert Interval objects and NA into an IntervalIndex'
        with pytest.raises(ValueError, match=msg):
            data.insert(1, 'foo')

        # invalid closed
        msg = 'inserted item must be closed on the same side as the index'
        for closed in {'left', 'right', 'both', 'neither'} - {item.closed}:
            with pytest.raises(ValueError, match=msg):
                bad_item = Interval(item.left, item.right, closed=closed)
                data.insert(1, bad_item)

        # GH 18295 (test missing)
        na_idx = IntervalIndex([np.nan], closed=data.closed)
        for na in (np.nan, pd.NaT, None):
            expected = data[:1].append(na_idx).append(data[1:])
            result = data.insert(1, na)
            tm.assert_index_equal(result, expected)

    def test_take(self, closed):
        index = self.create_index(closed=closed)

        result = index.take(range(10))
        tm.assert_index_equal(result, index)

        result = index.take([0, 0, 1])
        expected = IntervalIndex.from_arrays(
            [0, 0, 1], [1, 1, 2], closed=closed)
        tm.assert_index_equal(result, expected)

    def test_unique(self, closed):
        # unique non-overlapping
        idx = IntervalIndex.from_tuples(
            [(0, 1), (2, 3), (4, 5)], closed=closed)
        assert idx.is_unique is True

        # unique overlapping - distinct endpoints
        idx = IntervalIndex.from_tuples([(0, 1), (0.5, 1.5)], closed=closed)
        assert idx.is_unique is True

        # unique overlapping - shared endpoints
        idx = pd.IntervalIndex.from_tuples(
            [(1, 2), (1, 3), (2, 3)], closed=closed)
        assert idx.is_unique is True

        # unique nested
        idx = IntervalIndex.from_tuples([(-1, 1), (-2, 2)], closed=closed)
        assert idx.is_unique is True

        # duplicate
        idx = IntervalIndex.from_tuples(
            [(0, 1), (0, 1), (2, 3)], closed=closed)
        assert idx.is_unique is False

        # empty
        idx = IntervalIndex([], closed=closed)
        assert idx.is_unique is True

    def test_monotonic(self, closed):
        # increasing non-overlapping
        idx = IntervalIndex.from_tuples(
            [(0, 1), (2, 3), (4, 5)], closed=closed)
        assert idx.is_monotonic is True
        assert idx._is_strictly_monotonic_increasing is True
        assert idx.is_monotonic_decreasing is False
        assert idx._is_strictly_monotonic_decreasing is False

        # decreasing non-overlapping
        idx = IntervalIndex.from_tuples(
            [(4, 5), (2, 3), (1, 2)], closed=closed)
        assert idx.is_monotonic is False
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is True
        assert idx._is_strictly_monotonic_decreasing is True

        # unordered non-overlapping
        idx = IntervalIndex.from_tuples(
            [(0, 1), (4, 5), (2, 3)], closed=closed)
        assert idx.is_monotonic is False
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is False
        assert idx._is_strictly_monotonic_decreasing is False

        # increasing overlapping
        idx = IntervalIndex.from_tuples(
            [(0, 2), (0.5, 2.5), (1, 3)], closed=closed)
        assert idx.is_monotonic is True
        assert idx._is_strictly_monotonic_increasing is True
        assert idx.is_monotonic_decreasing is False
        assert idx._is_strictly_monotonic_decreasing is False

        # decreasing overlapping
        idx = IntervalIndex.from_tuples(
            [(1, 3), (0.5, 2.5), (0, 2)], closed=closed)
        assert idx.is_monotonic is False
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is True
        assert idx._is_strictly_monotonic_decreasing is True

        # unordered overlapping
        idx = IntervalIndex.from_tuples(
            [(0.5, 2.5), (0, 2), (1, 3)], closed=closed)
        assert idx.is_monotonic is False
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is False
        assert idx._is_strictly_monotonic_decreasing is False

        # increasing overlapping shared endpoints
        idx = pd.IntervalIndex.from_tuples(
            [(1, 2), (1, 3), (2, 3)], closed=closed)
        assert idx.is_monotonic is True
        assert idx._is_strictly_monotonic_increasing is True
        assert idx.is_monotonic_decreasing is False
        assert idx._is_strictly_monotonic_decreasing is False

        # decreasing overlapping shared endpoints
        idx = pd.IntervalIndex.from_tuples(
            [(2, 3), (1, 3), (1, 2)], closed=closed)
        assert idx.is_monotonic is False
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is True
        assert idx._is_strictly_monotonic_decreasing is True

        # stationary
        idx = IntervalIndex.from_tuples([(0, 1), (0, 1)], closed=closed)
        assert idx.is_monotonic is True
        assert idx._is_strictly_monotonic_increasing is False
        assert idx.is_monotonic_decreasing is True
        assert idx._is_strictly_monotonic_decreasing is False

        # empty
        idx = IntervalIndex([], closed=closed)
        assert idx.is_monotonic is True
        assert idx._is_strictly_monotonic_increasing is True
        assert idx.is_monotonic_decreasing is True
        assert idx._is_strictly_monotonic_decreasing is True

    @pytest.mark.skip(reason='not a valid repr as we use interval notation')
    def test_repr(self):
        i = IntervalIndex.from_tuples([(0, 1), (1, 2)], closed='right')
        expected = ("IntervalIndex(left=[0, 1],"
                    "\n              right=[1, 2],"
                    "\n              closed='right',"
                    "\n              dtype='interval[int64]')")
        assert repr(i) == expected

        i = IntervalIndex.from_tuples((Timestamp('20130101'),
                                       Timestamp('20130102')),
                                      (Timestamp('20130102'),
                                       Timestamp('20130103')),
                                      closed='right')
        expected = ("IntervalIndex(left=['2013-01-01', '2013-01-02'],"
                    "\n              right=['2013-01-02', '2013-01-03'],"
                    "\n              closed='right',"
                    "\n              dtype='interval[datetime64[ns]]')")
        assert repr(i) == expected

    @pytest.mark.skip(reason='not a valid repr as we use interval notation')
    def test_repr_max_seq_item_setting(self):
        super(TestIntervalIndex, self).test_repr_max_seq_item_setting()

    @pytest.mark.skip(reason='not a valid repr as we use interval notation')
    def test_repr_roundtrip(self):
        super(TestIntervalIndex, self).test_repr_roundtrip()

    def test_frame_repr(self):
        # https://github.com/pandas-dev/pandas/pull/24134/files
        df = pd.DataFrame({'A': [1, 2, 3, 4]},
                          index=pd.IntervalIndex.from_breaks([0, 1, 2, 3, 4]))
        result = repr(df)
        expected = (
            '        A\n'
            '(0, 1]  1\n'
            '(1, 2]  2\n'
            '(2, 3]  3\n'
            '(3, 4]  4'
        )
        assert result == expected

        # TODO: check this behavior is consistent with test_interval_new.py
    def test_get_item(self, closed):
        i = IntervalIndex.from_arrays((0, 1, np.nan), (1, 2, np.nan),
                                      closed=closed)
        assert i[0] == Interval(0.0, 1.0, closed=closed)
        assert i[1] == Interval(1.0, 2.0, closed=closed)
        assert isna(i[2])

        result = i[0:1]
        expected = IntervalIndex.from_arrays((0.,), (1.,), closed=closed)
        tm.assert_index_equal(result, expected)

        result = i[0:2]
        expected = IntervalIndex.from_arrays((0., 1), (1., 2.), closed=closed)
        tm.assert_index_equal(result, expected)

        result = i[1:3]
        expected = IntervalIndex.from_arrays((1., np.nan), (2., np.nan),
                                             closed=closed)
        tm.assert_index_equal(result, expected)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_get_loc_value(self):
        pytest.raises(KeyError, self.index.get_loc, 0)
        assert self.index.get_loc(0.5) == 0
        assert self.index.get_loc(1) == 0
        assert self.index.get_loc(1.5) == 1
        assert self.index.get_loc(2) == 1
        pytest.raises(KeyError, self.index.get_loc, -1)
        pytest.raises(KeyError, self.index.get_loc, 3)

        idx = IntervalIndex.from_tuples([(0, 2), (1, 3)])
        assert idx.get_loc(0.5) == 0
        assert idx.get_loc(1) == 0
        tm.assert_numpy_array_equal(idx.get_loc(1.5),
                                    np.array([0, 1], dtype='intp'))
        tm.assert_numpy_array_equal(np.sort(idx.get_loc(2)),
                                    np.array([0, 1], dtype='intp'))
        assert idx.get_loc(3) == 1
        pytest.raises(KeyError, idx.get_loc, 3.5)

        idx = IntervalIndex.from_arrays([0, 2], [1, 3])
        pytest.raises(KeyError, idx.get_loc, 1.5)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def slice_locs_cases(self, breaks):
        # TODO: same tests for more index types
        index = IntervalIndex.from_breaks([0, 1, 2], closed='right')
        assert index.slice_locs() == (0, 2)
        assert index.slice_locs(0, 1) == (0, 1)
        assert index.slice_locs(1, 1) == (0, 1)
        assert index.slice_locs(0, 2) == (0, 2)
        assert index.slice_locs(0.5, 1.5) == (0, 2)
        assert index.slice_locs(0, 0.5) == (0, 1)
        assert index.slice_locs(start=1) == (0, 2)
        assert index.slice_locs(start=1.2) == (1, 2)
        assert index.slice_locs(end=1) == (0, 1)
        assert index.slice_locs(end=1.1) == (0, 2)
        assert index.slice_locs(end=1.0) == (0, 1)
        assert index.slice_locs(-1, -1) == (0, 0)

        index = IntervalIndex.from_breaks([0, 1, 2], closed='neither')
        assert index.slice_locs(0, 1) == (0, 1)
        assert index.slice_locs(0, 2) == (0, 2)
        assert index.slice_locs(0.5, 1.5) == (0, 2)
        assert index.slice_locs(1, 1) == (1, 1)
        assert index.slice_locs(1, 2) == (1, 2)

        index = IntervalIndex.from_tuples([(0, 1), (2, 3), (4, 5)],
                                          closed='both')
        assert index.slice_locs(1, 1) == (0, 1)
        assert index.slice_locs(1, 2) == (0, 2)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_slice_locs_int64(self):
        self.slice_locs_cases([0, 1, 2])

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_slice_locs_float64(self):
        self.slice_locs_cases([0.0, 1.0, 2.0])

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def slice_locs_decreasing_cases(self, tuples):
        index = IntervalIndex.from_tuples(tuples)
        assert index.slice_locs(1.5, 0.5) == (1, 3)
        assert index.slice_locs(2, 0) == (1, 3)
        assert index.slice_locs(2, 1) == (1, 3)
        assert index.slice_locs(3, 1.1) == (0, 3)
        assert index.slice_locs(3, 3) == (0, 2)
        assert index.slice_locs(3.5, 3.3) == (0, 1)
        assert index.slice_locs(1, -3) == (2, 3)

        slice_locs = index.slice_locs(-1, -1)
        assert slice_locs[0] == slice_locs[1]

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_slice_locs_decreasing_int64(self):
        self.slice_locs_cases([(2, 4), (1, 3), (0, 2)])

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_slice_locs_decreasing_float64(self):
        self.slice_locs_cases([(2., 4.), (1., 3.), (0., 2.)])

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_slice_locs_fails(self):
        index = IntervalIndex.from_tuples([(1, 2), (0, 1), (2, 3)])
        with pytest.raises(KeyError):
            index.slice_locs(1, 2)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_get_loc_interval(self):
        assert self.index.get_loc(Interval(0, 1)) == 0
        assert self.index.get_loc(Interval(0, 0.5)) == 0
        assert self.index.get_loc(Interval(0, 1, 'left')) == 0
        pytest.raises(KeyError, self.index.get_loc, Interval(2, 3))
        pytest.raises(KeyError, self.index.get_loc,
                      Interval(-1, 0, 'left'))

    # Make consistent with test_interval_new.py (see #16316, #16386)
    @pytest.mark.parametrize('item', [3, Interval(1, 4)])
    def test_get_loc_length_one(self, item, closed):
        # GH 20921
        index = IntervalIndex.from_tuples([(0, 5)], closed=closed)
        result = index.get_loc(item)
        assert result == 0

    # Make consistent with test_interval_new.py (see #16316, #16386)
    @pytest.mark.parametrize('breaks', [
        date_range('20180101', periods=4),
        date_range('20180101', periods=4, tz='US/Eastern'),
        timedelta_range('0 days', periods=4)], ids=lambda x: str(x.dtype))
    def test_get_loc_datetimelike_nonoverlapping(self, breaks):
        # GH 20636
        # nonoverlapping = IntervalIndex method and no i8 conversion
        index = IntervalIndex.from_breaks(breaks)

        value = index[0].mid
        result = index.get_loc(value)
        expected = 0
        assert result == expected

        interval = Interval(index[0].left, index[1].right)
        result = index.get_loc(interval)
        expected = slice(0, 2)
        assert result == expected

    # Make consistent with test_interval_new.py (see #16316, #16386)
    @pytest.mark.parametrize('arrays', [
        (date_range('20180101', periods=4), date_range('20180103', periods=4)),
        (date_range('20180101', periods=4, tz='US/Eastern'),
         date_range('20180103', periods=4, tz='US/Eastern')),
        (timedelta_range('0 days', periods=4),
         timedelta_range('2 days', periods=4))], ids=lambda x: str(x[0].dtype))
    def test_get_loc_datetimelike_overlapping(self, arrays):
        # GH 20636
        # overlapping = IntervalTree method with i8 conversion
        index = IntervalIndex.from_arrays(*arrays)

        value = index[0].mid + Timedelta('12 hours')
        result = np.sort(index.get_loc(value))
        expected = np.array([0, 1], dtype='intp')
        assert tm.assert_numpy_array_equal(result, expected)

        interval = Interval(index[0].left, index[1].right)
        result = np.sort(index.get_loc(interval))
        expected = np.array([0, 1, 2], dtype='intp')
        assert tm.assert_numpy_array_equal(result, expected)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_get_indexer(self):
        actual = self.index.get_indexer([-1, 0, 0.5, 1, 1.5, 2, 3])
        expected = np.array([-1, -1, 0, 0, 1, 1, -1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(self.index)
        expected = np.array([0, 1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        index = IntervalIndex.from_breaks([0, 1, 2], closed='left')
        actual = index.get_indexer([-1, 0, 0.5, 1, 1.5, 2, 3])
        expected = np.array([-1, 0, 0, 1, 1, -1, -1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(index[:1])
        expected = np.array([0], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(index)
        expected = np.array([-1, 1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_get_indexer_subintervals(self):

        # TODO: is this right?
        # return indexers for wholly contained subintervals
        target = IntervalIndex.from_breaks(np.linspace(0, 2, 5))
        actual = self.index.get_indexer(target)
        expected = np.array([0, 0, 1, 1], dtype='p')
        tm.assert_numpy_array_equal(actual, expected)

        target = IntervalIndex.from_breaks([0, 0.67, 1.33, 2])
        actual = self.index.get_indexer(target)
        expected = np.array([0, 0, 1, 1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index.get_indexer(target[[0, -1]])
        expected = np.array([0, 1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        target = IntervalIndex.from_breaks([0, 0.33, 0.67, 1], closed='left')
        actual = self.index.get_indexer(target)
        expected = np.array([0, 0, 0], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

    # Make consistent with test_interval_new.py (see #16316, #16386)
    @pytest.mark.parametrize('item', [
        [3], np.arange(1, 5), [Interval(1, 4)], interval_range(1, 4)])
    def test_get_indexer_length_one(self, item, closed):
        # GH 17284
        index = IntervalIndex.from_tuples([(0, 5)], closed=closed)
        result = index.get_indexer(item)
        expected = np.array([0] * len(item), dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    # Make consistent with test_interval_new.py (see #16316, #16386)
    @pytest.mark.parametrize('arrays', [
        (date_range('20180101', periods=4), date_range('20180103', periods=4)),
        (date_range('20180101', periods=4, tz='US/Eastern'),
         date_range('20180103', periods=4, tz='US/Eastern')),
        (timedelta_range('0 days', periods=4),
         timedelta_range('2 days', periods=4))], ids=lambda x: str(x[0].dtype))
    def test_get_reindexer_datetimelike(self, arrays):
        # GH 20636
        index = IntervalIndex.from_arrays(*arrays)
        tuples = [(index[0].left, index[0].left + pd.Timedelta('12H')),
                  (index[-1].right - pd.Timedelta('12H'), index[-1].right)]
        target = IntervalIndex.from_tuples(tuples)

        result = index._get_reindexer(target)
        expected = np.array([0, 3], dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('breaks', [
        date_range('20180101', periods=4),
        date_range('20180101', periods=4, tz='US/Eastern'),
        timedelta_range('0 days', periods=4)], ids=lambda x: str(x.dtype))
    def test_maybe_convert_i8(self, breaks):
        # GH 20636
        index = IntervalIndex.from_breaks(breaks)

        # intervalindex
        result = index._maybe_convert_i8(index)
        expected = IntervalIndex.from_breaks(breaks.asi8)
        tm.assert_index_equal(result, expected)

        # interval
        interval = Interval(breaks[0], breaks[1])
        result = index._maybe_convert_i8(interval)
        expected = Interval(breaks[0].value, breaks[1].value)
        assert result == expected

        # datetimelike index
        result = index._maybe_convert_i8(breaks)
        expected = Index(breaks.asi8)
        tm.assert_index_equal(result, expected)

        # datetimelike scalar
        result = index._maybe_convert_i8(breaks[0])
        expected = breaks[0].value
        assert result == expected

        # list-like of datetimelike scalars
        result = index._maybe_convert_i8(list(breaks))
        expected = Index(breaks.asi8)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('breaks', [
        date_range('2018-01-01', periods=5),
        timedelta_range('0 days', periods=5)])
    def test_maybe_convert_i8_nat(self, breaks):
        # GH 20636
        index = IntervalIndex.from_breaks(breaks)

        to_convert = breaks._constructor([pd.NaT] * 3)
        expected = pd.Float64Index([np.nan] * 3)
        result = index._maybe_convert_i8(to_convert)
        tm.assert_index_equal(result, expected)

        to_convert = to_convert.insert(0, breaks[0])
        expected = expected.insert(0, float(breaks[0].value))
        result = index._maybe_convert_i8(to_convert)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('breaks', [
        np.arange(5, dtype='int64'),
        np.arange(5, dtype='float64')], ids=lambda x: str(x.dtype))
    @pytest.mark.parametrize('make_key', [
        IntervalIndex.from_breaks,
        lambda breaks: Interval(breaks[0], breaks[1]),
        lambda breaks: breaks,
        lambda breaks: breaks[0],
        list], ids=['IntervalIndex', 'Interval', 'Index', 'scalar', 'list'])
    def test_maybe_convert_i8_numeric(self, breaks, make_key):
        # GH 20636
        index = IntervalIndex.from_breaks(breaks)
        key = make_key(breaks)

        # no conversion occurs for numeric
        result = index._maybe_convert_i8(key)
        assert result is key

    @pytest.mark.parametrize('breaks1, breaks2', permutations([
        date_range('20180101', periods=4),
        date_range('20180101', periods=4, tz='US/Eastern'),
        timedelta_range('0 days', periods=4)], 2), ids=lambda x: str(x.dtype))
    @pytest.mark.parametrize('make_key', [
        IntervalIndex.from_breaks,
        lambda breaks: Interval(breaks[0], breaks[1]),
        lambda breaks: breaks,
        lambda breaks: breaks[0],
        list], ids=['IntervalIndex', 'Interval', 'Index', 'scalar', 'list'])
    def test_maybe_convert_i8_errors(self, breaks1, breaks2, make_key):
        # GH 20636
        index = IntervalIndex.from_breaks(breaks1)
        key = make_key(breaks2)

        msg = ('Cannot index an IntervalIndex of subtype {dtype1} with '
               'values of dtype {dtype2}')
        msg = re.escape(msg.format(dtype1=breaks1.dtype, dtype2=breaks2.dtype))
        with pytest.raises(ValueError, match=msg):
            index._maybe_convert_i8(key)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_contains(self):
        # Only endpoints are valid.
        i = IntervalIndex.from_arrays([0, 1], [1, 2])

        # Invalid
        assert 0 not in i
        assert 1 not in i
        assert 2 not in i

        # Valid
        assert Interval(0, 1) in i
        assert Interval(0, 2) in i
        assert Interval(0, 0.5) in i
        assert Interval(3, 5) not in i
        assert Interval(-1, 0, closed='left') not in i

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def testcontains(self):
        # can select values that are IN the range of a value
        i = IntervalIndex.from_arrays([0, 1], [1, 2])

        assert i.contains(0.1)
        assert i.contains(0.5)
        assert i.contains(1)
        assert i.contains(Interval(0, 1))
        assert i.contains(Interval(0, 2))

        # these overlaps completely
        assert i.contains(Interval(0, 3))
        assert i.contains(Interval(1, 3))

        assert not i.contains(20)
        assert not i.contains(-20)

    def test_dropna(self, closed):

        expected = IntervalIndex.from_tuples(
            [(0.0, 1.0), (1.0, 2.0)], closed=closed)

        ii = IntervalIndex.from_tuples([(0, 1), (1, 2), np.nan], closed=closed)
        result = ii.dropna()
        tm.assert_index_equal(result, expected)

        ii = IntervalIndex.from_arrays(
            [0, 1, np.nan], [1, 2, np.nan], closed=closed)
        result = ii.dropna()
        tm.assert_index_equal(result, expected)

    # TODO: check this behavior is consistent with test_interval_new.py
    def test_non_contiguous(self, closed):
        index = IntervalIndex.from_tuples([(0, 1), (2, 3)], closed=closed)
        target = [0.5, 1.5, 2.5]
        actual = index.get_indexer(target)
        expected = np.array([0, -1, 1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        assert 1.5 not in index

    @pytest.mark.parametrize("sort", [True, False])
    def test_union(self, closed, sort):
        index = self.create_index(closed=closed)
        other = IntervalIndex.from_breaks(range(5, 13), closed=closed)

        expected = IntervalIndex.from_breaks(range(13), closed=closed)
        result = index[::-1].union(other, sort=sort)
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        result = other[::-1].union(index, sort=sort)
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        tm.assert_index_equal(index.union(index, sort=sort), index)
        tm.assert_index_equal(index.union(index[:1], sort=sort), index)

        # GH 19101: empty result, same dtype
        index = IntervalIndex(np.array([], dtype='int64'), closed=closed)
        result = index.union(index, sort=sort)
        tm.assert_index_equal(result, index)

        # GH 19101: empty result, different dtypes
        other = IntervalIndex(np.array([], dtype='float64'), closed=closed)
        result = index.union(other, sort=sort)
        tm.assert_index_equal(result, index)

    @pytest.mark.parametrize("sort", [True, False])
    def test_intersection(self, closed, sort):
        index = self.create_index(closed=closed)
        other = IntervalIndex.from_breaks(range(5, 13), closed=closed)

        expected = IntervalIndex.from_breaks(range(5, 11), closed=closed)
        result = index[::-1].intersection(other, sort=sort)
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        result = other[::-1].intersection(index, sort=sort)
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        tm.assert_index_equal(index.intersection(index, sort=sort), index)

        # GH 19101: empty result, same dtype
        other = IntervalIndex.from_breaks(range(300, 314), closed=closed)
        expected = IntervalIndex(np.array([], dtype='int64'), closed=closed)
        result = index.intersection(other, sort=sort)
        tm.assert_index_equal(result, expected)

        # GH 19101: empty result, different dtypes
        breaks = np.arange(300, 314, dtype='float64')
        other = IntervalIndex.from_breaks(breaks, closed=closed)
        result = index.intersection(other, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("sort", [True, False])
    def test_difference(self, closed, sort):
        index = IntervalIndex.from_arrays([1, 0, 3, 2],
                                          [1, 2, 3, 4],
                                          closed=closed)
        result = index.difference(index[:1], sort=sort)
        expected = index[1:]
        if sort:
            expected = expected.sort_values()
        tm.assert_index_equal(result, expected)

        # GH 19101: empty result, same dtype
        result = index.difference(index, sort=sort)
        expected = IntervalIndex(np.array([], dtype='int64'), closed=closed)
        tm.assert_index_equal(result, expected)

        # GH 19101: empty result, different dtypes
        other = IntervalIndex.from_arrays(index.left.astype('float64'),
                                          index.right, closed=closed)
        result = index.difference(other, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("sort", [True, False])
    def test_symmetric_difference(self, closed, sort):
        index = self.create_index(closed=closed)
        result = index[1:].symmetric_difference(index[:-1], sort=sort)
        expected = IntervalIndex([index[0], index[-1]])
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        # GH 19101: empty result, same dtype
        result = index.symmetric_difference(index, sort=sort)
        expected = IntervalIndex(np.array([], dtype='int64'), closed=closed)
        if sort:
            tm.assert_index_equal(result, expected)
        assert tm.equalContents(result, expected)

        # GH 19101: empty result, different dtypes
        other = IntervalIndex.from_arrays(index.left.astype('float64'),
                                          index.right, closed=closed)
        result = index.symmetric_difference(other, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('op_name', [
        'union', 'intersection', 'difference', 'symmetric_difference'])
    @pytest.mark.parametrize("sort", [True, False])
    def test_set_operation_errors(self, closed, op_name, sort):
        index = self.create_index(closed=closed)
        set_op = getattr(index, op_name)

        # non-IntervalIndex
        msg = ('the other index needs to be an IntervalIndex too, but '
               'was type Int64Index')
        with pytest.raises(TypeError, match=msg):
            set_op(Index([1, 2, 3]), sort=sort)

        # mixed closed
        msg = ('can only do set operations between two IntervalIndex objects '
               'that are closed on the same side')
        for other_closed in {'right', 'left', 'both', 'neither'} - {closed}:
            other = self.create_index(closed=other_closed)
            with pytest.raises(ValueError, match=msg):
                set_op(other, sort=sort)

        # GH 19016: incompatible dtypes
        other = interval_range(Timestamp('20180101'), periods=9, closed=closed)
        msg = ('can only do {op} between two IntervalIndex objects that have '
               'compatible dtypes').format(op=op_name)
        with pytest.raises(TypeError, match=msg):
            set_op(other, sort=sort)

    def test_isin(self, closed):
        index = self.create_index(closed=closed)

        expected = np.array([True] + [False] * (len(index) - 1))
        result = index.isin(index[:1])
        tm.assert_numpy_array_equal(result, expected)

        result = index.isin([index[0]])
        tm.assert_numpy_array_equal(result, expected)

        other = IntervalIndex.from_breaks(np.arange(-2, 10), closed=closed)
        expected = np.array([True] * (len(index) - 1) + [False])
        result = index.isin(other)
        tm.assert_numpy_array_equal(result, expected)

        result = index.isin(other.tolist())
        tm.assert_numpy_array_equal(result, expected)

        for other_closed in {'right', 'left', 'both', 'neither'}:
            other = self.create_index(closed=other_closed)
            expected = np.repeat(closed == other_closed, len(index))
            result = index.isin(other)
            tm.assert_numpy_array_equal(result, expected)

            result = index.isin(other.tolist())
            tm.assert_numpy_array_equal(result, expected)

    def test_comparison(self):
        actual = Interval(0, 1) < self.index
        expected = np.array([False, True])
        tm.assert_numpy_array_equal(actual, expected)

        actual = Interval(0.5, 1.5) < self.index
        expected = np.array([False, True])
        tm.assert_numpy_array_equal(actual, expected)
        actual = self.index > Interval(0.5, 1.5)
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index == self.index
        expected = np.array([True, True])
        tm.assert_numpy_array_equal(actual, expected)
        actual = self.index <= self.index
        tm.assert_numpy_array_equal(actual, expected)
        actual = self.index >= self.index
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index < self.index
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(actual, expected)
        actual = self.index > self.index
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index == IntervalIndex.from_breaks([0, 1, 2], 'left')
        tm.assert_numpy_array_equal(actual, expected)

        actual = self.index == self.index.values
        tm.assert_numpy_array_equal(actual, np.array([True, True]))
        actual = self.index.values == self.index
        tm.assert_numpy_array_equal(actual, np.array([True, True]))
        actual = self.index <= self.index.values
        tm.assert_numpy_array_equal(actual, np.array([True, True]))
        actual = self.index != self.index.values
        tm.assert_numpy_array_equal(actual, np.array([False, False]))
        actual = self.index > self.index.values
        tm.assert_numpy_array_equal(actual, np.array([False, False]))
        actual = self.index.values > self.index
        tm.assert_numpy_array_equal(actual, np.array([False, False]))

        # invalid comparisons
        actual = self.index == 0
        tm.assert_numpy_array_equal(actual, np.array([False, False]))
        actual = self.index == self.index.left
        tm.assert_numpy_array_equal(actual, np.array([False, False]))

        with pytest.raises(TypeError, match='unorderable types'):
            self.index > 0
        with pytest.raises(TypeError, match='unorderable types'):
            self.index <= 0
        with pytest.raises(TypeError):
            self.index > np.arange(2)
        with pytest.raises(ValueError):
            self.index > np.arange(3)

    def test_missing_values(self, closed):
        idx = Index([np.nan, Interval(0, 1, closed=closed),
                     Interval(1, 2, closed=closed)])
        idx2 = IntervalIndex.from_arrays(
            [np.nan, 0, 1], [np.nan, 1, 2], closed=closed)
        assert idx.equals(idx2)

        with pytest.raises(ValueError):
            IntervalIndex.from_arrays(
                [np.nan, 0, 1], np.array([0, 1, 2]), closed=closed)

        tm.assert_numpy_array_equal(isna(idx),
                                    np.array([True, False, False]))

    def test_sort_values(self, closed):
        index = self.create_index(closed=closed)

        result = index.sort_values()
        tm.assert_index_equal(result, index)

        result = index.sort_values(ascending=False)
        tm.assert_index_equal(result, index[::-1])

        # with nan
        index = IntervalIndex([Interval(1, 2), np.nan, Interval(0, 1)])

        result = index.sort_values()
        expected = IntervalIndex([Interval(0, 1), Interval(1, 2), np.nan])
        tm.assert_index_equal(result, expected)

        result = index.sort_values(ascending=False)
        expected = IntervalIndex([np.nan, Interval(1, 2), Interval(0, 1)])
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('tz', [None, 'US/Eastern'])
    def test_datetime(self, tz):
        start = Timestamp('2000-01-01', tz=tz)
        dates = date_range(start=start, periods=10)
        index = IntervalIndex.from_breaks(dates)

        # test mid
        start = Timestamp('2000-01-01T12:00', tz=tz)
        expected = date_range(start=start, periods=9)
        tm.assert_index_equal(index.mid, expected)

        # __contains__ doesn't check individual points
        assert Timestamp('2000-01-01', tz=tz) not in index
        assert Timestamp('2000-01-01T12', tz=tz) not in index
        assert Timestamp('2000-01-02', tz=tz) not in index
        iv_true = Interval(Timestamp('2000-01-01T08', tz=tz),
                           Timestamp('2000-01-01T18', tz=tz))
        iv_false = Interval(Timestamp('1999-12-31', tz=tz),
                            Timestamp('2000-01-01', tz=tz))
        assert iv_true in index
        assert iv_false not in index

        # .contains does check individual points
        assert not index.contains(Timestamp('2000-01-01', tz=tz))
        assert index.contains(Timestamp('2000-01-01T12', tz=tz))
        assert index.contains(Timestamp('2000-01-02', tz=tz))
        assert index.contains(iv_true)
        assert not index.contains(iv_false)

        # test get_indexer
        start = Timestamp('1999-12-31T12:00', tz=tz)
        target = date_range(start=start, periods=7, freq='12H')
        actual = index.get_indexer(target)
        expected = np.array([-1, -1, 0, 0, 1, 1, 2], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

        start = Timestamp('2000-01-08T18:00', tz=tz)
        target = date_range(start=start, periods=7, freq='6H')
        actual = index.get_indexer(target)
        expected = np.array([7, 7, 8, 8, 8, 8, -1], dtype='intp')
        tm.assert_numpy_array_equal(actual, expected)

    def test_append(self, closed):

        index1 = IntervalIndex.from_arrays([0, 1], [1, 2], closed=closed)
        index2 = IntervalIndex.from_arrays([1, 2], [2, 3], closed=closed)

        result = index1.append(index2)
        expected = IntervalIndex.from_arrays(
            [0, 1, 1, 2], [1, 2, 2, 3], closed=closed)
        tm.assert_index_equal(result, expected)

        result = index1.append([index1, index2])
        expected = IntervalIndex.from_arrays(
            [0, 1, 0, 1, 1, 2], [1, 2, 1, 2, 2, 3], closed=closed)
        tm.assert_index_equal(result, expected)

        msg = ('can only append two IntervalIndex objects that are closed '
               'on the same side')
        for other_closed in {'left', 'right', 'both', 'neither'} - {closed}:
            index_other_closed = IntervalIndex.from_arrays(
                [0, 1], [1, 2], closed=other_closed)
            with pytest.raises(ValueError, match=msg):
                index1.append(index_other_closed)

    def test_is_non_overlapping_monotonic(self, closed):
        # Should be True in all cases
        tpls = [(0, 1), (2, 3), (4, 5), (6, 7)]
        idx = IntervalIndex.from_tuples(tpls, closed=closed)
        assert idx.is_non_overlapping_monotonic is True

        idx = IntervalIndex.from_tuples(tpls[::-1], closed=closed)
        assert idx.is_non_overlapping_monotonic is True

        # Should be False in all cases (overlapping)
        tpls = [(0, 2), (1, 3), (4, 5), (6, 7)]
        idx = IntervalIndex.from_tuples(tpls, closed=closed)
        assert idx.is_non_overlapping_monotonic is False

        idx = IntervalIndex.from_tuples(tpls[::-1], closed=closed)
        assert idx.is_non_overlapping_monotonic is False

        # Should be False in all cases (non-monotonic)
        tpls = [(0, 1), (2, 3), (6, 7), (4, 5)]
        idx = IntervalIndex.from_tuples(tpls, closed=closed)
        assert idx.is_non_overlapping_monotonic is False

        idx = IntervalIndex.from_tuples(tpls[::-1], closed=closed)
        assert idx.is_non_overlapping_monotonic is False

        # Should be False for closed='both', otherwise True (GH16560)
        if closed == 'both':
            idx = IntervalIndex.from_breaks(range(4), closed=closed)
            assert idx.is_non_overlapping_monotonic is False
        else:
            idx = IntervalIndex.from_breaks(range(4), closed=closed)
            assert idx.is_non_overlapping_monotonic is True

    @pytest.mark.parametrize('start, shift, na_value', [
        (0, 1, np.nan),
        (Timestamp('2018-01-01'), Timedelta('1 day'), pd.NaT),
        (Timedelta('0 days'), Timedelta('1 day'), pd.NaT)])
    def test_is_overlapping(self, start, shift, na_value, closed):
        # GH 23309
        # see test_interval_tree.py for extensive tests; interface tests here

        # non-overlapping
        tuples = [(start + n * shift, start + (n + 1) * shift)
                  for n in (0, 2, 4)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        assert index.is_overlapping is False

        # non-overlapping with NA
        tuples = [(na_value, na_value)] + tuples + [(na_value, na_value)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        assert index.is_overlapping is False

        # overlapping
        tuples = [(start + n * shift, start + (n + 2) * shift)
                  for n in range(3)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        assert index.is_overlapping is True

        # overlapping with NA
        tuples = [(na_value, na_value)] + tuples + [(na_value, na_value)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        assert index.is_overlapping is True

        # common endpoints
        tuples = [(start + n * shift, start + (n + 1) * shift)
                  for n in range(3)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        result = index.is_overlapping
        expected = closed == 'both'
        assert result is expected

        # common endpoints with NA
        tuples = [(na_value, na_value)] + tuples + [(na_value, na_value)]
        index = IntervalIndex.from_tuples(tuples, closed=closed)
        result = index.is_overlapping
        assert result is expected

    @pytest.mark.parametrize('tuples', [
        lzip(range(10), range(1, 11)),
        lzip(date_range('20170101', periods=10),
             date_range('20170101', periods=10)),
        lzip(timedelta_range('0 days', periods=10),
             timedelta_range('1 day', periods=10))])
    def test_to_tuples(self, tuples):
        # GH 18756
        idx = IntervalIndex.from_tuples(tuples)
        result = idx.to_tuples()
        expected = Index(com.asarray_tuplesafe(tuples))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('tuples', [
        lzip(range(10), range(1, 11)) + [np.nan],
        lzip(date_range('20170101', periods=10),
             date_range('20170101', periods=10)) + [np.nan],
        lzip(timedelta_range('0 days', periods=10),
             timedelta_range('1 day', periods=10)) + [np.nan]])
    @pytest.mark.parametrize('na_tuple', [True, False])
    def test_to_tuples_na(self, tuples, na_tuple):
        # GH 18756
        idx = IntervalIndex.from_tuples(tuples)
        result = idx.to_tuples(na_tuple=na_tuple)

        # check the non-NA portion
        expected_notna = Index(com.asarray_tuplesafe(tuples[:-1]))
        result_notna = result[:-1]
        tm.assert_index_equal(result_notna, expected_notna)

        # check the NA portion
        result_na = result[-1]
        if na_tuple:
            assert isinstance(result_na, tuple)
            assert len(result_na) == 2
            assert all(isna(x) for x in result_na)
        else:
            assert isna(result_na)

    def test_nbytes(self):
        # GH 19209
        left = np.arange(0, 4, dtype='i8')
        right = np.arange(1, 5, dtype='i8')

        result = IntervalIndex.from_arrays(left, right).nbytes
        expected = 64  # 4 * 8 * 2
        assert result == expected

    def test_itemsize(self):
        # GH 19209
        left = np.arange(0, 4, dtype='i8')
        right = np.arange(1, 5, dtype='i8')
        expected = 16  # 8 * 2

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = IntervalIndex.from_arrays(left, right).itemsize

        assert result == expected

    @pytest.mark.parametrize('new_closed', [
        'left', 'right', 'both', 'neither'])
    def test_set_closed(self, name, closed, new_closed):
        # GH 21670
        index = interval_range(0, 5, closed=closed, name=name)
        result = index.set_closed(new_closed)
        expected = interval_range(0, 5, closed=new_closed, name=name)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('bad_closed', ['foo', 10, 'LEFT', True, False])
    def test_set_closed_errors(self, bad_closed):
        # GH 21670
        index = interval_range(0, 5)
        msg = "invalid option for 'closed': {closed}".format(closed=bad_closed)
        with pytest.raises(ValueError, match=msg):
            index.set_closed(bad_closed)

    def test_is_all_dates(self):
        # GH 23576
        year_2017 = pd.Interval(pd.Timestamp('2017-01-01 00:00:00'),
                                pd.Timestamp('2018-01-01 00:00:00'))
        year_2017_index = pd.IntervalIndex([year_2017])
        assert not year_2017_index.is_all_dates
