import pytest
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, IntervalIndex, Interval
import pandas.util.testing as tm


class TestIntervalIndex(object):

    def setup_method(self, method):
        self.s = Series(np.arange(5), IntervalIndex.from_breaks(np.arange(6)))

    def test_loc_with_scalar(self):

        s = self.s
        expected = 0

        result = s.loc[0.5]
        assert result == expected

        result = s.loc[1]
        assert result == expected

        with pytest.raises(KeyError):
            s.loc[0]

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s.loc[:3])
        tm.assert_series_equal(expected, s.loc[:2.5])
        tm.assert_series_equal(expected, s.loc[0.1:2.5])
        tm.assert_series_equal(expected, s.loc[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s.loc[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s.loc[[2, 3, 4]])
        tm.assert_series_equal(expected, s.loc[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s.loc[s >= 2])

    def test_getitem_with_scalar(self):

        s = self.s
        expected = 0

        result = s[0.5]
        assert result == expected

        result = s[1]
        assert result == expected

        with pytest.raises(KeyError):
            s[0]

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s[:3])
        tm.assert_series_equal(expected, s[:2.5])
        tm.assert_series_equal(expected, s[0.1:2.5])
        tm.assert_series_equal(expected, s[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s[[2, 3, 4]])
        tm.assert_series_equal(expected, s[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s[s >= 2])

    def test_with_interval(self):

        s = self.s
        expected = 0

        result = s.loc[Interval(0, 1)]
        assert result == expected

        result = s[Interval(0, 1)]
        assert result == expected

        expected = s.iloc[3:5]
        result = s.loc[Interval(3, 6)]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[3:5]
        result = s.loc[[Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[3:5]
        result = s.loc[[Interval(3, 5)]]
        tm.assert_series_equal(expected, result)

        # missing
        with pytest.raises(KeyError):
            s.loc[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s.loc[Interval(5, 6)]

        with pytest.raises(KeyError):
            s[Interval(5, 6)]


    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_with_interval_updated_behavior(self):

        s = self.s
        expected = 0

        result = s.loc[Interval(0, 1)]
        assert result == expected

        result = s[Interval(0, 1)]
        assert result == expected

        # missing
        with pytest.raises(KeyError):
            s.loc[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s.loc[Interval(5, 6)]

        with pytest.raises(KeyError):
            s[Interval(5, 6)]

        # overlapping but not identical
        # right:
        with pytest.raises(KeyError):
            s.loc[Interval(4, 6)]

        with pytest.raises(KeyError):
            s.loc[[Interval(4, 6)]]

        # left:
        with pytest.raises(KeyError):
            s.loc[Interval(-1, 4)]

        with pytest.raises(KeyError):
            s.loc[[Interval(-1, 4)]]

        with pytest.raises(KeyError):
            s.loc[Interval(0, 1, closed='left')]

        with pytest.raises(KeyError):
            s.loc[[Interval(0, 1, closed='left')]]

        with pytest.raises(KeyError):
            s.loc[Interval(0, 1, closed='both')]

        with pytest.raises(KeyError):
            s.loc[[Interval(0, 1, closed='both')]]

        # inner
        with pytest.raises(KeyError):
            s.loc[Interval(2, 4)]

        with pytest.raises(KeyError):
            s.loc[[Interval(2, 4)]]

        with pytest.raises(KeyError):
            s.loc[Interval(4, 5, closed='left')]

        with pytest.raises(KeyError):
            s.loc[[Interval(4, 5, closed='left')]]

        with pytest.raises(KeyError):
            s.loc[Interval(4, 5, closed='both')]

        with pytest.raises(KeyError):
            s.loc[[Interval(4, 5, closed='both')]]

        # outer
        with pytest.raises(KeyError):
            s.loc[Interval(-1, 6)]

        with pytest.raises(KeyError):
            s.loc[[Interval(-1, 6)]]

        with pytest.raises(KeyError):
            s.loc[Interval(0, 6, closed='left')]

        with pytest.raises(KeyError):
            s.loc[[Interval(0, 6, closed='left')]]

        with pytest.raises(KeyError):
            s.loc[Interval(0, 5, closed='both')]

        with pytest.raises(KeyError):
            s.loc[[Interval(0, 5, closed='both')]]

    def test_with_slices(self):

        s = self.s

        # slice of interval
        with pytest.raises(NotImplementedError):
            result = s.loc[Interval(3, 6):]

        with pytest.raises(NotImplementedError):
            result = s[Interval(3, 6):]

        expected = s.iloc[3:5]
        result = s[[Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        # slice of scalar with step != 1
        with pytest.raises(ValueError):
            s[0:4:2]

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_with_slices_updated_behavior(self):

        s = self.s

        # slice of interval
        expected = s.iloc[4:]
        result = s.loc[Interval(3, 4):]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[4:]
        result = s[Interval(3, 4):]
        tm.assert_series_equal(expected, result)

        with pytest.raises(KeyError):
            s.loc[Interval(3, 6):]

        with pytest.raises(KeyError):
            s[Interval(3, 6):]

        with pytest.raises(KeyError):
            s.loc[Interval(3, 4, closed='left'):]

        with pytest.raises(KeyError):
            s[Interval(3, 4, closed='left'):]

        with pytest.raises(KeyError):
            s.loc[Interval(3, 4, closed='both'):]

        with pytest.raises(KeyError):
            s[Interval(3, 4, closed='both'):]

        # slice of scalar
        with pytest.raises(NotImplementedError):
            s[0:4] ## not sure what the behvaior should be here.

        # slice of scalar with step != 1
        with pytest.raises(ValueError):
            s[0:4:2] ## This should probably definitely fail I guess?

    def test_with_overlaps(self):

        s = self.s
        expected = s.iloc[[3, 4, 3, 4]]
        result = s.loc[[Interval(3, 6), Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        idx = IntervalIndex.from_tuples([(1, 5), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s[4]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s[[4]]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s.loc[[4]]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s[Interval(3, 5)]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s.loc[Interval(3, 5)]
        expected = s
        tm.assert_series_equal(expected, result)

        # doesn't intersect unique set of intervals
        with pytest.raises(KeyError):
            s[[Interval(3, 5)]]

        with pytest.raises(KeyError):
            s.loc[[Interval(3, 5)]]

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_with_overlaps_updated_behavior(self):

        idx = IntervalIndex.from_tuples([(1, 5), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        # scalar
        expected = s
        result = s[4]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s[[4]]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s.loc[[4]]
        tm.assert_series_equal(expected, result)

        # interval
        with pytest.raises(KeyError):
            s[Interval(3, 5)]

        with pytest.raises(KeyError):
            s[[Interval(3, 5)]]

        with pytest.raises(KeyError):
            s.loc[Interval(3, 5)]

        with pytest.raises(KeyError):
            s.loc[[Interval(3, 5)]]

    def test_non_unique(self):

        idx = IntervalIndex.from_tuples([(1, 3), (3, 7)])
        s = pd.Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        assert result == 0

        result = s.loc[[Interval(1, 3)]]
        expected = s.iloc[0:1]
        tm.assert_series_equal(expected, result)

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_non_unique_updated_behavior(self):

        # Actually I think we should remove this test? Not sure what exactly it's meant to gauge...
        pass

    def test_non_unique_moar(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        expected = s.iloc[[0, 1]]
        tm.assert_series_equal(expected, result)

        # non-unique index and slices not allowed
        with pytest.raises(ValueError):
            s.loc[Interval(1, 3):]

        with pytest.raises(ValueError):
            s[Interval(1, 3):]

        # non-unique
        with pytest.raises(ValueError):
            s[[Interval(1, 3)]]

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_non_unique_moar_updated_behavior(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        expected = s.iloc[[0, 1]]
        result = s.loc[Interval(1, 3)]
        tm.assert_series_equal(expected, result)

        # non-unique index and slices not allowed
        with pytest.raises(ValueError):
            s.loc[Interval(1, 3):]
        # this is confusing to me. I would have done:
        # expected = s
        # result = s.loc[Interval(1, 3):]
        # tm.assert_series_equal(expected, result)

        with pytest.raises(ValueError):
            s[Interval(1, 3):]
        # Same here:
        # expected = s
        # result = s[Interval(1, 3):]
        # tm.assert_series_equal(expected, result)

        # non-unique
        with pytest.raises(ValueError):
            s[[Interval(1, 3)]]
        # not sure why the behavior for [[]] is different than []...

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_non_unique_moar_updated_behavior(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        expected = s.iloc[[0, 1]]
        tm.assert_series_equal(expected, result)

        # non-unique index and slices not allowed
        with pytest.raises(ValueError):
            s.loc[Interval(1, 3):]

        with pytest.raises(ValueError):
            s[Interval(1, 3):]

        # non-unique
        with pytest.raises(ValueError):
            s[[Interval(1, 3)]]


    def test_non_matching(self):

        s = self.s

        # this is a departure from our current indexing scheme, but simpler
        with pytest.raises(KeyError):
            s.loc[[-1, 3, 4, 5]]

        with pytest.raises(KeyError):
            s.loc[[-1, 3]]

    def test_large_series(self):
        s = Series(np.arange(1000000),
                   index=IntervalIndex.from_breaks(np.arange(1000001)))

        result1 = s.loc[:80000]
        result2 = s.loc[0:80000]
        result3 = s.loc[0:80000:1]
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

    def test_loc_getitem_frame(self):

        df = DataFrame({'A': range(10)})
        s = pd.cut(df.A, 5)
        df['B'] = s
        df = df.set_index('B')

        result = df.loc[4]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError):
            df.loc[10]

        # single list-like
        result = df.loc[[4]]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        # non-unique
        result = df.loc[[4, 5]]
        expected = df.take([4, 5, 4, 5])
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError):
            df.loc[[10]]

        # partial missing
        with pytest.raises(KeyError):
            df.loc[[10, 4]]


    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_interval_covers(self):

        # class Interval:
        #     def covers(self, other: Interval) -> bool
        #     def covers(self, other: IntervalIndex) -> IntegerArray1D

        assert     Interval(1, 3).covers(Interval(1.5, 2.5))
        assert     Interval(1, 3).covers(Interval(1, 2))
        assert     Interval(1, 3).covers(Interval(2, 3))
        assert not Interval(1, 3).covers(Interval(0.5, 2.5))
        assert not Interval(1, 3).covers(Interval(1.5, 3.5))

        assert     Interval(1, 3, closed='right').covers(Interval(1, 3, closed='right'))
        assert not Interval(1, 3, closed='right').covers(Interval(1, 3, closed='left'))
        assert not Interval(1, 3, closed='right').covers(Interval(1, 3, closed='both'))

        assert not Interval(1, 3, closed='left').covers(Interval(1, 3, closed='right'))
        assert     Interval(1, 3, closed='left').covers(Interval(1, 3, closed='left'))
        assert not Interval(1, 3, closed='left').covers(Interval(1, 3, closed='both'))

        assert     Interval(1, 3, closed='both').covers(Interval(1, 3, closed='right'))
        assert     Interval(1, 3, closed='both').covers(Interval(1, 3, closed='left'))
        assert     Interval(1, 3, closed='both').covers(Interval(1, 3, closed='both'))

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)])

        assert     Interval(1, 3, closed='right').covers(idx) == np.array([1, 2])
        assert     Interval(0, 3, closed='right').covers(idx) == np.array([0, 1, 2])
        assert     Interval(0, 2, closed='right').covers(idx) == np.array([0])
        assert     Interval(2, 4, closed='right').covers(idx) == np.array([1])

        assert     Interval(1, 3, closed='left').covers(idx) == np.array([])
        assert     Interval(0, 3, closed='left').covers(idx) == np.array([0])
        assert     Interval(0, 2, closed='left').covers(idx) == np.array([0])
        assert     Interval(2, 4, closed='left').covers(idx) == np.array([1])

        assert     Interval(1, 3, closed='both').covers(idx) == np.array([1, 2])
        assert     Interval(0, 5, closed='both').covers(idx) == np.array([0, 1, 2])
        assert     Interval(0, 2, closed='both').covers(idx) == np.array([0])
        assert     Interval(2, 4, closed='both').covers(idx) == np.array([1])

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_interval_overlaps(self):

        # class Interval:
        #     def overlaps(self, other: Interval) -> bool
        #     def overlaps(self, other: IntervalIndex) -> IntegerArray1D

        assert     Interval(1, 3).overlaps(Interval(1.5, 2.5))
        assert     Interval(1, 3).overlaps(Interval(1, 2))
        assert     Interval(1, 3).overlaps(Interval(2, 3))
        assert     Interval(1, 3).overlaps(Interval(0.5, 2.5))
        assert     Interval(1, 3).overlaps(Interval(1.5, 3.5))

        assert not Interval(1, 3).overlaps(Interval(-1, 1))
        assert not Interval(1, 3).overlaps(Interval(3, 5))

        # right
        assert     Interval(1, 3, closed='right').overlaps(Interval(1, 3, closed='right'))
        assert     Interval(1, 3, closed='right').overlaps(Interval(1, 3, closed='left'))
        assert     Interval(1, 3, closed='right').overlaps(Interval(1, 3, closed='both'))

        assert not Interval(1, 3, closed='right').overlaps(Interval(-1, 1, closed='right'))
        assert not Interval(1, 3, closed='right').overlaps(Interval(-1, 1, closed='left'))
        assert not Interval(1, 3, closed='right').overlaps(Interval(-1, 1, closed='both'))

        assert not Interval(1, 3, closed='right').overlaps(Interval(3, 5, closed='right'))
        assert     Interval(1, 3, closed='right').overlaps(Interval(3, 5, closed='left'))
        assert     Interval(1, 3, closed='right').overlaps(Interval(3, 5, closed='both'))

        # left
        assert     Interval(1, 3, closed='left').overlaps(Interval(1, 3, closed='right'))
        assert     Interval(1, 3, closed='left').overlaps(Interval(1, 3, closed='left'))
        assert     Interval(1, 3, closed='left').overlaps(Interval(1, 3, closed='both'))

        assert not Interval(1, 3, closed='left').overlaps(Interval(-1, 1, closed='right'))
        assert not Interval(1, 3, closed='left').overlaps(Interval(-1, 1, closed='left'))
        assert not Interval(1, 3, closed='left').overlaps(Interval(-1, 1, closed='both'))

        assert not Interval(1, 3, closed='left').overlaps(Interval(3, 5, closed='right'))
        assert     Interval(1, 3, closed='left').overlaps(Interval(3, 5, closed='left'))
        assert     Interval(1, 3, closed='left').overlaps(Interval(3, 5, closed='both'))

        # both
        assert     Interval(1, 3, closed='both').overlaps(Interval(1, 3, closed='right'))
        assert     Interval(1, 3, closed='both').overlaps(Interval(1, 3, closed='left'))
        assert     Interval(1, 3, closed='both').overlaps(Interval(1, 3, closed='both'))

        assert     Interval(1, 3, closed='both').overlaps(Interval(-1, 1, closed='right'))
        assert not Interval(1, 3, closed='both').overlaps(Interval(-1, 1, closed='left'))
        assert     Interval(1, 3, closed='both').overlaps(Interval(-1, 1, closed='both'))

        assert not Interval(1, 3, closed='both').overlaps(Interval(3, 5, closed='right'))
        assert     Interval(1, 3, closed='both').overlaps(Interval(3, 5, closed='left'))
        assert     Interval(1, 3, closed='both').overlaps(Interval(3, 5, closed='both'))

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)])

        assert     Interval(1, 3, closed='right').overlaps(idx) == np.array([1, 2])
        assert     Interval(1, 2, closed='right').overlaps(idx) == np.array([2])
        assert     Interval(0, 2, closed='right').overlaps(idx) == np.array([0, 2])
        assert     Interval(3, 4, closed='right').overlaps(idx) == np.array([])

        assert     Interval(1, 3, closed='left').overlaps(idx) == np.array([0, 1, 2])
        assert     Interval(1, 2, closed='left').overlaps(idx) == np.array([0, 2])
        assert     Interval(0, 2, closed='left').overlaps(idx) == np.array([0, 2])
        assert     Interval(3, 4, closed='left').overlaps(idx) == np.array([3])

        assert     Interval(1, 3, closed='both').overlaps(idx) == np.array([0, 1, 2])
        assert     Interval(1, 2, closed='both').overlaps(idx) == np.array([0, 2])
        assert     Interval(0, 2, closed='both').overlaps(idx) == np.array([0, 2])
        assert     Interval(3, 4, closed='both').overlaps(idx) == np.array([3])

    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_intervalIndex_covers(self):

        # class IntervalIndex:
        #     def covers(self, other: Interval) -> IntegerArray1D
        #     def covers(self, other: IntervalIndex) -> Tuple[IntegerArray1D, IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)])

        assert     idx.covers(Interval(1, 3, closed='right')) == np.array([1, 2])
        assert     idx.covers(Interval(0, 3, closed='right')) == np.array([0, 1, 2])
        assert     idx.covers(Interval(0, 2, closed='right')) == np.array([0])
        assert     idx.covers(Interval(2, 4, closed='right')) == np.array([1])

        assert     idx.covers(Interval(1, 3, closed='left')) == np.array([])
        assert     idx.covers(Interval(0, 3, closed='left')) == np.array([0])
        assert     idx.covers(Interval(0, 2, closed='left')) == np.array([0])
        assert     idx.covers(Interval(2, 4, closed='left')) == np.array([1])

        assert     idx.covers(Interval(1, 3, closed='both')) == np.array([1, 2])
        assert     idx.covers(Interval(0, 5, closed='both')) == np.array([0, 1, 2])
        assert     idx.covers(Interval(0, 2, closed='both')) == np.array([0])
        assert     idx.covers(Interval(2, 4, closed='both')) == np.array([1])

        idx1 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='right')
        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='left')
        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='both')

        assert     idx.covers(idx1) == (np.array([0,1,2,2]), np.array([0,1,1,2]))  # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...
        assert     idx.covers(idx2) == (np.array([2]), np.array([1]))              # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...
        assert     idx.covers(idx3) == (np.array([0,1,2,2]), np.array([0,1,1,2]))  # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...


    @pytest.mark.xfail(reason="new indexing tests for issue 16316")
    def test_intervalIndex_overlaps(self):

        # class IntervalIndex:
        #     def overlaps(self, other: Interval) -> IntegerArray1D
        #     def overlaps(self, other: IntervalIndex) -> Tuple[IntegerArray1D, IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)])

        assert     idx.overlaps(Interval(1, 3, closed='right')) == np.array([1, 2])
        assert     idx.overlaps(Interval(1, 2, closed='right')) == np.array([2])
        assert     idx.overlaps(Interval(0, 2, closed='right')) == np.array([0, 2])
        assert     idx.overlaps(Interval(3, 4, closed='right')) == np.array([])

        assert     idx.overlaps(Interval(1, 3, closed='left')) == np.array([0, 1, 2])
        assert     idx.overlaps(Interval(1, 2, closed='left')) == np.array([0, 2])
        assert     idx.overlaps(Interval(0, 2, closed='left')) == np.array([0, 2])
        assert     idx.overlaps(Interval(3, 4, closed='left')) == np.array([3])

        assert     idx.overlaps(Interval(1, 3, closed='both')) == np.array([0, 1, 2])
        assert     idx.overlaps(Interval(1, 2, closed='both')) == np.array([0, 2])
        assert     idx.overlaps(Interval(0, 2, closed='both')) == np.array([0, 2])
        assert     idx.overlaps(Interval(3, 4, closed='both')) == np.array([3])

        idx1 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='right')
        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='left')
        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='both')

        assert     idx.overlaps(idx1) == (np.array([0,1,2,2]), np.array([0,1,1,2]))          # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...
        assert     idx.overlaps(idx2) == (np.array([0,0,1,1,2,2]), np.array([0,2,1,2,1,2]))  # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...
        assert     idx.overlaps(idx3) == (np.array([0,0,1,1,2,2]), np.array([0,2,1,2,1,2]))  # note: I had to choose an arbitrary ordering. If this test fails, double check the test too...

