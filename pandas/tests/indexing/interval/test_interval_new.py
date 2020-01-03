import re

import numpy as np
import pytest

from pandas import Interval, IntervalIndex, Series
import pandas._testing as tm


class TestIntervalIndex:
    def setup_method(self, method):
        self.s = Series(np.arange(5), IntervalIndex.from_breaks(np.arange(6)))

    def test_loc_with_interval(self):

        # loc with single label / list of labels:
        #   - Intervals: only exact matches
        #   - scalars: those that contain it

        s = self.s

        expected = 0
        result = s.loc[Interval(0, 1)]
        assert result == expected
        result = s[Interval(0, 1)]
        assert result == expected

        expected = s.iloc[3:5]
        result = s.loc[[Interval(3, 4), Interval(4, 5)]]
        tm.assert_series_equal(expected, result)
        result = s[[Interval(3, 4), Interval(4, 5)]]
        tm.assert_series_equal(expected, result)

        # missing or not exact
        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='left')")):
            s.loc[Interval(3, 5, closed="left")]

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='left')")):
            s[Interval(3, 5, closed="left")]

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='right')")):
            s[Interval(3, 5)]

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='right')")):
            s.loc[Interval(3, 5)]

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='right')")):
            s[Interval(3, 5)]

        with pytest.raises(
            KeyError, match=re.escape("Interval(-2, 0, closed='right')")
        ):
            s.loc[Interval(-2, 0)]

        with pytest.raises(
            KeyError, match=re.escape("Interval(-2, 0, closed='right')")
        ):
            s[Interval(-2, 0)]

        with pytest.raises(KeyError, match=re.escape("Interval(5, 6, closed='right')")):
            s.loc[Interval(5, 6)]

        with pytest.raises(KeyError, match=re.escape("Interval(5, 6, closed='right')")):
            s[Interval(5, 6)]

    def test_loc_with_scalar(self):

        # loc with single label / list of labels:
        #   - Intervals: only exact matches
        #   - scalars: those that contain it

        s = self.s

        assert s.loc[1] == 0
        assert s.loc[1.5] == 1
        assert s.loc[2] == 1

        assert s[1] == 0
        assert s[1.5] == 1
        assert s[2] == 1

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s.loc[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s.loc[[2, 3, 4]])
        tm.assert_series_equal(expected, s.loc[[1.5, 3, 4]])

        expected = s.iloc[[1, 1, 2, 1]]
        tm.assert_series_equal(expected, s.loc[[1.5, 2, 2.5, 1.5]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s.loc[s >= 2])

    def test_loc_with_slices(self):

        # loc with slices:
        #   - Interval objects: only works with exact matches
        #   - scalars: only works for non-overlapping, monotonic intervals,
        #     and start/stop select location based on the interval that
        #     contains them:
        #    (slice_loc(start, stop) == (idx.get_loc(start), idx.get_loc(stop))

        s = self.s

        # slice of interval

        expected = s.iloc[:3]
        result = s.loc[Interval(0, 1) : Interval(2, 3)]
        tm.assert_series_equal(expected, result)
        result = s[Interval(0, 1) : Interval(2, 3)]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[3:]
        result = s.loc[Interval(3, 4) :]
        tm.assert_series_equal(expected, result)
        result = s[Interval(3, 4) :]
        tm.assert_series_equal(expected, result)

        msg = "Interval objects are not currently supported"
        with pytest.raises(NotImplementedError, match=msg):
            s.loc[Interval(3, 6) :]

        with pytest.raises(NotImplementedError, match=msg):
            s[Interval(3, 6) :]

        with pytest.raises(NotImplementedError, match=msg):
            s.loc[Interval(3, 4, closed="left") :]

        with pytest.raises(NotImplementedError, match=msg):
            s[Interval(3, 4, closed="left") :]

        # TODO with non-existing intervals ?
        # s.loc[Interval(-1, 0):Interval(2, 3)]

        # slice of scalar

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s.loc[:3])
        tm.assert_series_equal(expected, s.loc[:2.5])
        tm.assert_series_equal(expected, s.loc[0.1:2.5])
        tm.assert_series_equal(expected, s.loc[-1:3])

        tm.assert_series_equal(expected, s[:3])
        tm.assert_series_equal(expected, s[:2.5])
        tm.assert_series_equal(expected, s[0.1:2.5])

        # slice of scalar with step != 1
        with pytest.raises(ValueError):
            s[0:4:2]

    def test_loc_with_overlap(self):

        idx = IntervalIndex.from_tuples([(1, 5), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        # scalar
        expected = s
        result = s.loc[4]
        tm.assert_series_equal(expected, result)

        result = s[4]
        tm.assert_series_equal(expected, result)

        result = s.loc[[4]]
        tm.assert_series_equal(expected, result)

        result = s[[4]]
        tm.assert_series_equal(expected, result)

        # interval
        expected = 0
        result = s.loc[Interval(1, 5)]
        result == expected

        result = s[Interval(1, 5)]
        result == expected

        expected = s
        result = s.loc[[Interval(1, 5), Interval(3, 7)]]
        tm.assert_series_equal(expected, result)

        result = s[[Interval(1, 5), Interval(3, 7)]]
        tm.assert_series_equal(expected, result)

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='right')")):
            s.loc[Interval(3, 5)]

        with pytest.raises(KeyError, match="^$"):
            s.loc[[Interval(3, 5)]]

        with pytest.raises(KeyError, match=re.escape("Interval(3, 5, closed='right')")):
            s[Interval(3, 5)]

        with pytest.raises(KeyError, match="^$"):
            s[[Interval(3, 5)]]

        # slices with interval (only exact matches)
        expected = s
        result = s.loc[Interval(1, 5) : Interval(3, 7)]
        tm.assert_series_equal(expected, result)

        result = s[Interval(1, 5) : Interval(3, 7)]
        tm.assert_series_equal(expected, result)

        msg = "'can only get slices from an IntervalIndex if bounds are"
        " non-overlapping and all monotonic increasing or decreasing'"
        with pytest.raises(KeyError, match=msg):
            s.loc[Interval(1, 6) : Interval(3, 8)]

        with pytest.raises(KeyError, match=msg):
            s[Interval(1, 6) : Interval(3, 8)]

        # slices with scalar raise for overlapping intervals
        # TODO KeyError is the appropriate error?
        with pytest.raises(KeyError, match=msg):
            s.loc[1:4]

    def test_non_unique(self):

        idx = IntervalIndex.from_tuples([(1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        assert result == 0

        result = s.loc[[Interval(1, 3)]]
        expected = s.iloc[0:1]
        tm.assert_series_equal(expected, result)

    def test_non_unique_moar(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        expected = s.iloc[[0, 1]]
        result = s.loc[Interval(1, 3)]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s.loc[Interval(1, 3) :]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s[Interval(1, 3) :]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[[0, 1]]
        result = s[[Interval(1, 3)]]
        tm.assert_series_equal(expected, result)
