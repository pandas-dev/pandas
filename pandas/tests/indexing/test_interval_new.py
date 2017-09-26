import pytest
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, IntervalIndex, Interval
from pandas.compat import product
import pandas.util.testing as tm


class TestIntervalIndex_new(object):

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

        with pytest.raises(NotImplementedError):
            s.loc[Interval(3, 6):]

        with pytest.raises(NotImplementedError):
            s[Interval(3, 6):]

        with pytest.raises(KeyError):
            s[Interval(3, 4, closed='both'):]

        # slice of scalar
        expected = s.iloc[:4]  # maybe [:5] ?
        result = s[0:4]
        tm.assert_series_equal(expected, result)

        # slice of scalar with step != 1
        with pytest.raises(NotImplementedError):
            s[0:4:2]

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
    def test_non_unique_moar_updated_behavior(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        expected = s.iloc[[0, 1]]
        result = s.loc[Interval(1, 3)]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s.loc[Interval(1, 3):]
        tm.assert_series_equal(expected, result)

        expected = s
        result = s[Interval(1, 3):]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[[0, 1]]
        result = s[[Interval(1, 3)]]
        tm.assert_series_equal(expected, result)
