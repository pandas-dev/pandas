import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, IntervalIndex, Series
import pandas.util.testing as tm


class TestIntervalIndex:
    def setup_method(self, method):
        self.s = Series(np.arange(5), IntervalIndex.from_breaks(np.arange(6)))

    def test_getitem_with_scalar(self):

        s = self.s

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s[:3])
        tm.assert_series_equal(expected, s[:2.5])
        tm.assert_series_equal(expected, s[0.1:2.5])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s[[2, 3, 4]])
        tm.assert_series_equal(expected, s[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s[s >= 2])

    @pytest.mark.parametrize("direction", ["increasing", "decreasing"])
    def test_nonoverlapping_monotonic(self, direction, closed):
        tpls = [(0, 1), (2, 3), (4, 5)]
        if direction == "decreasing":
            tpls = tpls[::-1]

        idx = IntervalIndex.from_tuples(tpls, closed=closed)
        s = Series(list("abc"), idx)

        for key, expected in zip(idx.left, s):
            if idx.closed_left:
                assert s[key] == expected
                assert s.loc[key] == expected
            else:
                with pytest.raises(KeyError, match=str(key)):
                    s[key]
                with pytest.raises(KeyError, match=str(key)):
                    s.loc[key]

        for key, expected in zip(idx.right, s):
            if idx.closed_right:
                assert s[key] == expected
                assert s.loc[key] == expected
            else:
                with pytest.raises(KeyError, match=str(key)):
                    s[key]
                with pytest.raises(KeyError, match=str(key)):
                    s.loc[key]

        for key, expected in zip(idx.mid, s):
            assert s[key] == expected
            assert s.loc[key] == expected

    def test_non_matching(self):
        s = self.s

        # this is a departure from our current
        # indexin scheme, but simpler
        with pytest.raises(KeyError, match="^$"):
            s.loc[[-1, 3, 4, 5]]

        with pytest.raises(KeyError, match="^$"):
            s.loc[[-1, 3]]

    def test_large_series(self):
        s = Series(
            np.arange(1000000), index=IntervalIndex.from_breaks(np.arange(1000001))
        )

        result1 = s.loc[:80000]
        result2 = s.loc[0:80000]
        result3 = s.loc[0:80000:1]
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

    def test_loc_getitem_frame(self):

        df = DataFrame({"A": range(10)})
        s = pd.cut(df.A, 5)
        df["B"] = s
        df = df.set_index("B")

        result = df.loc[4]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError, match="10"):
            df.loc[10]

        # single list-like
        result = df.loc[[4]]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        # non-unique
        result = df.loc[[4, 5]]
        expected = df.take([4, 5, 4, 5])
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError, match="^$"):
            df.loc[[10]]

        # partial missing
        with pytest.raises(KeyError, match="^$"):
            df.loc[[10, 4]]
