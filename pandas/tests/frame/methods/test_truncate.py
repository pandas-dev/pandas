import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm


class TestDataFrameTruncate:
    def test_truncate(self, datetime_frame):
        ts = datetime_frame[::3]

        start, end = datetime_frame.index[3], datetime_frame.index[6]

        start_missing = datetime_frame.index[2]
        end_missing = datetime_frame.index[7]

        # neither specified
        truncated = ts.truncate()
        tm.assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        tm.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        tm.assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        tm.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        tm.assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        tm.assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        tm.assert_frame_equal(truncated, expected)

        msg = "Truncate: 2000-01-06 00:00:00 must be after 2000-02-04 00:00:00"
        with pytest.raises(ValueError, match=msg):
            ts.truncate(
                before=ts.index[-1] - ts.index.freq, after=ts.index[0] + ts.index.freq
            )

    def test_truncate_copy(self, datetime_frame):
        index = datetime_frame.index
        truncated = datetime_frame.truncate(index[5], index[10])
        truncated.values[:] = 5.0
        assert not (datetime_frame.values[5:11] == 5).any()

    def test_truncate_nonsortedindex(self):
        # GH#17935

        df = pd.DataFrame({"A": ["a", "b", "c", "d", "e"]}, index=[5, 3, 2, 9, 0])
        msg = "truncate requires a sorted index"
        with pytest.raises(ValueError, match=msg):
            df.truncate(before=3, after=9)

        rng = pd.date_range("2011-01-01", "2012-01-01", freq="W")
        ts = pd.DataFrame(
            {"A": np.random.randn(len(rng)), "B": np.random.randn(len(rng))}, index=rng
        )
        msg = "truncate requires a sorted index"
        with pytest.raises(ValueError, match=msg):
            ts.sort_values("A", ascending=False).truncate(
                before="2011-11", after="2011-12"
            )

        df = pd.DataFrame(
            {
                3: np.random.randn(5),
                20: np.random.randn(5),
                2: np.random.randn(5),
                0: np.random.randn(5),
            },
            columns=[3, 20, 2, 0],
        )
        msg = "truncate requires a sorted index"
        with pytest.raises(ValueError, match=msg):
            df.truncate(before=2, after=20, axis=1)
