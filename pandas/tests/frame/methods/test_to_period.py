import numpy as np
import pytest

from pandas import DataFrame, date_range, period_range
import pandas._testing as tm


class TestToPeriod:
    def test_frame_to_period(self):
        K = 5

        dr = date_range("1/1/2000", "1/1/2001")
        pr = period_range("1/1/2000", "1/1/2001")
        df = DataFrame(np.random.randn(len(dr), K), index=dr)
        df["mix"] = "a"

        pts = df.to_period()
        exp = df.copy()
        exp.index = pr
        tm.assert_frame_equal(pts, exp)

        pts = df.to_period("M")
        tm.assert_index_equal(pts.index, exp.index.asfreq("M"))

        df = df.T
        pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = pr
        tm.assert_frame_equal(pts, exp)

        pts = df.to_period("M", axis=1)
        tm.assert_index_equal(pts.columns, exp.columns.asfreq("M"))

        msg = "No axis named 2 for object type <class 'pandas.core.frame.DataFrame'>"
        with pytest.raises(ValueError, match=msg):
            df.to_period(axis=2)
