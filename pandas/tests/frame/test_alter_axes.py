from datetime import datetime

import numpy as np
import pytz

from pandas import DataFrame, cut, date_range
import pandas._testing as tm


class TestDataFrameAlterAxes:
    def test_set_axis_setattr_index(self):
        # GH 6785
        # set the index manually

        df = DataFrame([{"ts": datetime(2014, 4, 1, tzinfo=pytz.utc), "foo": 1}])
        expected = df.set_index("ts")
        df.index = df["ts"]
        df.pop("ts")
        tm.assert_frame_equal(df, expected)

    def test_dti_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = date_range("2011/01/01", periods=6, freq="M", tz="US/Eastern")
        idx2 = date_range("2013", periods=6, freq="A", tz="Asia/Tokyo")

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.reindex(idx2)
        tm.assert_index_equal(df.index, idx2)

    def test_dti_set_index_reindex_with_tz(self):
        # GH 11314
        # with tz
        index = date_range(
            datetime(2015, 10, 1), datetime(2015, 10, 1, 23), freq="H", tz="US/Eastern"
        )
        df = DataFrame(np.random.randn(24, 1), columns=["a"], index=index)
        new_index = date_range(
            datetime(2015, 10, 2), datetime(2015, 10, 2, 23), freq="H", tz="US/Eastern"
        )

        result = df.set_index(new_index)
        assert result.index.freq == index.freq

    # Renaming

    def test_assign_columns(self, float_frame):
        float_frame["hi"] = "there"

        df = float_frame.copy()
        df.columns = ["foo", "bar", "baz", "quux", "foo2"]
        tm.assert_series_equal(float_frame["C"], df["baz"], check_names=False)
        tm.assert_series_equal(float_frame["hi"], df["foo2"], check_names=False)


class TestIntervalIndex:
    def test_set_reset_index(self):

        df = DataFrame({"A": range(10)})
        s = cut(df.A, 5)
        df["B"] = s
        df = df.set_index("B")

        df = df.reset_index()
