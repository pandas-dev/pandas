import numpy as np

from pandas import DataFrame, Index, PeriodIndex, period_range
import pandas._testing as tm


class TestPeriodIndex:
    def test_as_frame_columns(self):
        rng = period_range("1/1/2000", periods=5)
        df = DataFrame(np.random.randn(10, 5), columns=rng)

        ts = df[rng[0]]
        tm.assert_series_equal(ts, df.iloc[:, 0])

        # GH # 1211
        repr(df)

        ts = df["1/1/2000"]
        tm.assert_series_equal(ts, df.iloc[:, 0])

    def test_frame_setitem(self):
        rng = period_range("1/1/2000", periods=5, name="index")
        df = DataFrame(np.random.randn(5, 3), index=rng)

        df["Index"] = rng
        rs = Index(df["Index"])
        tm.assert_index_equal(rs, rng, check_names=False)
        assert rs.name == "Index"
        assert rng.name == "index"

        rs = df.reset_index().set_index("index")
        assert isinstance(rs.index, PeriodIndex)
        tm.assert_index_equal(rs.index, rng)

    def test_frame_index_to_string(self):
        index = PeriodIndex(["2011-1", "2011-2", "2011-3"], freq="M")
        frame = DataFrame(np.random.randn(3, 4), index=index)

        # it works!
        frame.to_string()
