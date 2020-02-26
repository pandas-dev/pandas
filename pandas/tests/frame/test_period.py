import numpy as np

from pandas import DataFrame, Index, PeriodIndex, period_range
import pandas._testing as tm


def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


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

    def test_align_frame(self):
        rng = period_range("1/1/2000", "1/1/2010", freq="A")
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected.values[1::2] = np.nan
        tm.assert_frame_equal(result, expected)

        result = ts + _permute(ts[::2])
        tm.assert_frame_equal(result, expected)
