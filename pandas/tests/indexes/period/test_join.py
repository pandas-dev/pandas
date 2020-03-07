import numpy as np
import pytest

from pandas._libs.tslibs import IncompatibleFrequency

from pandas import Index, PeriodIndex, period_range
import pandas._testing as tm


class TestJoin:
    def test_joins(self, join_type):
        index = period_range("1/1/2000", "1/20/2000", freq="D")

        joined = index.join(index[:-5], how=join_type)

        assert isinstance(joined, PeriodIndex)
        assert joined.freq == index.freq

    def test_join_self(self, join_type):
        index = period_range("1/1/2000", "1/20/2000", freq="D")

        res = index.join(index, how=join_type)
        assert index is res

    def test_join_does_not_recur(self):
        df = tm.makeCustomDataframe(
            3,
            2,
            data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type="p",
            r_idx_type="dt",
        )
        s = df.iloc[:2, 0]

        res = s.index.join(df.columns, how="outer")
        expected = Index([s.index[0], s.index[1], df.columns[0], df.columns[1]], object)
        tm.assert_index_equal(res, expected)

    def test_join_mismatched_freq_raises(self):
        index = period_range("1/1/2000", "1/20/2000", freq="D")
        index3 = period_range("1/1/2000", "1/20/2000", freq="2D")
        with pytest.raises(IncompatibleFrequency):
            index.join(index3)
