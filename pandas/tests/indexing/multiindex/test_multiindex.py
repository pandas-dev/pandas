import numpy as np
import pytest

import pandas._libs.index as _index
from pandas.errors import PerformanceWarning

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series
import pandas.util.testing as tm


class TestMultiIndexBasic:
    def test_multiindex_perf_warn(self):

        df = DataFrame(
            {
                "jim": [0, 0, 1, 1],
                "joe": ["x", "x", "z", "y"],
                "jolie": np.random.rand(4),
            }
        ).set_index(["jim", "joe"])

        with tm.assert_produces_warning(PerformanceWarning):
            df.loc[(1, "z")]

        df = df.iloc[[2, 1, 3, 0]]
        with tm.assert_produces_warning(PerformanceWarning):
            df.loc[(0,)]

    def test_multiindex_contains_dropped(self):
        # GH 19027
        # test that dropped MultiIndex levels are not in the MultiIndex
        # despite continuing to be in the MultiIndex's levels
        idx = MultiIndex.from_product([[1, 2], [3, 4]])
        assert 2 in idx
        idx = idx.drop(2)

        # drop implementation keeps 2 in the levels
        assert 2 in idx.levels[0]
        # but it should no longer be in the index itself
        assert 2 not in idx

        # also applies to strings
        idx = MultiIndex.from_product([["a", "b"], ["c", "d"]])
        assert "a" in idx
        idx = idx.drop("a")
        assert "a" in idx.levels[0]
        assert "a" not in idx

    @pytest.mark.parametrize(
        "data, expected",
        [
            (MultiIndex.from_product([(), ()]), True),
            (MultiIndex.from_product([(1, 2), (3, 4)]), True),
            (MultiIndex.from_product([("a", "b"), (1, 2)]), False),
        ],
    )
    def test_multiindex_is_homogeneous_type(self, data, expected):
        assert data._is_homogeneous_type is expected

    def test_indexing_over_hashtable_size_cutoff(self):
        n = 10000

        old_cutoff = _index._SIZE_CUTOFF
        _index._SIZE_CUTOFF = 20000

        s = Series(np.arange(n), MultiIndex.from_arrays((["a"] * n, np.arange(n))))

        # hai it works!
        assert s[("a", 5)] == 5
        assert s[("a", 6)] == 6
        assert s[("a", 7)] == 7

        _index._SIZE_CUTOFF = old_cutoff

    def test_multi_nan_indexing(self):

        # GH 3588
        df = DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            }
        )
        result = df.set_index(["a", "b"], drop=False)
        expected = DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            },
            index=[
                Index(["R1", "R2", np.nan, "R4"], name="a"),
                Index(["C1", "C2", "C3", "C4"], name="b"),
            ],
        )
        tm.assert_frame_equal(result, expected)

    def test_contains(self):
        # GH 24570
        tx = pd.timedelta_range("09:30:00", "16:00:00", freq="30 min")
        idx = MultiIndex.from_arrays([tx, np.arange(len(tx))])
        assert tx[0] in idx
        assert "element_not_exit" not in idx
        assert "0 day 09:30:00" in idx
