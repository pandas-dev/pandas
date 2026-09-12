import pytest

import pandas as pd
import pandas._testing as tm


class TestAstypeFromRangeIndex:
    def test_astype_interval_empty_rangeindex(self):
        # GH#68343
        rng = pd.RangeIndex(0)
        result = rng.astype("interval[int64, right]")
        expected = pd.IntervalIndex([])
        tm.assert_index_equal(result, expected)

    def test_astype_interval_nonempty_rangeindex_raises(self):
        # GH#68343: integer values are not intervals, matching Index([2, 3])
        rng = pd.RangeIndex(2, 4)
        with pytest.raises(TypeError, match="is not an interval"):
            rng.astype("interval[int64, right]")


class TestAstype:
    @pytest.mark.parametrize("ordered", [True, False])
    def test_astype_categorical_retains_ordered(self, ordered):
        index = pd.IntervalIndex.from_breaks(range(5))
        arr = index._data

        dtype = pd.CategoricalDtype(None, ordered=ordered)

        expected = pd.Categorical(list(arr), ordered=ordered)
        result = arr.astype(dtype)
        assert result.ordered is ordered
        tm.assert_categorical_equal(result, expected)

        # test IntervalIndex.astype while we're at it.
        result = index.astype(dtype)
        expected = pd.Index(expected)
        tm.assert_index_equal(result, expected)
