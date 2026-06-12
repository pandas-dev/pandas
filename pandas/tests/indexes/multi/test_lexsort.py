import numpy as np
import pytest

from pandas.errors import UnsortedIndexError

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
)
import pandas._testing as tm


class TestIsLexsorted:
    def test_is_lexsorted(self):
        levels = [[0, 1], [0, 1, 2]]

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
        )
        assert index._is_lexsorted()

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 2, 1]]
        )
        assert not index._is_lexsorted()

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 1, 0, 1, 1], [0, 1, 0, 2, 2, 1]]
        )
        assert not index._is_lexsorted()
        assert index._lexsort_depth == 0


class TestLexsortDepth:
    def test_lexsort_depth(self):
        # Test that lexsort_depth return the correct sortorder
        # when it was given to the MultiIndex const.
        # GH#28518

        levels = [[0, 1], [0, 1, 2]]

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]], sortorder=2
        )
        assert index._lexsort_depth == 2

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 2, 1]], sortorder=1
        )
        assert index._lexsort_depth == 1

        index = MultiIndex(
            levels=levels, codes=[[0, 0, 1, 0, 1, 1], [0, 1, 0, 2, 2, 1]], sortorder=0
        )
        assert index._lexsort_depth == 0

    def test_value_lexsort_depth_unsorted_levels(self):
        # GH#44380 sorted codes with unsorted levels should not
        # report a high value lexsort depth
        index = MultiIndex(
            levels=[["a", "c", "b"], [1, 2]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )
        # Codes are sorted [0,0,1,1,2,2] but levels ["a","c","b"] are not,
        # so the actual values (a,a,c,c,b,b) are NOT sorted
        assert index._lexsort_depth == 2
        assert index._value_lexsort_depth == 0

    def test_value_lexsort_depth_unused_level_entries(self):
        # GH#44380 unsorted entries that do not appear in the codes do not
        # affect the ordering of the values
        index = MultiIndex(
            levels=[["a", "c", "b"], [1, 2]],
            codes=[[0, 0, 2, 2], [0, 1, 0, 1]],
        )
        # values are (a,1),(a,2),(b,1),(b,2); the unused "c" entry is
        # irrelevant
        assert index.is_monotonic_increasing
        assert index._value_lexsort_depth == 2

        ser = Series(np.arange(4), index=index)
        result = ser.loc["a":"b"]
        tm.assert_series_equal(result, ser)

    def test_unsorted_levels_xs_slice_raises(self):
        # GH#44380 slices routed through get_loc_level (e.g. DataFrame.xs)
        # are gated like get_locs
        index = MultiIndex(
            levels=[["a", "c", "b"], [1, 2]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )
        df = DataFrame({"x": np.arange(6)}, index=index)
        with pytest.raises(UnsortedIndexError, match="lexsort depth"):
            df.xs(slice("a", "c"))
        with pytest.raises(UnsortedIndexError, match="lexsort depth"):
            index.get_loc_level(slice("a", "c"))

    def test_unsorted_levels_get_loc_uses_codes(self):
        # GH#44380 exact-label get_loc only needs grouped codes, not sorted
        # values, so it should return a slice without a PerformanceWarning
        index = MultiIndex(
            levels=[["a", "c", "b"], [1, 2]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )
        with tm.assert_produces_warning(None):
            assert index.get_loc("c") == slice(2, 4)

    def test_unsorted_levels_slice_raises(self):
        # GH#44380 slicing on a MultiIndex with unsorted levels and
        # sorted codes should raise UnsortedIndexError
        df = DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]}, index=["a", "c", "b"])
        stacked = df.stack()
        with pytest.raises(UnsortedIndexError, match="lexsort depth"):
            stacked.loc["a":"c"]

    def test_unsorted_levels_unstack_stack_roundtrip(self):
        # GH#44380 after unstack/stack/reindex cycle, slicing should
        # still raise UnsortedIndexError
        index = MultiIndex.from_product([["a", "c", "b"], [1, 2]])
        data = Series(np.arange(6), index=index)

        # Initial unsorted index correctly raises
        with pytest.raises(UnsortedIndexError, match="lexsort depth"):
            data.loc["a":"c"]

        # After unstack/stack cycle, should also raise
        data2 = data.unstack().reindex(["a", "c", "b"]).stack()
        with pytest.raises(UnsortedIndexError, match="lexsort depth"):
            data2.loc["a":"c"]
