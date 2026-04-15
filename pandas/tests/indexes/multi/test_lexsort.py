import numpy as np
import pytest

from pandas.errors import UnsortedIndexError

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
)


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

    def test_lexsort_depth_unsorted_levels(self):
        # GH#44380 sorted codes with unsorted levels should not
        # report a high lexsort depth
        index = MultiIndex(
            levels=[["a", "c", "b"], [1, 2]],
            codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
        )
        # Codes are sorted [0,0,1,1,2,2] but levels ["a","c","b"] are not,
        # so the actual values (a,a,c,c,b,b) are NOT sorted
        assert index._lexsort_depth == 0
        assert not index._is_lexsorted()

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
