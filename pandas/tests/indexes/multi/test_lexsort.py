import numpy as np
import pytest

import pandas as pd


class TestIsLexsorted:
    def test_is_lexsorted(self):
        levels = [[0, 1], [0, 1, 2]]

        index = pd.MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
        )
        assert index._is_lexsorted()

        index = pd.MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 2, 1]]
        )
        assert not index._is_lexsorted()

        index = pd.MultiIndex(
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

        index = pd.MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]], sortorder=2
        )
        assert index._lexsort_depth == 2

        index = pd.MultiIndex(
            levels=levels, codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 2, 1]], sortorder=1
        )
        assert index._lexsort_depth == 1

        index = pd.MultiIndex(
            levels=levels, codes=[[0, 0, 1, 0, 1, 1], [0, 1, 0, 2, 2, 1]], sortorder=0
        )
        assert index._lexsort_depth == 0


@pytest.mark.parametrize("na_position", ["first", "last"])
def test_lexsort_depth_sorted_with_nan(na_position):
    # GH#17931
    mi = pd.MultiIndex.from_arrays(
        [["a", None, "b", "a"], [np.nan, 1, 1, 2], [1, 2, None, 4]]
    )
    result = mi.sort_values(na_position=na_position)
    assert result._lexsort_depth == 3
