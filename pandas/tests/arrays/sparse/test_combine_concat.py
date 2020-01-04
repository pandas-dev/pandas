import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestSparseArrayConcat:
    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_basic(self, kind):
        a = pd.SparseArray([1, 0, 0, 2], kind=kind)
        b = pd.SparseArray([1, 0, 2, 2], kind=kind)

        result = pd.SparseArray._concat_same_type([a, b])
        # Can't make any assertions about the sparse index itself
        # since we aren't don't merge sparse blocs across arrays
        # in to_concat
        expected = np.array([1, 2, 1, 2, 2], dtype="int64")
        tm.assert_numpy_array_equal(result.sp_values, expected)
        assert result.kind == kind

    @pytest.mark.parametrize("kind", ["integer", "block"])
    def test_uses_first_kind(self, kind):
        other = "integer" if kind == "block" else "block"
        a = pd.SparseArray([1, 0, 0, 2], kind=kind)
        b = pd.SparseArray([1, 0, 2, 2], kind=other)

        result = pd.SparseArray._concat_same_type([a, b])
        expected = np.array([1, 2, 1, 2, 2], dtype="int64")
        tm.assert_numpy_array_equal(result.sp_values, expected)
        assert result.kind == kind
