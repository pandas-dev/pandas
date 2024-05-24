import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestRref:
    def test_rref_empty(self):
        df = pd.DataFrame([[]])
        msg = "Empty Dataframe has no row reduced echelon form"
        with pytest.raises(ValueError, match=msg):
            df.rref()

    def test_rref_invalid_values(self):
        df = pd.DataFrame([[1, 2], ["c", "d"]])
        msg = "The DataFrame must only contain numeric values"
        with pytest.raises(ValueError, match=msg):
            df.rref()

    def test_rref(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        result = df.rref().values
        expected = pd.DataFrame([[1, 0], [0, 1]]).values
        assert tm.shares_memory(result, df)
        assert np.array_equiv(result, expected)

    def test_copy_rref(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        result = df.rref(copy=True).values
        expected = pd.DataFrame([[1, 0], [0, 1]]).values
        assert not tm.shares_memory(result, df)
        assert np.array_equiv(result, expected)

    def test_rref_with_different_primitive_dtype(self):
        df = pd.DataFrame([[1, 2], [3, 4]], dtype=float)
        result = df.rref().values
        expected = pd.DataFrame([[1, 0], [0, 1]], dtype=float).values
        assert tm.shares_memory(result, df)
        assert np.array_equiv(result, expected)

    def test_rref_with_non_square_df(self):
        df = pd.DataFrame([[1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6]])
        result = df.rref().values
        expected = pd.DataFrame([[1, 0, -1], [0, 1, 2], [0, 0, 0], [0, 0, 0]]).values
        assert tm.shares_memory(result, df)
        assert np.array_equiv(result, expected)
