import numpy as np
import pandas as pd

import pandas.util._test_decorators as td
import pandas._testing as tm

from pandas.arrays import SparseArray
from pandas.core.arrays.sparse import SparseDtype

import pytest


class TestSparseDataFrameIndexing:
    def test_getitem_sparse_column(self):
        # https://github.com/pandas-dev/pandas/issues/23559
        data = SparseArray([0, 1])
        df = pd.DataFrame({"A": data})
        expected = pd.Series(data, name="A")
        result = df["A"]
        tm.assert_series_equal(result, expected)

        result = df.iloc[:, 0]
        tm.assert_series_equal(result, expected)

        result = df.loc[:, "A"]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("spmatrix_t", ["coo_matrix", "csc_matrix", "csr_matrix"])
    @td.skip_if_no_scipy
    def test_locindexer_from_spmatrix(self, spmatrix_t):
        import scipy.sparse

        spmatrix_t = getattr(scipy.sparse, spmatrix_t)

        spmatrix = spmatrix_t([[1.0, 0.0], [0.0, 0.0]], dtype=np.float64)
        df = pd.DataFrame.sparse.from_spmatrix(spmatrix)

        # regression test for #34526
        itr_idx = [1]
        result = df.loc[itr_idx].values
        expected = spmatrix.toarray()[itr_idx]
        tm.assert_numpy_array_equal(result, expected)

        # regression test for #34540
        result = df.loc[itr_idx].dtypes.values
        expected = np.full(2, SparseDtype(np.float64, fill_value=0))
        tm.assert_numpy_array_equal(result, expected)
