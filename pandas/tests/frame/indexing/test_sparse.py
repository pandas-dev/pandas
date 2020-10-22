import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm
from pandas.arrays import SparseArray
from pandas.core.arrays.sparse import SparseDtype


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
    @pytest.mark.parametrize("dtype", [np.int64, np.float64, complex])
    @td.skip_if_no_scipy
    def test_locindexer_from_spmatrix(self, spmatrix_t, dtype):
        import scipy.sparse

        spmatrix_t = getattr(scipy.sparse, spmatrix_t)

        # The bug is triggered by a sparse matrix with purely sparse columns.  So the
        # recipe below generates a rectangular matrix of dimension (5, 7) where all the
        # diagonal cells are ones, meaning the last two columns are purely sparse.
        rows, cols = 5, 7
        spmatrix = spmatrix_t(np.eye(rows, cols, dtype=dtype), dtype=dtype)
        df = pd.DataFrame.sparse.from_spmatrix(spmatrix)

        # regression test for #34526
        itr_idx = range(2, rows)
        result = df.loc[itr_idx].values
        expected = spmatrix.toarray()[itr_idx]
        tm.assert_numpy_array_equal(result, expected)

        # regression test for #34540
        result = df.loc[itr_idx].dtypes.values
        expected = np.full(cols, SparseDtype(dtype, fill_value=0))
        tm.assert_numpy_array_equal(result, expected)

    def test_reindex(self):
        # https://github.com/pandas-dev/pandas/issues/35286
        df = pd.DataFrame(
            {"A": [0, 1], "B": pd.array([0, 1], dtype=pd.SparseDtype("int64", 0))}
        )
        result = df.reindex([0, 2])
        expected = pd.DataFrame(
            {
                "A": [0.0, np.nan],
                "B": pd.array([0.0, np.nan], dtype=pd.SparseDtype("float64", 0.0)),
            },
            index=[0, 2],
        )
        tm.assert_frame_equal(result, expected)

    def test_all_sparse(self):
        df = pd.DataFrame({"A": pd.array([0, 0], dtype=pd.SparseDtype("int64"))})
        result = df.loc[[0, 1]]
        tm.assert_frame_equal(result, df)
