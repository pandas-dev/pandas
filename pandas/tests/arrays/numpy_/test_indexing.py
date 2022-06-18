import numpy as np

from pandas.core.dtypes.common import is_scalar

import pandas as pd
import pandas._testing as tm


class TestSearchsorted:
    def test_searchsorted_numeric_dtypes_scalar(self, any_real_numpy_dtype):
        arr = pd.array([1, 3, 90], dtype=any_real_numpy_dtype)
        result = arr.searchsorted(30)
        assert is_scalar(result)
        assert result == 2

        result = arr.searchsorted([30])
        expected = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_numeric_dtypes_vector(self, any_real_numpy_dtype):
        arr = pd.array([1, 3, 90], dtype=any_real_numpy_dtype)
        result = arr.searchsorted([2, 30])
        expected = np.array([1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_sorter(self, any_real_numpy_dtype):
        arr = pd.array([3, 1, 2], dtype=any_real_numpy_dtype)
        result = arr.searchsorted([0, 3], sorter=np.argsort(arr))
        expected = np.array([0, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
