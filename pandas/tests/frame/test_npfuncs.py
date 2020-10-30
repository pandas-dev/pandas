"""
Tests for np.foo applied to DataFrame, not necessarily ufuncs.
"""
import numpy as np

from pandas import Categorical, DataFrame
import pandas._testing as tm


class TestAsArray:
    def test_asarray_homogenous(self):
        df = DataFrame({"A": Categorical([1, 2]), "B": Categorical([1, 2])})
        result = np.asarray(df)
        # may change from object in the future
        expected = np.array([[1, 1], [2, 2]], dtype="object")
        tm.assert_numpy_array_equal(result, expected)
