"""
Tests for np.foo applied to DataFrame, not necessarily ufuncs.
"""
import numpy as np

from pandas import (
    Categorical,
    DataFrame,
)
import pandas._testing as tm


class TestAsArray:
    def test_asarray_homogenous(self):
        df = DataFrame({"A": Categorical([1, 2]), "B": Categorical([1, 2])})
        result = np.asarray(df)
        # may change from object in the future
        expected = np.array([[1, 1], [2, 2]], dtype="object")
        tm.assert_numpy_array_equal(result, expected)

    def test_np_sqrt(self, float_frame):
        with np.errstate(all="ignore"):
            result = np.sqrt(float_frame)
        assert isinstance(result, type(float_frame))
        assert result.index is float_frame.index
        assert result.columns is float_frame.columns

        tm.assert_frame_equal(result, float_frame.apply(np.sqrt))

    def test_np_ravel(self):
        # GH#26247 np.ravel() fails on list of DataFrames with column names
        x = np.zeros((10, 3))
        df = [DataFrame(batch.reshape(1, 3), columns=["1", "2", "3"]) for batch in x]
        result = np.ravel(df)
        expected = np.ravel([DataFrame(batch.reshape(1, 3)) for batch in x])
        assert all([a == b for a, b in zip(result, expected)])
