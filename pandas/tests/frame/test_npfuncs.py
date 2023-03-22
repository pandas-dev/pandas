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

    def test_sum_deprecated_axis_behavior(self):
        # GH#52042 deprecated behavior of df.sum(axis=None), which gets
        #  called when we do np.sum(df)

        arr = np.random.randn(4, 3)
        df = DataFrame(arr)

        msg = "The behavior of DataFrame.sum with axis=None is deprecated"
        with tm.assert_produces_warning(
            FutureWarning, match=msg, check_stacklevel=False
        ):
            res = np.sum(df)

        with tm.assert_produces_warning(FutureWarning, match=msg):
            expected = df.sum(axis=None)
        tm.assert_series_equal(res, expected)
