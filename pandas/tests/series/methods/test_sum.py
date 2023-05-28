import numpy as np

from pandas import Series
import pandas._testing as tm


class TestSeriesSum:
    def test_sum_uint64(self):
        # GH 53401
        s = Series([10000000000000000000], dtype="uint64")
        result = s.sum()
        expected = np.uint64(10000000000000000000)
        tm.assert_almost_equal(result, expected)
