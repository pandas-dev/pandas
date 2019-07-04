import pytest

from pandas import DataFrame
import pandas.util.testing as tm


class TestIndexingSlow:
    @pytest.mark.slow
    def test_large_dataframe_indexing(self):
        # GH10692
        result = DataFrame({"x": range(10 ** 6)}, dtype="int64")
        result.loc[len(result)] = len(result) + 1
        expected = DataFrame({"x": range(10 ** 6 + 1)}, dtype="int64")
        tm.assert_frame_equal(result, expected)
