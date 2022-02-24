from pandas.core.dtypes.common import is_scalar

import pandas as pd


class TestSearchsorted:
    def test_searchsorted(self, string_dtype):
        arr = pd.array(["a", "b", "c"], dtype=string_dtype)

        result = arr.searchsorted("a", side="left")
        assert is_scalar(result)
        assert result == 0

        result = arr.searchsorted("a", side="right")
        assert is_scalar(result)
        assert result == 1
