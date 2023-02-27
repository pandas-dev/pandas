import numpy as np

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


class TestDowncast:
    def test_downcast(self, using_copy_on_write):
        df = DataFrame({"a": [1.0, 2.0], "b": 1.5})
        df_orig = df.copy()
        result = df.downcast()

        assert not np.shares_memory(get_array(df, "a"), get_array(result, "a"))
        if using_copy_on_write:
            assert np.shares_memory(get_array(df, "b"), get_array(result, "b"))
        else:
            assert not np.shares_memory(get_array(df, "b"), get_array(result, "b"))

        result.iloc[0, 1] = 100.5
        tm.assert_frame_equal(df, df_orig)
