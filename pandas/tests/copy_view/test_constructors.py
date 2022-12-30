import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize("columns", [None, ["a"]])
@pytest.mark.parametrize("func", [lambda x: x, lambda x: x._mgr])
def test_dataframe_constructor_mgr(using_copy_on_write, func, columns):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()

    new_df = DataFrame(func(df))

    assert np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
    new_df.iloc[0] = 100

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
        tm.assert_frame_equal(df, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
        tm.assert_frame_equal(df, new_df)
