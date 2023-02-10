import numpy as np

from pandas import (
    Categorical,
    DataFrame,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_replace_categorical_inplace_reference(using_copy_on_write):
    df = DataFrame({"a": Categorical([1, 2, 3])})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]  # noqa
    df.replace(to_replace=[1], value=2, inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a").codes, arr_a.codes)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
        tm.assert_frame_equal(view, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a").codes, arr_a.codes)


def test_replace_inplace_reference(using_copy_on_write):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]  # noqa
    df.replace(to_replace=[1.5], value=15.5, inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), arr_a)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
    else:
        assert np.shares_memory(get_array(df, "a"), arr_a)
