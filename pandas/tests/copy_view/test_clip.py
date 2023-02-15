import numpy as np

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_clip_inplace_reference(using_copy_on_write):
    df = DataFrame({"a": [1.5, 2, 3]})
    df_copy = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]
    df.clip(lower=2, inplace=True)

    # Clip not actually inplace right now but could be
    assert not np.shares_memory(get_array(df, "a"), arr_a)

    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
        tm.assert_frame_equal(df_copy, view)
