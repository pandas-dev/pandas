import numpy as np

from pandas import Series
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_convert_dtypes_infer_objects(using_copy_on_write):
    ser = Series(["a", "b", "c"])
    ser_orig = ser.copy()
    result = ser.convert_dtypes(
        convert_integer=False,
        convert_boolean=False,
        convert_floating=False,
        convert_string=False,
    )

    if using_copy_on_write:
        assert np.shares_memory(get_array(ser), get_array(result))
    else:
        assert not np.shares_memory(get_array(ser), get_array(result))

    result.iloc[0] = "x"
    tm.assert_series_equal(ser, ser_orig)
