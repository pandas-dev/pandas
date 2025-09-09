import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


def test_divmod_bool_raises(box_with_array):
    # GH#46043 // raises, so divmod should too
    ser = Series([True, False])
    obj = tm.box_expected(ser, box_with_array)

    msg = "operator 'floordiv' not implemented for bool dtypes"
    with pytest.raises(NotImplementedError, match=msg):
        obj // obj

    if box_with_array is DataFrame:
        msg = "operator 'floordiv' not implemented for bool dtypes"
    else:
        msg = "operator 'divmod' not implemented for bool dtypes"
    with pytest.raises(NotImplementedError, match=msg):
        divmod(obj, obj)

    # go through __rdivmod__
    with pytest.raises(NotImplementedError, match=msg):
        divmod(True, obj)
