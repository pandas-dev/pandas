import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import TimedeltaArray


class TestTimedeltaArrayConstructor:
    def test_other_type_raises(self):
        msg = r"dtype bool cannot be converted to timedelta64\[ns\]"
        with pytest.raises(TypeError, match=msg):
            TimedeltaArray._from_sequence(np.array([1, 2, 3], dtype="bool"))

    def test_incorrect_dtype_raises(self):
        msg = "dtype 'category' is invalid, should be np.timedelta64 dtype"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence(
                np.array([1, 2, 3], dtype="i8"), dtype="category"
            )

        msg = "dtype 'int64' is invalid, should be np.timedelta64 dtype"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence(
                np.array([1, 2, 3], dtype="i8"), dtype=np.dtype("int64")
            )

        msg = r"dtype 'datetime64\[ns\]' is invalid, should be np.timedelta64 dtype"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence(
                np.array([1, 2, 3], dtype="i8"), dtype=np.dtype("M8[ns]")
            )

        msg = (
            r"dtype 'datetime64\[us, UTC\]' is invalid, should be np.timedelta64 dtype"
        )
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence(
                np.array([1, 2, 3], dtype="i8"), dtype="M8[us, UTC]"
            )

        msg = "Supported timedelta64 resolutions are 's', 'ms', 'us', 'ns'"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence(
                np.array([1, 2, 3], dtype="i8"), dtype=np.dtype("m8[Y]")
            )

    def test_copy(self):
        data = np.array([1, 2, 3], dtype="m8[ns]")
        arr = TimedeltaArray._from_sequence(data, copy=False)
        assert arr._ndarray is data

        arr = TimedeltaArray._from_sequence(data, copy=True)
        assert arr._ndarray is not data
        assert arr._ndarray.base is not data

    def test_from_sequence_dtype(self):
        msg = "dtype 'object' is invalid, should be np.timedelta64 dtype"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence([], dtype=object)


@pytest.mark.parametrize(
    "data",
    [
        np.array([1.5, np.nan, 90.0]),
        pd.array([1.5, None, 90.0], dtype="Float64"),
        [1.5, pd.NA, 90.0],
        [1.5, None, 90.0],
    ],
)
def test_from_sequence_numeric_honors_unit(data):
    # GH#63499 all-numeric data (including masked arrays and object input
    #  with missing values) is interpreted in the dtype's unit
    result = TimedeltaArray._from_sequence(data, dtype="m8[s]")
    expected = TimedeltaArray._from_sequence(np.array([1, "NaT", 90], dtype="m8[s]"))
    tm.assert_extension_array_equal(result, expected)
