import numpy as np
import pytest

from pandas.core.dtypes.cast import construct_1d_ndarray_preserving_na

import pandas._testing as tm


@pytest.mark.parametrize(
    "values, dtype, expected",
    [
        ([1, 2, 3], None, np.array([1, 2, 3])),
        (np.array([1, 2, 3]), None, np.array([1, 2, 3])),
        (["1", "2", None], None, np.array(["1", "2", None])),
        (["1", "2", None], np.dtype("str"), np.array(["1", "2", None])),
        ([1, 2, None], np.dtype("str"), np.array(["1", "2", None])),
    ],
)
def test_construct_1d_ndarray_preserving_na(values, dtype, expected):
    result = construct_1d_ndarray_preserving_na(values, dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)
