import re

import numpy as np
import pytest

from pandas.core.arrays.string_arrow import ArrowStringArray

pa = pytest.importorskip("pyarrow", minversion="1.0.0")


@pytest.mark.parametrize("chunked", [True, False])
@pytest.mark.parametrize("array", [np, pa])
def test_constructor_not_string_type_raises(array, chunked):
    arr = array.array([1, 2, 3])
    if chunked:
        if array is np:
            pytest.skip("chunked not applicable to numpy array")
        arr = pa.chunked_array(arr)
    if array is np:
        msg = "Unsupported type '<class 'numpy.ndarray'>' for ArrowStringArray"
    else:
        msg = re.escape(
            "ArrowStringArray requires a PyArrow (chunked) array of string type"
        )
    with pytest.raises(ValueError, match=msg):
        ArrowStringArray(arr)
