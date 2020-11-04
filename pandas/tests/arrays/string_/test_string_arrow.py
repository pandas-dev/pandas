import re

import pytest

from pandas.core.arrays.string_arrow import ArrowStringArray

pa = pytest.importorskip("pyarrow", minversion="1.0.0")


@pytest.mark.parametrize("chunked", [True, False])
def test_constructor_not_string_type_raises(chunked):
    arr = pa.array([1, 2, 3])
    if chunked:
        arr = pa.chunked_array(arr)
    msg = re.escape(
        "ArrowStringArray requires a PyArrow (chunked) array of string type"
    )
    with pytest.raises(ValueError, match=msg):
        ArrowStringArray(arr)
