import re

import numpy as np
import pytest

from pandas.compat import pa_version_under1p0

from pandas.core.arrays.string_arrow import (
    ArrowStringArray,
    ArrowStringDtype,
)


@pytest.mark.skipif(
    pa_version_under1p0,
    reason="pyarrow>=1.0.0 is required for PyArrow backed StringArray",
)
@pytest.mark.parametrize("chunked", [True, False])
@pytest.mark.parametrize("array", ["numpy", "pyarrow"])
def test_constructor_not_string_type_raises(array, chunked):
    import pyarrow as pa

    array = pa if array == "pyarrow" else np

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


@pytest.mark.skipif(
    not pa_version_under1p0,
    reason="pyarrow is installed",
)
def test_pyarrow_not_installed_raises():
    msg = re.escape("pyarrow>=1.0.0 is required for PyArrow backed StringArray")

    with pytest.raises(ImportError, match=msg):
        ArrowStringDtype()

    with pytest.raises(ImportError, match=msg):
        ArrowStringArray([])

    with pytest.raises(ImportError, match=msg):
        ArrowStringArray._from_sequence(["a", None, "b"])
