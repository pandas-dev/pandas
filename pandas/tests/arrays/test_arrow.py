import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.mark.parametrize(
    "pa_type", [pa.string(), pa.large_string(), pa.binary(), pa.large_binary()]
)
@pytest.mark.parametrize("box", [list, tuple])
def test_string_sequence_with_nan(pa_type, box, using_nan_is_na):
    # GH#64578
    is_binary = pa.types.is_binary(pa_type) or pa.types.is_large_binary(pa_type)
    text = b"nan" if is_binary else "nan"
    values = box([text, np.nan])
    dtype = pd.ArrowDtype(pa_type)

    if not using_nan_is_na:
        msg = "Expected bytes, got a 'float' object"
        with pytest.raises(pa.ArrowTypeError, match=msg):
            pd.array(values, dtype=dtype)
        return

    result = pd.array(values, dtype=dtype)
    expected = pd.array([text, None], dtype=dtype)
    tm.assert_extension_array_equal(result, expected)
