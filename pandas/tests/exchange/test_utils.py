import pandas as pd
import numpy as np
import pytest

from pandas.core.exchange.utils import dtype_to_arrow_c_fmt

# TODO: use ArrowSchema to get reference C-string.
# At the time, there is no way to access ArrowSchema holding a type format string from python.
# The only way to 'touch' it is to export the structure to a C-pointer:
# https://github.com/apache/arrow/blob/5680d209fd870f99134e2d7299b47acd90fabb8e/python/pyarrow/types.pxi#L230-L239
@pytest.mark.parametrize(
    "pandas_dtype, c_string",
    [
        (np.dtype("bool"), "b"),
        (np.dtype("int8"), "c"),
        (np.dtype("uint8"), "C"),
        (np.dtype("int16"), "s"),
        (np.dtype("uint16"), "S"),
        (np.dtype("int32"), "i"),
        (np.dtype("uint32"), "I"),
        (np.dtype("int64"), "l"),
        (np.dtype("uint64"), "L"),
        (np.dtype("float16"), "e"),
        (np.dtype("float32"), "f"),
        (np.dtype("float64"), "g"),
        (pd.Series(["a"]).dtype, "u"),
        (
            pd.Series([0]).astype("datetime64[ns]").dtype,
            "tsn:",
        ),
    ],
)
def test_dtype_to_arrow_c_fmt(pandas_dtype, c_string):  # noqa PR01
    """Test ``dtype_to_arrow_c_fmt`` utility function."""
    assert dtype_to_arrow_c_fmt(pandas_dtype) == c_string
