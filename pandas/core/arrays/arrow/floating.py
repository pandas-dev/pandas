from __future__ import annotations

import pyarrow as pa

from pandas.core.dtypes.dtypes import register_extension_dtype

from pandas.core.arrays.arrow.numeric import FloatingArrowDtype

_dtype_docstring = """
An ExtensionDtype for {dtype} data.

This dtype uses ``pa.null`` as missing value indicator.

Attributes
----------
None

Methods
-------
None
"""


@register_extension_dtype
class Float16ArrowDtype(FloatingArrowDtype):
    name = "float16"
    type = pa.float16()
    __doc__ = _dtype_docstring.format(dtype="float16")
    _dtype_checker = pa.is_float16()


@register_extension_dtype
class Float32ArrowDtype(FloatingArrowDtype):
    name = "float32"
    type = pa.float32()
    __doc__ = _dtype_docstring.format(dtype="float32")
    _dtype_checker = pa.is_float32()


@register_extension_dtype
class Float64ArrowDtype(FloatingArrowDtype):
    name = "float64"
    type = pa.float64()
    __doc__ = _dtype_docstring.format(dtype="float64")
    _dtype_checker = pa.is_float64()


INT_STR_TO_DTYPE: dict[str, FloatingArrowDtype] = {
    "float16": Float16ArrowDtype(),
    "float32": Float32ArrowDtype(),
    "float64": Float64ArrowDtype(),
}
