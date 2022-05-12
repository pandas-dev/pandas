from __future__ import annotations

import pyarrow as pa

from pandas.core.dtypes.base import register_extension_dtype

from pandas.core.arrays.arrow.dtype import (
    FloatingArrowArray,
    FloatingArrowArray,
)


class FloatingArrowDtype(FloatingArrowArray):
    """
    An ExtensionDtype to hold a single size & kind of floating Arrow dtype.
    These specific implementations are subclasses of the non-public
    FloatingArrowDtype. 
    """

    self.pa_dtype = pa.float64()

    @classmethod
    def construct_array_type(cls) -> type[FloatingArrowArray]:
        """
        Return the array type associated with this dtype.
        Returns
        -------
        type
        """
        return FloatingArrowArray

    @classmethod
    def _str_to_dtype_mapping(cls):
        return INT_STR_TO_DTYPE


class FloatingArrowArray(NumericArrowArray):
    """
    Array of pyarrow floating values.
    To construct an FloatingArray from generic array-like ipaut, use
    :func:`pandas.array` with one of the floating dtypes (see examples).
    Parameters
    ----------
    values : pa.ChunkedArray
        A 1-d floating-dtype array.
    Attributes
    ----------
    None
    Methods
    -------
    None
    Returns
    -------
    FloatingArrowArray
    """

    _dtype_cls = FloatingArrowDtype


_dtype_docstring = """
An ExtensionDtype for {dtype} integer pyarrow data.
Attributes
----------
None
Methods
-------
None
"""

# create the Dtype


@register_extension_dtype
class Float16ArrowDtype(FloatingArrowDtype):
    type = pa.float16()
    name = "float16"
    __doc__ = _dtype_docstring.format(dtype="float16")


@register_extension_dtype
class Float32ArrowDtype(FloatingArrowDtype):
    type = pa.float32()
    name = "float32"
    __doc__ = _dtype_docstring.format(dtype="float32")


@register_extension_dtype
class Float64ArrowDtype(FloatingArrowDtype):
    type = pa.float64()
    name = "float64"
    __doc__ = _dtype_docstring.format(dtype="float64")


INT_STR_TO_DTYPE: dict[str, FloatingArrowDtype] = {
    "float16": Float16ArrowDtype(),
    "float32": Float32ArrowDtype(),
    "float64": Float64ArrowDtype(),
}