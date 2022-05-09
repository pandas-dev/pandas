from __future__ import annotations

import pyarrow as pa

from pandas.core.dtypes.base import register_extension_dtype

from pandas.core.arrays.arrow.numeric import (
    NumericArrowArray,
    NumericArrowDtype,
)


class IntegerArrowDtype(NumericArrowDtype):
    """
    An ExtensionDtype to hold a single size & kind of integer Arrow dtype.

    These specific implementations are subclasses of the non-public
    IntegerArrowDtype. For example we have Int8ArrowDtype to represent signed int 8s.

    The attributes name & type are set when these subclasses are created.
    """

    _default_pa_dtype = pa.int64()
    _dtype_checker = pa.types.is_integer

    @classmethod
    def construct_array_type(cls) -> type[IntegerArrowArray]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return IntegerArrowArray

    @classmethod
    def _str_to_dtype_mapping(cls):
        return INT_STR_TO_DTYPE


class IntegerArrowArray(NumericArrowArray):
    """
    Array of pyarrow integer values.

    To construct an IntegerArray from generic array-like ipaut, use
    :func:`pandas.array` with one of the integer dtypes (see examples).

    Parameters
    ----------
    values : pa.ChunkedArray
        A 1-d integer-dtype array.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Returns
    -------
    IntegerArrowArray
    """

    _dtype_cls = IntegerArrowDtype


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
class Int8ArrowDtype(IntegerArrowDtype):
    type = pa.int8()
    name = "int8"
    __doc__ = _dtype_docstring.format(dtype="int8")


@register_extension_dtype
class Int16ArrowDtype(IntegerArrowDtype):
    type = pa.int16()
    name = "int16"
    __doc__ = _dtype_docstring.format(dtype="int16")


@register_extension_dtype
class Int32ArrowDtype(IntegerArrowDtype):
    type = pa.int32()
    name = "int32"
    __doc__ = _dtype_docstring.format(dtype="int32")


@register_extension_dtype
class Int64ArrowDtype(IntegerArrowDtype):
    type = pa.int64()
    name = "int64"
    __doc__ = _dtype_docstring.format(dtype="int64")


@register_extension_dtype
class UInt8ArrowDtype(IntegerArrowDtype):
    type = pa.uint8()
    name = "uint8"
    __doc__ = _dtype_docstring.format(dtype="uint8")


@register_extension_dtype
class UInt16ArrowDtype(IntegerArrowDtype):
    type = pa.uint16()
    name = "uint16"
    __doc__ = _dtype_docstring.format(dtype="uint16")


@register_extension_dtype
class UInt32ArrowDtype(IntegerArrowDtype):
    type = pa.uint32()
    name = "uint32"
    __doc__ = _dtype_docstring.format(dtype="uint32")


@register_extension_dtype
class UInt64ArrowDtype(IntegerArrowDtype):
    type = pa.uint64()
    name = "uint64"
    __doc__ = _dtype_docstring.format(dtype="uint64")


INT_STR_TO_DTYPE: dict[str, IntegerArrowDtype] = {
    "int8": Int8ArrowDtype(),
    "int16": Int16ArrowDtype(),
    "int32": Int32ArrowDtype(),
    "int64": Int64ArrowDtype(),
    "uint8": UInt8ArrowDtype(),
    "uint16": UInt16ArrowDtype(),
    "uint32": UInt32ArrowDtype(),
    "uint64": UInt64ArrowDtype(),
}
