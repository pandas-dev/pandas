from __future__ import annotations

from typing import (
    Any,
    Callable,
    TypeVar,
)

import pyarrow as pa

from pandas.errors import AbstractMethodError
from pandas.util._decorators import cache_readonly

from pandas.core.arrays.arrow.array import ArrowExtensionArray
from pandas.core.arrays.arrow.dtype import ArrowDtype

T = TypeVar("T", bound="FloatingArrowArray")


class FloatingArrowDtype(ArrowDtype):
    _default_pa_dtype: pa.null()
    _dtype_checker: Callable[[Any], bool]  # pa.types.is_<type>

    @property
    def _is_numeric(self) -> bool:
        return False

    @property
    def _is_float(self) -> bool:
        return True

    @classmethod
    def _str_to_dtype_mapping(cls):
        raise AbstractMethodError(cls)


class FloatingArrowArray(ArrowExtensionArray):
    """
    Base class for Floating dtypes.
    """

    _dtype_cls: type[FloatingArrowDtype]

    def __init__(self, values: pa.ChunkedArray) -> None:
        checker = self._dtype_cls._dtype_checker
        if not (isinstance(values, pa.ChunkedArray) and checker(values.type)):
            descr = (
                "floating"
            )
            raise TypeError(f"values should be {descr} arrow array.")
        super().__init__(values)

    @cache_readonly
    def dtype(self) -> FloatingArrowDtype:
        mapping = self._dtype_cls._str_to_dtype_mapping()
        return mapping[str(self._data.type)]

    @classmethod
    def _from_sequence(cls, scalars, *, dtype=None, copy: bool = False):
        if dtype is None:
            dtype = cls._dtype_cls._default_pa_dtype
        return cls(pa.chunked_array([scalars], type=dtype.type))

    @classmethod
    def _from_sequence_of_strings(cls, strings, *, dtype=None, copy: bool = False):
        from pandas.core.tools.numeric import to_numeric

        scalars = to_numeric(strings, errors="raise")
        return cls._from_sequence(scalars, dtype=dtype, copy=copy)