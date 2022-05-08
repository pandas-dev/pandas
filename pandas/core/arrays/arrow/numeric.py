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

T = TypeVar("T", bound="NumericArrowArray")


class NumericArrowDtype(ArrowDtype):
    _default_pa_dtype: pa.null()
    _dtype_checker: Callable[[Any], bool]  # pa.types.is_<type>

    @property
    def _is_numeric(self) -> bool:
        return True

    @cache_readonly
    def is_signed_integer(self) -> bool:
        return self.kind == "i"

    @cache_readonly
    def is_unsigned_integer(self) -> bool:
        return self.kind == "u"

    @classmethod
    def _str_to_dtype_mapping(cls):
        raise AbstractMethodError(cls)


class NumericArrowArray(ArrowExtensionArray):
    """
    Base class for Integer and Floating and Boolean dtypes.
    """

    _dtype_cls: type[NumericArrowDtype]

    def __init__(self, values: pa.ChunkedArray) -> None:
        checker = self._dtype_cls._dtype_checker
        if not (isinstance(values, pa.ChunkedArray) and checker(values.type)):
            descr = (
                "floating"
                if self._dtype_cls.kind == "f"  # type: ignore[comparison-overlap]
                else "integer"
            )
            raise TypeError(f"values should be {descr} arrow array.")
        super().__init__(values)

    @cache_readonly
    def dtype(self) -> NumericArrowDtype:
        mapping = self._dtype_cls._str_to_dtype_mapping()
        return mapping[str(self._data.type)]

    @classmethod
    def _from_sequence(cls, scalars, *, dtype=None, copy: bool = False):
        return cls(pa.chunked_array(scalars, type=cls._dtype_cls._default_pa_dtype))

    @classmethod
    def _from_sequence_of_strings(cls, strings, *, dtype=None, copy: bool = False):
        from pandas.core.tools.numeric import to_numeric

        scalars = to_numeric(strings, errors="raise")
        return cls._from_sequence(scalars, dtype=dtype, copy=copy)
