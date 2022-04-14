from __future__ import annotations

from typing import (
    Any,
    Callable,
)

import numpy as np
import pyarrow as pa

from pandas.errors import AbstractMethodError
from pandas.util._decorators import cache_readonly

from pandas.core.arrays.arrow.base import BaseArrowDtype
from pandas.core.arrays.masked import BaseMaskedArray


class NumericArrowDtype(BaseArrowDtype):

    _default_pa_dtype: pa.DataType
    _checker: Callable[[Any], bool]  # is_foo_dtype

    @cache_readonly
    def is_signed_integer(self) -> bool:
        return pa.types.is_signed_integer(self.type)

    @cache_readonly
    def is_unsigned_integer(self) -> bool:
        return pa.types.is_unsigned_integer(self.type)

    @property
    def _is_numeric(self) -> bool:
        return True

    def __from_arrow__(self, array: pa.Array | pa.ChunkedArray) -> BaseMaskedArray:
        """
        Construct IntegerArray/FloatingArray from pyarrow Array/ChunkedArray.
        """
        array_class = self.construct_array_type()
        return array_class(array)

    @classmethod
    def _str_to_dtype_mapping(cls):
        raise AbstractMethodError(cls)

    # TODO: Might not need this for pyarrow. Only used in _coerce_to_data_and_mask
    #  which is easily retrievable from pa.ChunkedArray
    @classmethod
    def _standardize_dtype(cls, dtype) -> NumericArrowDtype:
        """
        Convert a string representation or a pyarrow dtype to NumericArrowDtype.
        """
        if not issubclass(type(dtype), cls):
            mapping = cls._str_to_dtype_mapping()
            try:
                dtype = mapping[str(pa.type_for_alias(dtype))]
            except KeyError as err:
                raise ValueError(f"invalid dtype specified {dtype}") from err
        return dtype

    @classmethod
    def _safe_cast(cls, values: np.ndarray, dtype: np.dtype, copy: bool) -> np.ndarray:
        """
        Safely cast the values to the given dtype.

        "safe" in this context means the casting is lossless.
        """
        raise AbstractMethodError(cls)
