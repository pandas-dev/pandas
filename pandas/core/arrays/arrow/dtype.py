from __future__ import annotations

import numpy as np
import pyarrow as pa

from pandas._typing import DtypeObj
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.base import StorageExtensionDtype

from pandas.core.arrays.arrow import ArrowExtensionArray


class ArrowDtype(StorageExtensionDtype):
    """
    Base class for dtypes for BaseArrowArray subclasses.
    Modeled after BaseMaskedDtype
    """

    name: str
    base = None
    type: pa.DataType

    na_value = pa.NA

    def __init__(self, storage="pyarrow") -> None:
        super().__init__(storage)

    @cache_readonly
    def numpy_dtype(self) -> np.dtype:
        """Return an instance of the related numpy dtype"""
        return self.type.to_pandas_dtype()

    @cache_readonly
    def kind(self) -> str:
        return self.numpy_dtype.kind

    @cache_readonly
    def itemsize(self) -> int:
        """Return the number of bytes in this dtype"""
        return self.numpy_dtype.itemsize

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return ArrowExtensionArray

    @classmethod
    def construct_from_string(cls, string: str):
        """
        Construct this type from a string.

        Parameters
        ----------
        string : str
        """
        if not isinstance(string, str):
            raise TypeError(
                f"'construct_from_string' expects a string, got {type(string)}"
            )
        if string == f"{cls.name}[pyarrow]":
            return cls(storage="pyarrow")
        raise TypeError(f"Cannot construct a '{cls.__name__}' from '{string}'")

    @classmethod
    def from_numpy_dtype(cls, dtype: np.dtype) -> ArrowDtype:
        """
        Construct the ArrowDtype corresponding to the given numpy dtype.
        """
        # TODO: This may be incomplete
        pa_dtype = pa.from_numpy_dtype(dtype)
        if pa_dtype is cls.type:
            return cls()
        raise NotImplementedError(dtype)

    def _get_common_dtype(self, dtypes: list[DtypeObj]) -> DtypeObj | None:
        # We unwrap any masked dtypes, find the common dtype we would use
        #  for that, then re-mask the result.
        from pandas.core.dtypes.cast import find_common_type

        new_dtype = find_common_type(
            [
                dtype.numpy_dtype if isinstance(dtype, ArrowDtype) else dtype
                for dtype in dtypes
            ]
        )
        if not isinstance(new_dtype, np.dtype):
            # If we ever support e.g. Masked[DatetimeArray] then this will change
            return None
        try:
            return type(self).from_numpy_dtype(new_dtype)
        except (KeyError, NotImplementedError):
            return None

    def __from_arrow__(self, array: pa.Array | pa.ChunkedArray):
        """
        Construct IntegerArray/FloatingArray from pyarrow Array/ChunkedArray.
        """
        array_class = self.construct_array_type()
        return array_class(array)
