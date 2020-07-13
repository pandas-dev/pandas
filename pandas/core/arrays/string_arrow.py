from collections.abc import Iterable
from typing import Tuple, Type, Union

import numpy as np
import pyarrow as pa

from pandas._libs import missing as libmissing

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import register_extension_dtype

import pandas as pd
from pandas.api.types import is_integer
from pandas.core.arrays.base import ExtensionArray
from pandas.core.indexers import check_array_indexer


@register_extension_dtype
class ArrowStringDtype(ExtensionDtype):
    """
    Extension dtype for string data in a ``pyarrow.ChunkedArray``.

    .. versionadded:: 1.1.0

    .. warning::

       ArrowStringDtype is considered experimental. The implementation and
       parts of the API may change without warning.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Examples
    --------
    >>> pd.ArrowStringDtype()
    ArrowStringDtype
    """

    name = "arrow_string"

    #: StringDtype.na_value uses pandas.NA
    na_value = libmissing.NA

    @property
    def type(self) -> Type[str]:
        return str

    @classmethod
    def construct_array_type(cls) -> Type["ArrowStringArray"]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return ArrowStringArray

    def __hash__(self) -> int:
        return hash("ArrowStringDtype")

    def __repr__(self) -> str:
        return "ArrowStringDtype"

    def __from_arrow__(
        self, array: Union["pa.Array", "pa.ChunkedArray"]
    ) -> "ArrowStringArray":
        """
        Construct StringArray from pyarrow Array/ChunkedArray.
        """
        return ArrowStringArray(array)

    def __eq__(self, other) -> bool:
        """Check whether 'other' is equal to self.

        By default, 'other' is considered equal if
        * it's a string matching 'self.name'.
        * it's an instance of this type.

        Parameters
        ----------
        other : Any

        Returns
        -------
        bool
        """
        if isinstance(other, ArrowStringDtype):
            return True
        elif isinstance(other, str) and other == "arrow_string":
            return True
        else:
            return False


class ArrowStringArray(ExtensionArray):
    """
    Extension array for string data in a ``pyarrow.ChunkedArray``.

    .. versionadded:: 1.1.0

    .. warning::

       ArrowStringArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : pyarrow.Array or pyarrow.ChunkedArray
        The array of data.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    array
        The recommended function for creating a ArrowStringArray.
    Series.str
        The string methods are available on Series backed by
        a ArrowStringArray.

    Notes
    -----
    ArrowStringArray returns a BooleanArray for comparison methods.

    Examples
    --------
    >>> pd.array(['This is', 'some text', None, 'data.'], dtype="arrow_string")
    <ArrowStringArray>
    ['This is', 'some text', <NA>, 'data.']
    Length: 4, dtype: arrow_string
    """

    def __init__(self, values):
        if isinstance(values, pa.Array):
            self.data = pa.chunked_array([values])
        elif isinstance(values, pa.ChunkedArray):
            self.data = values
        else:
            raise ValueError(f"Unsupported type '{type(values)}' for ArrowStringArray")

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        # TODO(ARROW-9407): Accept pd.NA in Arrow
        scalars_corrected = [None if pd.isna(x) else x for x in scalars]
        return cls(pa.array(scalars_corrected, type=pa.string()))

    @property
    def dtype(self) -> ArrowStringDtype:
        """
        An instance of 'ArrowStringDtype'.
        """
        return ArrowStringDtype()

    def __array__(self, *args, **kwargs) -> "np.ndarray":
        """Correctly construct numpy arrays when passed to `np.asarray()`."""
        return self.data.__array__(*args, **kwargs)

    def __arrow_array__(self, type=None):
        """Convert myself to a pyarrow Array or ChunkedArray."""
        return self.data

    @property
    def size(self) -> int:
        """
        Return the number of elements in this array.

        Returns
        -------
        size : int
        """
        return len(self.data)

    @property
    def shape(self) -> Tuple[int]:
        """Return the shape of the data."""
        # This may be patched by pandas to support pseudo-2D operations.
        return (len(self.data),)

    @property
    def ndim(self) -> int:
        """Return the number of dimensions of the underlying data."""
        return 1

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self.data)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

    #     def _values_for_factorize(self):
    #         arr = self._ndarray.copy()
    #         mask = self.isna()
    #         arr[mask] = -1
    #         return arr, -1

    def __getitem__(self, item):
        # type (Any) -> Any
        """Select a subset of self.

        Parameters
        ----------
        item : int, slice, or ndarray
            * int: The position in 'self' to get.
            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None
            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : scalar or ExtensionArray

        Notes
        -----
        For scalar ``item``, return a scalar value suitable for the array's
        type. This should be an instance of ``self.dtype.type``.
        For slice ``key``, return an instance of ``ExtensionArray``, even
        if the slice is length 0 or 1.
        For a boolean mask, return an instance of ``ExtensionArray``, filtered
        to the values where ``item`` is True.
        """
        item = check_array_indexer(self, item)

        if isinstance(item, Iterable):
            raise NotImplementedError("Iterable")
            # if not is_array_like(item):
            #     item = np.array(item)
            # if is_integer_dtype(item) or (len(item) == 0):
            #     return self.take(item)
            # elif is_bool_dtype(item):
            #     indices = np.array(item)
            #     indices = np.argwhere(indices).flatten()
            #     return self.take(indices)
            # else:
            #     raise IndexError(
            #         """Only integers, slices and integer or
            #            boolean arrays are valid indices."""
            #     )
        elif is_integer(item):
            if item < 0:
                item += len(self)
            if item >= len(self):
                return None

        value = self.data[item]
        if isinstance(value, pa.ChunkedArray):
            return type(self)(value)
        else:
            return value.as_py()

    def __setitem__(self, key, value):
        raise NotImplementedError("__setitem__")

    def fillna(self, value=None, method=None, limit=None):
        raise NotImplementedError("fillna")

    #     def astype(self, dtype, copy=True):
    #         dtype = pandas_dtype(dtype)
    #         if isinstance(dtype, StringDtype):
    #             if copy:
    #                 return self.copy()
    #             return self
    #         elif isinstance(dtype, _IntegerDtype):
    #             arr = self._ndarray.copy()
    #             mask = self.isna()
    #             arr[mask] = 0
    #             values = arr.astype(dtype.numpy_dtype)
    #             return IntegerArray(values, mask, copy=False)

    #         return super().astype(dtype, copy)

    def _reduce(self, name, skipna=True, **kwargs):
        if name in ["min", "max"]:
            return getattr(self, name)(skipna=skipna)

        raise TypeError(f"Cannot perform reduction '{name}' with string dtype")

    #     def value_counts(self, dropna=False):
    #         from pandas import value_counts

    #         return value_counts(self._ndarray, dropna=dropna).astype("Int64")

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        size = 0
        for chunk in self.data.chunks:
            for buf in chunk.buffers():
                if buf is not None:
                    size += buf.size
        return size

    def isna(self) -> np.ndarray:
        """
        Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        return self.data.is_null()

    def copy(self):
        # type: () -> ExtensionArray
        """
        Return a copy of the array.

        Parameters
        ----------
        deep : bool, default False
            Also copy the underlying data backing this array.

        Returns
        -------
        ExtensionArray
        """
        return type(self)(self.data)
