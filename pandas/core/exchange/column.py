from functools import cached_property
from typing import (
    Any,
    Tuple,
)

import numpy as np

import pandas as pd
from pandas.api.types import (
    is_categorical_dtype,
    is_string_dtype,
)
from pandas.core.exchange.buffer import PandasBuffer
from pandas.core.exchange.dataframe_protocol import (
    Buffer,
    Column,
    ColumnNullType,
    DtypeKind,
)
from pandas.core.exchange.utils import (
    ArrowCTypes,
    Endianness,
    dtype_to_arrow_c_fmt,
)

_NP_KINDS = {
    "i": DtypeKind.INT,
    "u": DtypeKind.UINT,
    "f": DtypeKind.FLOAT,
    "b": DtypeKind.BOOL,
    "U": DtypeKind.STRING,
    "M": DtypeKind.DATETIME,
    "m": DtypeKind.DATETIME,
}

_NULL_DESCRIPTION = {
    DtypeKind.FLOAT: (ColumnNullType.USE_NAN, None),
    DtypeKind.DATETIME: (ColumnNullType.USE_NAN, None),
    DtypeKind.INT: (ColumnNullType.NON_NULLABLE, None),
    DtypeKind.UINT: (ColumnNullType.NON_NULLABLE, None),
    DtypeKind.BOOL: (ColumnNullType.NON_NULLABLE, None),
    # Null values for categoricals are stored as `-1` sentinel values
    # in the category date (e.g., `col.values.codes` is int8 np.ndarray)
    DtypeKind.CATEGORICAL: (ColumnNullType.USE_SENTINEL, -1),
    # follow Arrow in using 1 as valid value and 0 for missing/null value
    DtypeKind.STRING: (ColumnNullType.USE_BYTEMASK, 0),
}


class PandasColumn(Column):
    """
    A column object, with only the methods and properties required by the
    interchange protocol defined.
    A column can contain one or more chunks. Each chunk can contain up to three
    buffers - a data buffer, a mask buffer (depending on null representation),
    and an offsets buffer (if variable-size binary; e.g., variable-length
    strings).
    Note: this Column object can only be produced by ``__dataframe__``, so
          doesn't need its own version or ``__column__`` protocol.
    """

    def __init__(self, column: pd.Series, allow_copy: bool = True) -> None:
        """
        Note: doesn't deal with extension arrays yet, just assume a regular
        Series/ndarray for now.
        """
        if not isinstance(column, pd.Series):
            raise NotImplementedError(f"Columns of type {type(column)} not handled yet")

        # Store the column as a private attribute
        self._col = column
        self._allow_copy = allow_copy

    @property
    def size(self) -> int:
        """
        Size of the column, in elements.
        """
        return self._col.size

    @property
    def offset(self) -> int:
        """
        Offset of first element. Always zero.
        """
        # TODO: chunks are implemented now, probably this should return something
        return 0

    @cached_property
    def dtype(self):
        dtype = self._col.dtype

        if is_categorical_dtype(dtype):
            codes = self._col.values.codes
            (
                _,
                bitwidth,
                c_arrow_dtype_f_str,
                _,
            ) = self._dtype_from_pandasdtype(codes.dtype)
            return (
                DtypeKind.CATEGORICAL,
                bitwidth,
                c_arrow_dtype_f_str,
                "=",
            )
        elif is_string_dtype(dtype):
            return (DtypeKind.STRING, 8, dtype_to_arrow_c_fmt(dtype), "=")
        else:
            return self._dtype_from_pandasdtype(dtype)

    def _dtype_from_pandasdtype(self, dtype) -> Tuple[DtypeKind, int, str, str]:
        """
        See `self.dtype` for details.
        """
        # Note: 'c' (complex) not handled yet (not in array spec v1).
        #       'b', 'B' (bytes), 'S', 'a', (old-style string) 'V' (void) not handled
        #       datetime and timedelta both map to datetime (is timedelta handled?)

        kind = _NP_KINDS.get(dtype.kind, None)
        if kind is None:
            # Not a NumPy dtype. Check if it's a categorical maybe
            raise ValueError(f"Data type {dtype} not supported by exchange protocol")

        return (kind, dtype.itemsize * 8, dtype_to_arrow_c_fmt(dtype), dtype.byteorder)

    @property
    def describe_categorical(self):
        """
        If the dtype is categorical, there are two options:
        - There are only values in the data buffer.
        - There is a separate dictionary-style encoding for categorical values.
        Raises RuntimeError if the dtype is not categorical
        Content of returned dict:
            - "is_ordered" : bool, whether the ordering of dictionary indices is
                             semantically meaningful.
            - "is_dictionary" : bool, whether a dictionary-style mapping of
                                categorical values to other objects exists
            - "mapping" : dict, Python-level only (e.g. ``{int: str}``).
                          None if not a dictionary-style categorical.
        """
        if not self.dtype[0] == DtypeKind.CATEGORICAL:
            raise TypeError(
                "`describe_categorical only works on a column with "
                "categorical dtype!"
            )

        return {
            "is_ordered": self._col.cat.ordered,
            "is_dictionary": True,
            "mapping": dict(enumerate(self._col.cat.categories)),
        }

    @property
    def describe_null(self):
        kind = self.dtype[0]
        try:
            null, value = _NULL_DESCRIPTION[kind]
        except KeyError:
            raise NotImplementedError(f"Data type {kind} not yet supported")

        return null, value

    @cached_property
    def null_count(self) -> int:
        """
        Number of null elements. Should always be known.
        """
        return self._col.isna().sum()

    @property
    def metadata(self):
        """
        Store specific metadata of the column.
        """
        return {"index": self._col.index}

    def num_chunks(self) -> int:
        """
        Return the number of chunks the column consists of.
        """
        return 1

    def get_chunks(self, n_chunks=None):
        """
        Return an iterator yielding the chunks.
        See `DataFrame.get_chunks` for details on ``n_chunks``.
        """
        if n_chunks and n_chunks > 1:
            size = len(self._col)
            step = size // n_chunks
            if size % n_chunks != 0:
                step += 1
            for start in range(0, step * n_chunks, step):
                yield PandasColumn(
                    self._col.iloc[start : start + step], self._allow_copy
                )
        else:
            yield self

    def get_buffers(self):
        """
        Return a dictionary containing the underlying buffers.
        The returned dictionary has the following contents:
            - "data": a two-element tuple whose first element is a buffer
                      containing the data and whose second element is the data
                      buffer's associated dtype.
            - "validity": a two-element tuple whose first element is a buffer
                          containing mask values indicating missing data and
                          whose second element is the mask value buffer's
                          associated dtype. None if the null representation is
                          not a bit or byte mask.
            - "offsets": a two-element tuple whose first element is a buffer
                         containing the offset values for variable-size binary
                         data (e.g., variable-length strings) and whose second
                         element is the offsets buffer's associated dtype. None
                         if the data buffer does not have an associated offsets
                         buffer.
        """
        buffers = {}
        buffers["data"] = self._get_data_buffer()
        try:
            buffers["validity"] = self._get_validity_buffer()
        except:
            buffers["validity"] = None

        try:
            buffers["offsets"] = self._get_offsets_buffer()
        except:
            buffers["offsets"] = None

        return buffers

    def _get_data_buffer(
        self,
    ) -> Tuple[PandasBuffer, Any]:  # Any is for self.dtype tuple
        """
        Return the buffer containing the data and the buffer's associated dtype.
        """
        if self.dtype[0] in (
            DtypeKind.INT,
            DtypeKind.UINT,
            DtypeKind.FLOAT,
            DtypeKind.BOOL,
        ):
            buffer = PandasBuffer(self._col.to_numpy(), allow_copy=self._allow_copy)
            dtype = self.dtype
        elif self.dtype[0] == DtypeKind.CATEGORICAL:
            codes = self._col.values.codes
            buffer = PandasBuffer(codes, allow_copy=self._allow_copy)
            dtype = self._dtype_from_pandasdtype(codes.dtype)
        elif self.dtype[0] == DtypeKind.STRING:
            # Marshal the strings from a NumPy object array into a byte array
            buf = self._col.to_numpy()
            b = bytearray()

            # TODO: this for-loop is slow; can be implemented in Cython/C/C++ later
            for obj in buf:
                if isinstance(obj, str):
                    b.extend(obj.encode(encoding="utf-8"))

            # Convert the byte array to a Pandas "buffer" using a NumPy array as the backing store
            buffer = PandasBuffer(np.frombuffer(b, dtype="uint8"))

            # Define the dtype for the returned buffer
            dtype = (
                DtypeKind.STRING,
                8,
                ArrowCTypes.STRING,
                Endianness.NATIVE,
            )  # note: currently only support native endianness
        else:
            raise NotImplementedError(f"Data type {self._col.dtype} not handled yet")

        return buffer, dtype

    def _get_validity_buffer(self) -> Tuple[PandasBuffer, Any]:
        """
        Return the buffer containing the mask values indicating missing data and
        the buffer's associated dtype.
        Raises RuntimeError if null representation is not a bit or byte mask.
        """
        null, invalid = self.describe_null

        if self.dtype[0] == DtypeKind.STRING:
            # For now, have the mask array be comprised of bytes, rather than a bit array
            buf = self._col.to_numpy()
            mask = []

            # Determine the encoding for valid values
            valid = 1 if invalid == 0 else 0

            mask = [valid if isinstance(obj, str) else invalid for obj in buf]

            # Convert the mask array to a Pandas "buffer" using a NumPy array as the backing store
            buffer = PandasBuffer(np.asarray(mask, dtype="uint8"))

            # Define the dtype of the returned buffer
            dtype = (DtypeKind.UINT, 8, ArrowCTypes.UINT8, Endianness.NATIVE)

            return buffer, dtype

        if null == 0:
            msg = "This column is non-nullable so does not have a mask"
        elif null == 1:
            msg = "This column uses NaN as null so does not have a separate mask"
        else:
            raise NotImplementedError("See self.describe_null")

        raise RuntimeError(msg)

    def _get_offsets_buffer(self) -> Tuple[PandasBuffer, Any]:
        """
        Return the buffer containing the offset values for variable-size binary
        data (e.g., variable-length strings) and the buffer's associated dtype.
        Raises RuntimeError if the data buffer does not have an associated
        offsets buffer.
        """
        if self.dtype[0] == DtypeKind.STRING:
            # For each string, we need to manually determine the next offset
            values = self._col.to_numpy()
            ptr = 0
            offsets = [ptr] + [None] * len(values)
            for i, v in enumerate(values):
                # For missing values (in this case, `np.nan` values), we don't increment the pointer)
                if isinstance(v, str):
                    b = v.encode(encoding="utf-8")
                    ptr += len(b)

                offsets[i + 1] = ptr

            # Convert the list of offsets to a NumPy array of signed 64-bit integers (note: Arrow allows the offsets array to be either `int32` or `int64`; here, we default to the latter)
            buf = np.asarray(offsets, dtype="int64")

            # Convert the offsets to a Pandas "buffer" using the NumPy array as the backing store
            buffer = PandasBuffer(buf)

            # Assemble the buffer dtype info
            dtype = (
                DtypeKind.INT,
                64,
                ArrowCTypes.INT64,
                Endianness.NATIVE,
            )  # note: currently only support native endianness
        else:
            raise RuntimeError(
                "This column has a fixed-length dtype so does not have an offsets buffer"
            )

        return buffer, dtype
