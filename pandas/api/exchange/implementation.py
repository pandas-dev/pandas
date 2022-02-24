import collections.abc
import ctypes

from typing import Tuple, Any

from .dataframe_protocol import Buffer, Column, DataFrame as DataFrameXchg, DtypeKind, DlpackDeviceType

import pandas as pd
import numpy as np


def from_dataframe(df : DataFrameXchg,
                   allow_copy : bool = True) -> pd.DataFrame:
    """
    Construct a pandas DataFrame from ``df`` if it supports ``__dataframe__``
    """
    if isinstance(df, pd.DataFrame):
        return df

    if not hasattr(df, '__dataframe__'):
        raise ValueError("`df` does not support __dataframe__")

    return _from_dataframe(df.__dataframe__(allow_copy=allow_copy))


def _from_dataframe(df : DataFrameXchg) -> pd.DataFrame:
    """
    Note: not all cases are handled yet, only ones that can be implemented with
    only Pandas. Later, we need to implement/test support for categoricals,
    bit/byte masks, chunk handling, etc.
    """
    _buffers = []  # hold on to buffers, keeps memory alive
    result = []
    for chunk in df.get_chunks():
        # We need a dict of columns here, with each column being a numpy array (at
        # least for now, deal with non-numpy dtypes later).
        chunk_cols = {}
        for name in chunk.column_names():
            if not isinstance(name, str):
                raise ValueError(f"Column {name} is not a string")
            if name in chunk_cols:
                raise ValueError(f"Column {name} is not unique")
            col = chunk.get_column_by_name(name)
            if col.dtype[0] in (DtypeKind.INT, DtypeKind.UINT, DtypeKind.FLOAT, DtypeKind.BOOL):
                # Simple numerical or bool dtype, turn into numpy array
                chunk_cols[name], _buf = convert_column_to_ndarray(col)
            elif col.dtype[0] == DtypeKind.CATEGORICAL:
                chunk_cols[name], _buf = convert_categorical_column(col)
            elif col.dtype[0] == DtypeKind.STRING:
                chunk_cols[name], _buf = convert_string_column(col)
            else:
                raise NotImplementedError(f"Data type {col.dtype[0]} not handled yet")

            _buffers.append(_buf)

        df_new = pd.DataFrame(chunk_cols)
        result.append(df_new)

    df_new = pd.concat(result)
    df_new._buffers = _buffers
    return df_new


def convert_column_to_ndarray(col : Column) -> Tuple[np.ndarray, Buffer]:
    """
    Convert an int, uint, float or bool column to a numpy array.
    """
    if col.offset != 0:
        raise NotImplementedError("column.offset > 0 not handled yet")

    if col.describe_null[0] not in (0, 1):
        raise NotImplementedError("Null values represented as masks or "
                                  "sentinel values not handled yet")

    _buffer, _dtype = col.get_buffers()["data"]
    return buffer_to_ndarray(_buffer, _dtype), _buffer


def buffer_to_ndarray(_buffer : Buffer, _dtype) -> np.ndarray:
    # Handle the dtype
    kind = _dtype[0]
    bitwidth = _dtype[1]
    if _dtype[0] not in (DtypeKind.INT, DtypeKind.UINT, DtypeKind.FLOAT, DtypeKind.BOOL):
        raise RuntimeError("Not a boolean, integer or floating-point dtype")

    _ints = {8: np.int8, 16: np.int16, 32: np.int32, 64: np.int64}
    _uints = {8: np.uint8, 16: np.uint16, 32: np.uint32, 64: np.uint64}
    _floats = {32: np.float32, 64: np.float64}
    _np_dtypes = {0: _ints, 1: _uints, 2: _floats, 20: {8: bool}}
    column_dtype = _np_dtypes[kind][bitwidth]

    # No DLPack yet, so need to construct a new ndarray from the data pointer
    # and size in the buffer plus the dtype on the column
    ctypes_type = np.ctypeslib.as_ctypes_type(column_dtype)
    data_pointer = ctypes.cast(_buffer.ptr, ctypes.POINTER(ctypes_type))

    # NOTE: `x` does not own its memory, so the caller of this function must
    #       either make a copy or hold on to a reference of the column or
    #       buffer! (not done yet, this is pretty awful ...)
    x = np.ctypeslib.as_array(data_pointer,
                              shape=(_buffer.bufsize // (bitwidth//8),))

    return x


def convert_categorical_column(col : Column) -> Tuple[pd.Series, Buffer]:
    """
    Convert a categorical column to a Series instance.
    """
    ordered, is_dict, mapping = col.describe_categorical
    if not is_dict:
        raise NotImplementedError('Non-dictionary categoricals not supported yet')

    # If you want to cheat for testing (can't use `_col` in real-world code):
    #    categories = col._col.values.categories.values
    #    codes = col._col.values.codes
    categories = np.asarray(list(mapping.values()))
    codes_buffer, codes_dtype = col.get_buffers()["data"]
    codes = buffer_to_ndarray(codes_buffer, codes_dtype)
    values = categories[codes]

    # Seems like Pandas can only construct with non-null values, so need to
    # null out the nulls later
    cat = pd.Categorical(values, categories=categories, ordered=ordered)
    series = pd.Series(cat)
    null_kind = col.describe_null[0]
    if null_kind == 2:  # sentinel value
        sentinel = col.describe_null[1]
        series[codes == sentinel] = np.nan
    else:
        raise NotImplementedError("Only categorical columns with sentinel "
                                  "value supported at the moment")

    return series, codes_buffer


def convert_string_column(col : Column) -> Tuple[np.ndarray, dict]:
    """
    Convert a string column to a NumPy array.
    """
    # Retrieve the data buffers
    buffers = col.get_buffers()

    # Retrieve the data buffer containing the UTF-8 code units
    dbuffer, bdtype = buffers["data"]

    # Retrieve the offsets buffer containing the index offsets demarcating the beginning and end of each string
    obuffer, odtype = buffers["offsets"]

    # Retrieve the mask buffer indicating the presence of missing values
    mbuffer, mdtype = buffers["validity"]

    # Retrieve the missing value encoding
    null_kind, null_value = col.describe_null

    # Convert the buffers to NumPy arrays
    dt = (DtypeKind.UINT, 8, None, None)  # note: in order to go from STRING to an equivalent ndarray, we claim that the buffer is uint8 (i.e., a byte array)
    dbuf = buffer_to_ndarray(dbuffer, dt)

    obuf = buffer_to_ndarray(obuffer, odtype)
    mbuf = buffer_to_ndarray(mbuffer, mdtype)

    # Assemble the strings from the code units
    str_list = []
    for i in range(obuf.size-1):
        # Check for missing values
        if null_kind == 3:  # bit mask
            v = mbuf[i/8]
            if null_value == 1:
                v = ~v

            if v & (1<<(i%8)):
                str_list.append(np.nan)
                continue

        elif null_kind == 4 and mbuf[i] == null_value:  # byte mask
            str_list.append(np.nan)
            continue

        # Extract a range of code units
        units = dbuf[obuf[i]:obuf[i+1]]

        # Convert the list of code units to bytes
        b = bytes(units)

        # Create the string
        s = b.decode(encoding="utf-8")

        # Add to our list of strings
        str_list.append(s)

    # Convert the string list to a NumPy array
    return np.asarray(str_list, dtype="object"), buffers

# Implementation of interchange protocol
# --------------------------------------

class _PandasBuffer(Buffer):
    """
    Data in the buffer is guaranteed to be contiguous in memory.
    """

    def __init__(self, x : np.ndarray, allow_copy : bool = True) -> None:
        """
        Handle only regular columns (= numpy arrays) for now.
        """
        if not x.strides == (x.dtype.itemsize,):
            # The protocol does not support strided buffers, so a copy is
            # necessary. If that's not allowed, we need to raise an exception.
            if allow_copy:
                x = x.copy()
            else:
                raise RuntimeError("Exports cannot be zero-copy in the case "
                                   "of a non-contiguous buffer")

        # Store the numpy array in which the data resides as a private
        # attribute, so we can use it to retrieve the public attributes
        self._x = x

    @property
    def bufsize(self) -> int:
        """
        Buffer size in bytes.
        """
        return self._x.size * self._x.dtype.itemsize

    @property
    def ptr(self) -> int:
        """
        Pointer to start of the buffer as an integer.
        """
        return self._x.__array_interface__['data'][0]

    def __dlpack__(self):
        """
        DLPack not implemented in NumPy yet, so leave it out here.
        """
        raise NotImplementedError("__dlpack__")

    def __dlpack_device__(self) -> Tuple[DlpackDeviceType, int]:
        """
        Device type and device ID for where the data in the buffer resides.
        """
        return (DlpackDeviceType.CPU, None)

    def __repr__(self) -> str:
        return 'PandasBuffer(' + str({'bufsize': self.bufsize,
                                      'ptr': self.ptr,
                                      'device': self.__dlpack_device__()[0].name}
                                      ) + ')'

class _PandasColumn(Column):
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

    def __init__(self, column : pd.Series,
                 allow_copy : bool = True) -> None:
        """
        Note: doesn't deal with extension arrays yet, just assume a regular
        Series/ndarray for now.
        """
        if not isinstance(column, pd.Series):
            raise NotImplementedError("Columns of type {} not handled "
                                      "yet".format(type(column)))

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
        return 0

    @property
    def dtype(self):
        dtype = self._col.dtype

        # For now, assume that, if the column dtype is 'O' (i.e., `object`), then we have an array of strings
        if not isinstance(dtype, pd.CategoricalDtype) and dtype.kind == 'O':
            return (DtypeKind.STRING, 8, 'u', '=')

        return self._dtype_from_pandasdtype(dtype)

    def _dtype_from_pandasdtype(self, dtype) -> Tuple[DtypeKind, int, str, str]:
        """
        See `self.dtype` for details.
        """
        # Note: 'c' (complex) not handled yet (not in array spec v1).
        #       'b', 'B' (bytes), 'S', 'a', (old-style string) 'V' (void) not handled
        #       datetime and timedelta both map to datetime (is timedelta handled?)
        _np_kinds = {"i": DtypeKind.INT, "u": DtypeKind.UINT, "f": DtypeKind.FLOAT, "b": DtypeKind.BOOL,
                     "U": DtypeKind.STRING,
                     "M": DtypeKind.DATETIME, "m": DtypeKind.DATETIME}
        kind = _np_kinds.get(dtype.kind, None)
        if kind is None:
            # Not a NumPy dtype. Check if it's a categorical maybe
            if isinstance(dtype, pd.CategoricalDtype):
                kind = DtypeKind.CATEGORICAL
            else:
                raise ValueError(f"Data type {dtype} not supported by exchange"
                                 "protocol")

        if kind not in (DtypeKind.INT, DtypeKind.UINT, DtypeKind.FLOAT, DtypeKind.BOOL, DtypeKind.CATEGORICAL, DtypeKind.STRING):
            raise NotImplementedError(f"Data type {dtype} not handled yet")

        bitwidth = dtype.itemsize * 8
        format_str = dtype.str
        endianness = dtype.byteorder if not kind == DtypeKind.CATEGORICAL else '='
        return (kind, bitwidth, format_str, endianness)


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
            raise TypeError("`describe_categorical only works on a column with "
                            "categorical dtype!")

        ordered = self._col.dtype.ordered
        is_dictionary = True
        # NOTE: this shows the children approach is better, transforming
        # `categories` to a "mapping" dict is inefficient
        codes = self._col.values.codes  # ndarray, length `self.size`
        # categories.values is ndarray of length n_categories
        categories = self._col.values.categories.values
        mapping = {ix: val for ix, val in enumerate(categories)}
        return ordered, is_dictionary, mapping

    @property
    def describe_null(self):
        kind = self.dtype[0]
        value = None
        if kind == DtypeKind.FLOAT:
            null = 1  # np.nan
        elif kind == DtypeKind.DATETIME:
            null = 1  # np.datetime64('NaT')
        elif kind in (DtypeKind.INT, DtypeKind.UINT, DtypeKind.BOOL):
            # TODO: check if extension dtypes are used once support for them is
            #       implemented in this protocol code
            null = 0  # integer and boolean dtypes are non-nullable
        elif kind == DtypeKind.CATEGORICAL:
            # Null values for categoricals are stored as `-1` sentinel values
            # in the category date (e.g., `col.values.codes` is int8 np.ndarray)
            null = 2
            value = -1
        elif kind == DtypeKind.STRING:
            null = 4
            value = 0  # follow Arrow in using 1 as valid value and 0 for missing/null value
        else:
            raise NotImplementedError(f"Data type {self.dtype} not yet supported")

        return null, value

    @property
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
        return {}

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
                step +=1
            for start in range(0, step * n_chunks, step):
                yield _PandasColumn(self._col.iloc[start:start + step], self._allow_copy)
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

    def _get_data_buffer(self) -> Tuple[_PandasBuffer, Any]:  # Any is for self.dtype tuple
        """
        Return the buffer containing the data and the buffer's associated dtype.
        """
        if self.dtype[0] in (DtypeKind.INT, DtypeKind.UINT, DtypeKind.FLOAT, DtypeKind.BOOL):
            buffer = _PandasBuffer(
                self._col.to_numpy(), allow_copy=self._allow_copy)
            dtype = self.dtype
        elif self.dtype[0] == DtypeKind.CATEGORICAL:
            codes = self._col.values.codes
            buffer = _PandasBuffer(
                codes, allow_copy=self._allow_copy)
            dtype = self._dtype_from_pandasdtype(codes.dtype)
        elif self.dtype[0] == DtypeKind.STRING:
            # Marshal the strings from a NumPy object array into a byte array
            buf = self._col.to_numpy()
            b = bytearray()

            # TODO: this for-loop is slow; can be implemented in Cython/C/C++ later
            for i in range(buf.size):
                if type(buf[i]) == str:
                    b.extend(buf[i].encode(encoding="utf-8"))

            # Convert the byte array to a Pandas "buffer" using a NumPy array as the backing store
            buffer = _PandasBuffer(np.frombuffer(b, dtype="uint8"))

            # Define the dtype for the returned buffer
            dtype = (DtypeKind.STRING, 8, "u", "=")  # note: currently only support native endianness
        else:
            raise NotImplementedError(f"Data type {self._col.dtype} not handled yet")

        return buffer, dtype

    def _get_validity_buffer(self) -> Tuple[_PandasBuffer, Any]:
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
            if invalid == 0:
                valid = 1
            else:
                valid = 0

            for i in range(buf.size):
                if type(buf[i]) == str:
                    v = valid
                else:
                    v = invalid

                mask.append(v)

            # Convert the mask array to a Pandas "buffer" using a NumPy array as the backing store
            buffer = _PandasBuffer(np.asarray(mask, dtype="uint8"))

            # Define the dtype of the returned buffer
            dtype = (DtypeKind.UINT, 8, "C", "=")

            return buffer, dtype

        if null == 0:
            msg = "This column is non-nullable so does not have a mask"
        elif null == 1:
            msg = "This column uses NaN as null so does not have a separate mask"
        else:
            raise NotImplementedError("See self.describe_null")

        raise RuntimeError(msg)

    def _get_offsets_buffer(self) -> Tuple[_PandasBuffer, Any]:
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
            offsets = [ptr]
            for v in values:
                # For missing values (in this case, `np.nan` values), we don't increment the pointer)
                if type(v) == str:
                    b = v.encode(encoding="utf-8")
                    ptr += len(b)

                offsets.append(ptr)

            # Convert the list of offsets to a NumPy array of signed 64-bit integers (note: Arrow allows the offsets array to be either `int32` or `int64`; here, we default to the latter)
            buf = np.asarray(offsets, dtype="int64")

            # Convert the offsets to a Pandas "buffer" using the NumPy array as the backing store
            buffer = _PandasBuffer(buf)

            # Assemble the buffer dtype info
            dtype = (DtypeKind.INT, 64, 'l', "=")  # note: currently only support native endianness
        else:
            raise RuntimeError("This column has a fixed-length dtype so does not have an offsets buffer")

        return buffer, dtype


class _PandasDataFrameXchg(DataFrameXchg):
    """
    A data frame class, with only the methods required by the interchange
    protocol defined.
    Instances of this (private) class are returned from
    ``pd.DataFrame.__dataframe__`` as objects with the methods and
    attributes defined on this class.
    """
    def __init__(self, df : pd.DataFrame, nan_as_null : bool = False,
                 allow_copy : bool = True) -> None:
        """
        Constructor - an instance of this (private) class is returned from
        `pd.DataFrame.__dataframe__`.
        """
        self._df = df
        # ``nan_as_null`` is a keyword intended for the consumer to tell the
        # producer to overwrite null values in the data with ``NaN`` (or ``NaT``).
        # This currently has no effect; once support for nullable extension
        # dtypes is added, this value should be propagated to columns.
        self._nan_as_null = nan_as_null
        self._allow_copy = allow_copy

    @property
    def metadata(self):
        # `index` isn't a regular column, and the protocol doesn't support row
        # labels - so we export it as Pandas-specific metadata here.
        return {"pandas.index": self._df.index}

    def num_columns(self) -> int:
        return len(self._df.columns)

    def num_rows(self) -> int:
        return len(self._df)

    def num_chunks(self) -> int:
        return 1

    def column_names(self):
        return self._df.columns.tolist()

    def get_column(self, i: int) -> _PandasColumn:
        return _PandasColumn(
            self._df.iloc[:, i], allow_copy=self._allow_copy)

    def get_column_by_name(self, name: str) -> _PandasColumn:
        return _PandasColumn(
            self._df[name], allow_copy=self._allow_copy)

    def get_columns(self):
        return [_PandasColumn(self._df[name], allow_copy=self._allow_copy)
                for name in self._df.columns]

    def select_columns(self, indices):
        if not isinstance(indices, collections.abc.Sequence):
            raise ValueError("`indices` is not a sequence")
        if not isinstance(indices, list):
            indices = list(indices)

        return _PandasDataFrameXchg(self._df.iloc[:, indices], self._nan_as_null, self._allow_copy)

    def select_columns_by_name(self, names):
        if not isinstance(names, collections.abc.Sequence):
            raise ValueError("`names` is not a sequence")
        if not isinstance(names, list):
            names = list(names)

        return _PandasDataFrameXchg(self._df.loc[:, names], self._nan_as_null, self._allow_copy)

    def get_chunks(self, n_chunks=None):
        """
        Return an iterator yielding the chunks.
        """
        if n_chunks and n_chunks > 1:
            size = len(self._df)
            step = size // n_chunks
            if size % n_chunks != 0:
                step +=1
            for start in range(0, step * n_chunks, step):
                yield _PandasDataFrameXchg(self._df.iloc[start:start + step, :], self._nan_as_null, self._allow_copy)
        else:
            yield self
