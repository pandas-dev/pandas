"""Coders for strings."""

from __future__ import annotations

from functools import partial

import numpy as np

from xarray.coding.variables import (
    VariableCoder,
    lazy_elemwise_func,
    pop_to,
    safe_setitem,
    unpack_for_decoding,
    unpack_for_encoding,
)
from xarray.core import indexing
from xarray.core.utils import module_available
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array

HAS_NUMPY_2_0 = module_available("numpy", minversion="2.0.0.dev0")


def create_vlen_dtype(element_type):
    if element_type not in (str, bytes):
        raise TypeError(f"unsupported type for vlen_dtype: {element_type!r}")
    # based on h5py.special_dtype
    return np.dtype("O", metadata={"element_type": element_type})


def check_vlen_dtype(dtype):
    if dtype.kind != "O" or dtype.metadata is None:
        return None
    else:
        # check xarray (element_type) as well as h5py (vlen)
        return dtype.metadata.get("element_type", dtype.metadata.get("vlen"))


def is_unicode_dtype(dtype):
    return dtype.kind == "U" or check_vlen_dtype(dtype) is str


def is_bytes_dtype(dtype):
    return dtype.kind == "S" or check_vlen_dtype(dtype) is bytes


class EncodedStringCoder(VariableCoder):
    """Transforms between unicode strings and fixed-width UTF-8 bytes."""

    def __init__(self, allows_unicode=True):
        self.allows_unicode = allows_unicode

    def encode(self, variable: Variable, name=None) -> Variable:
        dims, data, attrs, encoding = unpack_for_encoding(variable)

        contains_unicode = is_unicode_dtype(data.dtype)
        encode_as_char = encoding.get("dtype") == "S1"
        if encode_as_char:
            del encoding["dtype"]  # no longer relevant

        if contains_unicode and (encode_as_char or not self.allows_unicode):
            if "_FillValue" in attrs:
                raise NotImplementedError(
                    f"variable {name!r} has a _FillValue specified, but "
                    "_FillValue is not yet supported on unicode strings: "
                    "https://github.com/pydata/xarray/issues/1647"
                )

            string_encoding = encoding.pop("_Encoding", "utf-8")
            safe_setitem(attrs, "_Encoding", string_encoding, name=name)
            # TODO: figure out how to handle this in a lazy way with dask
            data = encode_string_array(data, string_encoding)

            return Variable(dims, data, attrs, encoding)
        else:
            variable.encoding = encoding
            return variable

    def decode(self, variable: Variable, name=None) -> Variable:
        dims, data, attrs, encoding = unpack_for_decoding(variable)

        if "_Encoding" in attrs:
            string_encoding = pop_to(attrs, encoding, "_Encoding")
            func = partial(decode_bytes_array, encoding=string_encoding)
            data = lazy_elemwise_func(data, func, np.dtype(object))

        return Variable(dims, data, attrs, encoding)


def decode_bytes_array(bytes_array, encoding="utf-8"):
    # This is faster than using np.char.decode() or np.vectorize()
    bytes_array = np.asarray(bytes_array)
    decoded = [x.decode(encoding) for x in bytes_array.ravel()]
    return np.array(decoded, dtype=object).reshape(bytes_array.shape)


def encode_string_array(string_array, encoding="utf-8"):
    string_array = np.asarray(string_array)
    encoded = [x.encode(encoding) for x in string_array.ravel()]
    return np.array(encoded, dtype=bytes).reshape(string_array.shape)


def ensure_fixed_length_bytes(var: Variable) -> Variable:
    """Ensure that a variable with vlen bytes is converted to fixed width."""
    if check_vlen_dtype(var.dtype) is bytes:
        dims, data, attrs, encoding = unpack_for_encoding(var)
        # TODO: figure out how to handle this with dask
        data = np.asarray(data, dtype=np.bytes_)
        return Variable(dims, data, attrs, encoding)
    else:
        return var


class CharacterArrayCoder(VariableCoder):
    """Transforms between arrays containing bytes and character arrays."""

    def encode(self, variable, name=None):
        variable = ensure_fixed_length_bytes(variable)

        dims, data, attrs, encoding = unpack_for_encoding(variable)
        if data.dtype.kind == "S" and encoding.get("dtype") is not str:
            data = bytes_to_char(data)
            if "char_dim_name" in encoding.keys():
                char_dim_name = encoding.pop("char_dim_name")
            else:
                char_dim_name = f"string{data.shape[-1]}"
            dims = dims + (char_dim_name,)
        return Variable(dims, data, attrs, encoding)

    def decode(self, variable, name=None):
        dims, data, attrs, encoding = unpack_for_decoding(variable)

        if data.dtype == "S1" and dims:
            encoding["char_dim_name"] = dims[-1]
            dims = dims[:-1]
            data = char_to_bytes(data)
        return Variable(dims, data, attrs, encoding)


def bytes_to_char(arr):
    """Convert numpy/dask arrays from fixed width bytes to characters."""
    if arr.dtype.kind != "S":
        raise ValueError("argument must have a fixed-width bytes dtype")

    if is_chunked_array(arr):
        chunkmanager = get_chunked_array_type(arr)

        return chunkmanager.map_blocks(
            _numpy_bytes_to_char,
            arr,
            dtype="S1",
            chunks=arr.chunks + ((arr.dtype.itemsize,)),
            new_axis=[arr.ndim],
        )
    return _numpy_bytes_to_char(arr)


def _numpy_bytes_to_char(arr):
    """Like netCDF4.stringtochar, but faster and more flexible."""
    # adapt handling of copy-kwarg to numpy 2.0
    # see https://github.com/numpy/numpy/issues/25916
    # and https://github.com/numpy/numpy/pull/25922
    copy = None if HAS_NUMPY_2_0 else False
    # ensure the array is contiguous
    arr = np.array(arr, copy=copy, order="C", dtype=np.bytes_)
    return arr.reshape(arr.shape + (1,)).view("S1")


def char_to_bytes(arr):
    """Convert numpy/dask arrays from characters to fixed width bytes."""
    if arr.dtype != "S1":
        raise ValueError("argument must have dtype='S1'")

    if not arr.ndim:
        # no dimension to concatenate along
        return arr

    size = arr.shape[-1]

    if not size:
        # can't make an S0 dtype
        return np.zeros(arr.shape[:-1], dtype=np.bytes_)

    if is_chunked_array(arr):
        chunkmanager = get_chunked_array_type(arr)

        if len(arr.chunks[-1]) > 1:
            raise ValueError(
                "cannot stacked dask character array with "
                f"multiple chunks in the last dimension: {arr}"
            )

        dtype = np.dtype("S" + str(arr.shape[-1]))
        return chunkmanager.map_blocks(
            _numpy_char_to_bytes,
            arr,
            dtype=dtype,
            chunks=arr.chunks[:-1],
            drop_axis=[arr.ndim - 1],
        )
    else:
        return StackedBytesArray(arr)


def _numpy_char_to_bytes(arr):
    """Like netCDF4.chartostring, but faster and more flexible."""
    # adapt handling of copy-kwarg to numpy 2.0
    # see https://github.com/numpy/numpy/issues/25916
    # and https://github.com/numpy/numpy/pull/25922
    copy = None if HAS_NUMPY_2_0 else False
    # based on: http://stackoverflow.com/a/10984878/809705
    arr = np.array(arr, copy=copy, order="C")
    dtype = "S" + str(arr.shape[-1])
    return arr.view(dtype).reshape(arr.shape[:-1])


class StackedBytesArray(indexing.ExplicitlyIndexedNDArrayMixin):
    """Wrapper around array-like objects to create a new indexable object where
    values, when accessed, are automatically stacked along the last dimension.

    >>> indexer = indexing.BasicIndexer((slice(None),))
    >>> StackedBytesArray(np.array(["a", "b", "c"], dtype="S1"))[indexer]
    array(b'abc', dtype='|S3')
    """

    def __init__(self, array):
        """
        Parameters
        ----------
        array : array-like
            Original array of values to wrap.
        """
        if array.dtype != "S1":
            raise ValueError(
                "can only use StackedBytesArray if argument has dtype='S1'"
            )
        self.array = indexing.as_indexable(array)

    @property
    def dtype(self):
        return np.dtype("S" + str(self.array.shape[-1]))

    @property
    def shape(self) -> tuple[int, ...]:
        return self.array.shape[:-1]

    def __repr__(self):
        return f"{type(self).__name__}({self.array!r})"

    def _vindex_get(self, key):
        return _numpy_char_to_bytes(self.array.vindex[key])

    def _oindex_get(self, key):
        return _numpy_char_to_bytes(self.array.oindex[key])

    def __getitem__(self, key):
        # require slicing the last dimension completely
        key = type(key)(indexing.expanded_indexer(key.tuple, self.array.ndim))
        if key.tuple[-1] != slice(None):
            raise IndexError("too many indices")
        return _numpy_char_to_bytes(self.array[key])
