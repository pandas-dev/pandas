from __future__ import annotations

from collections.abc import Callable, Hashable, MutableMapping
from typing import TYPE_CHECKING, Any, Union

import numpy as np

from xarray.core import indexing
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array

if TYPE_CHECKING:
    T_VarTuple = tuple[tuple[Hashable, ...], Any, dict, dict]
    T_Name = Union[Hashable, None]


class SerializationWarning(RuntimeWarning):
    """Warnings about encoding/decoding issues in serialization."""


class VariableCoder:
    """Base class for encoding and decoding transformations on variables.

    We use coders for transforming variables between xarray's data model and
    a format suitable for serialization. For example, coders apply CF
    conventions for how data should be represented in netCDF files.

    Subclasses should implement encode() and decode(), which should satisfy
    the identity ``coder.decode(coder.encode(variable)) == variable``. If any
    options are necessary, they should be implemented as arguments to the
    __init__ method.

    The optional name argument to encode() and decode() exists solely for the
    sake of better error messages, and should correspond to the name of
    variables in the underlying store.
    """

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        """Convert an encoded variable to a decoded variable"""
        raise NotImplementedError()

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        """Convert a decoded variable to an encoded variable"""
        raise NotImplementedError()


class _ElementwiseFunctionArray(indexing.ExplicitlyIndexedNDArrayMixin):
    """Lazily computed array holding values of elemwise-function.

    Do not construct this object directly: call lazy_elemwise_func instead.

    Values are computed upon indexing or coercion to a NumPy array.
    """

    def __init__(self, array, func: Callable, dtype: np.typing.DTypeLike | None):
        assert not is_chunked_array(array)
        self.array = indexing.as_indexable(array)
        self.func = func
        self._dtype = dtype

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(self._dtype)

    def transpose(self, order):
        # For elementwise functions, we can compose transpose and function application
        return type(self)(self.array.transpose(order), self.func, self.dtype)

    def _oindex_get(self, key):
        return type(self)(self.array.oindex[key], self.func, self.dtype)

    def _vindex_get(self, key):
        return type(self)(self.array.vindex[key], self.func, self.dtype)

    def __getitem__(self, key):
        return type(self)(self.array[key], self.func, self.dtype)

    def get_duck_array(self):
        return self.func(self.array.get_duck_array())

    async def async_get_duck_array(self):
        return self.func(await self.array.async_get_duck_array())

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.array!r}, func={self.func!r}, dtype={self.dtype!r})"


def lazy_elemwise_func(array, func: Callable, dtype: np.typing.DTypeLike | None):
    """Lazily apply an element-wise function to an array.
    Parameters
    ----------
    array : any valid value of Variable._data
    func : callable
        Function to apply to indexed slices of an array. For use with dask,
        this should be a pickle-able object.
    dtype : coercible to np.dtype
        Dtype for the result of this function.

    Returns
    -------
    Either a dask.array.Array or _ElementwiseFunctionArray.
    """
    if is_chunked_array(array):
        chunkmanager = get_chunked_array_type(array)

        return chunkmanager.map_blocks(func, array, dtype=dtype)  # type: ignore[arg-type]
    else:
        return _ElementwiseFunctionArray(array, func, dtype)


def safe_setitem(dest, key: Hashable, value, name: T_Name = None):
    if key in dest:
        var_str = f" on variable {name!r}" if name else ""
        raise ValueError(
            f"Key '{key}' already exists in attrs{var_str}, and will not be overwritten. "
            "This is probably an encoding field used by xarray to describe "
            "how a variable is serialized. To proceed, remove this key from "
            "the variable's attributes manually."
        )
    dest[key] = value


def pop_to(
    source: MutableMapping, dest: MutableMapping, key: Hashable, name: T_Name = None
) -> Any:
    """
    A convenience function which pops a key k from source to dest.
    None values are not passed on.  If k already exists in dest an
    error is raised.
    """
    value = source.pop(key, None)
    if value is not None:
        safe_setitem(dest, key, value, name=name)
    return value


def unpack_for_encoding(var: Variable) -> T_VarTuple:
    return var.dims, var.data, var.attrs.copy(), var.encoding.copy()


def unpack_for_decoding(var: Variable) -> T_VarTuple:
    return var.dims, var._data, var.attrs.copy(), var.encoding.copy()
