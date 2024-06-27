"""Coders for individual Variable objects."""

from __future__ import annotations

import warnings
from collections.abc import Hashable, MutableMapping
from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Union

import numpy as np
import pandas as pd

from xarray.core import dtypes, duck_array_ops, indexing
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
        """Convert an decoded variable to a encoded variable"""
        raise NotImplementedError()


class _ElementwiseFunctionArray(indexing.ExplicitlyIndexedNDArrayMixin):
    """Lazily computed array holding values of elemwise-function.

    Do not construct this object directly: call lazy_elemwise_func instead.

    Values are computed upon indexing or coercion to a NumPy array.
    """

    def __init__(self, array, func: Callable, dtype: np.typing.DTypeLike):
        assert not is_chunked_array(array)
        self.array = indexing.as_indexable(array)
        self.func = func
        self._dtype = dtype

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(self._dtype)

    def _oindex_get(self, key):
        return type(self)(self.array.oindex[key], self.func, self.dtype)

    def _vindex_get(self, key):
        return type(self)(self.array.vindex[key], self.func, self.dtype)

    def __getitem__(self, key):
        return type(self)(self.array[key], self.func, self.dtype)

    def get_duck_array(self):
        return self.func(self.array.get_duck_array())

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.array!r}, func={self.func!r}, dtype={self.dtype!r})"


class NativeEndiannessArray(indexing.ExplicitlyIndexedNDArrayMixin):
    """Decode arrays on the fly from non-native to native endianness

    This is useful for decoding arrays from netCDF3 files (which are all
    big endian) into native endianness, so they can be used with Cython
    functions, such as those found in bottleneck and pandas.

    >>> x = np.arange(5, dtype=">i2")

    >>> x.dtype
    dtype('>i2')

    >>> NativeEndiannessArray(x).dtype
    dtype('int16')

    >>> indexer = indexing.BasicIndexer((slice(None),))
    >>> NativeEndiannessArray(x)[indexer].dtype
    dtype('int16')
    """

    __slots__ = ("array",)

    def __init__(self, array) -> None:
        self.array = indexing.as_indexable(array)

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(self.array.dtype.kind + str(self.array.dtype.itemsize))

    def _oindex_get(self, key):
        return np.asarray(self.array.oindex[key], dtype=self.dtype)

    def _vindex_get(self, key):
        return np.asarray(self.array.vindex[key], dtype=self.dtype)

    def __getitem__(self, key) -> np.ndarray:
        return np.asarray(self.array[key], dtype=self.dtype)


class BoolTypeArray(indexing.ExplicitlyIndexedNDArrayMixin):
    """Decode arrays on the fly from integer to boolean datatype

    This is useful for decoding boolean arrays from integer typed netCDF
    variables.

    >>> x = np.array([1, 0, 1, 1, 0], dtype="i1")

    >>> x.dtype
    dtype('int8')

    >>> BoolTypeArray(x).dtype
    dtype('bool')

    >>> indexer = indexing.BasicIndexer((slice(None),))
    >>> BoolTypeArray(x)[indexer].dtype
    dtype('bool')
    """

    __slots__ = ("array",)

    def __init__(self, array) -> None:
        self.array = indexing.as_indexable(array)

    @property
    def dtype(self) -> np.dtype:
        return np.dtype("bool")

    def _oindex_get(self, key):
        return np.asarray(self.array.oindex[key], dtype=self.dtype)

    def _vindex_get(self, key):
        return np.asarray(self.array.vindex[key], dtype=self.dtype)

    def __getitem__(self, key) -> np.ndarray:
        return np.asarray(self.array[key], dtype=self.dtype)


def lazy_elemwise_func(array, func: Callable, dtype: np.typing.DTypeLike):
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


def unpack_for_encoding(var: Variable) -> T_VarTuple:
    return var.dims, var.data, var.attrs.copy(), var.encoding.copy()


def unpack_for_decoding(var: Variable) -> T_VarTuple:
    return var.dims, var._data, var.attrs.copy(), var.encoding.copy()


def safe_setitem(dest, key: Hashable, value, name: T_Name = None):
    if key in dest:
        var_str = f" on variable {name!r}" if name else ""
        raise ValueError(
            f"failed to prevent overwriting existing key {key} in attrs{var_str}. "
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


def _apply_mask(
    data: np.ndarray,
    encoded_fill_values: list,
    decoded_fill_value: Any,
    dtype: np.typing.DTypeLike,
) -> np.ndarray:
    """Mask all matching values in a NumPy arrays."""
    data = np.asarray(data, dtype=dtype)
    condition = False
    for fv in encoded_fill_values:
        condition |= data == fv
    return np.where(condition, decoded_fill_value, data)


def _is_time_like(units):
    # test for time-like
    if units is None:
        return False
    time_strings = [
        "days",
        "hours",
        "minutes",
        "seconds",
        "milliseconds",
        "microseconds",
        "nanoseconds",
    ]
    units = str(units)
    # to prevent detecting units like `days accumulated` as time-like
    # special casing for datetime-units and timedelta-units (GH-8269)
    if "since" in units:
        from xarray.coding.times import _unpack_netcdf_time_units

        try:
            _unpack_netcdf_time_units(units)
        except ValueError:
            return False
        return True
    else:
        return any(tstr == units for tstr in time_strings)


def _check_fill_values(attrs, name, dtype):
    """ "Check _FillValue and missing_value if available.

    Return dictionary with raw fill values and set with encoded fill values.

    Issue SerializationWarning if appropriate.
    """
    raw_fill_dict = {}
    [
        pop_to(attrs, raw_fill_dict, attr, name=name)
        for attr in ("missing_value", "_FillValue")
    ]
    encoded_fill_values = set()
    for k in list(raw_fill_dict):
        v = raw_fill_dict[k]
        kfill = {fv for fv in np.ravel(v) if not pd.isnull(fv)}
        if not kfill and np.issubdtype(dtype, np.integer):
            warnings.warn(
                f"variable {name!r} has non-conforming {k!r} "
                f"{v!r} defined, dropping {k!r} entirely.",
                SerializationWarning,
                stacklevel=3,
            )
            del raw_fill_dict[k]
        else:
            encoded_fill_values |= kfill

        if len(encoded_fill_values) > 1:
            warnings.warn(
                f"variable {name!r} has multiple fill values "
                f"{encoded_fill_values} defined, decoding all values to NaN.",
                SerializationWarning,
                stacklevel=3,
            )

    return raw_fill_dict, encoded_fill_values


class CFMaskCoder(VariableCoder):
    """Mask or unmask fill values according to CF conventions."""

    def encode(self, variable: Variable, name: T_Name = None):
        dims, data, attrs, encoding = unpack_for_encoding(variable)

        dtype = np.dtype(encoding.get("dtype", data.dtype))
        fv = encoding.get("_FillValue")
        mv = encoding.get("missing_value")
        # to properly handle _FillValue/missing_value below [a], [b]
        # we need to check if unsigned data is written as signed data
        unsigned = encoding.get("_Unsigned") is not None

        fv_exists = fv is not None
        mv_exists = mv is not None

        if not fv_exists and not mv_exists:
            return variable

        if fv_exists and mv_exists and not duck_array_ops.allclose_or_equiv(fv, mv):
            raise ValueError(
                f"Variable {name!r} has conflicting _FillValue ({fv}) and missing_value ({mv}). Cannot encode data."
            )

        if fv_exists:
            # Ensure _FillValue is cast to same dtype as data's
            # [a] need to skip this if _Unsigned is available
            if not unsigned:
                encoding["_FillValue"] = dtype.type(fv)
            fill_value = pop_to(encoding, attrs, "_FillValue", name=name)

        if mv_exists:
            # try to use _FillValue, if it exists to align both values
            # or use missing_value and ensure it's cast to same dtype as data's
            # [b] need to provide mv verbatim if _Unsigned is available
            encoding["missing_value"] = attrs.get(
                "_FillValue",
                (dtype.type(mv) if not unsigned else mv),
            )
            fill_value = pop_to(encoding, attrs, "missing_value", name=name)

        # apply fillna
        if not pd.isnull(fill_value):
            # special case DateTime to properly handle NaT
            if _is_time_like(attrs.get("units")) and data.dtype.kind in "iu":
                data = duck_array_ops.where(
                    data != np.iinfo(np.int64).min, data, fill_value
                )
            else:
                data = duck_array_ops.fillna(data, fill_value)

        return Variable(dims, data, attrs, encoding, fastpath=True)

    def decode(self, variable: Variable, name: T_Name = None):
        raw_fill_dict, encoded_fill_values = _check_fill_values(
            variable.attrs, name, variable.dtype
        )

        if raw_fill_dict:
            dims, data, attrs, encoding = unpack_for_decoding(variable)
            [
                safe_setitem(encoding, attr, value, name=name)
                for attr, value in raw_fill_dict.items()
            ]

            if encoded_fill_values:
                # special case DateTime to properly handle NaT
                dtype: np.typing.DTypeLike
                decoded_fill_value: Any
                if _is_time_like(attrs.get("units")) and data.dtype.kind in "iu":
                    dtype, decoded_fill_value = np.int64, np.iinfo(np.int64).min
                else:
                    if "scale_factor" not in attrs and "add_offset" not in attrs:
                        dtype, decoded_fill_value = dtypes.maybe_promote(data.dtype)
                    else:
                        dtype, decoded_fill_value = (
                            _choose_float_dtype(data.dtype, attrs),
                            np.nan,
                        )

                transform = partial(
                    _apply_mask,
                    encoded_fill_values=encoded_fill_values,
                    decoded_fill_value=decoded_fill_value,
                    dtype=dtype,
                )
                data = lazy_elemwise_func(data, transform, dtype)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


def _scale_offset_decoding(data, scale_factor, add_offset, dtype: np.typing.DTypeLike):
    data = data.astype(dtype=dtype, copy=True)
    if scale_factor is not None:
        data *= scale_factor
    if add_offset is not None:
        data += add_offset
    return data


def _choose_float_dtype(
    dtype: np.dtype, mapping: MutableMapping
) -> type[np.floating[Any]]:
    """Return a float dtype that can losslessly represent `dtype` values."""
    # check scale/offset first to derive wanted float dtype
    # see https://github.com/pydata/xarray/issues/5597#issuecomment-879561954
    scale_factor = mapping.get("scale_factor")
    add_offset = mapping.get("add_offset")
    if scale_factor is not None or add_offset is not None:
        # get the type from scale_factor/add_offset to determine
        # the needed floating point type
        if scale_factor is not None:
            scale_type = np.dtype(type(scale_factor))
        if add_offset is not None:
            offset_type = np.dtype(type(add_offset))
        # CF conforming, both scale_factor and add-offset are given and
        # of same floating point type (float32/64)
        if (
            add_offset is not None
            and scale_factor is not None
            and offset_type == scale_type
            and scale_type in [np.float32, np.float64]
        ):
            # in case of int32 -> we need upcast to float64
            # due to precision issues
            if dtype.itemsize == 4 and np.issubdtype(dtype, np.integer):
                return np.float64
            return scale_type.type
        # Not CF conforming and add_offset given:
        # A scale factor is entirely safe (vanishing into the mantissa),
        # but a large integer offset could lead to loss of precision.
        # Sensitivity analysis can be tricky, so we just use a float64
        # if there's any offset at all - better unoptimised than wrong!
        if add_offset is not None:
            return np.float64
        # return dtype depending on given scale_factor
        return scale_type.type
    # If no scale_factor or add_offset is given, use some general rules.
    # Keep float32 as-is. Upcast half-precision to single-precision,
    # because float16 is "intended for storage but not computation"
    if dtype.itemsize <= 4 and np.issubdtype(dtype, np.floating):
        return np.float32
    # float32 can exactly represent all integers up to 24 bits
    if dtype.itemsize <= 2 and np.issubdtype(dtype, np.integer):
        return np.float32
    # For all other types and circumstances, we just use float64.
    # (safe because eg. complex numbers are not supported in NetCDF)
    return np.float64


class CFScaleOffsetCoder(VariableCoder):
    """Scale and offset variables according to CF conventions.

    Follows the formula:
        decode_values = encoded_values * scale_factor + add_offset
    """

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        dims, data, attrs, encoding = unpack_for_encoding(variable)

        if "scale_factor" in encoding or "add_offset" in encoding:
            # if we have a _FillValue/masked_value we do not want to cast now
            # but leave that to CFMaskCoder
            dtype = data.dtype
            if "_FillValue" not in encoding and "missing_value" not in encoding:
                dtype = _choose_float_dtype(data.dtype, encoding)
            # but still we need a copy prevent changing original data
            data = duck_array_ops.astype(data, dtype=dtype, copy=True)
        if "add_offset" in encoding:
            data -= pop_to(encoding, attrs, "add_offset", name=name)
        if "scale_factor" in encoding:
            data /= pop_to(encoding, attrs, "scale_factor", name=name)

        return Variable(dims, data, attrs, encoding, fastpath=True)

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        _attrs = variable.attrs
        if "scale_factor" in _attrs or "add_offset" in _attrs:
            dims, data, attrs, encoding = unpack_for_decoding(variable)

            scale_factor = pop_to(attrs, encoding, "scale_factor", name=name)
            add_offset = pop_to(attrs, encoding, "add_offset", name=name)
            if np.ndim(scale_factor) > 0:
                scale_factor = np.asarray(scale_factor).item()
            if np.ndim(add_offset) > 0:
                add_offset = np.asarray(add_offset).item()
            # if we have a _FillValue/masked_value we already have the wanted
            # floating point dtype here (via CFMaskCoder), so no check is necessary
            # only check in other cases
            dtype = data.dtype
            if "_FillValue" not in encoding and "missing_value" not in encoding:
                dtype = _choose_float_dtype(dtype, encoding)

            transform = partial(
                _scale_offset_decoding,
                scale_factor=scale_factor,
                add_offset=add_offset,
                dtype=dtype,
            )
            data = lazy_elemwise_func(data, transform, dtype)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


class UnsignedIntegerCoder(VariableCoder):
    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        # from netCDF best practices
        # https://docs.unidata.ucar.edu/nug/current/best_practices.html#bp_Unsigned-Data
        #     "_Unsigned = "true" to indicate that
        #      integer data should be treated as unsigned"
        if variable.encoding.get("_Unsigned", "false") == "true":
            dims, data, attrs, encoding = unpack_for_encoding(variable)

            pop_to(encoding, attrs, "_Unsigned")
            signed_dtype = np.dtype(f"i{data.dtype.itemsize}")
            if "_FillValue" in attrs:
                new_fill = signed_dtype.type(attrs["_FillValue"])
                attrs["_FillValue"] = new_fill
            data = duck_array_ops.astype(duck_array_ops.around(data), signed_dtype)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        if "_Unsigned" in variable.attrs:
            dims, data, attrs, encoding = unpack_for_decoding(variable)
            unsigned = pop_to(attrs, encoding, "_Unsigned")

            if data.dtype.kind == "i":
                if unsigned == "true":
                    unsigned_dtype = np.dtype(f"u{data.dtype.itemsize}")
                    transform = partial(np.asarray, dtype=unsigned_dtype)
                    data = lazy_elemwise_func(data, transform, unsigned_dtype)
                    if "_FillValue" in attrs:
                        new_fill = unsigned_dtype.type(attrs["_FillValue"])
                        attrs["_FillValue"] = new_fill
            elif data.dtype.kind == "u":
                if unsigned == "false":
                    signed_dtype = np.dtype(f"i{data.dtype.itemsize}")
                    transform = partial(np.asarray, dtype=signed_dtype)
                    data = lazy_elemwise_func(data, transform, signed_dtype)
                    if "_FillValue" in attrs:
                        new_fill = signed_dtype.type(attrs["_FillValue"])
                        attrs["_FillValue"] = new_fill
            else:
                warnings.warn(
                    f"variable {name!r} has _Unsigned attribute but is not "
                    "of integer type. Ignoring attribute.",
                    SerializationWarning,
                    stacklevel=3,
                )

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


class DefaultFillvalueCoder(VariableCoder):
    """Encode default _FillValue if needed."""

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        dims, data, attrs, encoding = unpack_for_encoding(variable)
        # make NaN the fill value for float types
        if (
            "_FillValue" not in attrs
            and "_FillValue" not in encoding
            and np.issubdtype(variable.dtype, np.floating)
        ):
            attrs["_FillValue"] = variable.dtype.type(np.nan)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        raise NotImplementedError()


class BooleanCoder(VariableCoder):
    """Code boolean values."""

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        if (
            (variable.dtype == bool)
            and ("dtype" not in variable.encoding)
            and ("dtype" not in variable.attrs)
        ):
            dims, data, attrs, encoding = unpack_for_encoding(variable)
            attrs["dtype"] = "bool"
            data = duck_array_ops.astype(data, dtype="i1", copy=True)

            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        if variable.attrs.get("dtype", False) == "bool":
            dims, data, attrs, encoding = unpack_for_decoding(variable)
            # overwrite (!) dtype in encoding, and remove from attrs
            # needed for correct subsequent encoding
            encoding["dtype"] = attrs.pop("dtype")
            data = BoolTypeArray(data)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


class EndianCoder(VariableCoder):
    """Decode Endianness to native."""

    def encode(self):
        raise NotImplementedError()

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        dims, data, attrs, encoding = unpack_for_decoding(variable)
        if not data.dtype.isnative:
            data = NativeEndiannessArray(data)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable


class NonStringCoder(VariableCoder):
    """Encode NonString variables if dtypes differ."""

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        if "dtype" in variable.encoding and variable.encoding["dtype"] not in (
            "S1",
            str,
        ):
            dims, data, attrs, encoding = unpack_for_encoding(variable)
            dtype = np.dtype(encoding.pop("dtype"))
            if dtype != variable.dtype:
                if np.issubdtype(dtype, np.integer):
                    if (
                        np.issubdtype(variable.dtype, np.floating)
                        and "_FillValue" not in variable.attrs
                        and "missing_value" not in variable.attrs
                    ):
                        warnings.warn(
                            f"saving variable {name} with floating "
                            "point data as an integer dtype without "
                            "any _FillValue to use for NaNs",
                            SerializationWarning,
                            stacklevel=10,
                        )
                    data = np.around(data)
                data = data.astype(dtype=dtype)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self):
        raise NotImplementedError()


class ObjectVLenStringCoder(VariableCoder):
    def encode(self):
        raise NotImplementedError

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        if variable.dtype == object and variable.encoding.get("dtype", False) == str:
            variable = variable.astype(variable.encoding["dtype"])
            return variable
        else:
            return variable


class NativeEnumCoder(VariableCoder):
    """Encode Enum into variable dtype metadata."""

    def encode(self, variable: Variable, name: T_Name = None) -> Variable:
        if (
            "dtype" in variable.encoding
            and np.dtype(variable.encoding["dtype"]).metadata
            and "enum" in variable.encoding["dtype"].metadata
        ):
            dims, data, attrs, encoding = unpack_for_encoding(variable)
            data = data.astype(dtype=variable.encoding.pop("dtype"))
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        raise NotImplementedError()
