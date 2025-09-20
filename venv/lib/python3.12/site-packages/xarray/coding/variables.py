"""Coders for individual Variable objects."""

from __future__ import annotations

import warnings
from collections.abc import Hashable, MutableMapping
from functools import partial
from typing import TYPE_CHECKING, Any, Union

import numpy as np
import pandas as pd

from xarray.coding.common import (
    SerializationWarning,
    VariableCoder,
    lazy_elemwise_func,
    pop_to,
    safe_setitem,
    unpack_for_decoding,
    unpack_for_encoding,
)
from xarray.coding.times import CFDatetimeCoder, CFTimedeltaCoder
from xarray.core import dtypes, duck_array_ops, indexing
from xarray.core.types import Self
from xarray.core.variable import Variable

if TYPE_CHECKING:
    T_VarTuple = tuple[tuple[Hashable, ...], Any, dict, dict]
    T_Name = Union[Hashable, None]


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
        return type(self)(self.array.oindex[key])

    def _vindex_get(self, key):
        return type(self)(self.array.vindex[key])

    def __getitem__(self, key) -> Self:
        return type(self)(self.array[key])

    def get_duck_array(self):
        return duck_array_ops.astype(self.array.get_duck_array(), dtype=self.dtype)

    def transpose(self, order):
        return type(self)(self.array.transpose(order))


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
        return type(self)(self.array.oindex[key])

    def _vindex_get(self, key):
        return type(self)(self.array.vindex[key])

    def __getitem__(self, key) -> Self:
        return type(self)(self.array[key])

    def get_duck_array(self):
        return duck_array_ops.astype(self.array.get_duck_array(), dtype=self.dtype)

    def transpose(self, order):
        return type(self)(self.array.transpose(order))


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
    # return "datetime" for datetime-like
    # return "timedelta" for timedelta-like
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
        return "datetime"
    else:
        return "timedelta" if any(tstr == units for tstr in time_strings) else False


def _check_fill_values(attrs, name, dtype):
    """Check _FillValue and missing_value if available.

    Return dictionary with raw fill values and set with encoded fill values.

    Issue SerializationWarning if appropriate.
    """
    raw_fill_dict = {}
    for attr in ("missing_value", "_FillValue"):
        pop_to(attrs, raw_fill_dict, attr, name=name)
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


def _convert_unsigned_fill_value(
    name: T_Name,
    data: Any,
    unsigned: str,
    raw_fill_value: Any,
    encoded_fill_values: set,
) -> Any:
    if data.dtype.kind == "i":
        if unsigned == "true":
            unsigned_dtype = np.dtype(f"u{data.dtype.itemsize}")
            transform = partial(np.asarray, dtype=unsigned_dtype)
            if raw_fill_value is not None:
                new_fill = np.array(raw_fill_value, dtype=data.dtype)
                encoded_fill_values.remove(raw_fill_value)
                # use view here to prevent OverflowError
                encoded_fill_values.add(new_fill.view(unsigned_dtype).item())
            data = lazy_elemwise_func(data, transform, unsigned_dtype)
    elif data.dtype.kind == "u":
        if unsigned == "false":
            signed_dtype = np.dtype(f"i{data.dtype.itemsize}")
            transform = partial(np.asarray, dtype=signed_dtype)
            data = lazy_elemwise_func(data, transform, signed_dtype)
            if raw_fill_value is not None:
                new_fill = signed_dtype.type(raw_fill_value)
                encoded_fill_values.remove(raw_fill_value)
                encoded_fill_values.add(new_fill)
    else:
        warnings.warn(
            f"variable {name!r} has _Unsigned attribute but is not "
            "of integer type. Ignoring attribute.",
            SerializationWarning,
            stacklevel=3,
        )
    return data


def _encode_unsigned_fill_value(
    name: T_Name,
    fill_value: Any,
    encoded_dtype: np.dtype,
) -> Any:
    try:
        if hasattr(fill_value, "item"):
            # if numpy type, convert to python native integer to determine overflow
            # otherwise numpy unsigned ints will silently cast to the signed counterpart
            fill_value = fill_value.item()
        # passes if provided fill value fits in encoded on-disk type
        new_fill = encoded_dtype.type(fill_value)
    except OverflowError:
        encoded_kind_str = "signed" if encoded_dtype.kind == "i" else "unsigned"
        warnings.warn(
            f"variable {name!r} will be stored as {encoded_kind_str} integers "
            f"but _FillValue attribute can't be represented as a "
            f"{encoded_kind_str} integer.",
            SerializationWarning,
            stacklevel=3,
        )
        # user probably provided the fill as the in-memory dtype,
        # convert to on-disk type to match CF standard
        orig_kind = "u" if encoded_dtype.kind == "i" else "i"
        orig_dtype = np.dtype(f"{orig_kind}{encoded_dtype.itemsize}")
        # use view here to prevent OverflowError
        new_fill = np.array(fill_value, dtype=orig_dtype).view(encoded_dtype).item()
    return new_fill


class CFMaskCoder(VariableCoder):
    """Mask or unmask fill values according to CF conventions."""

    def __init__(
        self,
        decode_times: bool | CFDatetimeCoder = False,
        decode_timedelta: bool | CFTimedeltaCoder = False,
    ) -> None:
        self.decode_times = decode_times
        self.decode_timedelta = decode_timedelta

    def encode(self, variable: Variable, name: T_Name = None):
        dims, data, attrs, encoding = unpack_for_encoding(variable)

        dtype = np.dtype(encoding.get("dtype", data.dtype))
        # from netCDF best practices
        # https://docs.unidata.ucar.edu/nug/current/best_practices.html#bp_Unsigned-Data
        #     "_Unsigned = "true" to indicate that
        #      integer data should be treated as unsigned"
        has_unsigned = encoding.get("_Unsigned") is not None
        fv = encoding.get("_FillValue")
        mv = encoding.get("missing_value")
        fill_value = None

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
            # but not for packed data
            if has_unsigned:
                encoding["_FillValue"] = _encode_unsigned_fill_value(name, fv, dtype)
            elif "add_offset" not in encoding and "scale_factor" not in encoding:
                encoding["_FillValue"] = dtype.type(fv)
            else:
                encoding["_FillValue"] = fv
            fill_value = pop_to(encoding, attrs, "_FillValue", name=name)

        if mv_exists:
            # try to use _FillValue, if it exists to align both values
            # or use missing_value and ensure it's cast to same dtype as data's
            # but not for packed data
            encoding["missing_value"] = attrs.get(
                "_FillValue",
                (
                    _encode_unsigned_fill_value(name, mv, dtype)
                    if has_unsigned
                    else (
                        dtype.type(mv)
                        if "add_offset" not in encoding
                        and "scale_factor" not in encoding
                        else mv
                    )
                ),
            )
            fill_value = pop_to(encoding, attrs, "missing_value", name=name)

        # apply fillna
        if fill_value is not None and not pd.isnull(fill_value):
            # special case DateTime to properly handle NaT
            if _is_time_like(attrs.get("units")):
                if data.dtype.kind in "iu":
                    data = duck_array_ops.where(
                        data != np.iinfo(np.int64).min, data, fill_value
                    )
                else:
                    # if we have float data (data was packed prior masking)
                    # we just fillna
                    data = duck_array_ops.fillna(data, fill_value)
                    # but if the fill_value is of integer type
                    # we need to round and cast
                    if np.array(fill_value).dtype.kind in "iu":
                        data = duck_array_ops.astype(
                            duck_array_ops.around(data), type(fill_value)
                        )
            else:
                data = duck_array_ops.fillna(data, fill_value)

        if fill_value is not None and has_unsigned:
            pop_to(encoding, attrs, "_Unsigned")
            # XXX: Is this actually needed? Doesn't the backend handle this?
            # two-stage casting to prevent undefined cast from float to unsigned int
            # first float -> int with corresponding itemsize
            # second int -> int/uint to final itemsize
            signed_dtype = np.dtype(f"i{data.itemsize}")
            data = duck_array_ops.astype(
                duck_array_ops.astype(
                    duck_array_ops.around(data), signed_dtype, copy=False
                ),
                dtype,
                copy=False,
            )
            attrs["_FillValue"] = fill_value

        return Variable(dims, data, attrs, encoding, fastpath=True)

    def decode(self, variable: Variable, name: T_Name = None):
        raw_fill_dict, encoded_fill_values = _check_fill_values(
            variable.attrs, name, variable.dtype
        )
        if "_Unsigned" not in variable.attrs and not raw_fill_dict:
            return variable

        dims, data, attrs, encoding = unpack_for_decoding(variable)

        # Even if _Unsigned is used, retain on-disk _FillValue
        for attr, value in raw_fill_dict.items():
            safe_setitem(encoding, attr, value, name=name)

        if "_Unsigned" in attrs:
            unsigned = pop_to(attrs, encoding, "_Unsigned")
            data = _convert_unsigned_fill_value(
                name,
                data,
                unsigned,
                raw_fill_dict.get("_FillValue"),
                encoded_fill_values,
            )

        if encoded_fill_values:
            dtype: np.typing.DTypeLike
            decoded_fill_value: Any
            # in case of packed data we have to decode into float
            # in any case
            if "scale_factor" in attrs or "add_offset" in attrs:
                dtype, decoded_fill_value = (
                    _choose_float_dtype(data.dtype, attrs),
                    np.nan,
                )
            else:
                # in case of no-packing special case DateTime/Timedelta to properly
                # handle NaT, we need to check if time-like will be decoded
                # or not in further processing
                is_time_like = _is_time_like(attrs.get("units"))
                if (
                    (is_time_like == "datetime" and self.decode_times)
                    or (is_time_like == "timedelta" and self.decode_timedelta)
                ) and data.dtype.kind in "iu":
                    dtype = np.int64
                    decoded_fill_value = np.iinfo(np.int64).min
                else:
                    dtype, decoded_fill_value = dtypes.maybe_promote(data.dtype)

            transform = partial(
                _apply_mask,
                encoded_fill_values=encoded_fill_values,
                decoded_fill_value=decoded_fill_value,
                dtype=dtype,
            )
            data = lazy_elemwise_func(data, transform, dtype)

        return Variable(dims, data, attrs, encoding, fastpath=True)


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
    # Todo: with nc-complex from netcdf4-python >= 1.7.0 this is available
    # (safe because eg. complex numbers are not supported in NetCDF)
    return np.float64


class CFScaleOffsetCoder(VariableCoder):
    """Scale and offset variables according to CF conventions.

    Follows the formula:
        decode_values = encoded_values * scale_factor + add_offset
    """

    def __init__(
        self,
        decode_times: bool | CFDatetimeCoder = False,
        decode_timedelta: bool | CFTimedeltaCoder = False,
    ) -> None:
        self.decode_times = decode_times
        self.decode_timedelta = decode_timedelta

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
            if duck_array_ops.ndim(scale_factor) > 0:
                scale_factor = np.asarray(scale_factor).item()
            if duck_array_ops.ndim(add_offset) > 0:
                add_offset = np.asarray(add_offset).item()
            # if we have a _FillValue/masked_value in encoding we already have the wanted
            # floating point dtype here (via CFMaskCoder), so no check is necessary
            # only check in other cases and for time-like
            dtype = data.dtype
            is_time_like = _is_time_like(attrs.get("units"))
            if (
                ("_FillValue" not in encoding and "missing_value" not in encoding)
                or (is_time_like == "datetime" and self.decode_times)
                or (is_time_like == "timedelta" and self.decode_timedelta)
            ):
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
                    data = duck_array_ops.round(data)
                data = duck_array_ops.astype(data, dtype=dtype)
            return Variable(dims, data, attrs, encoding, fastpath=True)
        else:
            return variable

    def decode(self):
        raise NotImplementedError()


class ObjectVLenStringCoder(VariableCoder):
    def encode(self):
        raise NotImplementedError

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        if variable.dtype.kind == "O" and variable.encoding.get("dtype", False) is str:
            variable = variable.astype(variable.encoding["dtype"])
            return variable
        else:
            return variable


class Numpy2StringDTypeCoder(VariableCoder):
    # Convert Numpy 2 StringDType arrays to object arrays for backwards compatibility
    # TODO: remove this if / when we decide to allow StringDType arrays in Xarray
    def encode(self):
        raise NotImplementedError

    def decode(self, variable: Variable, name: T_Name = None) -> Variable:
        if variable.dtype.kind == "T":
            return variable.astype(object)
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
