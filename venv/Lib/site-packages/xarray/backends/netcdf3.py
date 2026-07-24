from __future__ import annotations

import unicodedata

import numpy as np

from xarray import coding
from xarray.core.variable import Variable

# Special characters that are permitted in netCDF names except in the
# 0th position of the string
_specialchars = '_.@+- !"#$%&\\()*,:;<=>?[]^`{|}~'

# The following are reserved names in CDL and may not be used as names of
# variables, dimension, attributes
_reserved_names = {
    "byte",
    "char",
    "short",
    "ushort",
    "int",
    "uint",
    "int64",
    "uint64",
    "float",
    "real",
    "double",
    "bool",
    "string",
}

# These data-types aren't supported by netCDF3, so they are automatically
# coerced instead as indicated by the "coerce_nc3_dtype" function
_nc3_dtype_coercions = {
    "int64": "int32",
    "uint64": "int32",
    "uint32": "int32",
    "uint16": "int16",
    "uint8": "int8",
    "bool": "int8",
}

# encode all strings as UTF-8
STRING_ENCODING = "utf-8"
COERCION_VALUE_ERROR = (
    "could not safely cast array from {dtype} to {new_dtype}. While it is not "
    "always the case, a common reason for this is that xarray has deemed it "
    "safest to encode np.datetime64[ns] or np.timedelta64[ns] values with "
    "int64 values representing units of 'nanoseconds'. This is either due to "
    "the fact that the times are known to require nanosecond precision for an "
    "accurate round trip, or that the times are unknown prior to writing due "
    "to being contained in a chunked array. Ways to work around this are "
    "either to use a backend that supports writing int64 values, or to "
    "manually specify the encoding['units'] and encoding['dtype'] (e.g. "
    "'seconds since 1970-01-01' and np.dtype('int32')) on the time "
    "variable(s) such that the times can be serialized in a netCDF3 file "
    "(note that depending on the situation, however, this latter option may "
    "result in an inaccurate round trip)."
)


def coerce_nc3_dtype(arr):
    """Coerce an array to a data type that can be stored in a netCDF-3 file

    This function performs the dtype conversions as specified by the
    ``_nc3_dtype_coercions`` mapping:
        int64  -> int32
        uint64 -> int32
        uint32 -> int32
        uint16 -> int16
        uint8  -> int8
        bool   -> int8

    Data is checked for equality, or equivalence (non-NaN values) using the
    ``(cast_array == original_array).all()``.
    """
    dtype = str(arr.dtype)
    if dtype in _nc3_dtype_coercions:
        new_dtype = _nc3_dtype_coercions[dtype]
        # TODO: raise a warning whenever casting the data-type instead?
        cast_arr = arr.astype(new_dtype)
        if not (cast_arr == arr).all():
            raise ValueError(
                COERCION_VALUE_ERROR.format(dtype=dtype, new_dtype=new_dtype)
            )
        arr = cast_arr
    return arr


def encode_nc3_attr_value(value):
    if isinstance(value, bytes):
        pass
    elif isinstance(value, str):
        value = value.encode(STRING_ENCODING)
    else:
        value = coerce_nc3_dtype(np.atleast_1d(value))
        if value.ndim > 1:
            raise ValueError("netCDF attributes must be 1-dimensional")
    return value


def encode_nc3_attrs(attrs):
    return {k: encode_nc3_attr_value(v) for k, v in attrs.items()}


def _maybe_prepare_times(var):
    # checks for integer-based time-like and
    # replaces np.iinfo(np.int64).min with _FillValue or np.nan
    # this keeps backwards compatibility

    data = var.data
    if data.dtype.kind in "iu":
        units = var.attrs.get("units", None)
        if units is not None and coding.variables._is_time_like(units):
            mask = data == np.iinfo(np.int64).min
            if mask.any():
                data = np.where(mask, var.attrs.get("_FillValue", np.nan), data)
    return data


def encode_nc3_variable(var, name=None):
    for coder in [
        coding.strings.EncodedStringCoder(allows_unicode=False),
        coding.strings.CharacterArrayCoder(),
    ]:
        var = coder.encode(var, name=name)
    data = _maybe_prepare_times(var)
    data = coerce_nc3_dtype(data)
    attrs = encode_nc3_attrs(var.attrs)
    return Variable(var.dims, data, attrs, var.encoding)


def _isalnumMUTF8(c):
    """Return True if the given UTF-8 encoded character is alphanumeric
    or multibyte.

    Input is not checked!
    """
    return c.isalnum() or (len(c.encode("utf-8")) > 1)


def is_valid_nc3_name(s):
    """Test whether an object can be validly converted to a netCDF-3
    dimension, variable or attribute name

    Earlier versions of the netCDF C-library reference implementation
    enforced a more restricted set of characters in creating new names,
    but permitted reading names containing arbitrary bytes. This
    specification extends the permitted characters in names to include
    multi-byte UTF-8 encoded Unicode and additional printing characters
    from the US-ASCII alphabet. The first character of a name must be
    alphanumeric, a multi-byte UTF-8 character, or '_' (reserved for
    special names with meaning to implementations, such as the
    "_FillValue" attribute). Subsequent characters may also include
    printing special characters, except for '/' which is not allowed in
    names. Names that have trailing space characters are also not
    permitted.
    """
    if not isinstance(s, str):
        return False
    num_bytes = len(s.encode("utf-8"))
    return (
        (unicodedata.normalize("NFC", s) == s)
        and (s not in _reserved_names)
        and (num_bytes >= 0)
        and ("/" not in s)
        and (s[-1] != " ")
        and (_isalnumMUTF8(s[0]) or (s[0] == "_"))
        and all(_isalnumMUTF8(c) or c in _specialchars for c in s)
    )
