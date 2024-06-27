# -#- coding: utf-8 -#-
"""
Deal with parquet logical types (aka converted types), higher-order things built from primitive types.

The implementations in this class are pure python for the widest compatibility,
but they're not necessarily the most performant.
"""

import logging
import numpy as np
import pandas as pd

from fastparquet import parquet_thrift
from fastparquet.cencoding import time_shift
from fastparquet.json import json_decoder

logger = logging.getLogger('parquet')  # pylint: disable=invalid-name

try:
    from bson import BSON
    unbson = BSON.decode
    tobson = BSON.encode
except ImportError:  # pragma: no cover
    try:
        import bson
        unbson = bson.loads
        tobson = bson.dumps
    except:
        def unbson(x):
            raise ImportError("BSON not found")
        def tobson(x):
            raise ImportError("BSON not found")

# Explicitly use numpy type in order to avoid promotion errors due to NEP 50 in numpy >= 2
DAYS_TO_NANOS = np.int64(86400000000000)
"""Number of nanoseconds in a day. Used to convert a Date to a date"""
nat = np.datetime64('NaT').view('int64')

simple = {
    parquet_thrift.Type.INT32: np.dtype('int32'),
    parquet_thrift.Type.INT64: np.dtype('int64'),
    parquet_thrift.Type.FLOAT: np.dtype('float32'),
    parquet_thrift.Type.DOUBLE: np.dtype('float64'),
    parquet_thrift.Type.BOOLEAN: np.dtype('bool'),
    parquet_thrift.Type.INT96: np.dtype('S12'),
    parquet_thrift.Type.BYTE_ARRAY: np.dtype("O"),
    parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY: np.dtype("O")
}
complex = {
    parquet_thrift.ConvertedType.UTF8: np.dtype("O"),
    parquet_thrift.ConvertedType.DECIMAL: np.dtype('float64'),
    parquet_thrift.ConvertedType.UINT_8: np.dtype("uint8"),
    parquet_thrift.ConvertedType.UINT_16: np.dtype("uint16"),
    parquet_thrift.ConvertedType.UINT_32: np.dtype('uint32'),
    parquet_thrift.ConvertedType.UINT_64: np.dtype('uint64'),
    parquet_thrift.ConvertedType.INT_8: np.dtype("int8"),
    parquet_thrift.ConvertedType.INT_16: np.dtype("int16"),
    parquet_thrift.ConvertedType.INT_32: np.dtype('int32'),
    parquet_thrift.ConvertedType.INT_64: np.dtype('int64'),
    parquet_thrift.ConvertedType.TIME_MILLIS: np.dtype('<m8[ms]'),
    parquet_thrift.ConvertedType.DATE: np.dtype('<M8[ns]'),
    parquet_thrift.ConvertedType.TIMESTAMP_MILLIS: np.dtype('<M8[ms]'),
    parquet_thrift.ConvertedType.TIME_MICROS: np.dtype('<m8[us]'),
    parquet_thrift.ConvertedType.TIMESTAMP_MICROS: np.dtype('<M8[us]')
}
nullable = {
    np.dtype('int8'): pd.Int8Dtype(),
    np.dtype('int16'): pd.Int16Dtype(),
    np.dtype('int32'): pd.Int32Dtype(),
    np.dtype('int64'): pd.Int64Dtype(),
    np.dtype('uint8'): pd.UInt8Dtype(),
    np.dtype('uint16'): pd.UInt16Dtype(),
    np.dtype('uint32'): pd.UInt32Dtype(),
    np.dtype('uint64'): pd.UInt64Dtype(),
    np.dtype('bool'): pd.BooleanDtype()
}
pandas_nullable = {
    "Int8": pd.Int8Dtype(),
    "Int16": pd.Int16Dtype(),
    "Int32": pd.Int32Dtype(),
    "Int64": pd.Int64Dtype(),
    "UInt8": pd.UInt8Dtype(),
    "UInt16": pd.UInt16Dtype(),
    "UInt32": pd.UInt32Dtype(),
    "UInt64": pd.UInt64Dtype(),
    "boolean": pd.BooleanDtype()
}


def _logical_to_time_dtype(logical_timestamp_type):
    if getattr(logical_timestamp_type.unit, "NANOS", None) is not None:
        unit = "ns"
    elif getattr(logical_timestamp_type.unit, "MICROS", None) is not None:
        unit = "us"
    elif getattr(logical_timestamp_type.unit, "MILLIS", None) is not None:
        unit = "ms"
    else:
        raise ValueError("Timestamp ")

    return np.dtype(f"<M8[{unit}]")


def typemap(se, md=None):
    """Get the final dtype - no actual conversion"""
    md = md or {}
    md = md.get(se.name, {})
    if md and ("Int" in md["numpy_type"] or md["numpy_type"] == "boolean"):
        # arrow has numpy and pandas types swapped
        return pandas_nullable[md["numpy_type"]]
    if md and ("Int" in md["pandas_type"] or md["pandas_type"] == "boolean"):
        return pandas_nullable[md["pandas_type"]]
    if se.logicalType is not None and se.logicalType.TIMESTAMP is not None:
        return _logical_to_time_dtype(se.logicalType.TIMESTAMP)
    if se.converted_type is None:
        if se.type in simple:
            return simple[se.type]
        else:
            return np.dtype("S%i" % se.type_length)
    if md and "time" in md.get("numpy_type", ""):
        return np.dtype(md["numpy_type"])
    if se.converted_type in complex:
        return complex[se.converted_type]
    return np.dtype("O")


def converts_inplace(se):
    """when converting, reuses input array"""
    if se.type == parquet_thrift.Type.BOOLEAN:
        return False  # always needs unpacking
    ctype = se.converted_type
    if ctype is None:
        return True
    if se.type == parquet_thrift.Type.BYTE_ARRAY:
        return ctype == parquet_thrift.ConvertedType.UTF8
    if ctype in [
        parquet_thrift.ConvertedType.DATE,
        parquet_thrift.ConvertedType.TIME_MILLIS,
        parquet_thrift.ConvertedType.TIMESTAMP_MILLIS,
        parquet_thrift.ConvertedType.TIME_MICROS,
        parquet_thrift.ConvertedType.TIMESTAMP_MICROS
    ]:
        return True
    if getattr(se.logicalType, "TIMESTAMP", None) is not None:
        # this will be nanos, since micro and milli hit block above
        return True
    return False


def convert(data, se, timestamp96=True, dtype=None):
    """Convert known types from primitive to rich.

    Parameters
    ----------
    data: pandas series of primitive type
    se: a schema element.
    timestamp96: convert int96 as if it were written by mr-parquet
    """
    ctype = se.converted_type
    if se.type == parquet_thrift.Type.INT96 and timestamp96:
        data2 = data.view([('ns', 'i8'), ('day', 'i4')])
        # TODO: this should be ms unit, now that we can?
        return ((data2['day'] - np.int64(2440588)) * DAYS_TO_NANOS +
                data2['ns']).view('M8[ns]')
    if se.logicalType is not None and se.logicalType.TIMESTAMP is not None:
        dt = _logical_to_time_dtype(se.logicalType.TIMESTAMP)
        return data.view(dt)
    if ctype is None:
        return data
    if ctype == parquet_thrift.ConvertedType.UTF8:
        if data.dtype != "O" or (len(data) == 1 and not isinstance(data[0], str)):
            # fixed string
            import pandas as pd
            return pd.Series(data).str.decode("utf8").values
        # already converted in speedups.unpack_byte_array
        return data
    if ctype == parquet_thrift.ConvertedType.DECIMAL:
        scale_factor = 10**-se.scale
        if data.dtype.kind in ['i', 'f']:
            return data * scale_factor
        else:  # byte-string
            # NB: general but slow method
            # could optimize when data.dtype.itemsize <= 8
            # TODO: easy cythonize (but rare)
            # TODO: extension point for pandas-decimal (no conversion needed)
            return np.array([
                int.from_bytes(
                    data.data[i:i + 1], byteorder='big', signed=True
                ) * scale_factor
                for i in range(len(data))
            ])
    elif ctype == parquet_thrift.ConvertedType.DATE:
        data = data * DAYS_TO_NANOS
        return data.view('datetime64[ns]')
    elif ctype == parquet_thrift.ConvertedType.TIME_MILLIS:
        # this was not covered by new pandas time units
        data = data.astype('int64', copy=False)
        time_shift(data, 1000000)
        return data.view('timedelta64[ns]')
    elif ctype == parquet_thrift.ConvertedType.TIMESTAMP_MILLIS:
        return data.view('datetime64[ms]')
    elif ctype == parquet_thrift.ConvertedType.TIME_MICROS:
        return data.view('timedelta64[us]')
    elif ctype == parquet_thrift.ConvertedType.TIMESTAMP_MICROS:
        return data.view('datetime64[us]')
    elif ctype == parquet_thrift.ConvertedType.UINT_8:
        # TODO: return strided views?
        #  data.view('uint8')[::data.itemsize].view(out_dtype)
        return data.astype(np.uint8)
    elif ctype == parquet_thrift.ConvertedType.UINT_16:
        return data.astype(np.uint16)
    elif ctype == parquet_thrift.ConvertedType.UINT_32:
        return data.astype(np.uint32)
    elif ctype == parquet_thrift.ConvertedType.UINT_64:
        return data.astype(np.uint64)
    elif ctype == parquet_thrift.ConvertedType.INT_8:
        return data.astype(np.int8)
    elif ctype == parquet_thrift.ConvertedType.INT_16:
        return data.astype(np.int16)
    elif ctype == parquet_thrift.ConvertedType.INT_32:
        return data.astype(np.int32)
    elif ctype == parquet_thrift.ConvertedType.INT_64:
        return data.astype(np.int64)
    elif ctype == parquet_thrift.ConvertedType.JSON:
        if isinstance(data, list) or data.dtype != "O":
            out = np.empty(len(data), dtype="O")
        else:
            out = data
        # TODO: unnecessary list - loop would save memory, and can cythonize
        decoder = json_decoder()
        out[:] = [decoder(d) for d in data]
        return out
    elif ctype == parquet_thrift.ConvertedType.BSON:
        if isinstance(data, list) or data.dtype != "O":
            out = np.empty(len(data), dtype="O")
        else:
            out = data
        # TODO: unnecessary list - loop would save memory, and can cythonize
        #  and could use better BSON lib (bson-numpy, python-bsonjs)?
        out[:] = [unbson(d) for d in data]
        return out
    elif ctype == parquet_thrift.ConvertedType.INTERVAL:
        # for those that understand, output is month, day, ms
        # maybe should convert to timedelta
        return data.view('<u4').reshape((len(data), -1))
    else:
        logger.info("Converted type '%s'' not handled",
                    parquet_thrift.ConvertedType._VALUES_TO_NAMES[ctype])  # pylint:disable=protected-access
    return data
