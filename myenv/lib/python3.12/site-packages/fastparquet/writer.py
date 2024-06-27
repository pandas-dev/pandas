import ast
from copy import copy
import itertools
import json
import os
import struct

import numpy as np
import pandas as pd
from pandas.core.arrays.masked import BaseMaskedDtype

from fastparquet.util import join_path

from fastparquet import parquet_thrift, __version__, cencoding
from fastparquet.api import ParquetFile, partitions, part_ids
from fastparquet.compression import compress_data
from fastparquet.converted_types import tobson
from fastparquet.json import json_encoder
from fastparquet.util import (default_open, default_mkdirs, check_column_names,
                   created_by, get_column_metadata,
                   norm_col_name, path_string, reset_row_idx, get_fs,
                   update_custom_metadata)
from fastparquet.speedups import array_encode_utf8, pack_byte_array
from fastparquet.cencoding import NumpyIO, ThriftObject, from_buffer
from decimal import Decimal

MARKER = b'PAR1'
ROW_GROUP_SIZE = 50_000_000
NaT = np.timedelta64(None).tobytes()  # require numpy version >= 1.7
nat = np.datetime64('NaT').view('int64')

typemap = {  # primitive type, converted type, bit width
    'boolean': (parquet_thrift.Type.BOOLEAN, None, 1),
    'Int32': (parquet_thrift.Type.INT32, None, 32),
    'Int64': (parquet_thrift.Type.INT64, None, 64),
    'Int8': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.INT_8, 8),
    'Int16': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.INT_16, 16),
    'UInt8': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_8, 8),
    'UInt16': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_16, 16),
    'UInt32': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_32, 32),
    'UInt64': (parquet_thrift.Type.INT64, parquet_thrift.ConvertedType.UINT_64, 64),
    'bool': (parquet_thrift.Type.BOOLEAN, None, 1),
    'int32': (parquet_thrift.Type.INT32, None, 32),
    'int64': (parquet_thrift.Type.INT64, None, 64),
    'int8': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.INT_8, 8),
    'int16': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.INT_16, 16),
    'uint8': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_8, 8),
    'uint16': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_16, 16),
    'uint32': (parquet_thrift.Type.INT32, parquet_thrift.ConvertedType.UINT_32, 32),
    'uint64': (parquet_thrift.Type.INT64, parquet_thrift.ConvertedType.UINT_64, 64),
    'float32': (parquet_thrift.Type.FLOAT, None, 32),
    'float64': (parquet_thrift.Type.DOUBLE, None, 64),
    'float16': (parquet_thrift.Type.FLOAT, None, 16),
    'Float32': (parquet_thrift.Type.FLOAT, None, 32),
    'Float64': (parquet_thrift.Type.DOUBLE, None, 64),
    'Float16': (parquet_thrift.Type.FLOAT, None, 16),
}

revmap = {parquet_thrift.Type.INT32: np.int32,
          parquet_thrift.Type.INT64: np.int64,
          parquet_thrift.Type.FLOAT: np.float32,
          parquet_thrift.Type.DOUBLE: np.float64}

pdoptional_to_numpy_typemap = {
    pd.Int8Dtype(): np.int8,
    pd.Int16Dtype(): np.int16,
    pd.Int32Dtype(): np.int32,
    pd.Int64Dtype(): np.int64,
    pd.UInt8Dtype(): np.uint8,
    pd.UInt16Dtype(): np.uint16,
    pd.UInt32Dtype(): np.uint32,
    pd.UInt64Dtype(): np.uint64,
    pd.BooleanDtype(): bool
}


def find_type(data, fixed_text=None, object_encoding=None, times='int64',
              is_index:bool=None):
    """ Get appropriate typecodes for column dtype

    Data conversion do not happen here, see convert().

    The user is expected to transform their data into the appropriate dtype
    before saving to parquet, we will not make any assumptions for them.

    Known types that cannot be represented (must be first converted another
    type or to raw binary): float128, complex

    Parameters
    ----------
    data: pd.Series
    fixed_text: int or None
        For str and bytes, the fixed-string length to use. If None, object
        column will remain variable length.
    object_encoding: None or infer|bytes|utf8|json|bson|bool|int|int32|float
        How to encode object type into bytes. If None, bytes is assumed;
        if 'infer', type is guessed from 10 first non-null values.
    times: 'int64'|'int96'
        Normal integers or 12-byte encoding for timestamps.
    is_index: bool, optional
        Set `True` if column storing a row index, `False` otherwise. Required
        if column name is a tuple (when dataframe managed has a column
        multi-index). In this case, with this flag set `True`, name of columns
        used to store a row index are reset from tuple to simple string.

    Returns
    -------
    - a thrift schema element
    - a thrift typecode to be passed to the column chunk writer
    - converted data (None if convert is False)

    """
    dtype = data.dtype
    logical_type = None
    if dtype.name in typemap:
        type, converted_type, width = typemap[dtype.name]
    elif "S" in str(dtype)[:2] or "U" in str(dtype)[:2]:
        type, converted_type, width = (parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY,
                                       None, dtype.itemsize)
    elif dtype == "O":
        if object_encoding == 'infer':
            object_encoding = infer_object_encoding(data)

        if object_encoding == 'utf8':
            type, converted_type, width = (parquet_thrift.Type.BYTE_ARRAY,
                                           parquet_thrift.ConvertedType.UTF8,
                                           None)
        elif object_encoding in ['bytes', None]:
            type, converted_type, width = (parquet_thrift.Type.BYTE_ARRAY, None,
                                           None)
        elif object_encoding == 'json':
            type, converted_type, width = (parquet_thrift.Type.BYTE_ARRAY,
                                           parquet_thrift.ConvertedType.JSON,
                                           None)
        elif object_encoding == 'bson':
            type, converted_type, width = (parquet_thrift.Type.BYTE_ARRAY,
                                           parquet_thrift.ConvertedType.BSON,
                                           None)
        elif object_encoding == 'bool':
            type, converted_type, width = (parquet_thrift.Type.BOOLEAN, None,
                                           1)
        elif object_encoding == 'int':
            type, converted_type, width = (parquet_thrift.Type.INT64, None,
                                           64)
        elif object_encoding == 'int32':
            type, converted_type, width = (parquet_thrift.Type.INT32, None,
                                           32)
        elif object_encoding == 'float':
            type, converted_type, width = (parquet_thrift.Type.DOUBLE, None,
                                           64)
        elif object_encoding == 'decimal':
            type, converted_type, width = (parquet_thrift.Type.DOUBLE, None,
                                           64)
        else:
            raise ValueError('Object encoding (%s) not one of '
                             'infer|utf8|bytes|json|bson|bool|int|int32|float|decimal' %
                             object_encoding)
        if fixed_text:
            width = fixed_text
            type = parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY
    elif dtype.kind == "M":
        if times == 'int64':
            # output will have the same resolution as original data, for resolution <= ms
            tz = getattr(dtype, "tz", None) is not None
            if "ns" in dtype.str:
                type = parquet_thrift.Type.INT64
                converted_type = None
                logical_type = parquet_thrift.LogicalType(
                    TIMESTAMP=parquet_thrift.TimestampType(
                        isAdjustedToUTC=tz,
                        unit=parquet_thrift.TimeUnit(NANOS=parquet_thrift.NanoSeconds())
                    )
                )
                width = None
            elif "us" in dtype.str:
                type, converted_type, width = (
                    parquet_thrift.Type.INT64,
                    parquet_thrift.ConvertedType.TIMESTAMP_MICROS, None
                )
                logical_type = ThriftObject.from_fields(
                    "LogicalType",
                    TIMESTAMP=ThriftObject.from_fields(
                        "TimestampType",
                        isAdjustedToUTC=True,
                        unit=ThriftObject.from_fields("TimeUnit", MICROS={})
                    )
                )

            else:
                type, converted_type, width = (
                    parquet_thrift.Type.INT64,
                    parquet_thrift.ConvertedType.TIMESTAMP_MILLIS, None
                )
                logical_type = ThriftObject.from_fields(
                    "LogicalType",
                    TIMESTAMP=ThriftObject.from_fields(
                        "TimestampType",
                        isAdjustedToUTC=True,
                        unit=ThriftObject.from_fields("TimeUnit", MILLIS={})
                    )
                )
        elif times == 'int96':
            type, converted_type, width = (parquet_thrift.Type.INT96, None,
                                           None)
        else:
            raise ValueError(
                    "Parameter times must be [int64|int96], not %s" % times)
        # warning removed as irrelevant for most users
        # if hasattr(dtype, 'tz') and str(dtype.tz) != 'UTC':
        #     warnings.warn(
        #         'Coercing datetimes to UTC before writing the parquet file, the timezone is stored in the metadata. '
        #         'Reading back with fastparquet/pyarrow will restore the timezone properly.'
        #     )
    elif dtype.kind == "m":
        type, converted_type, width = (parquet_thrift.Type.INT64,
                                       parquet_thrift.ConvertedType.TIME_MICROS, None)
    elif "string" in str(dtype):
        type, converted_type, width = (parquet_thrift.Type.BYTE_ARRAY,
                                       parquet_thrift.ConvertedType.UTF8,
                                       None)
    else:
        raise ValueError("Don't know how to convert data type: %s" % dtype)
    se = parquet_thrift.SchemaElement(
        name=norm_col_name(data.name, is_index), type_length=width,
        converted_type=converted_type, type=type,
        repetition_type=parquet_thrift.FieldRepetitionType.REQUIRED,
        logicalType=logical_type,
        i32=True
    )
    return se, type


def convert(data, se):
    """Convert data according to the schema encoding"""
    dtype = data.dtype
    type = se.type
    converted_type = se.converted_type
    if dtype.name in typemap:
        if type in revmap:
            out = data.values.astype(revmap[type], copy=False)
        elif type == parquet_thrift.Type.BOOLEAN:
            # TODO: with our own bitpack writer, no need to copy for
            #  the padding
            padded = np.pad(data.values, (0, 8 - (len(data) % 8)),
                                'constant', constant_values=(0, 0))
            out = np.packbits(padded.reshape(-1, 8)[:, ::-1].ravel())
        elif dtype.name in typemap:
            out = data.values
    elif "S" in str(dtype)[:2] or "U" in str(dtype)[:2]:
        out = data.values
    elif dtype == "O":
        # TODO: nullable types
        try:
            if converted_type == parquet_thrift.ConvertedType.UTF8:
                # getattr for new pandas StringArray
                # TODO: to bytes in one step
                out = array_encode_utf8(data)
            elif converted_type == parquet_thrift.ConvertedType.DECIMAL:
                out = data.values.astype(np.float64, copy=False)
            elif converted_type is None:
                if type in revmap:
                    out = data.values.astype(revmap[type], copy=False)
                elif type == parquet_thrift.Type.BOOLEAN:
                    # TODO: with our own bitpack writer, no need to copy for
                    #  the padding
                    padded = np.pad(data.values, (0, 8 - (len(data) % 8)),
                                        'constant', constant_values=(0, 0))
                    out = np.packbits(padded.reshape(-1, 8)[:, ::-1].ravel())
                else:
                    out = data.values
            elif converted_type == parquet_thrift.ConvertedType.JSON:
                encoder = json_encoder()
                # TODO: avoid list. np.fromiter can be used with numpy >= 1.23.0,
                #  but older versions don't support object arrays.
                out = np.array([encoder(x) for x in data], dtype="O")
            elif converted_type == parquet_thrift.ConvertedType.BSON:
                out = data.map(tobson).values
            if type == parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY:
                out = out.astype('S%i' % se.type_length)
        except Exception as e:
            ct = parquet_thrift.ConvertedType._VALUES_TO_NAMES[
                converted_type] if converted_type is not None else None
            raise ValueError('Error converting column "%s" to bytes using '
                             'encoding %s. Original error: '
                             '%s' % (data.name, ct, e))
    elif str(dtype) == "string":
        try:
            if converted_type == parquet_thrift.ConvertedType.UTF8:
                # TODO: into bytes in one step
                out = array_encode_utf8(data)
            elif converted_type is None:
                out = data.values
            if type == parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY:
                out = out.astype('S%i' % se.type_length)
        except Exception as e:  # pragma: no cover
            ct = parquet_thrift.ConvertedType._VALUES_TO_NAMES[
                converted_type] if converted_type is not None else None
            raise ValueError('Error converting column "%s" to bytes using '
                             'encoding %s. Original error: '
                             '%s' % (data.name, ct, e))

    elif converted_type == parquet_thrift.ConvertedType.TIME_MICROS:
        # TODO: shift inplace
        if data.dtype == "m8[ns]":
            out = np.empty(len(data), 'int64')
            time_shift(data.values.view('int64'), out)
        else:
            # assuming ms or us
            out = data.values
    elif type == parquet_thrift.Type.INT96 and dtype.kind == 'M':
        ns_per_day = (24 * 3600 * 1000000000)
        day = data.values.view('int64') // ns_per_day + 2440588
        ns = (data.values.view('int64') % ns_per_day)  # - ns_per_day // 2
        out = np.empty(len(data), dtype=[('ns', 'i8'), ('day', 'i4')])
        out['ns'] = ns
        out['day'] = day
    elif dtype.kind == "M":
        out = data.values.view("int64")
    else:
        raise ValueError("Don't know how to convert data type: %s" % dtype)
    return out


def infer_object_encoding(data):
    """Guess object type from first 10 non-na values by iteration"""
    if data.empty:
        return "utf8"
    t = None
    s = 0
    encs = {
        str: "utf8",
        bytes: "bytes",
        list: "json",
        dict: "json",
        bool: "bool",
        Decimal: "decimal",
        int: "int",
        float: "float",
        np.floating: "float",
        np.str_: "utf8"
    }
    for i in data:
        try:
            if i is None or i is pd.NA or i is pd.NaT or i is np.nan or pd.isna(i):
                continue
        except (ValueError, TypeError):
            pass
        tt = type(i)
        if tt in encs:
            tt = encs[tt]
            if t is None:
                t = tt
            elif t != tt:
                raise ValueError("Can't infer object conversion type: %s" % data)
            s += 1
        else:
            raise ValueError("Can't infer object conversion type: %s" % data)
        if s > 10:
            break
    return t


def time_shift(indata, outdata, factor=1000):
    outdata.view("int64")[:] = np.where(
        indata.view('int64') == nat,
        nat,
        indata.view('int64') // factor
    )


def encode_plain(data, se):
    """PLAIN encoding; returns byte representation"""
    out = convert(data, se)
    if se.type == parquet_thrift.Type.BYTE_ARRAY:
        return pack_byte_array(list(out))
    else:
        return out.tobytes()


def encode_dict(data, _):
    """ The data part of dictionary encoding is always int8/16, with RLE/bitpack
    """
    width = data.values.dtype.itemsize * 8
    buf = np.empty(10, dtype=np.uint8)
    o = NumpyIO(buf)
    o.write_byte(width)
    bit_packed_count = (len(data) + 7) // 8
    cencoding.encode_unsigned_varint(bit_packed_count << 1 | 1, o)  # write run header
    # TODO: `bytes`, `tobytes` makes copy, and adding bytes also makes copy
    return bytes(o.so_far()) + data.values.tobytes()


encode = {
    'PLAIN': encode_plain,
    'RLE_DICTIONARY': encode_dict,
}


def make_definitions(data, no_nulls, datapage_version=1):
    """For data that can contain NULLs, produce definition levels binary
    data: either bitpacked bools, or (if number of nulls == 0), single RLE
    block."""
    buf = np.empty(10, dtype=np.uint8)
    temp = NumpyIO(buf)

    if no_nulls:
        # no nulls at all
        l = len(data)
        cencoding.encode_unsigned_varint(l << 1, temp)
        temp.write_byte(1)
        if datapage_version == 1:
            # TODO: adding bytes causes copy
            block = struct.pack('<I', temp.tell()) + temp.so_far()
        else:
            block = bytes(temp.so_far())
        out = data
    else:
        se = parquet_thrift.SchemaElement(type=parquet_thrift.Type.BOOLEAN)
        dnn = data.notnull()
        out = encode_plain(dnn, se)

        cencoding.encode_unsigned_varint(len(out) << 1 | 1, temp)
        head = temp.so_far()

        # TODO: adding bytes causes copy
        if datapage_version == 1:
            block = struct.pack('<I', len(head) + len(out)) + head + out
        else:
            # no need to write length, it's in the header
            # head.write(out)?
            block = bytes(head) + out
        out = data[dnn]
    return block, out


DATAPAGE_VERSION = 2 if os.environ.get("FASTPARQUET_DATAPAGE_V2", False) else 1
MAX_PAGE_SIZE = 500 * 2**20


def _rows_per_page(data, selement, has_nulls=True, page_size=None):
    page_size = page_size or MAX_PAGE_SIZE
    if isinstance(data.dtype, pd.CategoricalDtype):
        bytes_per_element = data.cat.codes.dtype.itemsize
    elif selement.type == parquet_thrift.Type.BOOLEAN:
        bytes_per_element = 0.125
    elif selement.type == parquet_thrift.Type.INT64:
        bytes_per_element = 8
    elif selement.type == parquet_thrift.Type.INT32:
        bytes_per_element = 4
    elif isinstance(data.dtype, BaseMaskedDtype) and data.dtype in pdoptional_to_numpy_typemap:
        bytes_per_element = np.dtype(pdoptional_to_numpy_typemap[data.dtype]).itemsize
    elif data.dtype == "object" or str(data.dtype) == "string":
        dd = data.iloc[:1000]
        d2 = dd[dd.notnull()]
        try:
            sample = d2.str.len()
            chrs = sample.sum()
            bytes_per_element = chrs / (len(sample) or 4) + 4
        except AttributeError:
            # could not apply str to this type of object
            # this estimate is probably grossly wrong
            bytes_per_element = (dd.memory_usage(index=False, deep=True) / (len(dd) or 1)) or 16
    else:
        bytes_per_element = data.dtype.itemsize

    return int(page_size // (bytes_per_element + 0.125 * has_nulls))


def write_column(f, data0, selement, compression=None, datapage_version=None,
                 stats=True):
    """
    Write a single column of data to an open Parquet file

    Parameters
    ----------
    f: open binary file
    data0: pandas Series
    selement: thrift SchemaElement
        produced by ``find_type``
    compression: str, dict, or None
        if ``str``, must be one of the keys in ``compression.compress``
        if ``dict``, must have key ``"type"`` which specifies the compression
        type to use, which must be one of the keys in ``compression.compress``,
        and may optionally have key ``"args`` which should be a dictionary of
        options to pass to the underlying compression engine.
    datapage_version: None or int
        Uses data-page version 1. If 2, uses v2. If None (default), given by values
        of global DATAPAGE_VERSION (set by environment variable FASTPARQUET_DATAPAGE_V2
        at import time).
    stats: bool
        Whether to calculate and write summary statistics

    Returns
    -------
    chunk: ColumnChunk structure

    """
    datapage_version = datapage_version or DATAPAGE_VERSION
    has_nulls = selement.repetition_type == parquet_thrift.FieldRepetitionType.OPTIONAL
    tot_rows = len(data0)
    encoding = "PLAIN"
    first_page = True
    cats = False
    name = data0.name
    diff = 0
    max, min = None, None
    column_chunk_start = f.tell()
    rows_per_page = _rows_per_page(data0, selement, has_nulls)
    row_offsets = list(range(0, len(data0), rows_per_page)) + [len(data0)]
    global_num_nulls = 0
    dict_page_offset = None
    data_page_offset = column_chunk_start

    # column global stats
    if isinstance(data0.dtype, pd.CategoricalDtype) and stats:
        try:
            dnnu = data0.unique().as_ordered()
            max, min = dnnu.max(), dnnu.min()
            if pd.isna(max):
                stats = False
            else:
                if selement.type == parquet_thrift.Type.BYTE_ARRAY:
                    if selement.converted_type is not None:
                        max = encode['PLAIN'](pd.Series([max]), selement)[4:]
                        min = encode['PLAIN'](pd.Series([min]), selement)[4:]
                else:
                    max = encode['PLAIN'](pd.Series([max]), selement)
                    min = encode['PLAIN'](pd.Series([min]), selement)
        except (TypeError, ValueError):
            stats = False
    elif stats:
        try:
            max, min = data0.max(), data0.min()
            if pd.isna(max):
                stats = False
            else:
                if selement.type == parquet_thrift.Type.BYTE_ARRAY:
                    if selement.converted_type is not None:
                        # max = max.encode("utf8") ?
                        max = encode['PLAIN'](pd.Series([max], name=name), selement)[4:]
                        min = encode['PLAIN'](pd.Series([min], name=name), selement)[4:]
                else:
                    max = encode['PLAIN'](pd.Series([max], name=name, dtype=data0.dtype), selement)
                    min = encode['PLAIN'](pd.Series([min], name=name, dtype=data0.dtype), selement)
        except (TypeError, ValueError):
            stats = False

    for row_start, row_end in zip(row_offsets[:-1], row_offsets[1:]):
        data = data0.iloc[row_start:row_end]
        if has_nulls:
            if isinstance(data.dtype, pd.CategoricalDtype):
                num_nulls = (data.cat.codes == -1).sum()
            else:
                num_nulls = len(data) - data.count()
            definition_data, data = make_definitions(data, num_nulls == 0, datapage_version=datapage_version)
            # make_definitions returns `data` with all nulls dropped
            # the null-stripped `data` can be converted from Optional Types to
            # their numpy counterparts
            if isinstance(data.dtype, BaseMaskedDtype) and data.dtype in pdoptional_to_numpy_typemap:
                data = data.astype(pdoptional_to_numpy_typemap[data.dtype], copy=False)
            if data.dtype.kind == "O" and not isinstance(data.dtype, pd.CategoricalDtype):
                try:
                    if selement.type == parquet_thrift.Type.INT64:
                        data = data.astype("int64", copy=False)
                    elif selement.type == parquet_thrift.Type.INT32:
                        data = data.astype("int32", copy=False)
                    elif selement.type == parquet_thrift.Type.BOOLEAN:
                        data = data.astype(bool, copy=False)
                except ValueError as e:
                    t = parquet_thrift.Type._VALUES_TO_NAMES[selement.type]
                    raise ValueError('Error converting column "%s" to primitive '
                                     'type %s. Original error: '
                                     '%s' % (data.name, t, e))
        else:
            definition_data = b""
            num_nulls = 0
        num_nulls = int(num_nulls)
        global_num_nulls += num_nulls

        # No nested field handling (encode those as J/BSON)
        repetition_data = b""

        if isinstance(data.dtype, pd.CategoricalDtype):
            if first_page:
                # make "index page"
                dict_page_offset = column_chunk_start
                dph = parquet_thrift.DictionaryPageHeader(
                    num_values=check_32(len(data.cat.categories)),
                    encoding=parquet_thrift.Encoding.PLAIN,
                    i32=1
                )
                bdata = encode['PLAIN'](pd.Series(data.cat.categories), selement)
                l0 = len(bdata)
                if compression and compression.upper() != "UNCOMPRESSED":
                    bdata = compress_data(bdata, compression)
                    l1 = len(bdata)
                else:
                    l1 = l0
                diff += l0 - l1
                ph = parquet_thrift.PageHeader(
                        type=parquet_thrift.PageType.DICTIONARY_PAGE,
                        uncompressed_page_size=check_32(l0),
                        compressed_page_size=check_32(l1),
                        dictionary_page_header=dph, crc=None, i32=1)

                write_thrift(f, ph)
                f.write(bdata)
                data_page_offset = f.tell()
                ncats = len(data.cat.categories)
                dcat = data.cat.categories.dtype
                cats = True
                encoding = "RLE_DICTIONARY"
            data = data.cat.codes
        if str(data0.dtype) in ['int8', 'int16', 'uint8', 'uint16']:
            # PLAIN encoding must be upcast to parquet primitive
            data = data.astype('int32')

        if datapage_version == 1:
            bdata = b"".join([
                repetition_data, definition_data, encode[encoding](data, selement), 8 * b'\x00'
            ])
            dph = parquet_thrift.DataPageHeader(
                num_values=check_32(row_end - row_start),
                encoding=getattr(parquet_thrift.Encoding, encoding),
                definition_level_encoding=parquet_thrift.Encoding.RLE,
                repetition_level_encoding=parquet_thrift.Encoding.BIT_PACKED,
                i32=1
            )
            l0 = len(bdata)

            if compression:
                bdata = compress_data(bdata, compression)
                l1 = len(bdata)
            else:
                l1 = l0
            diff += l0 - l1

            ph = parquet_thrift.PageHeader(type=parquet_thrift.PageType.DATA_PAGE,
                                           uncompressed_page_size=check_32(l0),
                                           compressed_page_size=check_32(l1),
                                           data_page_header=dph, i32=1)
            write_thrift(f, ph)
            f.write(bdata)
        elif datapage_version == 2:
            is_compressed = isinstance(compression, dict) or (
                compression is not None and compression.upper() != "UNCOMPRESSED")
            dph = parquet_thrift.DataPageHeaderV2(
                num_values=check_32(row_end - row_start),
                num_nulls=check_32(num_nulls),
                num_rows=check_32(row_end - row_start),
                encoding=getattr(parquet_thrift.Encoding, encoding),
                definition_levels_byte_length=len(definition_data),
                repetition_levels_byte_length=0,  # len(repetition_data),
                is_compressed=is_compressed,
                statistics=None,
                i32=1
            )
            bdata = encode[encoding](data, selement)
            lb = len(bdata)
            if is_compressed:
                bdata = compress_data(bdata, compression)
                diff += lb - len(bdata)
            else:
                diff += 0
            ph = parquet_thrift.PageHeader(
                type=parquet_thrift.PageType.DATA_PAGE_V2,
                uncompressed_page_size=check_32(lb + len(definition_data)),
                compressed_page_size=check_32(len(bdata) + len(definition_data)),
                data_page_header_v2=dph, i32=1)
            write_thrift(f, ph)
            # f.write(repetition_data)  # no-op
            f.write(definition_data)
            f.write(bdata)
        first_page = False

    compressed_size = f.tell() - column_chunk_start
    uncompressed_size = compressed_size + diff

    # encoding stats for thrift metadata
    if cats:
        p = [
            parquet_thrift.PageEncodingStats(
                page_type=parquet_thrift.PageType.DICTIONARY_PAGE,
                encoding=parquet_thrift.Encoding.PLAIN, count=1, i32=1),
            parquet_thrift.PageEncodingStats(
                page_type=parquet_thrift.PageType.DATA_PAGE,
                encoding=parquet_thrift.Encoding.RLE_DICTIONARY,
                count=len(row_offsets) - 1, i32=1),
        ]
        encodings = [parquet_thrift.Encoding.PLAIN,
                     parquet_thrift.Encoding.RLE_DICTIONARY]

    else:
        p = [parquet_thrift.PageEncodingStats(
             page_type=parquet_thrift.PageType.DATA_PAGE,
             encoding=parquet_thrift.Encoding.PLAIN,
             count=len(row_offsets) - 1, i32=1)]
        encodings = [parquet_thrift.Encoding.PLAIN]

    if isinstance(compression, dict):
        algorithm = compression.get("type", None)
    else:
        algorithm = compression

    # output thrift metadata
    if stats:
        s = parquet_thrift.Statistics(max=max, min=min, null_count=global_num_nulls)
    else:
        s = parquet_thrift.Statistics(null_count=global_num_nulls)

    kvm = []
    if isinstance(name, (list, tuple)):
        name = str(tuple(name))
    cmd = ThriftObject.from_fields(
        "ColumnMetaData",
        type=selement.type, path_in_schema=[name],
        encodings=encodings,
        codec=(getattr(parquet_thrift.CompressionCodec, algorithm.upper())
               if algorithm else 0),
        num_values=tot_rows,
        statistics=s,
        data_page_offset=data_page_offset,
        dictionary_page_offset=dict_page_offset,
        encoding_stats=p,
        key_value_metadata=kvm,
        total_uncompressed_size=uncompressed_size,
        total_compressed_size=compressed_size,
        i32list=[1, 4]
    )
    if cats:
        kvm.append(
            parquet_thrift.KeyValue(key='num_categories', value=str(ncats)))
        kvm.append(
            parquet_thrift.KeyValue(key='numpy_dtype', value=str(data.dtype)))
        kvm.append(
            parquet_thrift.KeyValue(key='label_dtype', value=str(dcat)))
    chunk = parquet_thrift.ColumnChunk(file_offset=column_chunk_start,
                                       meta_data=cmd,
                                       file_path=None)
    return chunk


class DataFrameSizeWarning(UserWarning):
    pass


def make_row_group(f, data, schema, compression=None, stats=True):
    """ Make a single row group of a Parquet file """
    rows = len(data)
    if rows == 0:
        return
    if isinstance(data.columns, pd.MultiIndex):
        if any(not isinstance(c, (bytes, str)) for c in itertools.chain(*data.columns.values)):
            raise ValueError('Column names must be multi-index, str or bytes:',
                             {c: type(c) for c in data.columns
                              if not isinstance(c, (bytes, str))})

    else:
        if any(not isinstance(c, (bytes, str)) for c in data):
            raise ValueError('Column names must be multi-index, str or bytes:',
                             {c: type(c) for c in data.columns
                              if not isinstance(c, (bytes, str))})

    cols = []
    for column in schema:
        if column.type is not None:
            if isinstance(compression, dict):
                comp = compression.get(column.name, None)
                if comp is None:
                    comp = compression.get('_default', None)
            else:
                comp = compression
            if isinstance(data.columns, pd.MultiIndex):
                try:
                    name = ast.literal_eval(column.name)
                except ValueError:
                    name = column.name
                coldata = data[name]
            else:
                coldata = data[column.name]
            if isinstance(stats, int):
                st = stats
            elif stats == "auto":
                st = coldata.dtype.kind in ["i", "u", "f", "M"]
            else:
                st = column.name in stats
            chunk = write_column(f, coldata, column,
                                 compression=comp, stats=st)
            cols.append(chunk)
    rg = ThriftObject.from_fields(
        "RowGroup", num_rows=rows, columns=cols,
        total_byte_size=sum([c.meta_data.total_uncompressed_size for c in cols]))
    return rg


def make_part_file(f, data, schema, compression=None, fmd=None,
                   stats=True):
    if len(data) == 0:
        return
    with f as f:
        f.write(MARKER)
        rg = make_row_group(f, data, schema, compression=compression,
                            stats=stats)
        if fmd is None:
            fmd = parquet_thrift.FileMetaData(num_rows=rg.num_rows,
                                              schema=schema,
                                              version=1,
                                              created_by=created_by,
                                              row_groups=[rg],
                                              i32list=[1])
            foot_size = write_thrift(f, fmd)
            f.write(struct.pack(b"<I", foot_size))
        else:
            fmd = copy(fmd)
            fmd.row_groups = [rg]
            fmd.num_rows = rg.num_rows
            foot_size = write_thrift(f, fmd)
            f.write(struct.pack(b"<I", foot_size))
        f.write(MARKER)
    return rg


def make_metadata(data, has_nulls=True, ignore_columns=None, fixed_text=None,
                  object_encoding=None, times='int64', index_cols=None, partition_cols=None,
                  cols_dtype="object"):
    if ignore_columns is None:
        ignore_columns = []
    if index_cols is None:
        index_cols = []
    if partition_cols is None:
        partition_cols = []
    if not data.columns.is_unique:
        raise ValueError('Cannot create parquet dataset with duplicate'
                         ' column names (%s)' % data.columns)
    index_cols_orig = None
    if isinstance(data.columns, pd.MultiIndex):
        if isinstance(index_cols, list) and index_cols != []:
            index_cols_orig = copy(index_cols)
            # TODO: for loop required to manage row multi-index.
            name = index_cols[0][0]
            index_cols = [{'field_name': name,
                           'metadata': None,
                           'name': name,
                           'numpy_type': 'object',
                           'pandas_type': 'mixed-integer'}]
        ci = [
            get_column_metadata(ser, n)
            for ser, n
            in zip(data.columns.levels, data.columns.names)
        ]
    else:
        ci = [{'name': data.columns.name,
               'field_name': data.columns.name,
               'pandas_type': 'mixed-integer',
               'numpy_type': str(cols_dtype),
               'metadata': None}]
    if not isinstance(index_cols, list):
        start = index_cols.start
        stop = index_cols.stop
        step = index_cols.step

        index_cols = [{'name': index_cols.name,
                       'start': start,
                       'stop': stop,
                       'step': step,
                       'kind': 'range'}]
    pandas_metadata = {'index_columns': index_cols,
                       'partition_columns': [],
                       'columns': [],
                       'column_indexes': ci,
                       'creator': {'library': 'fastparquet',
                                   'version': __version__},
                       'pandas_version': pd.__version__}
    root = parquet_thrift.SchemaElement(name=b'schema',
                                        num_children=0,
                                        i32=True)
    meta = parquet_thrift.KeyValue(key=b"pandas", value=None)
    fmd = ThriftObject.from_fields("FileMetaData", num_rows=len(data),
                                   schema=None,
                                   version=1,
                                   created_by=created_by.encode(),
                                   row_groups=[],
                                   key_value_metadata=[meta],
                                   i32list=[1])

    object_encoding = object_encoding or {}
    for column in partition_cols:
        pandas_metadata['partition_columns'].append(get_column_metadata(data[column], column))
    schema = [root]
    for column in data.columns:
        if column in ignore_columns:
            continue
        oencoding = (object_encoding if isinstance(object_encoding, str)
                     else object_encoding.get(column, None))
        pandas_metadata['columns'].append(
            get_column_metadata(data[column], column, object_dtype=oencoding))
        fixed = None if fixed_text is None else fixed_text.get(column, None)
        is_index = (column in index_cols_orig) if index_cols_orig else None
        if isinstance(data[column].dtype, pd.CategoricalDtype):
            se, type = find_type(data[column].cat.categories, fixed_text=fixed,
                                 object_encoding=oencoding, is_index=is_index)
            se.name = column
        else:
            se, type = find_type(data[column], fixed_text=fixed,
                                 object_encoding=oencoding, times=times,
                                 is_index=is_index)
        col_has_nulls = has_nulls
        if has_nulls is None:
            se.repetition_type = data[column].dtype == "O"
        elif has_nulls is not True and has_nulls is not False:
            col_has_nulls = column in has_nulls
        if col_has_nulls:
            se.repetition_type = parquet_thrift.FieldRepetitionType.OPTIONAL
        schema.append(se)
        root[5] += 1
    fmd.schema = schema
    meta.value = json.dumps(pandas_metadata, sort_keys=True).encode()
    return fmd


def write_simple(fn, data, fmd, row_group_offsets=None, compression=None,
                 open_with=default_open, has_nulls=None, append=False,
                 stats=True):
    """
    Write to one single file (for file_scheme='simple')

    Parameters
    ----------
    fn: string
        Parquet collection to write to, gathered a single file.
    data: pandas dataframe or iterable of pandas dataframe.
        The table to write. Index of the dataframe is not written.
        If an iterable of dataframe, each one is written as a row group.
    fmd: thrift object
        Parquet file metadata.
    row_group_offsets: int or list of ints,
        If int, row-groups will be approximately this many rows, rounded down
        to make row groups about the same size;
        If a list, the explicit index values to start new row groups;
        If `None`, set to 50_000_000.
    compression:
        Compression to apply to each column, e.g. ``GZIP`` or ``SNAPPY`` or a
        ``dict`` like ``{"col1": "SNAPPY", "col2": None}`` to specify per
        column compression types.
    open_with: function
        When called with a f(path, mode), returns an open file-like object.
    append: bool (False)
        If False, construct data-set from scratch; if True, add new row-group(s)
        to existing data-set. In the latter case, the data-set must exist,
        and the schema must match the input data.
    stats: True|False|list(str)
        Whether to calculate and write summary statistics.
        If True (default), do it for every column;
        if False, never do; and if a list of str, do it only for those
        specified columns.
    """
    if isinstance(data, pd.DataFrame):
        data = iter_dataframe(data, row_group_offsets)
    mode = 'rb+' if append else 'wb'
    if hasattr(fn, "write"):
        of = fn
    else:
        of = open_with(fn, mode)
    with of as f:
        if append:
            f.seek(-8, 2)
            head_size = struct.unpack('<I', f.read(4))[0]
            f.seek(-(head_size+8), 2)
        else:
            f.write(MARKER)
        rgs = fmd.row_groups
        for i, row_group in enumerate(data):
            rg = make_row_group(f, row_group, fmd.schema,
                                compression=compression, stats=stats)
            if rg is not None:
                rgs.append(rg)

        fmd.row_groups = rgs
        fmd.num_rows = sum(rg.num_rows for rg in rgs)
        foot_size = write_thrift(f, fmd)
        f.write(struct.pack(b"<I", foot_size))
        f.write(MARKER)


def write_multi(dn, data, fmd, row_group_offsets=None, compression=None,
                file_scheme='hive', write_fmd=True, open_with=default_open,
                mkdirs=None, partition_on=[], append=False, stats=True):
    """Write to separate parquet files.
    
    Write data following `file_scheme='hive'`, `'drill'` or `'flat'`.
    
    Parameters
    ----------
    dn: string
        Directory path containing the parquet collection to write to.
    data: pandas dataframe or iterable of pandas dataframe.
        The table to write. Index of the dataframe is not written.
        If an iterable of dataframe, each one is written as a row group.
    fmd: thrift object
        Parquet file metadata. `fmd` is modified inplace.
    row_group_offsets: int or list of ints,
        If int, row-groups will be approximately this many rows, rounded down
        to make row groups about the same size;
        If a list, the explicit index values to start new row groups;
        If `None`, set to 50_000_000.
    compression:
        Compression to apply to each column, e.g. ``GZIP`` or ``SNAPPY`` or a
        ``dict`` like ``{"col1": "SNAPPY", "col2": None}`` to specify per
        column compression types.
        By default, do not compress.
        Please, review full description of this parameter in `write` docstring.
    file_scheme: 'hive'|'drill', default 'hive'
        If hive or drill: each row group is in a separate file, and a separate
        file (called "_metadata") contains the metadata.
    write_fmd: bool, default True
        Write updated common metadata to disk.
    open_with: function
        When called with a f(path, mode), returns an open file-like object.
    mkdirs: function
        When called with a path/URL, creates any necessary dictionaries to
        make that location writable, e.g., ``os.makedirs``. This is not
        necessary if using the simple file scheme.
        If not provided, set to 'default_mkdirs'.
    partition_on: list of column names
        Passed to groupby in order to split data within each row-group,
        producing a structured directory tree. Note: as with pandas, null
        values will be dropped. Ignored if `file_scheme` is simple.
    append: bool, default False
        If False, construct dataset from scratch;
        If True, add new row-group(s) to existing dataset.The data-set must
        exist, and the schema must match the input data.
    stats: True|False|list of str
        Whether to calculate and write summary statistics.
        If True (default), do it for every column;
        If False, never do;
        If a list of str, do it only for those specified columns.
    """
    if mkdirs is None:
        mkdirs = default_mkdirs
    if not append:
        # New dataset.
        i_offset = 0
        mkdirs(dn)
    else:
        i_offset = find_max_part(fmd.row_groups)
    if isinstance(data, pd.DataFrame):
        data = iter_dataframe(data, row_group_offsets)
    rg_list = fmd.row_groups
    for i, row_group in enumerate(data):
        part = 'part.%i.parquet' % (i + i_offset)
        if partition_on:
            rgs = partition_on_columns(row_group, partition_on, dn, part,
                                       fmd,  compression, open_with, mkdirs,
                                       with_field=file_scheme == 'hive',
                                       stats=stats)
            rg_list.extend(rgs)
        else:
            partname = join_path(dn, part)
            with open_with(partname, 'wb') as f2:
                rg = make_part_file(f2, row_group, fmd.schema,
                                    compression=compression, fmd=fmd,
                                    stats=stats)
            for chunk in rg.columns:
                chunk.file_path = part
            rg_list.append(rg)
        fmd.row_groups = rg_list
    fmd.num_rows = sum(rg.num_rows for rg in fmd.row_groups)
    if write_fmd:
        write_common_metadata(join_path(dn, '_metadata'), fmd, open_with,
                              no_row_groups=False)
        write_common_metadata(join_path(dn, '_common_metadata'), fmd,
                              open_with)


def iter_dataframe(data, row_group_offsets=None):
    """Yield data in chunk.
    
    Parameters
    ----------
    data : dataframe
        Pandas dataframe to slice.
    row_group_offsets: int or list of ints
        If int, row-groups will be approximately this many rows, rounded down
        to make row groups about the same size;
        If a list, the explicit index values to start new row groups;
        If `None`, set to 50_000_000.

    Yields
    ------
    dataframe
        Chunk of data.
    """
    if row_group_offsets is None:
        row_group_offsets = ROW_GROUP_SIZE
    # TODO
    # Could be extended to accept target size in memory for a row group (MB or
    # GB), instead of target number of rows.
    if isinstance(row_group_offsets, int):
        n_rows = len(data)
        if not row_group_offsets:     # if row group is 0.
            row_group_offsets = [0]
        else:
            nparts = max((n_rows - 1) // row_group_offsets + 1, 1)
            chunksize = max(min((n_rows - 1) // nparts + 1, n_rows), 1)
            row_group_offsets = list(range(0, n_rows, chunksize))
    for i, start in enumerate(row_group_offsets):
        end = (row_group_offsets[i+1] if i < (len(row_group_offsets) - 1)
               else None)
        yield data.iloc[start:end]


def write(filename, data, row_group_offsets=None,
          compression=None, file_scheme='simple', open_with=default_open,
          mkdirs=None, has_nulls=True, write_index=None,
          partition_on=[], fixed_text=None, append=False,
          object_encoding='infer', times='int64',
          custom_metadata=None, stats="auto"):
    """Write pandas dataframe to filename with parquet format.

    Parameters
    ----------
    filename: str or pathlib.Path
        Parquet collection to write to, either a single file (if file_scheme
        is simple) or a directory containing the metadata and data-files.
    data: pandas dataframe
        The table to write.
    row_group_offsets: int or list of int
        If int, row-groups will be approximately this many rows, rounded down
        to make row groups about the same size;
        If a list, the explicit index values to start new row groups;
        If `None`, set to 50_000_000.
        In case of partitioning the data, final row-groups size can be reduced
        significantly further by the partitioning, occuring as a subsequent
        step.
    compression: str, dict
        compression to apply to each column, e.g. ``GZIP`` or ``SNAPPY`` or a
        ``dict`` like ``{"col1": "SNAPPY", "col2": None}`` to specify per
        column compression types.
        In both cases, the compressor settings would be the underlying
        compressor defaults. To pass arguments to the underlying compressor,
        each ``dict`` entry should itself be a dictionary::

            {
                col1: {
                    "type": "LZ4",
                    "args": {
                        "mode": "high_compression",
                        "compression": 9
                     }
                },
                col2: {
                    "type": "SNAPPY",
                    "args": None
                }
                "_default": {
                    "type": "GZIP",
                    "args": None
                }
            }

        where ``"type"`` specifies the compression type to use, and ``"args"``
        specifies a ``dict`` that will be turned into keyword arguments for
        the compressor.
        If the dictionary contains a "_default" entry, this will be used for any
        columns not explicitly specified in the dictionary.
    file_scheme: 'simple'|'hive'|'drill'
        If simple: all goes in a single file
        If hive or drill: each row group is in a separate file, and a separate
        file (called "_metadata") contains the metadata.
    open_with: function
        When called with a f(path, mode), returns an open file-like object
    mkdirs: function
        When called with a path/URL, creates any necessary dictionaries to
        make that location writable, e.g., ``os.makedirs``. This is not
        necessary if using the simple file scheme
    has_nulls: bool, 'infer' or list of strings
        Whether columns can have nulls. If a list of strings, those given
        columns will be marked as "optional" in the metadata, and include
        null definition blocks on disk. Some data types (floats and times)
        can instead use the sentinel values NaN and NaT, which are not the same
        as NULL in parquet, but functionally act the same in many cases,
        particularly if converting back to pandas later. A value of 'infer'
        will assume nulls for object columns and not otherwise.
        Ignored if appending to an existing parquet data-set.
    write_index: boolean
        Whether or not to write the index to a separate column.  By default we
        write the index *if* it is not 0, 1, ..., n.
        Ignored if appending to an existing parquet data-set.
    partition_on: string or list of string
        Column names passed to groupby in order to split data within each
        row-group, producing a structured directory tree. Note: as with pandas,
        null values will be dropped. Ignored if file_scheme is simple.
        Checked when appending to an existing parquet dataset that requested
        partition column names match those of existing parquet data-set.
    fixed_text: {column: int length} or None
        For bytes or str columns, values will be converted
        to fixed-length strings of the given length for the given columns
        before writing, potentially providing a large speed
        boost. The length applies to the binary representation *after*
        conversion for utf8, json or bson.
        Ignored if appending to an existing parquet dataset.
    append: bool (False) or 'overwrite'
        If False, construct data-set from scratch; if True, add new row-group(s)
        to existing data-set. In the latter case, the data-set must exist,
        and the schema must match the input data.

        If 'overwrite', existing partitions will be replaced in-place, where
        the given data has any rows within a given partition. To use this,
        the existing dataset had to be written with these other parameters set
        to specific values, or will raise ValueError:

           *  ``file_scheme='hive'``
           *  ``partition_on`` set to at least one column name.

    object_encoding: str or {col: type}
        For object columns, this gives the data type, so that the values can
        be encoded to bytes.
        Possible values are bytes|utf8|json|bson|bool|int|int32|decimal,
        where bytes is assumed if not specified (i.e., no conversion). The
        special value 'infer' will cause the type to be guessed from the first
        ten non-null values. The decimal.Decimal type is a valid choice, but will
        result in float encoding with possible loss of accuracy.
        Ignored if appending to an existing parquet data-set.
    times: 'int64' (default), or 'int96':
        In "int64" mode, datetimes are written as 8-byte integers, us
        resolution; in "int96" mode, they are written as 12-byte blocks, with
        the first 8 bytes as ns within the day, the next 4 bytes the julian day.
        'int96' mode is included only for compatibility.
        Ignored if appending to an existing parquet data-set.
    custom_metadata: dict
        Key-value metadata to write
        Ignored if appending to an existing parquet data-set.
    stats: True|False|list(str)|"auto"
        Whether to calculate and write summary statistics.
        If True, do it for every column;
        If False, never do;
        And if a list of str, do it only for those specified columns.
        "auto" (default) means True for any int/float or timestamp column

    Examples
    --------
    >>> fastparquet.write('myfile.parquet', df)  # doctest: +SKIP
    """
    custom_metadata = custom_metadata or {}
    if getattr(data, "attrs", None):
        custom_metadata["PANDAS_ATTRS"] = json.dumps(data.attrs)
    if file_scheme not in ('simple', 'hive', 'drill'):
        raise ValueError( 'File scheme should be simple|hive|drill, not '
                         f'{file_scheme}.')
    fs, filename, open_with, mkdirs = get_fs(filename, open_with, mkdirs)

    if append == 'overwrite':
        overwrite(dirpath=filename, data=data,
                  row_group_offsets=row_group_offsets, compression=compression,
                  open_with=open_with, mkdirs=mkdirs, remove_with=None,
                  stats=stats)
        return
    if isinstance(partition_on, str):
        partition_on = [partition_on]
    if append:
        pf = ParquetFile(filename, open_with=open_with)
        if pf._get_index():
            # Format dataframe (manage row index).
            data = reset_row_idx(data)
        if file_scheme == 'simple':
            # Case 'simple'
            if pf.file_scheme not in ['simple', 'empty']:
                raise ValueError( 'File scheme requested is simple, but '
                                 f'existing file scheme is {pf.file_scheme}.')
        else:
            # Case 'hive', 'drill'
            if pf.file_scheme not in ['hive', 'empty', 'flat']:
                raise ValueError(f'Requested file scheme is {file_scheme}, but'
                                  ' existing file scheme is not.')
            if tuple(partition_on) != tuple(pf.cats):
                raise ValueError('When appending, partitioning columns must '
                                 'match existing data')
        pf.write_row_groups(data, row_group_offsets, sort_key=None,
                            sort_pnames=False, compression=compression,
                            write_fmd=True, open_with=open_with,
                            mkdirs=mkdirs, stats=stats)
    else:
        # Case 'append=False'.
        # Define 'index_cols' to be recorded in metadata.
        cols_dtype = data.columns.dtype
        if (write_index or write_index is None
                and not isinstance(data.index, pd.RangeIndex)):
            # Keep name(s) of index to metadata.
            cols = set(data)
            data = reset_row_idx(data)
            index_cols = [c for c in data if c not in cols]
        elif write_index is None and isinstance(data.index, pd.RangeIndex):
            # write_index=None, range to metadata
            index_cols = data.index
        else:
            # write_index=False
            index_cols = []
        # Initialize common metadata.
        if str(has_nulls) == 'infer':
            has_nulls = None
        check_column_names(data.columns, partition_on, fixed_text,
                           object_encoding, has_nulls)
        ignore = partition_on if file_scheme != 'simple' else []
        fmd = make_metadata(data, has_nulls=has_nulls, ignore_columns=ignore,
                            fixed_text=fixed_text,
                            object_encoding=object_encoding,
                            times=times, index_cols=index_cols,
                            partition_cols=partition_on, cols_dtype=cols_dtype)
        if custom_metadata:
            kvm = fmd.key_value_metadata or []
            kvm.extend(
                [
                    parquet_thrift.KeyValue(key=key, value=value)
                    for key, value in custom_metadata.items()
                ]
            )
            fmd.key_value_metadata = kvm

        if file_scheme == 'simple':
            # Case 'simple'
            write_simple(filename, data, fmd,
                         row_group_offsets=row_group_offsets,
                         compression=compression, open_with=open_with,
                         has_nulls=None, append=False, stats=stats)
        else:
            # Case 'hive', 'drill'
            write_multi(filename, data, fmd,
                        row_group_offsets=row_group_offsets,
                        compression=compression, file_scheme=file_scheme,
                        write_fmd=True, open_with=open_with,
                        mkdirs=mkdirs, partition_on=partition_on,
                        append=False, stats=stats)


def find_max_part(row_groups):
    """
    Find the highest integer matching "**part.*.parquet" in referenced paths.
    """
    pids = part_ids(row_groups)
    if pids:
        return max(pids) + 1
    else:
        return 0


def partition_on_columns(data, columns, root_path, partname, fmd,
                         compression, open_with, mkdirs, with_field=True,
                         stats=True):
    """
    Split each row-group by the given columns

    Each combination of column values (determined by pandas groupby) will
    be written in structured directories.
    """
    # Pandas groupby has by default 'sort=True' meaning groups are sorted
    # between them on key.
    gb = data.groupby(columns if len(columns) > 1 else columns[0], observed=False)
    remaining = list(data)
    for column in columns:
        remaining.remove(column)
    if not remaining:
        raise ValueError("Cannot include all columns in partition_on")
    rgs = []
    for key, group in sorted(gb):
        if group.empty:
            continue
        df = group[remaining]
        if not isinstance(key, tuple):
            key = (key,)
        if with_field:
            path = join_path(*(
                "%s=%s" % (name, path_string(val))
                for name, val in zip(columns, key)
            ))
        else:
            path = join_path(*("%s" % val for val in key))
        relname = join_path(path, partname)
        mkdirs(join_path(root_path, path))
        fullname = join_path(root_path, path, partname)
        with open_with(fullname, 'wb') as f2:
            rg = make_part_file(f2, df, fmd.schema,
                                compression=compression, fmd=fmd, stats=stats)
        if rg is not None:
            for chunk in rg.columns:
                chunk.file_path = relname
            rgs.append(rg)
    return rgs


def write_common_metadata(fn, fmd, open_with=default_open,
                          no_row_groups=True):
    """
    For hive-style parquet, write schema in special shared file

    Parameters
    ----------
    fn: str
        Filename to write to
    fmd: thrift FileMetaData
        Information to write
    open_with: func
        To use to create writable file as f(path, mode)
    no_row_groups: bool (True)
        Strip out row groups from metadata before writing - used for "common
        metadata" files, containing only the schema.
    """
    consolidate_categories(fmd)
    with open_with(fn, 'wb') as f:
        f.write(MARKER)
        if no_row_groups:
            fmd = copy(fmd)
            fmd.row_groups = []
            foot_size = write_thrift(f, fmd)
        else:
            foot_size = write_thrift(f, fmd)
        f.write(struct.pack(b"<I", foot_size))
        f.write(MARKER)


def consolidate_categories(fmd):
    key_value = [k for k in fmd.key_value_metadata or []
                 if k.key == b'pandas']
    if not key_value:
        # no pandas categories
        return
    key_value = key_value[0]
    meta = json.loads(key_value.value)
    cats = [c for c in meta['columns']
            if 'num_categories' in (c['metadata'] or [])]
    for cat in cats:
        for rg in fmd.row_groups:
            for col in rg.columns:
                if ".".join(col.meta_data.path_in_schema) == cat['name']:
                    ncats = [k.value for k in (col.meta_data.key_value_metadata or [])
                             if k.key == b'num_categories']
                    if ncats and int(ncats[0]) > cat['metadata'][
                            'num_categories']:
                        cat['metadata']['num_categories'] = int(ncats[0])
    key_value[2] = json.dumps(meta, sort_keys=True).encode()


def merge(file_list, verify_schema=True, open_with=default_open,
          root=False):
    """
    Create a logical data-set out of multiple parquet files.

    The files referenced in file_list must either be in the same directory,
    or at the same level within a structured directory, where the directories
    give partitioning information. The schemas of the files should also be
    consistent.

    Parameters
    ----------
    file_list: list of paths or ParquetFile instances
    verify_schema: bool (True)
        If True, will first check that all the schemas in the input files are
        identical.
    open_with: func
        Used for opening a file for writing as f(path, mode). If input list
        is ParquetFile instances, will be inferred from the first one of these.
    root: str
        If passing a list of files, the top directory of the data-set may
        be ambiguous for partitioning where the upmost field has only one
        value. Use this to specify the data'set root directory, if required.

    Returns
    -------
    ParquetFile instance corresponding to the merged data.
    """
    out = ParquetFile(file_list, verify_schema, open_with, root)
    out._write_common_metadata(open_with)
    return out


def overwrite(dirpath, data, row_group_offsets=None, sort_pnames:bool=True,
              compression=None, open_with=default_open, mkdirs=None,
              remove_with=None, stats=True):
    """Merge new data to existing parquet dataset.

    This function requires existing data on disk, written with 'hive' format.

    This function is a work-in-progress. Several update modes can be envisaged
    and in the mid term, this function will provide a skeleton for achieving
    update of an existing dataset with new data.

    With current version, the only *update mode* supported is
    ``overwrite_partitioned``. With this mode, row-groups on disk that have
    partition values overlapping with those of new data are removed first
    before new data is added.

    Parameters
    ----------
    dirpath : str
        Directory path containing a parquet dataset, written with hive format,
        and with defined partitions.
    data : pandas dataframe
        The table to write.
    row_group_offsets : int or list of int, optional
        If int, row-groups will be approximately this many rows, rounded down
        to make row groups about the same size;
        If a list, the explicit index values to start new row groups;
        If `None`, set to 50_000_000.
        In case of partitioning the data, final row-groups size can be reduced
        significantly further by the partitioning, occuring as a subsequent
        step.
    sort_pnames: bool, default True
        Align name of part files with position of the 1st row group they
        contain.
    compression : str or dict, optional
        Compression to apply to each column, e.g. ``GZIP`` or ``SNAPPY`` or a
        ``dict`` like ``{"col1": "SNAPPY", "col2": None}`` to specify per
        column compression types.
        By default, do not compress.
        Please, review full description of this parameter in `write` docstring.
    open_with : function, optional
        When called with a f(path, mode), returns an open file-like object.
    mkdirs : function, optional
        When called with a path/URL, creates any necessary dictionaries to
        make that location writable, e.g., ``os.makedirs``.
    remove_with : function, optional
        When called with f(path), removes file or directory specified by
        `path` (and any contained files).
    stats: True|False|list of str
        Whether to calculate and write summary statistics.
        If True (default), do it for every column;
        If False, never do;
        And if a list of ``str``, do it only for those specified columns.
    """
    pf = ParquetFile(dirpath, open_with=open_with)
    if (pf.file_scheme == 'simple'
        or (pf.file_scheme == 'empty' and pf.fn[-9:] != '_metadata')):
        raise ValueError('Not possible to overwrite with simple file '
                         'scheme.')
    defined_partitions = list(pf.cats)
    if not defined_partitions:
        raise ValueError('No partitioning column has been set in existing '
                         'dataset. Overwrite of partitions is not possible.')
    # 1st step (from existing data).
    # Define 'sort_key' function to be used to sort all row groups once those
    # of new data will have been added.   
    # 'partitions_starts' is a `dict` that keeps index of 1st row group for
    # each partition in existing data.
    n_rgs = len(pf.row_groups)
    max_idx = n_rgs-1
    partitions_starts = {partitions(rg): (max_idx-i)
                         for i, rg in enumerate(reversed(pf.row_groups))}
    def sort_key(row_group) -> int:
        """Return 1st row-group index with same partition.
        
        If no partition matching, returns an index larger than the 1st
        row-group indexes of any existing partitions.
        """
        # Taking n_rgs (=len(pf.row_groups)) as index for row-groups without
        # matching partition among existing ones is overkill but works.
        rg_partition = partitions(row_group)
        return (partitions_starts[rg_partition]
                if (rg_partition in partitions_starts) else n_rgs)
    # 2nd step (from new and existing data).
    # Identify row groups from existing data with same partition values as
    # those in new data.
    partition_values_in_new = pd.unique(data.loc[:,defined_partitions]
                                            .astype(str).agg('/'.join, axis=1))
    rgs_to_remove = filter(lambda rg : (partitions(rg, True)
                                        in partition_values_in_new),
                           pf.row_groups)
    # 3rd step (on new data).
    # Format new data so that it can be written to disk.
    if pf._get_index():
        # Reset index of pandas dataframe.
        data = reset_row_idx(data)
    # 4th step: write new data, remove previously existing row groups,
    # sort row groups and write updated metadata.
    pf.write_row_groups(data, row_group_offsets=row_group_offsets,
                        sort_key=sort_key, compression=compression,
                        write_fmd=False, open_with=open_with, mkdirs=mkdirs,
                        stats=stats)
    pf.remove_row_groups(rgs_to_remove, sort_pnames=sort_pnames,
                         write_fmd=True, open_with=open_with,
                         remove_with=remove_with)

def write_thrift(f, obj):
    # TODO inline this
    if obj.thrift_name == "FileMetaData":
        for kv in obj.key_value_metadata:
            if not isinstance(kv.key, (bytes, str)):
                raise TypeError(f"KeyValue key expected `str` or `bytes`, got: {kv.key!r}")
            if not isinstance(kv.value, (bytes, str)):
                raise TypeError(f"KeyValue value expected `str` or `bytes`, got: {kv.value!r}")
    return f.write(obj.to_bytes())

def update_file_custom_metadata(path: str, custom_metadata: dict,
                                is_metadata_file: bool = None):
    """Update metadata in file without rewriting data portion if a data file.

    This function updates only the user key-values metadata, not the whole
    metadata of a parquet file.
    Update strategy depends if key found in new custom metadata is also found
    in already existing custom metadata within thrift object, as well as its
    value.
        
      - If not found in existing, it is added.
      - If found in existing, it is updated.
      - If its value is `None`, it is not added, and if found in existing,
        it is removed from existing.

    Parameters
    ----------
    path : str
        Local path to file.
    custom_metadata : dict
        Key-value metadata to update in thrift object.
        The values must be strings or binary. To pass a dictionary, serialize it as json string then encode it in binary.
    is_metadata_file : bool, default None
        Define if target file is a pure metadata file, or is a parquet data
        file. If `None`, is set depending file name.

          - if ending with '_metadata', it assumes file is a metadata file.
          - otherwise, it assumes it is a parquet data file.

    Notes
    -----
    This method does not work for remote files.
    """
    if is_metadata_file is None:
        if path[-9:] == '_metadata':
            is_metadata_file = True
        else:
            is_metadata_file = False
    with open(path, "rb+") as f:
        if is_metadata_file:
            # For pure metadata file, metadata starts just four bytes in.
            loc = 4
        else:
            loc0 = f.seek(-8, 2)
            size = int.from_bytes(f.read(4), "little")
            loc = loc0 - size
        f.seek(loc)
        data = f.read()
        fmd = from_buffer(data, "FileMetaData")
        update_custom_metadata(fmd, custom_metadata)
        f.seek(loc)
        foot_size = write_thrift(f, fmd)
        f.write(struct.pack(b"<I", foot_size))
        f.write(b"PAR1")


def check_32(x):
    if x > 2**31:
        raise OverflowError
    return x
