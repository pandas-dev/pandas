# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.


import ast
from collections.abc import Sequence
from concurrent import futures
# import threading submodule upfront to avoid partially initialized
# module bug (ARROW-11983)
import concurrent.futures.thread  # noqa
from copy import deepcopy
import decimal
from itertools import zip_longest
import json
import operator
import re
import warnings

try:
    import numpy as np
except ImportError:
    np = None
import pyarrow as pa
from pyarrow.lib import _pandas_api, frombytes, is_threading_enabled  # noqa


_logical_type_map = {}
_numpy_logical_type_map = {}
_pandas_logical_type_map = {}


def get_logical_type_map():
    global _logical_type_map

    if not _logical_type_map:
        _logical_type_map.update({
            pa.lib.Type_NA: 'empty',
            pa.lib.Type_BOOL: 'bool',
            pa.lib.Type_INT8: 'int8',
            pa.lib.Type_INT16: 'int16',
            pa.lib.Type_INT32: 'int32',
            pa.lib.Type_INT64: 'int64',
            pa.lib.Type_UINT8: 'uint8',
            pa.lib.Type_UINT16: 'uint16',
            pa.lib.Type_UINT32: 'uint32',
            pa.lib.Type_UINT64: 'uint64',
            pa.lib.Type_HALF_FLOAT: 'float16',
            pa.lib.Type_FLOAT: 'float32',
            pa.lib.Type_DOUBLE: 'float64',
            pa.lib.Type_DATE32: 'date',
            pa.lib.Type_DATE64: 'date',
            pa.lib.Type_TIME32: 'time',
            pa.lib.Type_TIME64: 'time',
            pa.lib.Type_BINARY: 'bytes',
            pa.lib.Type_FIXED_SIZE_BINARY: 'bytes',
            pa.lib.Type_STRING: 'unicode',
        })
    return _logical_type_map


def get_logical_type(arrow_type):
    logical_type_map = get_logical_type_map()

    try:
        return logical_type_map[arrow_type.id]
    except KeyError:
        if isinstance(arrow_type, pa.lib.DictionaryType):
            return 'categorical'
        elif isinstance(arrow_type, pa.lib.ListType):
            return 'list[{}]'.format(get_logical_type(arrow_type.value_type))
        elif isinstance(arrow_type, pa.lib.TimestampType):
            return 'datetimetz' if arrow_type.tz is not None else 'datetime'
        elif pa.types.is_decimal(arrow_type):
            return 'decimal'
        return 'object'


def get_numpy_logical_type_map():
    global _numpy_logical_type_map
    if not _numpy_logical_type_map:
        _numpy_logical_type_map.update({
            np.bool_: 'bool',
            np.int8: 'int8',
            np.int16: 'int16',
            np.int32: 'int32',
            np.int64: 'int64',
            np.uint8: 'uint8',
            np.uint16: 'uint16',
            np.uint32: 'uint32',
            np.uint64: 'uint64',
            np.float32: 'float32',
            np.float64: 'float64',
            'datetime64[D]': 'date',
            np.str_: 'string',
            np.bytes_: 'bytes',
        })
    return _numpy_logical_type_map


def get_logical_type_from_numpy(pandas_collection):
    numpy_logical_type_map = get_numpy_logical_type_map()
    try:
        return numpy_logical_type_map[pandas_collection.dtype.type]
    except KeyError:
        if hasattr(pandas_collection.dtype, 'tz'):
            return 'datetimetz'
        # See https://github.com/pandas-dev/pandas/issues/24739 (infer_dtype will
        # result in "datetime64" without unit, while pandas astype requires a unit)
        if str(pandas_collection.dtype).startswith('datetime64'):
            return str(pandas_collection.dtype)
        result = _pandas_api.infer_dtype(pandas_collection)
        if result == 'string':
            return 'unicode'
        return result


def get_extension_dtype_info(column):
    dtype = column.dtype
    if str(dtype) == 'category':
        cats = getattr(column, 'cat', column)
        assert cats is not None
        metadata = {
            'num_categories': len(cats.categories),
            'ordered': cats.ordered,
        }
        physical_dtype = str(cats.codes.dtype)
    elif hasattr(dtype, 'tz'):
        metadata = {'timezone': pa.lib.tzinfo_to_string(dtype.tz)}
        physical_dtype = 'datetime64[ns]'
    else:
        metadata = None
        physical_dtype = str(dtype)
    return physical_dtype, metadata


def get_column_metadata(column, name, arrow_type, field_name):
    """Construct the metadata for a given column

    Parameters
    ----------
    column : pandas.Series or pandas.Index
    name : str
    arrow_type : pyarrow.DataType
    field_name : str
        Equivalent to `name` when `column` is a `Series`, otherwise if `column`
        is a pandas Index then `field_name` will not be the same as `name`.
        This is the name of the field in the arrow Table's schema.

    Returns
    -------
    dict
    """
    logical_type = get_logical_type(arrow_type)

    string_dtype, extra_metadata = get_extension_dtype_info(column)
    if logical_type == 'decimal':
        extra_metadata = {
            'precision': arrow_type.precision,
            'scale': arrow_type.scale,
        }
        string_dtype = 'object'

    if (
        name is not None
        and not (isinstance(name, float) and np.isnan(name))
        and not isinstance(name, str)
    ):
        raise TypeError(
            'Column name must be a string. Got column {} of type {}'.format(
                name, type(name).__name__
            )
        )

    assert isinstance(field_name, str), str(type(field_name))
    return {
        'name': name,
        'field_name': field_name,
        'pandas_type': logical_type,
        'numpy_type': string_dtype,
        'metadata': extra_metadata,
    }


def construct_metadata(columns_to_convert, df, column_names, index_levels,
                       index_descriptors, preserve_index, types,
                       column_field_names=None):
    """Returns a dictionary containing enough metadata to reconstruct a pandas
    DataFrame as an Arrow Table, including index columns.

    Parameters
    ----------
    columns_to_convert : list[pd.Series]
    df : pandas.DataFrame
    column_names : list[str | None]
    column_field_names: list[str]
    index_levels : List[pd.Index]
    index_descriptors : List[Dict]
    preserve_index : bool
    types : List[pyarrow.DataType]

    Returns
    -------
    dict
    """
    if column_field_names is None:
        # backwards compatibility for external projects that are using
        # `construct_metadata` such as cudf
        # see https://github.com/apache/arrow/pull/44963#discussion_r1875771953
        column_field_names = [str(name) for name in column_names]

    num_serialized_index_levels = len([descr for descr in index_descriptors
                                       if not isinstance(descr, dict)])
    # Use ntypes instead of Python shorthand notation [:-len(x)] as [:-0]
    # behaves differently to what we want.
    ntypes = len(types)
    df_types = types[:ntypes - num_serialized_index_levels]
    index_types = types[ntypes - num_serialized_index_levels:]

    column_metadata = []
    for col, name, field_name, arrow_type in zip(columns_to_convert, column_names,
                                                 column_field_names, df_types):
        metadata = get_column_metadata(col, name=name,
                                       arrow_type=arrow_type,
                                       field_name=field_name)
        column_metadata.append(metadata)

    index_column_metadata = []
    if preserve_index is not False:
        non_str_index_names = []
        for level, arrow_type, descriptor in zip(index_levels, index_types,
                                                 index_descriptors):
            if isinstance(descriptor, dict):
                # The index is represented in a non-serialized fashion,
                # e.g. RangeIndex
                continue

            if level.name is not None and not isinstance(level.name, str):
                non_str_index_names.append(level.name)

            metadata = get_column_metadata(
                level,
                name=_column_name_to_strings(level.name),
                arrow_type=arrow_type,
                field_name=descriptor,
            )
            index_column_metadata.append(metadata)

        if len(non_str_index_names) > 0:
            warnings.warn(
                f"The DataFrame has non-str index name `{non_str_index_names}`"
                " which will be converted to string"
                " and not roundtrip correctly.",
                UserWarning, stacklevel=4)

        column_indexes = []

        levels = getattr(df.columns, 'levels', [df.columns])
        names = getattr(df.columns, 'names', [df.columns.name])
        for level, name in zip(levels, names):
            metadata = _get_simple_index_descriptor(level, name)
            column_indexes.append(metadata)
    else:
        index_descriptors = index_column_metadata = column_indexes = []

    return {
        b'pandas': json.dumps({
            'index_columns': index_descriptors,
            'column_indexes': column_indexes,
            'columns': column_metadata + index_column_metadata,
            'creator': {
                'library': 'pyarrow',
                'version': pa.__version__
            },
            'pandas_version': _pandas_api.version
        }).encode('utf8')
    }


def _get_simple_index_descriptor(level, name):
    string_dtype, extra_metadata = get_extension_dtype_info(level)
    pandas_type = get_logical_type_from_numpy(level)
    if 'mixed' in pandas_type:
        warnings.warn(
            "The DataFrame has column names of mixed type. They will be "
            "converted to strings and not roundtrip correctly.",
            UserWarning, stacklevel=4)
    if pandas_type == 'unicode':
        assert not extra_metadata
        extra_metadata = {'encoding': 'UTF-8'}
    return {
        'name': name,
        'field_name': name,
        'pandas_type': pandas_type,
        'numpy_type': string_dtype,
        'metadata': extra_metadata,
    }


def _column_name_to_strings(name):
    """Convert a column name (or level) to either a string or a recursive
    collection of strings.

    Parameters
    ----------
    name : str or tuple

    Returns
    -------
    value : str or tuple

    Examples
    --------
    >>> name = 'foo'
    >>> _column_name_to_strings(name)
    'foo'
    >>> name = ('foo', 'bar')
    >>> _column_name_to_strings(name)
    "('foo', 'bar')"
    >>> import pandas as pd
    >>> name = (1, pd.Timestamp('2017-02-01 00:00:00'))
    >>> _column_name_to_strings(name)
    "('1', '2017-02-01 00:00:00')"
    """
    if isinstance(name, str):
        return name
    elif isinstance(name, bytes):
        # XXX: should we assume that bytes in Python 3 are UTF-8?
        return name.decode('utf8')
    elif isinstance(name, tuple):
        return str(tuple(map(_column_name_to_strings, name)))
    elif isinstance(name, Sequence):
        raise TypeError("Unsupported type for MultiIndex level")
    elif name is None or (isinstance(name, float) and np.isnan(name)):
        return name
    return str(name)


def _index_level_name(index, i, column_names):
    """Return the name of an index level or a default name if `index.name` is
    None or is already a column name.

    Parameters
    ----------
    index : pandas.Index
    i : int

    Returns
    -------
    name : str
    """
    if index.name is not None and index.name not in column_names:
        return _column_name_to_strings(index.name)
    else:
        return '__index_level_{:d}__'.format(i)


def _get_columns_to_convert(df, schema, preserve_index, columns):
    columns = _resolve_columns_of_interest(df, schema, columns)

    if not df.columns.is_unique:
        raise ValueError(
            'Duplicate column names found: {}'.format(list(df.columns))
        )

    if schema is not None:
        return _get_columns_to_convert_given_schema(df, schema, preserve_index)

    column_names = []
    column_field_names = []

    index_levels = (
        _get_index_level_values(df.index) if preserve_index is not False
        else []
    )

    columns_to_convert = []
    convert_fields = []

    for name in columns:
        col = df[name]
        name = _column_name_to_strings(name)

        if _pandas_api.is_sparse(col):
            raise TypeError(
                "Sparse pandas data (column {}) not supported.".format(name))

        columns_to_convert.append(col)
        convert_fields.append(None)
        column_names.append(name)
        column_field_names.append(str(name))

    index_descriptors = []
    index_column_names = []
    for i, index_level in enumerate(index_levels):
        name = _index_level_name(index_level, i, column_names)
        if (isinstance(index_level, _pandas_api.pd.RangeIndex) and
                preserve_index is None):
            descr = _get_range_index_descriptor(index_level)
        else:
            columns_to_convert.append(index_level)
            convert_fields.append(None)
            descr = name
            index_column_names.append(name)
        index_descriptors.append(descr)

    all_names = column_field_names + index_column_names

    # all_names : all of the columns in the resulting table including the data
    # columns and serialized index columns
    # column_names : the names of the data columns
    # index_column_names : the names of the serialized index columns
    # index_descriptors : descriptions of each index to be used for
    # reconstruction
    # index_levels : the extracted index level values
    # columns_to_convert : assembled raw data (both data columns and indexes)
    # to be converted to Arrow format
    # columns_fields : specified column to use for coercion / casting
    # during serialization, if a Schema was provided
    return (all_names, column_names, column_field_names, index_column_names,
            index_descriptors, index_levels, columns_to_convert, convert_fields)


def _get_columns_to_convert_given_schema(df, schema, preserve_index):
    """
    Specialized version of _get_columns_to_convert in case a Schema is
    specified.
    In that case, the Schema is used as the single point of truth for the
    table structure (types, which columns are included, order of columns, ...).
    """
    column_names = []
    columns_to_convert = []
    convert_fields = []
    index_descriptors = []
    index_column_names = []
    index_levels = []

    for name in schema.names:
        try:
            col = df[name]
            is_index = False
        except KeyError:
            try:
                col = _get_index_level(df, name)
            except (KeyError, IndexError):
                # name not found as index level
                raise KeyError(
                    "name '{}' present in the specified schema is not found "
                    "in the columns or index".format(name))
            if preserve_index is False:
                raise ValueError(
                    "name '{}' present in the specified schema corresponds "
                    "to the index, but 'preserve_index=False' was "
                    "specified".format(name))
            elif (preserve_index is None and
                    isinstance(col, _pandas_api.pd.RangeIndex)):
                raise ValueError(
                    "name '{}' is present in the schema, but it is a "
                    "RangeIndex which will not be converted as a column "
                    "in the Table, but saved as metadata-only not in "
                    "columns. Specify 'preserve_index=True' to force it "
                    "being added as a column, or remove it from the "
                    "specified schema".format(name))
            is_index = True

        if _pandas_api.is_sparse(col):
            raise TypeError(
                "Sparse pandas data (column {}) not supported.".format(name))

        field = schema.field(name)
        columns_to_convert.append(col)
        convert_fields.append(field)
        column_names.append(name)

        if is_index:
            index_column_names.append(name)
            index_descriptors.append(name)
            index_levels.append(col)

    all_names = column_names + index_column_names

    return (all_names, column_names, column_names, index_column_names,
            index_descriptors, index_levels, columns_to_convert, convert_fields)


def _get_index_level(df, name):
    """
    Get the index level of a DataFrame given 'name' (column name in an arrow
    Schema).
    """
    key = name
    if name not in df.index.names and _is_generated_index_name(name):
        # we know we have an autogenerated name => extract number and get
        # the index level positionally
        key = int(name[len("__index_level_"):-2])
    return df.index.get_level_values(key)


def _level_name(name):
    # preserve type when default serializable, otherwise str it
    try:
        json.dumps(name)
        return name
    except TypeError:
        return str(name)


def _get_range_index_descriptor(level):
    # public start/stop/step attributes added in pandas 0.25.0
    return {
        'kind': 'range',
        'name': _level_name(level.name),
        'start': _pandas_api.get_rangeindex_attribute(level, 'start'),
        'stop': _pandas_api.get_rangeindex_attribute(level, 'stop'),
        'step': _pandas_api.get_rangeindex_attribute(level, 'step')
    }


def _get_index_level_values(index):
    n = len(getattr(index, 'levels', [index]))
    return [index.get_level_values(i) for i in range(n)]


def _resolve_columns_of_interest(df, schema, columns):
    if schema is not None and columns is not None:
        raise ValueError('Schema and columns arguments are mutually '
                         'exclusive, pass only one of them')
    elif schema is not None:
        columns = schema.names
    elif columns is not None:
        columns = [c for c in columns if c in df.columns]
    else:
        columns = df.columns

    return columns


def dataframe_to_types(df, preserve_index, columns=None):
    (all_names,
     column_names,
     column_field_names,
     _,
     index_descriptors,
     index_columns,
     columns_to_convert,
     _) = _get_columns_to_convert(df, None, preserve_index, columns)

    types = []
    # If pandas knows type, skip conversion
    for c in columns_to_convert:
        values = c.values
        if _pandas_api.is_categorical(values):
            type_ = pa.array(c, from_pandas=True).type
        elif _pandas_api.is_extension_array_dtype(values):
            empty = c.head(0) if isinstance(
                c, _pandas_api.pd.Series) else c[:0]
            type_ = pa.array(empty, from_pandas=True).type
        else:
            values, type_ = get_datetimetz_type(values, c.dtype, None)
            type_ = pa.lib._ndarray_to_arrow_type(values, type_)
            if type_ is None:
                type_ = pa.array(c, from_pandas=True).type
        types.append(type_)

    metadata = construct_metadata(
        columns_to_convert, df, column_names, index_columns, index_descriptors,
        preserve_index, types, column_field_names=column_field_names
    )

    return all_names, types, metadata


def dataframe_to_arrays(df, schema, preserve_index, nthreads=1, columns=None,
                        safe=True):
    (all_names,
     column_names,
     column_field_names,
     index_column_names,
     index_descriptors,
     index_columns,
     columns_to_convert,
     convert_fields) = _get_columns_to_convert(df, schema, preserve_index,
                                               columns)

    # NOTE(wesm): If nthreads=None, then we use a heuristic to decide whether
    # using a thread pool is worth it. Currently the heuristic is whether the
    # nrows > 100 * ncols and ncols > 1.
    if nthreads is None:
        nrows, ncols = len(df), len(df.columns)
        if nrows > ncols * 100 and ncols > 1:
            nthreads = pa.cpu_count()
        else:
            nthreads = 1
    # if we don't have threading in libarrow, don't use threading here either
    if not is_threading_enabled():
        nthreads = 1

    def convert_column(col, field):
        if field is None:
            field_nullable = True
            type_ = None
        else:
            field_nullable = field.nullable
            type_ = field.type

        try:
            result = pa.array(col, type=type_, from_pandas=True, safe=safe)
        except (pa.ArrowInvalid,
                pa.ArrowNotImplementedError,
                pa.ArrowTypeError) as e:
            e.args += ("Conversion failed for column {!s} with type {!s}"
                       .format(col.name, col.dtype),)
            raise e
        if not field_nullable and result.null_count > 0:
            raise ValueError("Field {} was non-nullable but pandas column "
                             "had {} null values".format(str(field),
                                                         result.null_count))
        return result

    def _can_definitely_zero_copy(arr):
        return (isinstance(arr, np.ndarray) and
                arr.flags.contiguous and
                issubclass(arr.dtype.type, np.integer))

    if nthreads == 1:
        arrays = [convert_column(c, f)
                  for c, f in zip(columns_to_convert, convert_fields)]
    else:
        arrays = []
        with futures.ThreadPoolExecutor(nthreads) as executor:
            for c, f in zip(columns_to_convert, convert_fields):
                if _can_definitely_zero_copy(c.values):
                    arrays.append(convert_column(c, f))
                else:
                    arrays.append(executor.submit(convert_column, c, f))

        for i, maybe_fut in enumerate(arrays):
            if isinstance(maybe_fut, futures.Future):
                arrays[i] = maybe_fut.result()

    types = [x.type for x in arrays]

    if schema is None:
        fields = []
        for name, type_ in zip(all_names, types):
            fields.append(pa.field(name, type_))
        schema = pa.schema(fields)

    pandas_metadata = construct_metadata(
        columns_to_convert, df, column_names, index_columns, index_descriptors,
        preserve_index, types, column_field_names=column_field_names
    )
    metadata = deepcopy(schema.metadata) if schema.metadata else dict()
    metadata.update(pandas_metadata)
    schema = schema.with_metadata(metadata)

    # If dataframe is empty but with RangeIndex ->
    # remember the length of the indexes
    n_rows = None
    if len(arrays) == 0:
        try:
            kind = index_descriptors[0]["kind"]
            if kind == "range":
                start = index_descriptors[0]["start"]
                stop = index_descriptors[0]["stop"]
                step = index_descriptors[0]["step"]
                n_rows = len(range(start, stop, step))
        except IndexError:
            pass

    return arrays, schema, n_rows


def get_datetimetz_type(values, dtype, type_):
    if values.dtype.type != np.datetime64:
        return values, type_

    if _pandas_api.is_datetimetz(dtype) and type_ is None:
        # If no user type passed, construct a tz-aware timestamp type
        tz = dtype.tz
        unit = dtype.unit
        type_ = pa.timestamp(unit, tz)
    elif type_ is None:
        # Trust the NumPy dtype
        type_ = pa.from_numpy_dtype(values.dtype)

    return values, type_

# ----------------------------------------------------------------------
# Converting pyarrow.Table efficiently to pandas.DataFrame


def _reconstruct_block(item, columns=None, extension_columns=None, return_block=True):
    """
    Construct a pandas Block from the `item` dictionary coming from pyarrow's
    serialization or returned by arrow::python::ConvertTableToPandas.

    This function takes care of converting dictionary types to pandas
    categorical, Timestamp-with-timezones to the proper pandas Block, and
    conversion to pandas ExtensionBlock

    Parameters
    ----------
    item : dict
        For basic types, this is a dictionary in the form of
        {'block': np.ndarray of values, 'placement': pandas block placement}.
        Additional keys are present for other types (dictionary, timezone,
        object).
    columns :
        Column names of the table being constructed, used for extension types
    extension_columns : dict
        Dictionary of {column_name: pandas_dtype} that includes all columns
        and corresponding dtypes that will be converted to a pandas
        ExtensionBlock.

    Returns
    -------
    pandas Block

    """
    import pandas.core.internals as _int

    block_arr = item.get('block', None)
    placement = item['placement']
    if 'dictionary' in item:
        arr = _pandas_api.categorical_type.from_codes(
            block_arr, categories=item['dictionary'],
            ordered=item['ordered'])
    elif 'timezone' in item:
        unit, _ = np.datetime_data(block_arr.dtype)
        dtype = make_datetimetz(unit, item['timezone'])
        if _pandas_api.is_ge_v21():
            arr = _pandas_api.pd.array(
                block_arr.view("int64"), dtype=dtype, copy=False
            )
        else:
            arr = block_arr
            if return_block:
                block = _int.make_block(block_arr, placement=placement,
                                        klass=_int.DatetimeTZBlock,
                                        dtype=dtype)
                return block
    elif 'py_array' in item:
        # create ExtensionBlock
        arr = item['py_array']
        assert len(placement) == 1
        name = columns[placement[0]]
        pandas_dtype = extension_columns[name]
        if not hasattr(pandas_dtype, '__from_arrow__'):
            raise ValueError("This column does not support to be converted "
                             "to a pandas ExtensionArray")
        arr = pandas_dtype.__from_arrow__(arr)
    else:
        arr = block_arr

    if return_block:
        return _int.make_block(arr, placement=placement)
    else:
        return arr, placement


def make_datetimetz(unit, tz):
    if _pandas_api.is_v1():
        unit = 'ns'  # ARROW-3789: Coerce date/timestamp types to datetime64[ns]
    tz = pa.lib.string_to_tzinfo(tz)
    return _pandas_api.datetimetz_type(unit, tz=tz)


def table_to_dataframe(
    options, table, categories=None, ignore_metadata=False, types_mapper=None
):
    all_columns = []
    column_indexes = []
    pandas_metadata = table.schema.pandas_metadata

    if not ignore_metadata and pandas_metadata is not None:
        all_columns = pandas_metadata['columns']
        column_indexes = pandas_metadata.get('column_indexes', [])
        index_descriptors = pandas_metadata['index_columns']
        table = _add_any_metadata(table, pandas_metadata)
        table, index = _reconstruct_index(table, index_descriptors,
                                          all_columns, types_mapper)
        ext_columns_dtypes = _get_extension_dtypes(
            table, all_columns, types_mapper, options, categories)
    else:
        index = _pandas_api.pd.RangeIndex(table.num_rows)
        ext_columns_dtypes = _get_extension_dtypes(
            table, [], types_mapper, options, categories
        )

    _check_data_column_metadata_consistency(all_columns)
    columns = _deserialize_column_index(table, all_columns, column_indexes)

    column_names = table.column_names
    result = pa.lib.table_to_blocks(options, table, categories,
                                    list(ext_columns_dtypes.keys()))
    if _pandas_api.is_ge_v3():
        from pandas.api.internals import create_dataframe_from_blocks

        blocks = [
            _reconstruct_block(
                item, column_names, ext_columns_dtypes, return_block=False)
            for item in result
        ]
        df = create_dataframe_from_blocks(blocks, index=index, columns=columns)
        return df
    else:
        from pandas.core.internals import BlockManager
        from pandas import DataFrame

        blocks = [
            _reconstruct_block(item, column_names, ext_columns_dtypes)
            for item in result
        ]
        axes = [columns, index]
        mgr = BlockManager(blocks, axes)
        if _pandas_api.is_ge_v21():
            df = DataFrame._from_mgr(mgr, mgr.axes)
        else:
            df = DataFrame(mgr)
        return df


# Set of the string repr of all numpy dtypes that can be stored in a pandas
# dataframe (complex not included since not supported by Arrow)
_pandas_supported_numpy_types = {
    "int8", "int16", "int32", "int64",
    "uint8", "uint16", "uint32", "uint64",
    "float16", "float32", "float64",
    "object", "bool"
}


def _get_extension_dtypes(table, columns_metadata, types_mapper, options, categories):
    """
    Based on the stored column pandas metadata and the extension types
    in the arrow schema, infer which columns should be converted to a
    pandas extension dtype.

    The 'numpy_type' field in the column metadata stores the string
    representation of the original pandas dtype (and, despite its name,
    not the 'pandas_type' field).
    Based on this string representation, a pandas/numpy dtype is constructed
    and then we can check if this dtype supports conversion from arrow.

    """
    strings_to_categorical = options["strings_to_categorical"]
    categories = categories or []

    ext_columns = {}

    # older pandas version that does not yet support extension dtypes
    if _pandas_api.extension_dtype is None:
        return ext_columns

    # use the specified mapping of built-in arrow types to pandas dtypes
    if types_mapper:
        for field in table.schema:
            typ = field.type
            pandas_dtype = types_mapper(typ)
            if pandas_dtype is not None:
                ext_columns[field.name] = pandas_dtype

    # infer from extension type in the schema
    for field in table.schema:
        typ = field.type
        if field.name not in ext_columns and isinstance(typ, pa.BaseExtensionType):
            try:
                pandas_dtype = typ.to_pandas_dtype()
            except NotImplementedError:
                pass
            else:
                ext_columns[field.name] = pandas_dtype

    # infer the extension columns from the pandas metadata
    for col_meta in columns_metadata:
        try:
            name = col_meta['field_name']
        except KeyError:
            name = col_meta['name']
        dtype = col_meta['numpy_type']

        if name not in ext_columns and dtype not in _pandas_supported_numpy_types:
            # pandas_dtype is expensive, so avoid doing this for types
            # that are certainly numpy dtypes
            pandas_dtype = _pandas_api.pandas_dtype(dtype)
            if isinstance(pandas_dtype, _pandas_api.extension_dtype):
                if isinstance(pandas_dtype, _pandas_api.pd.StringDtype):
                    # when the metadata indicate to use the string dtype,
                    # ignore this in case:
                    # - it is specified to convert strings / this column to categorical
                    # - the column itself is dictionary encoded and would otherwise be
                    #   converted to categorical
                    if strings_to_categorical or name in categories:
                        continue
                    try:
                        if pa.types.is_dictionary(table.schema.field(name).type):
                            continue
                    except KeyError:
                        pass
                if hasattr(pandas_dtype, "__from_arrow__"):
                    ext_columns[name] = pandas_dtype

    # for pandas 3.0+, use pandas' new default string dtype
    if _pandas_api.uses_string_dtype() and not strings_to_categorical:
        for field in table.schema:
            if field.name not in ext_columns and (
                pa.types.is_string(field.type)
                or pa.types.is_large_string(field.type)
                or pa.types.is_string_view(field.type)
            ) and field.name not in categories:
                ext_columns[field.name] = _pandas_api.pd.StringDtype(na_value=np.nan)

    return ext_columns


def _check_data_column_metadata_consistency(all_columns):
    # It can never be the case in a released version of pyarrow that
    # c['name'] is None *and* 'field_name' is not a key in the column metadata,
    # because the change to allow c['name'] to be None and the change to add
    # 'field_name' are in the same release (0.8.0)
    assert all(
        (c['name'] is None and 'field_name' in c) or c['name'] is not None
        for c in all_columns
    )


def _deserialize_column_index(block_table, all_columns, column_indexes):
    if all_columns:
        columns_name_dict = {
            c.get('field_name', _column_name_to_strings(c['name'])): c['name']
            for c in all_columns
        }
        columns_values = [
            columns_name_dict.get(name, name) for name in block_table.column_names
        ]
    else:
        columns_values = block_table.column_names

    # Construct the base index
    if len(column_indexes) > 1:
        # If we're passed multiple column indexes then evaluate with
        # ast.literal_eval, since the column index values show up as a list of
        # tuples
        columns = _pandas_api.pd.MultiIndex.from_tuples(
            list(map(ast.literal_eval, columns_values)),
            names=[col_index['name'] for col_index in column_indexes],
        )
    else:
        columns = _pandas_api.pd.Index(
            columns_values, name=column_indexes[0]["name"] if column_indexes else None
        )

    # if we're reconstructing the index
    if len(column_indexes) > 0:
        columns = _reconstruct_columns_from_metadata(columns, column_indexes)

    return columns


def _reconstruct_index(table, index_descriptors, all_columns, types_mapper=None):
    # 0. 'field_name' is the name of the column in the arrow Table
    # 1. 'name' is the user-facing name of the column, that is, it came from
    #    pandas
    # 2. 'field_name' and 'name' differ for index columns
    # 3. We fall back on c['name'] for backwards compatibility
    field_name_to_metadata = {
        c.get('field_name', c['name']): c
        for c in all_columns
    }

    # Build up a list of index columns and names while removing those columns
    # from the original table
    index_arrays = []
    index_names = []
    result_table = table
    for descr in index_descriptors:
        if isinstance(descr, str):
            result_table, index_level, index_name = _extract_index_level(
                table, result_table, descr, field_name_to_metadata, types_mapper)
            if index_level is None:
                # ARROW-1883: the serialized index column was not found
                continue
        elif descr['kind'] == 'range':
            index_name = descr['name']
            index_level = _pandas_api.pd.RangeIndex(descr['start'],
                                                    descr['stop'],
                                                    step=descr['step'],
                                                    name=index_name)
            if len(index_level) != len(table):
                # Possibly the result of munged metadata
                continue
        else:
            raise ValueError("Unrecognized index kind: {}"
                             .format(descr['kind']))
        index_arrays.append(index_level)
        index_names.append(index_name)

    pd = _pandas_api.pd

    # Reconstruct the row index
    if len(index_arrays) > 1:
        index = pd.MultiIndex.from_arrays(index_arrays, names=index_names)
    elif len(index_arrays) == 1:
        index = index_arrays[0]
        if not isinstance(index, pd.Index):
            # Box anything that wasn't boxed above
            index = pd.Index(index, name=index_names[0])
    else:
        index = pd.RangeIndex(table.num_rows)

    return result_table, index


def _extract_index_level(table, result_table, field_name,
                         field_name_to_metadata, types_mapper=None):
    logical_name = field_name_to_metadata[field_name]['name']
    index_name = _backwards_compatible_index_name(field_name, logical_name)
    i = table.schema.get_field_index(field_name)

    if i == -1:
        # The serialized index column was removed by the user
        return result_table, None, None

    col = table.column(i)
    index_level = col.to_pandas(types_mapper=types_mapper)
    index_level.name = None
    result_table = result_table.remove_column(
        result_table.schema.get_field_index(field_name)
    )
    return result_table, index_level, index_name


def _backwards_compatible_index_name(raw_name, logical_name):
    """Compute the name of an index column that is compatible with older
    versions of :mod:`pyarrow`.

    Parameters
    ----------
    raw_name : str
    logical_name : str

    Returns
    -------
    result : str

    Notes
    -----
    * Part of :func:`~pyarrow.pandas_compat.table_to_blockmanager`
    """
    # Part of table_to_blockmanager
    if raw_name == logical_name and _is_generated_index_name(raw_name):
        return None
    else:
        return logical_name


def _is_generated_index_name(name):
    pattern = r'^__index_level_\d+__$'
    return re.match(pattern, name) is not None


def get_pandas_logical_type_map():
    global _pandas_logical_type_map

    if not _pandas_logical_type_map:
        _pandas_logical_type_map.update({
            'date': 'datetime64[D]',
            'datetime': 'datetime64[ns]',
            'datetimetz': 'datetime64[ns]',
            'unicode': 'str',
            'bytes': np.bytes_,
            'string': 'str',
            'integer': np.int64,
            'floating': np.float64,
            'decimal': np.object_,
            'empty': np.object_,
        })
    return _pandas_logical_type_map


def _pandas_type_to_numpy_type(pandas_type):
    """Get the numpy dtype that corresponds to a pandas type.

    Parameters
    ----------
    pandas_type : str
        The result of a call to pandas.lib.infer_dtype.

    Returns
    -------
    dtype : np.dtype
        The dtype that corresponds to `pandas_type`.
    """
    pandas_logical_type_map = get_pandas_logical_type_map()
    try:
        return pandas_logical_type_map[pandas_type]
    except KeyError:
        if 'mixed' in pandas_type:
            # catching 'mixed', 'mixed-integer' and 'mixed-integer-float'
            return np.object_
        return np.dtype(pandas_type)


def _reconstruct_columns_from_metadata(columns, column_indexes):
    """Construct a pandas MultiIndex from `columns` and column index metadata
    in `column_indexes`.

    Parameters
    ----------
    columns : List[pd.Index]
        The columns coming from a pyarrow.Table
    column_indexes : List[Dict[str, str]]
        The column index metadata deserialized from the JSON schema metadata
        in a :class:`~pyarrow.Table`.

    Returns
    -------
    result : MultiIndex
        The index reconstructed using `column_indexes` metadata with levels of
        the correct type.

    Notes
    -----
    * Part of :func:`~pyarrow.pandas_compat.table_to_blockmanager`
    """
    pd = _pandas_api.pd
    # Get levels and labels, and provide sane defaults if the index has a
    # single level to avoid if/else spaghetti.
    levels = getattr(columns, 'levels', None) or [columns]
    labels = getattr(columns, 'codes', None) or [None]

    # Convert each level to the dtype provided in the metadata
    levels_dtypes = [
        (level, col_index.get('pandas_type', str(level.dtype)),
         col_index.get('numpy_type', None))
        for level, col_index in zip_longest(
            levels, column_indexes, fillvalue={}
        )
    ]

    new_levels = []
    encoder = operator.methodcaller('encode', 'UTF-8')

    for level, pandas_dtype, numpy_dtype in levels_dtypes:
        dtype = _pandas_type_to_numpy_type(pandas_dtype)
        # Since our metadata is UTF-8 encoded, Python turns things that were
        # bytes into unicode strings when json.loads-ing them. We need to
        # convert them back to bytes to preserve metadata.
        if dtype == np.bytes_:
            level = level.map(encoder)
        # ARROW-13756: if index is timezone aware DataTimeIndex
        elif pandas_dtype == "datetimetz":
            tz = pa.lib.string_to_tzinfo(
                column_indexes[0]['metadata']['timezone'])
            level = pd.to_datetime(level, utc=True).tz_convert(tz)
            if _pandas_api.is_ge_v3():
                # with pandas 3+, to_datetime returns a unit depending on the string
                # data, so we restore it to the original unit from the metadata
                level = level.as_unit(np.datetime_data(dtype)[0])
        # GH-41503: if the column index was decimal, restore to decimal
        elif pandas_dtype == "decimal":
            level = _pandas_api.pd.Index([decimal.Decimal(i) for i in level])
        elif (
            level.dtype == "str" and numpy_dtype == "object"
            and ("mixed" in pandas_dtype or pandas_dtype in ["unicode", "string"])
        ):
            # the metadata indicate that the original dataframe used object dtype,
            # but ignore this and keep string dtype if:
            # - the original columns used mixed types -> we don't attempt to faithfully
            #   roundtrip in this case, but keep the column names as strings
            # - the original columns were inferred to be strings but stored in object
            #   dtype -> we don't restore the object dtype because all metadata
            #   generated using pandas < 3 will have this case by default, and
            #   for pandas >= 3 we want to use the default string dtype for .columns
            new_levels.append(level)
            continue
        elif level.dtype != dtype:
            level = level.astype(dtype)
        # ARROW-9096: if original DataFrame was upcast we keep that
        if level.dtype != numpy_dtype and pandas_dtype != "datetimetz":
            level = level.astype(numpy_dtype)

        new_levels.append(level)

    if len(new_levels) > 1:
        return pd.MultiIndex(new_levels, labels, names=columns.names)
    else:
        return pd.Index(new_levels[0], dtype=new_levels[0].dtype, name=columns.name)


def _add_any_metadata(table, pandas_metadata):
    modified_columns = {}
    modified_fields = {}

    schema = table.schema

    index_columns = pandas_metadata['index_columns']
    # only take index columns into account if they are an actual table column
    index_columns = [idx_col for idx_col in index_columns
                     if isinstance(idx_col, str)]
    n_index_levels = len(index_columns)
    n_columns = len(pandas_metadata['columns']) - n_index_levels

    # Add time zones
    for i, col_meta in enumerate(pandas_metadata['columns']):

        raw_name = col_meta.get('field_name')
        if not raw_name:
            # deal with metadata written with arrow < 0.8 or fastparquet
            raw_name = col_meta['name']
            if i >= n_columns:
                # index columns
                raw_name = index_columns[i - n_columns]
            if raw_name is None:
                raw_name = 'None'

        idx = schema.get_field_index(raw_name)
        if idx != -1:
            if col_meta['pandas_type'] == 'datetimetz':
                col = table[idx]
                if not isinstance(col.type, pa.lib.TimestampType):
                    continue
                metadata = col_meta['metadata']
                if not metadata:
                    continue
                metadata_tz = metadata.get('timezone')
                if metadata_tz and metadata_tz != col.type.tz:
                    converted = col.to_pandas()
                    tz_aware_type = pa.timestamp('ns', tz=metadata_tz)
                    with_metadata = pa.Array.from_pandas(converted,
                                                         type=tz_aware_type)

                    modified_fields[idx] = pa.field(schema[idx].name,
                                                    tz_aware_type)
                    modified_columns[idx] = with_metadata

    if len(modified_columns) > 0:
        columns = []
        fields = []
        for i in range(len(table.schema)):
            if i in modified_columns:
                columns.append(modified_columns[i])
                fields.append(modified_fields[i])
            else:
                columns.append(table[i])
                fields.append(table.schema[i])
        return pa.Table.from_arrays(columns, schema=pa.schema(fields))
    else:
        return table


# ----------------------------------------------------------------------
# Helper functions used in lib


def make_tz_aware(series, tz):
    """
    Make a datetime64 Series timezone-aware for the given tz
    """
    tz = pa.lib.string_to_tzinfo(tz)
    series = (series.dt.tz_localize('utc')
                    .dt.tz_convert(tz))
    return series
