"""
Functions for preparing various inputs passed to the DataFrame or Series
constructors before passing them to a BlockManager.
"""
from collections import OrderedDict

import numpy as np
import numpy.ma as ma

from pandas._libs import lib
from pandas._libs.tslibs import IncompatibleFrequency
import pandas.compat as compat
from pandas.compat import lmap, lrange, raise_with_traceback

from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar, construct_1d_ndarray_preserving_na,
    construct_1d_object_array_from_listlike, infer_dtype_from_scalar,
    maybe_cast_to_datetime, maybe_cast_to_integer_array, maybe_castable,
    maybe_convert_platform, maybe_infer_to_datetimelike, maybe_upcast)
from pandas.core.dtypes.common import (
    is_categorical_dtype, is_datetime64tz_dtype, is_dtype_equal,
    is_extension_array_dtype, is_extension_type, is_float_dtype,
    is_integer_dtype, is_iterator, is_list_like, is_object_dtype, pandas_dtype)
from pandas.core.dtypes.generic import (
    ABCDataFrame, ABCDatetimeIndex, ABCIndexClass, ABCPandasArray,
    ABCPeriodIndex, ABCSeries, ABCTimedeltaIndex)
from pandas.core.dtypes.missing import isna

from pandas.core import algorithms, common as com
from pandas.core.arrays import Categorical, ExtensionArray, period_array
from pandas.core.index import (
    Index, _get_objs_combined_axis, _union_indexes, ensure_index)
from pandas.core.indexes import base as ibase
from pandas.core.internals import (
    create_block_manager_from_arrays, create_block_manager_from_blocks)
from pandas.core.internals.arrays import extract_array

# ---------------------------------------------------------------------
# BlockManager Interface


def arrays_to_mgr(arrays, arr_names, index, columns, dtype=None):
    """
    Segregate Series based on type and coerce into matrices.

    Needs to handle a lot of exceptional cases.
    """
    # figure out the index, if necessary
    if index is None:
        index = extract_index(arrays)
    else:
        index = ensure_index(index)

    # don't force copy because getting jammed in an ndarray anyway
    arrays = _homogenize(arrays, index, dtype)

    # from BlockManager perspective
    axes = [ensure_index(columns), index]

    return create_block_manager_from_arrays(arrays, arr_names, axes)


def masked_rec_array_to_mgr(data, index, columns, dtype, copy):
    """
    Extract from a masked rec array and create the manager.
    """

    # essentially process a record array then fill it
    fill_value = data.fill_value
    fdata = ma.getdata(data)
    if index is None:
        index = get_names_from_index(fdata)
        if index is None:
            index = ibase.default_index(len(data))
    index = ensure_index(index)

    if columns is not None:
        columns = ensure_index(columns)
    arrays, arr_columns = to_arrays(fdata, columns)

    # fill if needed
    new_arrays = []
    for fv, arr, col in zip(fill_value, arrays, arr_columns):
        mask = ma.getmaskarray(data[col])
        if mask.any():
            arr, fv = maybe_upcast(arr, fill_value=fv, copy=True)
            arr[mask] = fv
        new_arrays.append(arr)

    # create the manager
    arrays, arr_columns = reorder_arrays(new_arrays, arr_columns, columns)
    if columns is None:
        columns = arr_columns

    mgr = arrays_to_mgr(arrays, arr_columns, index, columns, dtype)

    if copy:
        mgr = mgr.copy()
    return mgr


# ---------------------------------------------------------------------
# DataFrame Constructor Interface

def init_ndarray(values, index, columns, dtype=None, copy=False):
    # input must be a ndarray, list, Series, index

    if isinstance(values, ABCSeries):
        if columns is None:
            if values.name is not None:
                columns = [values.name]
        if index is None:
            index = values.index
        else:
            values = values.reindex(index)

        # zero len case (GH #2234)
        if not len(values) and columns is not None and len(columns):
            values = np.empty((0, 1), dtype=object)

    # we could have a categorical type passed or coerced to 'category'
    # recast this to an arrays_to_mgr
    if (is_categorical_dtype(getattr(values, 'dtype', None)) or
            is_categorical_dtype(dtype)):

        if not hasattr(values, 'dtype'):
            values = prep_ndarray(values, copy=copy)
            values = values.ravel()
        elif copy:
            values = values.copy()

        index, columns = _get_axes(len(values), 1, index, columns)
        return arrays_to_mgr([values], columns, index, columns,
                             dtype=dtype)
    elif (is_datetime64tz_dtype(values) or
          is_extension_array_dtype(values)):
        # GH#19157
        if columns is None:
            columns = [0]
        return arrays_to_mgr([values], columns, index, columns,
                             dtype=dtype)

    # by definition an array here
    # the dtypes will be coerced to a single dtype
    values = prep_ndarray(values, copy=copy)

    if dtype is not None:
        if not is_dtype_equal(values.dtype, dtype):
            try:
                values = values.astype(dtype)
            except Exception as orig:
                e = ValueError("failed to cast to '{dtype}' (Exception "
                               "was: {orig})".format(dtype=dtype,
                                                     orig=orig))
                raise_with_traceback(e)

    index, columns = _get_axes(*values.shape, index=index, columns=columns)
    values = values.T

    # if we don't have a dtype specified, then try to convert objects
    # on the entire block; this is to convert if we have datetimelike's
    # embedded in an object type
    if dtype is None and is_object_dtype(values):
        values = maybe_infer_to_datetimelike(values)

    return create_block_manager_from_blocks([values], [columns, index])


def init_dict(data, index, columns, dtype=None):
    """
    Segregate Series based on type and coerce into matrices.
    Needs to handle a lot of exceptional cases.
    """
    if columns is not None:
        from pandas.core.series import Series
        arrays = Series(data, index=columns, dtype=object)
        data_names = arrays.index

        missing = arrays.isnull()
        if index is None:
            # GH10856
            # raise ValueError if only scalars in dict
            index = extract_index(arrays[~missing])
        else:
            index = ensure_index(index)

        # no obvious "empty" int column
        if missing.any() and not is_integer_dtype(dtype):
            if dtype is None or np.issubdtype(dtype, np.flexible):
                # GH#1783
                nan_dtype = object
            else:
                nan_dtype = dtype
            val = construct_1d_arraylike_from_scalar(np.nan, len(index),
                                                     nan_dtype)
            arrays.loc[missing] = [val] * missing.sum()

    else:
        keys = com.dict_keys_to_ordered_list(data)
        columns = data_names = Index(keys)
        # GH#24096 need copy to be deep for datetime64tz case
        # TODO: See if we can avoid these copies
        arrays = [data[k] if not is_datetime64tz_dtype(data[k]) else
                  data[k].copy(deep=True) for k in keys]
    return arrays_to_mgr(arrays, data_names, index, columns, dtype=dtype)


# ---------------------------------------------------------------------

def prep_ndarray(values, copy=True):
    if not isinstance(values, (np.ndarray, ABCSeries, Index)):
        if len(values) == 0:
            return np.empty((0, 0), dtype=object)

        def convert(v):
            return maybe_convert_platform(v)

        # we could have a 1-dim or 2-dim list here
        # this is equiv of np.asarray, but does object conversion
        # and platform dtype preservation
        try:
            if is_list_like(values[0]) or hasattr(values[0], 'len'):
                values = np.array([convert(v) for v in values])
            elif isinstance(values[0], np.ndarray) and values[0].ndim == 0:
                # GH#21861
                values = np.array([convert(v) for v in values])
            else:
                values = convert(values)
        except (ValueError, TypeError):
            values = convert(values)

    else:

        # drop subclass info, do not copy data
        values = np.asarray(values)
        if copy:
            values = values.copy()

    if values.ndim == 1:
        values = values.reshape((values.shape[0], 1))
    elif values.ndim != 2:
        raise ValueError('Must pass 2-d input')

    return values


def _homogenize(data, index, dtype=None):
    oindex = None
    homogenized = []

    for val in data:
        if isinstance(val, ABCSeries):
            if dtype is not None:
                val = val.astype(dtype)
            if val.index is not index:
                # Forces alignment. No need to copy data since we
                # are putting it into an ndarray later
                val = val.reindex(index, copy=False)
        else:
            if isinstance(val, dict):
                if oindex is None:
                    oindex = index.astype('O')

                if isinstance(index, (ABCDatetimeIndex, ABCTimedeltaIndex)):
                    val = com.dict_compat(val)
                else:
                    val = dict(val)
                val = lib.fast_multiget(val, oindex.values, default=np.nan)
            val = sanitize_array(val, index, dtype=dtype, copy=False,
                                 raise_cast_failure=False)

        homogenized.append(val)

    return homogenized


def extract_index(data):
    index = None
    if len(data) == 0:
        index = Index([])
    elif len(data) > 0:
        raw_lengths = []
        indexes = []

        have_raw_arrays = False
        have_series = False
        have_dicts = False

        for val in data:
            if isinstance(val, ABCSeries):
                have_series = True
                indexes.append(val.index)
            elif isinstance(val, dict):
                have_dicts = True
                indexes.append(list(val.keys()))
            elif is_list_like(val) and getattr(val, 'ndim', 1) == 1:
                have_raw_arrays = True
                raw_lengths.append(len(val))

        if not indexes and not raw_lengths:
            raise ValueError('If using all scalar values, you must pass'
                             ' an index')

        if have_series or have_dicts:
            index = _union_indexes(indexes)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError('arrays must all be same length')

            if have_dicts:
                raise ValueError('Mixing dicts with non-Series may lead to '
                                 'ambiguous ordering.')

            if have_series:
                if lengths[0] != len(index):
                    msg = ('array length {length} does not match index '
                           'length {idx_len}'
                           .format(length=lengths[0], idx_len=len(index)))
                    raise ValueError(msg)
            else:
                index = ibase.default_index(lengths[0])

    return ensure_index(index)


def reorder_arrays(arrays, arr_columns, columns):
    # reorder according to the columns
    if (columns is not None and len(columns) and arr_columns is not None and
            len(arr_columns)):
        indexer = ensure_index(arr_columns).get_indexer(columns)
        arr_columns = ensure_index([arr_columns[i] for i in indexer])
        arrays = [arrays[i] for i in indexer]
    return arrays, arr_columns


def get_names_from_index(data):
    has_some_name = any(getattr(s, 'name', None) is not None for s in data)
    if not has_some_name:
        return ibase.default_index(len(data))

    index = lrange(len(data))
    count = 0
    for i, s in enumerate(data):
        n = getattr(s, 'name', None)
        if n is not None:
            index[i] = n
        else:
            index[i] = 'Unnamed {count}'.format(count=count)
            count += 1

    return index


def _get_axes(N, K, index, columns):
    # helper to create the axes as indexes
    # return axes or defaults

    if index is None:
        index = ibase.default_index(N)
    else:
        index = ensure_index(index)

    if columns is None:
        columns = ibase.default_index(K)
    else:
        columns = ensure_index(columns)
    return index, columns


# ---------------------------------------------------------------------
# Conversion of Inputs to Arrays

def to_arrays(data, columns, coerce_float=False, dtype=None):
    """
    Return list of arrays, columns.
    """
    if isinstance(data, ABCDataFrame):
        if columns is not None:
            arrays = [data._ixs(i, axis=1).values
                      for i, col in enumerate(data.columns) if col in columns]
        else:
            columns = data.columns
            arrays = [data._ixs(i, axis=1).values for i in range(len(columns))]

        return arrays, columns

    if not len(data):
        if isinstance(data, np.ndarray):
            columns = data.dtype.names
            if columns is not None:
                return [[]] * len(columns), columns
        return [], []  # columns if columns is not None else []
    if isinstance(data[0], (list, tuple)):
        return _list_to_arrays(data, columns, coerce_float=coerce_float,
                               dtype=dtype)
    elif isinstance(data[0], compat.Mapping):
        return _list_of_dict_to_arrays(data, columns,
                                       coerce_float=coerce_float, dtype=dtype)
    elif isinstance(data[0], ABCSeries):
        return _list_of_series_to_arrays(data, columns,
                                         coerce_float=coerce_float,
                                         dtype=dtype)
    elif isinstance(data[0], Categorical):
        if columns is None:
            columns = ibase.default_index(len(data))
        return data, columns
    elif (isinstance(data, (np.ndarray, ABCSeries, Index)) and
          data.dtype.names is not None):

        columns = list(data.dtype.names)
        arrays = [data[k] for k in columns]
        return arrays, columns
    else:
        # last ditch effort
        data = lmap(tuple, data)
        return _list_to_arrays(data, columns, coerce_float=coerce_float,
                               dtype=dtype)


def _list_to_arrays(data, columns, coerce_float=False, dtype=None):
    if len(data) > 0 and isinstance(data[0], tuple):
        content = list(lib.to_object_array_tuples(data).T)
    else:
        # list of lists
        content = list(lib.to_object_array(data).T)
    return _convert_object_array(content, columns, dtype=dtype,
                                 coerce_float=coerce_float)


def _list_of_series_to_arrays(data, columns, coerce_float=False, dtype=None):
    if columns is None:
        columns = _get_objs_combined_axis(data, sort=False)

    indexer_cache = {}

    aligned_values = []
    for s in data:
        index = getattr(s, 'index', None)
        if index is None:
            index = ibase.default_index(len(s))

        if id(index) in indexer_cache:
            indexer = indexer_cache[id(index)]
        else:
            indexer = indexer_cache[id(index)] = index.get_indexer(columns)

        values = com.values_from_object(s)
        aligned_values.append(algorithms.take_1d(values, indexer))

    values = np.vstack(aligned_values)

    if values.dtype == np.object_:
        content = list(values.T)
        return _convert_object_array(content, columns, dtype=dtype,
                                     coerce_float=coerce_float)
    else:
        return values.T, columns


def _list_of_dict_to_arrays(data, columns, coerce_float=False, dtype=None):
    if columns is None:
        gen = (list(x.keys()) for x in data)
        sort = not any(isinstance(d, OrderedDict) for d in data)
        columns = lib.fast_unique_multiple_list_gen(gen, sort=sort)

    # assure that they are of the base dict class and not of derived
    # classes
    data = [(type(d) is dict) and d or dict(d) for d in data]

    content = list(lib.dicts_to_array(data, list(columns)).T)
    return _convert_object_array(content, columns, dtype=dtype,
                                 coerce_float=coerce_float)


def _convert_object_array(content, columns, coerce_float=False, dtype=None):
    if columns is None:
        columns = ibase.default_index(len(content))
    else:
        if len(columns) != len(content):  # pragma: no cover
            # caller's responsibility to check for this...
            raise AssertionError('{col:d} columns passed, passed data had '
                                 '{con} columns'.format(col=len(columns),
                                                        con=len(content)))

    # provide soft conversion of object dtypes
    def convert(arr):
        if dtype != object and dtype != np.object:
            arr = lib.maybe_convert_objects(arr, try_float=coerce_float)
            arr = maybe_cast_to_datetime(arr, dtype)
        return arr

    arrays = [convert(arr) for arr in content]

    return arrays, columns


# ---------------------------------------------------------------------
# Series-Based

def sanitize_index(data, index, copy=False):
    """
    Sanitize an index type to return an ndarray of the underlying, pass
    through a non-Index.
    """

    if index is None:
        return data

    if len(data) != len(index):
        raise ValueError('Length of values does not match length of index')

    if isinstance(data, ABCIndexClass) and not copy:
        pass
    elif isinstance(data, (ABCPeriodIndex, ABCDatetimeIndex)):
        data = data._values
        if copy:
            data = data.copy()

    elif isinstance(data, np.ndarray):

        # coerce datetimelike types
        if data.dtype.kind in ['M', 'm']:
            data = sanitize_array(data, index, copy=copy)

    return data


def sanitize_array(data, index, dtype=None, copy=False,
                   raise_cast_failure=False):
    """
    Sanitize input data to an ndarray, copy if specified, coerce to the
    dtype if specified.
    """
    if dtype is not None:
        dtype = pandas_dtype(dtype)

    if isinstance(data, ma.MaskedArray):
        mask = ma.getmaskarray(data)
        if mask.any():
            data, fill_value = maybe_upcast(data, copy=True)
            data.soften_mask()  # set hardmask False if it was True
            data[mask] = fill_value
        else:
            data = data.copy()

    data = extract_array(data, extract_numpy=True)

    # GH#846
    if isinstance(data, np.ndarray):

        if dtype is not None:
            subarr = np.array(data, copy=False)

            # possibility of nan -> garbage
            if is_float_dtype(data.dtype) and is_integer_dtype(dtype):
                try:
                    subarr = _try_cast(data, True, dtype, copy,
                                       True)
                except ValueError:
                    if copy:
                        subarr = data.copy()
            else:
                subarr = _try_cast(data, True, dtype, copy, raise_cast_failure)
        elif isinstance(data, Index):
            # don't coerce Index types
            # e.g. indexes can have different conversions (so don't fast path
            # them)
            # GH#6140
            subarr = sanitize_index(data, index, copy=copy)
        else:

            # we will try to copy be-definition here
            subarr = _try_cast(data, True, dtype, copy, raise_cast_failure)

    elif isinstance(data, ExtensionArray):
        if isinstance(data, ABCPandasArray):
            # We don't want to let people put our PandasArray wrapper
            # (the output of Series/Index.array), into a Series. So
            # we explicitly unwrap it here.
            subarr = data.to_numpy()
        else:
            subarr = data

        # everything else in this block must also handle ndarray's,
        # becuase we've unwrapped PandasArray into an ndarray.

        if dtype is not None:
            subarr = data.astype(dtype)

        if copy:
            subarr = data.copy()
        return subarr

    elif isinstance(data, (list, tuple)) and len(data) > 0:
        if dtype is not None:
            try:
                subarr = _try_cast(data, False, dtype, copy,
                                   raise_cast_failure)
            except Exception:
                if raise_cast_failure:  # pragma: no cover
                    raise
                subarr = np.array(data, dtype=object, copy=copy)
                subarr = lib.maybe_convert_objects(subarr)

        else:
            subarr = maybe_convert_platform(data)

        subarr = maybe_cast_to_datetime(subarr, dtype)

    elif isinstance(data, range):
        # GH#16804
        arr = np.arange(data.start, data.stop, data.step, dtype='int64')
        subarr = _try_cast(arr, False, dtype, copy, raise_cast_failure)
    else:
        subarr = _try_cast(data, False, dtype, copy, raise_cast_failure)

    # scalar like, GH
    if getattr(subarr, 'ndim', 0) == 0:
        if isinstance(data, list):  # pragma: no cover
            subarr = np.array(data, dtype=object)
        elif index is not None:
            value = data

            # figure out the dtype from the value (upcast if necessary)
            if dtype is None:
                dtype, value = infer_dtype_from_scalar(value)
            else:
                # need to possibly convert the value here
                value = maybe_cast_to_datetime(value, dtype)

            subarr = construct_1d_arraylike_from_scalar(
                value, len(index), dtype)

        else:
            return subarr.item()

    # the result that we want
    elif subarr.ndim == 1:
        if index is not None:

            # a 1-element ndarray
            if len(subarr) != len(index) and len(subarr) == 1:
                subarr = construct_1d_arraylike_from_scalar(
                    subarr[0], len(index), subarr.dtype)

    elif subarr.ndim > 1:
        if isinstance(data, np.ndarray):
            raise Exception('Data must be 1-dimensional')
        else:
            subarr = com.asarray_tuplesafe(data, dtype=dtype)

    # This is to prevent mixed-type Series getting all casted to
    # NumPy string type, e.g. NaN --> '-1#IND'.
    if issubclass(subarr.dtype.type, str):
        # GH#16605
        # If not empty convert the data to dtype
        # GH#19853: If data is a scalar, subarr has already the result
        if not lib.is_scalar(data):
            if not np.all(isna(data)):
                data = np.array(data, dtype=dtype, copy=False)
            subarr = np.array(data, dtype=object, copy=copy)

    if is_object_dtype(subarr.dtype) and dtype != 'object':
        inferred = lib.infer_dtype(subarr, skipna=False)
        if inferred == 'period':
            try:
                subarr = period_array(subarr)
            except IncompatibleFrequency:
                pass

    return subarr


def _try_cast(arr, take_fast_path, dtype, copy, raise_cast_failure):

    # perf shortcut as this is the most common case
    if take_fast_path:
        if maybe_castable(arr) and not copy and dtype is None:
            return arr

    try:
        # GH#15832: Check if we are requesting a numeric dype and
        # that we can convert the data to the requested dtype.
        if is_integer_dtype(dtype):
            subarr = maybe_cast_to_integer_array(arr, dtype)

        subarr = maybe_cast_to_datetime(arr, dtype)
        # Take care in creating object arrays (but iterators are not
        # supported):
        if is_object_dtype(dtype) and (is_list_like(subarr) and
                                       not (is_iterator(subarr) or
                                       isinstance(subarr, np.ndarray))):
            subarr = construct_1d_object_array_from_listlike(subarr)
        elif not is_extension_type(subarr):
            subarr = construct_1d_ndarray_preserving_na(subarr, dtype,
                                                        copy=copy)
    except (ValueError, TypeError):
        if is_categorical_dtype(dtype):
            # We *do* allow casting to categorical, since we know
            # that Categorical is the only array type for 'category'.
            subarr = Categorical(arr, dtype.categories,
                                 ordered=dtype.ordered)
        elif is_extension_array_dtype(dtype):
            # create an extension array from its dtype
            array_type = dtype.construct_array_type()._from_sequence
            subarr = array_type(arr, dtype=dtype, copy=copy)
        elif dtype is not None and raise_cast_failure:
            raise
        else:
            subarr = np.array(arr, dtype=object, copy=copy)
    return subarr
