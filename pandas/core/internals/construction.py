"""
Functions for preparing various inputs passed to the DataFrame or Series
constructors before passing them to a BlockManager.
"""
from collections import abc

import numpy as np
import numpy.ma as ma

from pandas._libs import lib

from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar,
    maybe_cast_to_datetime,
    maybe_convert_platform,
    maybe_infer_to_datetimelike,
    maybe_upcast,
)
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeIndex,
    ABCIndexClass,
    ABCPeriodIndex,
    ABCSeries,
    ABCTimedeltaIndex,
)

from pandas.core import algorithms, common as com
from pandas.core.arrays import Categorical
from pandas.core.construction import extract_array, sanitize_array
from pandas.core.indexes import base as ibase
from pandas.core.indexes.api import (
    Index,
    ensure_index,
    get_objs_combined_axis,
    union_indexes,
)
from pandas.core.internals import (
    create_block_manager_from_arrays,
    create_block_manager_from_blocks,
)

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


def masked_rec_array_to_mgr(data, index, columns, dtype, copy: bool):
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
        # TODO: numpy docs suggest fv must be scalar, but could it be
        #  non-scalar for object dtype?
        assert lib.is_scalar(fv), fv
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
    if is_categorical_dtype(getattr(values, "dtype", None)) or is_categorical_dtype(
        dtype
    ):

        if not hasattr(values, "dtype"):
            values = _prep_ndarray(values, copy=copy)
            values = values.ravel()
        elif copy:
            values = values.copy()

        index, columns = _get_axes(len(values), 1, index, columns)
        return arrays_to_mgr([values], columns, index, columns, dtype=dtype)
    elif is_extension_array_dtype(values) or is_extension_array_dtype(dtype):
        # GH#19157

        if isinstance(values, np.ndarray) and values.ndim > 1:
            # GH#12513 a EA dtype passed with a 2D array, split into
            #  multiple EAs that view the values
            values = [values[:, n] for n in range(values.shape[1])]
        else:
            values = [values]

        if columns is None:
            columns = list(range(len(values)))
        return arrays_to_mgr(values, columns, index, columns, dtype=dtype)

    # by definition an array here
    # the dtypes will be coerced to a single dtype
    values = _prep_ndarray(values, copy=copy)

    if dtype is not None:
        if not is_dtype_equal(values.dtype, dtype):
            try:
                values = values.astype(dtype)
            except Exception as orig:
                # e.g. ValueError when trying to cast object dtype to float64
                raise ValueError(
                    f"failed to cast to '{dtype}' (Exception was: {orig})"
                ) from orig

    index, columns = _get_axes(*values.shape, index=index, columns=columns)
    values = values.T

    # if we don't have a dtype specified, then try to convert objects
    # on the entire block; this is to convert if we have datetimelike's
    # embedded in an object type
    if dtype is None and is_object_dtype(values):

        if values.ndim == 2 and values.shape[0] != 1:
            # transpose and separate blocks

            dvals_list = [maybe_infer_to_datetimelike(row) for row in values]
            for n in range(len(dvals_list)):
                if isinstance(dvals_list[n], np.ndarray):
                    dvals_list[n] = dvals_list[n].reshape(1, -1)

            from pandas.core.internals.blocks import make_block

            # TODO: What about re-joining object columns?
            block_values = [
                make_block(dvals_list[n], placement=[n]) for n in range(len(dvals_list))
            ]

        else:
            datelike_vals = maybe_infer_to_datetimelike(values)
            block_values = [datelike_vals]
    else:
        block_values = [values]

    return create_block_manager_from_blocks(block_values, [columns, index])


def init_dict(data, index, columns, dtype=None):
    """
    Segregate Series based on type and coerce into matrices.
    Needs to handle a lot of exceptional cases.
    """
    if columns is not None:
        from pandas.core.series import Series

        arrays = Series(data, index=columns, dtype=object)
        data_names = arrays.index

        missing = arrays.isna()
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
            val = construct_1d_arraylike_from_scalar(np.nan, len(index), nan_dtype)
            arrays.loc[missing] = [val] * missing.sum()

    else:
        keys = list(data.keys())
        columns = data_names = Index(keys)
        arrays = (com.maybe_iterable_to_list(data[k]) for k in keys)
        # GH#24096 need copy to be deep for datetime64tz case
        # TODO: See if we can avoid these copies
        arrays = [
            arr if not isinstance(arr, ABCIndexClass) else arr._data for arr in arrays
        ]
        arrays = [
            arr if not is_datetime64tz_dtype(arr) else arr.copy() for arr in arrays
        ]
    return arrays_to_mgr(arrays, data_names, index, columns, dtype=dtype)


# ---------------------------------------------------------------------


def _prep_ndarray(values, copy: bool = True) -> np.ndarray:
    if not isinstance(values, (np.ndarray, ABCSeries, Index)):
        if len(values) == 0:
            return np.empty((0, 0), dtype=object)
        elif isinstance(values, range):
            arr = np.arange(values.start, values.stop, values.step, dtype="int64")
            return arr[..., np.newaxis]

        def convert(v):
            return maybe_convert_platform(v)

        # we could have a 1-dim or 2-dim list here
        # this is equiv of np.asarray, but does object conversion
        # and platform dtype preservation
        try:
            if is_list_like(values[0]) or hasattr(values[0], "len"):
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
        raise ValueError("Must pass 2-d input")

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
                    oindex = index.astype("O")

                if isinstance(index, (ABCDatetimeIndex, ABCTimedeltaIndex)):
                    val = com.dict_compat(val)
                else:
                    val = dict(val)
                val = lib.fast_multiget(val, oindex.values, default=np.nan)
            val = sanitize_array(
                val, index, dtype=dtype, copy=False, raise_cast_failure=False
            )

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
            elif is_list_like(val) and getattr(val, "ndim", 1) == 1:
                have_raw_arrays = True
                raw_lengths.append(len(val))

        if not indexes and not raw_lengths:
            raise ValueError("If using all scalar values, you must pass an index")

        if have_series:
            index = union_indexes(indexes)
        elif have_dicts:
            index = union_indexes(indexes, sort=False)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError("arrays must all be same length")

            if have_dicts:
                raise ValueError(
                    "Mixing dicts with non-Series may lead to ambiguous ordering."
                )

            if have_series:
                if lengths[0] != len(index):
                    msg = (
                        f"array length {lengths[0]} does not match index "
                        f"length {len(index)}"
                    )
                    raise ValueError(msg)
            else:
                index = ibase.default_index(lengths[0])

    return ensure_index(index)


def reorder_arrays(arrays, arr_columns, columns):
    # reorder according to the columns
    if (
        columns is not None
        and len(columns)
        and arr_columns is not None
        and len(arr_columns)
    ):
        indexer = ensure_index(arr_columns).get_indexer(columns)
        arr_columns = ensure_index([arr_columns[i] for i in indexer])
        arrays = [arrays[i] for i in indexer]
    return arrays, arr_columns


def get_names_from_index(data):
    has_some_name = any(getattr(s, "name", None) is not None for s in data)
    if not has_some_name:
        return ibase.default_index(len(data))

    index = list(range(len(data)))
    count = 0
    for i, s in enumerate(data):
        n = getattr(s, "name", None)
        if n is not None:
            index[i] = n
        else:
            index[i] = f"Unnamed {count}"
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
            arrays = [
                data._ixs(i, axis=1).values
                for i, col in enumerate(data.columns)
                if col in columns
            ]
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
        return _list_to_arrays(data, columns, coerce_float=coerce_float, dtype=dtype)
    elif isinstance(data[0], abc.Mapping):
        return _list_of_dict_to_arrays(
            data, columns, coerce_float=coerce_float, dtype=dtype
        )
    elif isinstance(data[0], ABCSeries):
        return _list_of_series_to_arrays(
            data, columns, coerce_float=coerce_float, dtype=dtype
        )
    elif isinstance(data[0], Categorical):
        if columns is None:
            columns = ibase.default_index(len(data))
        return data, columns
    elif (
        isinstance(data, (np.ndarray, ABCSeries, Index))
        and data.dtype.names is not None
    ):

        columns = list(data.dtype.names)
        arrays = [data[k] for k in columns]
        return arrays, columns
    else:
        # last ditch effort
        data = [tuple(x) for x in data]
        return _list_to_arrays(data, columns, coerce_float=coerce_float, dtype=dtype)


def _list_to_arrays(data, columns, coerce_float=False, dtype=None):
    if len(data) > 0 and isinstance(data[0], tuple):
        content = list(lib.to_object_array_tuples(data).T)
    else:
        # list of lists
        content = list(lib.to_object_array(data).T)
    # gh-26429 do not raise user-facing AssertionError
    try:
        result = _convert_object_array(
            content, columns, dtype=dtype, coerce_float=coerce_float
        )
    except AssertionError as e:
        raise ValueError(e) from e
    return result


def _list_of_series_to_arrays(data, columns, coerce_float=False, dtype=None):
    if columns is None:
        # We know pass_data is non-empty because data[0] is a Series
        pass_data = [x for x in data if isinstance(x, (ABCSeries, ABCDataFrame))]
        columns = get_objs_combined_axis(pass_data, sort=False)

    indexer_cache = {}

    aligned_values = []
    for s in data:
        index = getattr(s, "index", None)
        if index is None:
            index = ibase.default_index(len(s))

        if id(index) in indexer_cache:
            indexer = indexer_cache[id(index)]
        else:
            indexer = indexer_cache[id(index)] = index.get_indexer(columns)

        values = extract_array(s, extract_numpy=True)
        aligned_values.append(algorithms.take_1d(values, indexer))

    values = np.vstack(aligned_values)

    if values.dtype == np.object_:
        content = list(values.T)
        return _convert_object_array(
            content, columns, dtype=dtype, coerce_float=coerce_float
        )
    else:
        return values.T, columns


def _list_of_dict_to_arrays(data, columns, coerce_float=False, dtype=None):
    """
    Convert list of dicts to numpy arrays

    if `columns` is not passed, column names are inferred from the records
    - for OrderedDict and dicts, the column names match
      the key insertion-order from the first record to the last.
    - For other kinds of dict-likes, the keys are lexically sorted.

    Parameters
    ----------
    data : iterable
        collection of records (OrderedDict, dict)
    columns: iterables or None
    coerce_float : bool
    dtype : np.dtype

    Returns
    -------
    tuple
        arrays, columns
    """
    if columns is None:
        gen = (list(x.keys()) for x in data)
        sort = not any(isinstance(d, dict) for d in data)
        columns = lib.fast_unique_multiple_list_gen(gen, sort=sort)

    # assure that they are of the base dict class and not of derived
    # classes
    data = [(type(d) is dict) and d or dict(d) for d in data]

    content = list(lib.dicts_to_array(data, list(columns)).T)
    return _convert_object_array(
        content, columns, dtype=dtype, coerce_float=coerce_float
    )


def _convert_object_array(content, columns, coerce_float=False, dtype=None):
    if columns is None:
        columns = ibase.default_index(len(content))
    else:
        if len(columns) != len(content):  # pragma: no cover
            # caller's responsibility to check for this...
            raise AssertionError(
                f"{len(columns)} columns passed, passed data had "
                f"{len(content)} columns"
            )

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


def sanitize_index(data, index: Index):
    """
    Sanitize an index type to return an ndarray of the underlying, pass
    through a non-Index.
    """
    if len(data) != len(index):
        raise ValueError("Length of values does not match length of index")

    if isinstance(data, ABCIndexClass):
        pass
    elif isinstance(data, (ABCPeriodIndex, ABCDatetimeIndex)):
        data = data._values

    elif isinstance(data, np.ndarray):

        # coerce datetimelike types
        if data.dtype.kind in ["M", "m"]:
            data = sanitize_array(data, index, copy=False)

    return data
