"""
Functions for preparing various inputs passed to the DataFrame or Series
constructors before passing them to a BlockManager.
"""
from __future__ import annotations

from collections import abc
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Hashable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
import numpy.ma as ma

from pandas._libs import lib
from pandas._typing import (
    Axis,
    DtypeObj,
    Manager,
    Scalar,
)

from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar,
    construct_1d_ndarray_preserving_na,
    dict_compat,
    maybe_cast_to_datetime,
    maybe_convert_platform,
    maybe_infer_to_datetimelike,
    maybe_upcast,
)
from pandas.core.dtypes.common import (
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_ea_dtype,
    is_extension_array_dtype,
    is_integer_dtype,
    is_list_like,
    is_named_tuple,
    is_object_dtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeIndex,
    ABCIndex,
    ABCSeries,
    ABCTimedeltaIndex,
)

from pandas.core import (
    algorithms,
    common as com,
)
from pandas.core.arrays import (
    Categorical,
    ExtensionArray,
)
from pandas.core.construction import (
    extract_array,
    sanitize_array,
)
from pandas.core.indexes import base as ibase
from pandas.core.indexes.api import (
    Index,
    ensure_index,
    get_objs_combined_axis,
    union_indexes,
)
from pandas.core.internals.blocks import block_shape
from pandas.core.internals.managers import (
    create_block_manager_from_arrays,
    create_block_manager_from_blocks,
)

if TYPE_CHECKING:
    from numpy.ma.mrecords import MaskedRecords

# ---------------------------------------------------------------------
# BlockManager Interface


def arrays_to_mgr(
    arrays,
    arr_names,
    index,
    columns,
    dtype: Optional[DtypeObj] = None,
    verify_integrity: bool = True,
):
    """
    Segregate Series based on type and coerce into matrices.

    Needs to handle a lot of exceptional cases.
    """
    arr_names = ensure_index(arr_names)

    if verify_integrity:
        # figure out the index, if necessary
        if index is None:
            index = extract_index(arrays)
        else:
            index = ensure_index(index)

        # don't force copy because getting jammed in an ndarray anyway
        arrays = _homogenize(arrays, index, dtype)

        columns = ensure_index(columns)
    else:
        columns = ensure_index(columns)
        index = ensure_index(index)

    # from BlockManager perspective
    axes = [columns, index]

    return create_block_manager_from_arrays(arrays, arr_names, axes)


def masked_rec_array_to_mgr(
    data: MaskedRecords, index, columns, dtype: Optional[DtypeObj], copy: bool
):
    """
    Extract from a masked rec array and create the manager.
    """
    # essentially process a record array then fill it
    fdata = ma.getdata(data)
    if index is None:
        index = _get_names_from_index(fdata)
        if index is None:
            index = ibase.default_index(len(data))
    index = ensure_index(index)

    if columns is not None:
        columns = ensure_index(columns)
    arrays, arr_columns = to_arrays(fdata, columns)

    # fill if needed
    new_arrays = []
    for col in arr_columns:
        arr = data[col]
        fv = arr.fill_value

        mask = ma.getmaskarray(arr)
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


def mgr_to_mgr(mgr, typ: str):
    """
    Convert to specific type of Manager. Does not copy if the type is already
    correct. Does not guarantee a copy otherwise.
    """
    from pandas.core.internals import (
        ArrayManager,
        BlockManager,
    )

    new_mgr: Manager

    if typ == "block":
        if isinstance(mgr, BlockManager):
            new_mgr = mgr
        else:
            new_mgr = arrays_to_mgr(
                mgr.arrays, mgr.axes[0], mgr.axes[1], mgr.axes[0], dtype=None
            )
    elif typ == "array":
        if isinstance(mgr, ArrayManager):
            new_mgr = mgr
        else:
            arrays = [mgr.iget_values(i).copy() for i in range(len(mgr.axes[0]))]
            new_mgr = ArrayManager(arrays, [mgr.axes[1], mgr.axes[0]])
    else:
        raise ValueError(f"'typ' needs to be one of {{'block', 'array'}}, got '{type}'")
    return new_mgr


# ---------------------------------------------------------------------
# DataFrame Constructor Interface


def init_ndarray(values, index, columns, dtype: Optional[DtypeObj], copy: bool):
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

    if is_ea_dtype(values) or is_extension_array_dtype(dtype):
        # GH#19157

        if isinstance(values, np.ndarray) and values.ndim > 1:
            # GH#12513 a EA dtype passed with a 2D array, split into
            #  multiple EAs that view the values
            values = [values[:, n] for n in range(values.shape[1])]
        else:
            values = [values]

        if columns is None:
            columns = Index(range(len(values)))

        return arrays_to_mgr(values, columns, index, columns, dtype=dtype)

    if is_datetime64tz_dtype(values):
        # TODO: combine into _prep_ndarray?
        values = extract_array(values, extract_numpy=True)
        if copy:
            values = values.copy()
        if values.ndim == 1:
            values = values.reshape(-1, 1)

    else:
        # by definition an array here
        # the dtypes will be coerced to a single dtype
        values = _prep_ndarray(values, copy=copy)

    if dtype is not None and not is_dtype_equal(values.dtype, dtype):
        try:
            values = construct_1d_ndarray_preserving_na(
                values.ravel(), dtype=dtype, copy=False
            ).reshape(values.shape)
        except Exception as orig:
            # e.g. ValueError when trying to cast object dtype to float64
            raise ValueError(
                f"failed to cast to '{dtype}' (Exception was: {orig})"
            ) from orig

    # _prep_ndarray ensures that values.ndim == 2 at this point
    index, columns = _get_axes(
        values.shape[0], values.shape[1], index=index, columns=columns
    )
    values = values.T

    # if we don't have a dtype specified, then try to convert objects
    # on the entire block; this is to convert if we have datetimelike's
    # embedded in an object type
    if dtype is None and is_object_dtype(values.dtype):

        if values.ndim == 2 and values.shape[0] != 1:
            # transpose and separate blocks

            # TODO: do this in one go
            dvals_list = [maybe_infer_to_datetimelike(row) for row in values]
            dvals_list = [extract_array(x, extract_numpy=True) for x in dvals_list]
            # TODO: unpack DatetimeIndex directly in maybe_infer_to_datetimelike
            for n in range(len(dvals_list)):
                dvals_list[n] = dvals_list[n].reshape(1, -1)

            from pandas.core.internals.blocks import make_block

            # TODO: What about re-joining object columns?
            block_values = [
                make_block(dvals_list[n], placement=[n], ndim=2)
                for n in range(len(dvals_list))
            ]

        else:
            datelike_vals = maybe_infer_to_datetimelike(values)
            block_values = [datelike_vals]
            block_values = [extract_array(x, extract_numpy=True) for x in block_values]
            if values.ndim == 2:
                block_values = [block_shape(x, 2) for x in block_values]

    else:
        block_values = [values]

    return create_block_manager_from_blocks(block_values, [columns, index])


def init_dict(data: Dict, index, columns, dtype: Optional[DtypeObj] = None):
    """
    Segregate Series based on type and coerce into matrices.
    Needs to handle a lot of exceptional cases.
    """
    arrays: Union[Sequence[Any], Series]

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
            if dtype is None or (
                not is_extension_array_dtype(dtype)
                and np.issubdtype(dtype, np.flexible)
            ):
                # GH#1783
                nan_dtype = np.dtype(object)
            else:
                nan_dtype = dtype
            val = construct_1d_arraylike_from_scalar(np.nan, len(index), nan_dtype)
            arrays.loc[missing] = [val] * missing.sum()

    else:
        keys = list(data.keys())
        columns = data_names = Index(keys)
        arrays = [com.maybe_iterable_to_list(data[k]) for k in keys]
        # GH#24096 need copy to be deep for datetime64tz case
        # TODO: See if we can avoid these copies
        arrays = [arr if not isinstance(arr, ABCIndex) else arr._data for arr in arrays]
        arrays = [
            arr if not is_datetime64tz_dtype(arr) else arr.copy() for arr in arrays
        ]
    return arrays_to_mgr(arrays, data_names, index, columns, dtype=dtype)


def nested_data_to_arrays(
    data: Sequence,
    columns: Optional[Index],
    index: Optional[Index],
    dtype: Optional[DtypeObj],
):
    """
    Convert a single sequence of arrays to multiple arrays.
    """
    # By the time we get here we have already checked treat_as_nested(data)

    if is_named_tuple(data[0]) and columns is None:
        columns = data[0]._fields

    arrays, columns = to_arrays(data, columns, dtype=dtype)
    columns = ensure_index(columns)

    if index is None:
        if isinstance(data[0], ABCSeries):
            index = _get_names_from_index(data)
        elif isinstance(data[0], Categorical):
            # GH#38845 hit in test_constructor_categorical
            index = ibase.default_index(len(data[0]))
        else:
            index = ibase.default_index(len(data))

    return arrays, columns, index


def treat_as_nested(data) -> bool:
    """
    Check if we should use nested_data_to_arrays.
    """
    return (
        len(data) > 0
        and is_list_like(data[0])
        and getattr(data[0], "ndim", 1) == 1
        and not (isinstance(data, ExtensionArray) and data.ndim == 2)
    )


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
        raise ValueError(f"Must pass 2-d input. shape={values.shape}")

    return values


def _homogenize(data, index, dtype: Optional[DtypeObj]):
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
                    val = dict_compat(val)
                else:
                    val = dict(val)
                val = lib.fast_multiget(val, oindex._values, default=np.nan)
            val = sanitize_array(
                val, index, dtype=dtype, copy=False, raise_cast_failure=False
            )

        homogenized.append(val)

    return homogenized


def extract_index(data) -> Index:
    """
    Try to infer an Index from the passed data, raise ValueError on failure.
    """
    index = None
    if len(data) == 0:
        index = Index([])
    elif len(data) > 0:
        raw_lengths = []
        indexes: List[Union[List[Hashable], Index]] = []

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
                raise ValueError("All arrays must be of the same length")

            if have_dicts:
                raise ValueError(
                    "Mixing dicts with non-Series may lead to ambiguous ordering."
                )

            if have_series:
                assert index is not None  # for mypy
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


def _get_names_from_index(data):
    has_some_name = any(getattr(s, "name", None) is not None for s in data)
    if not has_some_name:
        return ibase.default_index(len(data))

    index: List[Hashable] = list(range(len(data)))
    count = 0
    for i, s in enumerate(data):
        n = getattr(s, "name", None)
        if n is not None:
            index[i] = n
        else:
            index[i] = f"Unnamed {count}"
            count += 1

    return index


def _get_axes(
    N: int, K: int, index: Optional[Index], columns: Optional[Index]
) -> Tuple[Index, Index]:
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


def dataclasses_to_dicts(data):
    """
    Converts a list of dataclass instances to a list of dictionaries.

    Parameters
    ----------
    data : List[Type[dataclass]]

    Returns
    --------
    list_dict : List[dict]

    Examples
    --------
    >>> @dataclass
    >>> class Point:
    ...     x: int
    ...     y: int

    >>> dataclasses_to_dicts([Point(1,2), Point(2,3)])
    [{"x":1,"y":2},{"x":2,"y":3}]

    """
    from dataclasses import asdict

    return list(map(asdict, data))


# ---------------------------------------------------------------------
# Conversion of Inputs to Arrays


def to_arrays(data, columns, dtype: Optional[DtypeObj] = None):
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

    elif isinstance(data[0], Categorical):
        if columns is None:
            columns = ibase.default_index(len(data))
        return data, columns

    elif isinstance(data, np.ndarray) and data.dtype.names is not None:
        # e.g. recarray
        columns = list(data.dtype.names)
        arrays = [data[k] for k in columns]
        return arrays, columns

    if isinstance(data[0], (list, tuple)):
        content, columns = _list_to_arrays(data, columns)
    elif isinstance(data[0], abc.Mapping):
        content, columns = _list_of_dict_to_arrays(data, columns)
    elif isinstance(data[0], ABCSeries):
        content, columns = _list_of_series_to_arrays(data, columns)
    else:
        # last ditch effort
        data = [tuple(x) for x in data]
        content, columns = _list_to_arrays(data, columns)

    content, columns = _finalize_columns_and_data(content, columns, dtype)
    return content, columns


def _list_to_arrays(
    data: List[Scalar],
    columns: Union[Index, List],
) -> Tuple[List[Scalar], Union[Index, List[Axis]]]:
    # Note: we already check len(data) > 0 before getting hre
    if isinstance(data[0], tuple):
        content = lib.to_object_array_tuples(data)
    else:
        # list of lists
        content = lib.to_object_array(data)
    return content, columns


def _list_of_series_to_arrays(
    data: List,
    columns: Union[Index, List],
) -> Tuple[List[Scalar], Union[Index, List[Axis]]]:
    if columns is None:
        # We know pass_data is non-empty because data[0] is a Series
        pass_data = [x for x in data if isinstance(x, (ABCSeries, ABCDataFrame))]
        columns = get_objs_combined_axis(pass_data, sort=False)

    indexer_cache: Dict[int, Scalar] = {}

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
        aligned_values.append(algorithms.take_nd(values, indexer))

    content = np.vstack(aligned_values)

    return content, columns


def _list_of_dict_to_arrays(
    data: List[Dict],
    columns: Union[Index, List],
) -> Tuple[List[Scalar], Union[Index, List[Axis]]]:
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

    content = lib.dicts_to_array(data, list(columns))
    return content, columns


def _finalize_columns_and_data(
    content: np.ndarray,
    columns: Optional[Union[Index, List]],
    dtype: Optional[DtypeObj],
) -> Tuple[List[np.ndarray], Union[Index, List[Axis]]]:
    """
    Ensure we have valid columns, cast object dtypes if possible.
    """
    content = list(content.T)

    try:
        columns = _validate_or_indexify_columns(content, columns)
    except AssertionError as err:
        # GH#26429 do not raise user-facing AssertionError
        raise ValueError(err) from err

    if len(content) and content[0].dtype == np.object_:
        content = _convert_object_array(content, dtype=dtype)
    return content, columns


def _validate_or_indexify_columns(
    content: List, columns: Optional[Union[Index, List]]
) -> Union[Index, List[Axis]]:
    """
    If columns is None, make numbers as column names; Otherwise, validate that
    columns have valid length.

    Parameters
    ----------
    content: list of data
    columns: Iterable or None

    Returns
    -------
    columns: If columns is Iterable, return as is; If columns is None, assign
    positional column index value as columns.

    Raises
    ------
    1. AssertionError when content is not composed of list of lists, and if
        length of columns is not equal to length of content.
    2. ValueError when content is list of lists, but length of each sub-list
        is not equal
    3. ValueError when content is list of lists, but length of sub-list is
        not equal to length of content
    """
    if columns is None:
        columns = ibase.default_index(len(content))
    else:

        # Add mask for data which is composed of list of lists
        is_mi_list = isinstance(columns, list) and all(
            isinstance(col, list) for col in columns
        )

        if not is_mi_list and len(columns) != len(content):  # pragma: no cover
            # caller's responsibility to check for this...
            raise AssertionError(
                f"{len(columns)} columns passed, passed data had "
                f"{len(content)} columns"
            )
        elif is_mi_list:

            # check if nested list column, length of each sub-list should be equal
            if len({len(col) for col in columns}) > 1:
                raise ValueError(
                    "Length of columns passed for MultiIndex columns is different"
                )

            # if columns is not empty and length of sublist is not equal to content
            elif columns and len(columns[0]) != len(content):
                raise ValueError(
                    f"{len(columns[0])} columns passed, passed data had "
                    f"{len(content)} columns"
                )
    return columns


def _convert_object_array(
    content: List[Scalar], dtype: Optional[DtypeObj] = None
) -> List[Scalar]:
    """
    Internal function to convert object array.

    Parameters
    ----------
    content: list of processed data records
    dtype: np.dtype, default is None

    Returns
    -------
    arrays: casted content if not object dtype, otherwise return as is in list.
    """
    # provide soft conversion of object dtypes
    def convert(arr):
        if dtype != np.dtype("O"):
            arr = lib.maybe_convert_objects(arr)
            arr = maybe_cast_to_datetime(arr, dtype)
        return arr

    arrays = [convert(arr) for arr in content]

    return arrays


# ---------------------------------------------------------------------
# Series-Based


def sanitize_index(data, index: Index):
    """
    Sanitize an index type to return an ndarray of the underlying, pass
    through a non-Index.
    """
    if len(data) != len(index):
        raise ValueError(
            "Length of values "
            f"({len(data)}) "
            "does not match length of index "
            f"({len(index)})"
        )

    if isinstance(data, np.ndarray):

        # coerce datetimelike types
        if data.dtype.kind in ["M", "m"]:
            data = sanitize_array(data, index, copy=False)

    return data
