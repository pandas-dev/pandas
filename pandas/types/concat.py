"""
Utility functions related to concat
"""

import numpy as np
import pandas.core.common as com
import pandas.tslib as tslib
from pandas import compat
from pandas.compat import map


def get_dtype_kinds(l):
    """
    Parameters
    ----------
    l : list of arrays

    Returns
    -------
    a set of kinds that exist in this list of arrays
    """

    typs = set()
    for arr in l:

        dtype = arr.dtype
        if com.is_categorical_dtype(dtype):
            typ = 'category'
        elif com.is_sparse(arr):
            typ = 'sparse'
        elif com.is_datetimetz(arr):
            typ = 'datetimetz'
        elif com.is_datetime64_dtype(dtype):
            typ = 'datetime'
        elif com.is_timedelta64_dtype(dtype):
            typ = 'timedelta'
        elif com.is_object_dtype(dtype):
            typ = 'object'
        elif com.is_bool_dtype(dtype):
            typ = 'bool'
        else:
            typ = dtype.kind
        typs.add(typ)
    return typs


def _get_series_result_type(result):
    """
    return appropriate class of Series concat
    input is either dict or array-like
    """
    if isinstance(result, dict):
        # concat Series with axis 1
        if all(com.is_sparse(c) for c in compat.itervalues(result)):
            from pandas.sparse.api import SparseDataFrame
            return SparseDataFrame
        else:
            from pandas.core.frame import DataFrame
            return DataFrame

    elif com.is_sparse(result):
        # concat Series with axis 1
        from pandas.sparse.api import SparseSeries
        return SparseSeries
    else:
        from pandas.core.series import Series
        return Series


def _get_frame_result_type(result, objs):
    """
    return appropriate class of DataFrame-like concat
    if any block is SparseBlock, return SparseDataFrame
    otherwise, return 1st obj
    """
    if any(b.is_sparse for b in result.blocks):
        from pandas.sparse.api import SparseDataFrame
        return SparseDataFrame
    else:
        return objs[0]


def _concat_compat(to_concat, axis=0):
    """
    provide concatenation of an array of arrays each of which is a single
    'normalized' dtypes (in that for example, if it's object, then it is a
    non-datetimelike and provide a combined dtype for the resulting array that
    preserves the overall dtype if possible)

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation

    Returns
    -------
    a single array, preserving the combined dtypes
    """

    # filter empty arrays
    # 1-d dtypes always are included here
    def is_nonempty(x):
        try:
            return x.shape[axis] > 0
        except Exception:
            return True

    nonempty = [x for x in to_concat if is_nonempty(x)]

    # If all arrays are empty, there's nothing to convert, just short-cut to
    # the concatenation, #3121.
    #
    # Creating an empty array directly is tempting, but the winnings would be
    # marginal given that it would still require shape & dtype calculation and
    # np.concatenate which has them both implemented is compiled.

    typs = get_dtype_kinds(to_concat)

    # these are mandated to handle empties as well
    if 'datetime' in typs or 'datetimetz' in typs or 'timedelta' in typs:
        return _concat_datetime(to_concat, axis=axis, typs=typs)

    elif 'sparse' in typs:
        return _concat_sparse(to_concat, axis=axis, typs=typs)

    elif 'category' in typs:
        return _concat_categorical(to_concat, axis=axis)

    if not nonempty:
        # we have all empties, but may need to coerce the result dtype to
        # object if we have non-numeric type operands (numpy would otherwise
        # cast this to float)
        typs = get_dtype_kinds(to_concat)
        if len(typs) != 1:

            if (not len(typs - set(['i', 'u', 'f'])) or
                    not len(typs - set(['bool', 'i', 'u']))):
                # let numpy coerce
                pass
            else:
                # coerce to object
                to_concat = [x.astype('object') for x in to_concat]

    return np.concatenate(to_concat, axis=axis)


def _concat_categorical(to_concat, axis=0):
    """Concatenate an object/categorical array of arrays, each of which is a
    single dtype

    Parameters
    ----------
    to_concat : array of arrays
    axis : int
        Axis to provide concatenation in the current implementation this is
        always 0, e.g. we only have 1D categoricals

    Returns
    -------
    Categorical
        A single array, preserving the combined dtypes
    """

    from pandas.core.categorical import Categorical

    def convert_categorical(x):
        # coerce to object dtype
        if com.is_categorical_dtype(x.dtype):
            return x.get_values()
        return x.ravel()

    if get_dtype_kinds(to_concat) - set(['object', 'category']):
        # convert to object type and perform a regular concat
        return _concat_compat([np.array(x, copy=False, dtype=object)
                               for x in to_concat], axis=0)

    # we could have object blocks and categoricals here
    # if we only have a single categoricals then combine everything
    # else its a non-compat categorical
    categoricals = [x for x in to_concat if com.is_categorical_dtype(x.dtype)]

    # validate the categories
    categories = categoricals[0]
    rawcats = categories.categories
    for x in categoricals[1:]:
        if not categories.is_dtype_equal(x):
            raise ValueError("incompatible categories in categorical concat")

    # we've already checked that all categoricals are the same, so if their
    # length is equal to the input then we have all the same categories
    if len(categoricals) == len(to_concat):
        # concating numeric types is much faster than concating object types
        # and fastpath takes a shorter path through the constructor
        return Categorical(np.concatenate([x.codes for x in to_concat],
                                          axis=0),
                           rawcats, ordered=categoricals[0].ordered,
                           fastpath=True)
    else:
        concatted = np.concatenate(list(map(convert_categorical, to_concat)),
                                   axis=0)
        return Categorical(concatted, rawcats)


def _concat_datetime(to_concat, axis=0, typs=None):
    """
    provide concatenation of an datetimelike array of arrays each of which is a
    single M8[ns], datetimet64[ns, tz] or m8[ns] dtype

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation
    typs : set of to_concat dtypes

    Returns
    -------
    a single array, preserving the combined dtypes
    """

    def convert_to_pydatetime(x, axis):
        # coerce to an object dtype

        # if dtype is of datetimetz or timezone
        if x.dtype.kind == com._NS_DTYPE.kind:
            if getattr(x, 'tz', None) is not None:
                x = x.asobject.values
            else:
                shape = x.shape
                x = tslib.ints_to_pydatetime(x.view(np.int64).ravel())
                x = x.reshape(shape)

        elif x.dtype == com._TD_DTYPE:
            shape = x.shape
            x = tslib.ints_to_pytimedelta(x.view(np.int64).ravel())
            x = x.reshape(shape)

        if axis == 1:
            x = np.atleast_2d(x)
        return x

    if typs is None:
        typs = get_dtype_kinds(to_concat)

    # must be single dtype
    if len(typs) == 1:

        if 'datetimetz' in typs:
            # datetime with no tz should be stored as "datetime" in typs,
            # thus no need to care

            # we require ALL of the same tz for datetimetz
            tzs = set([x.tz for x in to_concat])
            if len(tzs) == 1:
                from pandas.tseries.index import DatetimeIndex
                new_values = np.concatenate([x.tz_localize(None).asi8
                                             for x in to_concat])
                return DatetimeIndex(new_values, tz=list(tzs)[0])

        elif 'datetime' in typs:
            new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                        axis=axis)
            return new_values.view(com._NS_DTYPE)

        elif 'timedelta' in typs:
            new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                        axis=axis)
            return new_values.view(com._TD_DTYPE)

    # need to coerce to object
    to_concat = [convert_to_pydatetime(x, axis) for x in to_concat]
    return np.concatenate(to_concat, axis=axis)


def _concat_sparse(to_concat, axis=0, typs=None):
    """
    provide concatenation of an sparse/dense array of arrays each of which is a
    single dtype

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation
    typs : set of to_concat dtypes

    Returns
    -------
    a single array, preserving the combined dtypes
    """

    from pandas.sparse.array import SparseArray, _make_index

    def convert_sparse(x, axis):
        # coerce to native type
        if isinstance(x, SparseArray):
            x = x.get_values()
        x = x.ravel()
        if axis > 0:
            x = np.atleast_2d(x)
        return x

    if typs is None:
        typs = com.get_dtype_kinds(to_concat)

    if len(typs) == 1:
        # concat input as it is if all inputs are sparse
        # and have the same fill_value
        fill_values = set(c.fill_value for c in to_concat)
        if len(fill_values) == 1:
            sp_values = [c.sp_values for c in to_concat]
            indexes = [c.sp_index.to_int_index() for c in to_concat]

            indices = []
            loc = 0
            for idx in indexes:
                indices.append(idx.indices + loc)
                loc += idx.length
            sp_values = np.concatenate(sp_values)
            indices = np.concatenate(indices)
            sp_index = _make_index(loc, indices, kind=to_concat[0].sp_index)

            return SparseArray(sp_values, sparse_index=sp_index,
                               fill_value=to_concat[0].fill_value)

    # input may be sparse / dense mixed and may have different fill_value
    # input must contain sparse at least 1
    sparses = [c for c in to_concat if com.is_sparse(c)]
    fill_values = [c.fill_value for c in sparses]
    sp_indexes = [c.sp_index for c in sparses]

    # densify and regular concat
    to_concat = [convert_sparse(x, axis) for x in to_concat]
    result = np.concatenate(to_concat, axis=axis)

    if not len(typs - set(['sparse', 'f', 'i'])):
        # sparsify if inputs are sparse and dense numerics
        # first sparse input's fill_value and SparseIndex is used
        result = SparseArray(result.ravel(), fill_value=fill_values[0],
                             kind=sp_indexes[0])
    else:
        # coerce to object if needed
        result = result.astype('object')
    return result
