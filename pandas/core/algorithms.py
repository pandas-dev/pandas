"""
Generic data algorithms. This module is experimental at the moment and not
intended for public consumption
"""
from __future__ import division
from warnings import warn
import numpy as np

from pandas import compat, lib, tslib, _np_version_under1p8
from pandas.types.cast import _maybe_promote
from pandas.types.generic import ABCSeries, ABCIndex
from pandas.types.common import (is_integer_dtype,
                                 is_int64_dtype,
                                 is_categorical_dtype,
                                 is_extension_type,
                                 is_datetimetz,
                                 is_period_dtype,
                                 is_period_arraylike,
                                 is_float_dtype,
                                 needs_i8_conversion,
                                 is_categorical,
                                 is_datetime64_dtype,
                                 is_timedelta64_dtype,
                                 is_scalar,
                                 _ensure_platform_int,
                                 _ensure_object,
                                 _ensure_float64,
                                 _ensure_int64,
                                 is_list_like)
from pandas.types.missing import isnull

import pandas.core.common as com
import pandas.algos as algos
import pandas.hashtable as htable
from pandas.compat import string_types
from pandas.tslib import iNaT


# --------------- #
# top-level algos #
# --------------- #

def match(to_match, values, na_sentinel=-1):
    """
    Compute locations of to_match into values

    Parameters
    ----------
    to_match : array-like
        values to find positions of
    values : array-like
        Unique set of values
    na_sentinel : int, default -1
        Value to mark "not found"

    Examples
    --------

    Returns
    -------
    match : ndarray of integers
    """
    values = com._asarray_tuplesafe(values)
    if issubclass(values.dtype.type, string_types):
        values = np.array(values, dtype='O')

    f = lambda htype, caster: _match_generic(to_match, values, htype, caster)
    result = _hashtable_algo(f, values.dtype, np.int64)

    if na_sentinel != -1:

        # replace but return a numpy array
        # use a Series because it handles dtype conversions properly
        from pandas.core.series import Series
        result = Series(result.ravel()).replace(-1, na_sentinel).values.\
            reshape(result.shape)

    return result


def _match_generic(values, index, table_type, type_caster):
    values = type_caster(values)
    index = type_caster(index)
    table = table_type(min(len(index), 1000000))
    table.map_locations(index)
    return table.lookup(values)


def unique(values):
    """
    Compute unique values (not necessarily sorted) efficiently from input array
    of values

    Parameters
    ----------
    values : array-like

    Returns
    -------
    uniques
    """
    values = com._asarray_tuplesafe(values)

    f = lambda htype, caster: _unique_generic(values, htype, caster)
    return _hashtable_algo(f, values.dtype)


def _unique_generic(values, table_type, type_caster):
    values = type_caster(values)
    table = table_type(min(len(values), 1000000))
    uniques = table.unique(values)
    return type_caster(uniques)


def isin(comps, values):
    """
    Compute the isin boolean array

    Parameters
    ----------
    comps: array-like
    values: array-like

    Returns
    -------
    boolean array same length as comps
    """

    if not is_list_like(comps):
        raise TypeError("only list-like objects are allowed to be passed"
                        " to isin(), you passed a "
                        "[{0}]".format(type(comps).__name__))
    comps = np.asarray(comps)
    if not is_list_like(values):
        raise TypeError("only list-like objects are allowed to be passed"
                        " to isin(), you passed a "
                        "[{0}]".format(type(values).__name__))
    if not isinstance(values, np.ndarray):
        values = list(values)

    # GH11232
    # work-around for numpy < 1.8 and comparisions on py3
    # faster for larger cases to use np.in1d
    if (_np_version_under1p8 and compat.PY3) or len(comps) > 1000000:
        f = lambda x, y: np.in1d(x, np.asarray(list(y)))
    else:
        f = lambda x, y: lib.ismember_int64(x, set(y))

    # may need i8 conversion for proper membership testing
    if is_datetime64_dtype(comps):
        from pandas.tseries.tools import to_datetime
        values = to_datetime(values)._values.view('i8')
        comps = comps.view('i8')
    elif is_timedelta64_dtype(comps):
        from pandas.tseries.timedeltas import to_timedelta
        values = to_timedelta(values)._values.view('i8')
        comps = comps.view('i8')
    elif is_int64_dtype(comps):
        pass
    else:
        f = lambda x, y: lib.ismember(x, set(values))

    return f(comps, values)


def safe_sort(values, labels=None, na_sentinel=-1, assume_unique=False):
    """
    Sort ``values`` and reorder corresponding ``labels``.
    ``values`` should be unique if ``labels`` is not None.
    Safe for use with mixed types (int, str), orders ints before strs.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    values : list-like
        Sequence; must be unique if ``labels`` is not None.
    labels : list_like
        Indices to ``values``. All out of bound indices are treated as
        "not found" and will be masked with ``na_sentinel``.
    na_sentinel : int, default -1
        Value in ``labels`` to mark "not found".
        Ignored when ``labels`` is None.
    assume_unique : bool, default False
        When True, ``values`` are assumed to be unique, which can speed up
        the calculation. Ignored when ``labels`` is None.

    Returns
    -------
    ordered : ndarray
        Sorted ``values``
    new_labels : ndarray
        Reordered ``labels``; returned when ``labels`` is not None.

    Raises
    ------
    TypeError
        * If ``values`` is not list-like or if ``labels`` is neither None
        nor list-like
        * If ``values`` cannot be sorted
    ValueError
        * If ``labels`` is not None and ``values`` contain duplicates.
    """
    if not is_list_like(values):
        raise TypeError("Only list-like objects are allowed to be passed to"
                        "safe_sort as values")
    values = np.array(values, copy=False)

    def sort_mixed(values):
        # order ints before strings, safe in py3
        str_pos = np.array([isinstance(x, string_types) for x in values],
                           dtype=bool)
        nums = np.sort(values[~str_pos])
        strs = np.sort(values[str_pos])
        return _ensure_object(np.concatenate([nums, strs]))

    sorter = None
    if compat.PY3 and lib.infer_dtype(values) == 'mixed-integer':
        # unorderable in py3 if mixed str/int
        ordered = sort_mixed(values)
    else:
        try:
            sorter = values.argsort()
            ordered = values.take(sorter)
        except TypeError:
            # try this anyway
            ordered = sort_mixed(values)

    # labels:

    if labels is None:
        return ordered

    if not is_list_like(labels):
        raise TypeError("Only list-like objects or None are allowed to be"
                        "passed to safe_sort as labels")
    labels = _ensure_platform_int(np.asarray(labels))

    from pandas import Index
    if not assume_unique and not Index(values).is_unique:
        raise ValueError("values should be unique if labels is not None")

    if sorter is None:
        # mixed types
        (hash_klass, _), values = _get_data_algo(values, _hashtables)
        t = hash_klass(len(values))
        t.map_locations(values)
        sorter = _ensure_platform_int(t.lookup(ordered))

    reverse_indexer = np.empty(len(sorter), dtype=np.int_)
    reverse_indexer.put(sorter, np.arange(len(sorter)))

    mask = (labels < -len(values)) | (labels >= len(values)) | \
        (labels == na_sentinel)

    # (Out of bound indices will be masked with `na_sentinel` next, so we may
    # deal with them here without performance loss using `mode='wrap'`.)
    new_labels = reverse_indexer.take(labels, mode='wrap')
    np.putmask(new_labels, mask, na_sentinel)

    return ordered, _ensure_platform_int(new_labels)


def factorize(values, sort=False, order=None, na_sentinel=-1, size_hint=None):
    """
    Encode input values as an enumerated type or categorical variable

    Parameters
    ----------
    values : ndarray (1-d)
        Sequence
    sort : boolean, default False
        Sort by values
    na_sentinel : int, default -1
        Value to mark "not found"
    size_hint : hint to the hashtable sizer

    Returns
    -------
    labels : the indexer to the original array
    uniques : ndarray (1-d) or Index
        the unique values. Index is returned when passed values is Index or
        Series

    note: an array of Periods will ignore sort as it returns an always sorted
    PeriodIndex
    """
    from pandas import Index, Series, DatetimeIndex, PeriodIndex

    # handling two possibilities here
    # - for a numpy datetimelike simply view as i8 then cast back
    # - for an extension datetimelike view as i8 then
    #   reconstruct from boxed values to transfer metadata
    dtype = None
    if needs_i8_conversion(values):
        if is_period_dtype(values):
            values = PeriodIndex(values)
            vals = values.asi8
        elif is_datetimetz(values):
            values = DatetimeIndex(values)
            vals = values.asi8
        else:
            # numpy dtype
            dtype = values.dtype
            vals = values.view(np.int64)
    else:
        vals = np.asarray(values)

    (hash_klass, vec_klass), vals = _get_data_algo(vals, _hashtables)

    table = hash_klass(size_hint or len(vals))
    uniques = vec_klass()
    labels = table.get_labels(vals, uniques, 0, na_sentinel, True)

    labels = _ensure_platform_int(labels)

    uniques = uniques.to_array()

    if sort and len(uniques) > 0:
        uniques, labels = safe_sort(uniques, labels, na_sentinel=na_sentinel,
                                    assume_unique=True)

    if dtype is not None:
        uniques = uniques.astype(dtype)

    if isinstance(values, Index):
        uniques = values._shallow_copy(uniques, name=None)
    elif isinstance(values, Series):
        uniques = Index(uniques)
    return labels, uniques


def value_counts(values, sort=True, ascending=False, normalize=False,
                 bins=None, dropna=True):
    """
    Compute a histogram of the counts of non-null values.

    Parameters
    ----------
    values : ndarray (1-d)
    sort : boolean, default True
        Sort by values
    ascending : boolean, default False
        Sort in ascending order
    normalize: boolean, default False
        If True then compute a relative histogram
    bins : integer, optional
        Rather than count values, group them into half-open bins,
        convenience for pd.cut, only works with numeric data
    dropna : boolean, default True
        Don't include counts of NaN

    Returns
    -------
    value_counts : Series

    """
    from pandas.core.series import Series
    name = getattr(values, 'name', None)

    if bins is not None:
        try:
            from pandas.tools.tile import cut
            values = Series(values).values
            cat, bins = cut(values, bins, retbins=True)
        except TypeError:
            raise TypeError("bins argument only works with numeric data.")
        values = cat.codes

    if is_extension_type(values) and not is_datetimetz(values):
        # handle Categorical and sparse,
        # datetime tz can be handeled in ndarray path
        result = Series(values).values.value_counts(dropna=dropna)
        result.name = name
        counts = result.values
    else:
        # ndarray path. pass original to handle DatetimeTzBlock
        keys, counts = _value_counts_arraylike(values, dropna=dropna)

        from pandas import Index, Series
        if not isinstance(keys, Index):
            keys = Index(keys)
        result = Series(counts, index=keys, name=name)

    if bins is not None:
        # TODO: This next line should be more efficient
        result = result.reindex(np.arange(len(cat.categories)),
                                fill_value=0)
        result.index = bins[:-1]

    if sort:
        result = result.sort_values(ascending=ascending)

    if normalize:
        result = result / float(counts.sum())

    return result


def _value_counts_arraylike(values, dropna=True):
    is_datetimetz_type = is_datetimetz(values)
    is_period_type = (is_period_dtype(values) or
                      is_period_arraylike(values))

    orig = values

    from pandas.core.series import Series
    values = Series(values).values
    dtype = values.dtype

    if needs_i8_conversion(dtype) or is_period_type:

        from pandas.tseries.index import DatetimeIndex
        from pandas.tseries.period import PeriodIndex

        if is_period_type:
            # values may be an object
            values = PeriodIndex(values)
            freq = values.freq

        values = values.view(np.int64)
        keys, counts = htable.value_count_int64(values, dropna)

        if dropna:
            msk = keys != iNaT
            keys, counts = keys[msk], counts[msk]

        # convert the keys back to the dtype we came in
        keys = keys.astype(dtype)

        # dtype handling
        if is_datetimetz_type:
            keys = DatetimeIndex._simple_new(keys, tz=orig.dtype.tz)
        if is_period_type:
            keys = PeriodIndex._simple_new(keys, freq=freq)

    elif is_integer_dtype(dtype):
        values = _ensure_int64(values)
        keys, counts = htable.value_count_int64(values, dropna)
    elif is_float_dtype(dtype):
        values = _ensure_float64(values)
        keys, counts = htable.value_count_float64(values, dropna)
    else:
        values = _ensure_object(values)
        mask = isnull(values)
        keys, counts = htable.value_count_object(values, mask)
        if not dropna and mask.any():
            keys = np.insert(keys, 0, np.NaN)
            counts = np.insert(counts, 0, mask.sum())

    return keys, counts


def duplicated(values, keep='first'):
    """
    Return boolean ndarray denoting duplicate values

    .. versionadded:: 0.19.0

    Parameters
    ----------
    keep : {'first', 'last', False}, default 'first'
        - ``first`` : Mark duplicates as ``True`` except for the first
          occurrence.
        - ``last`` : Mark duplicates as ``True`` except for the last
          occurrence.
        - False : Mark all duplicates as ``True``.

    Returns
    -------
    duplicated : ndarray
    """

    dtype = values.dtype

    # no need to revert to original type
    if needs_i8_conversion(dtype):
        values = values.view(np.int64)
    elif is_period_arraylike(values):
        from pandas.tseries.period import PeriodIndex
        values = PeriodIndex(values).asi8
    elif is_categorical_dtype(dtype):
        values = values.values.codes
    elif isinstance(values, (ABCSeries, ABCIndex)):
        values = values.values

    if is_integer_dtype(dtype):
        values = _ensure_int64(values)
        duplicated = htable.duplicated_int64(values, keep=keep)
    elif is_float_dtype(dtype):
        values = _ensure_float64(values)
        duplicated = htable.duplicated_float64(values, keep=keep)
    else:
        values = _ensure_object(values)
        duplicated = htable.duplicated_object(values, keep=keep)

    return duplicated


def mode(values):
    """Returns the mode or mode(s) of the passed Series or ndarray (sorted)"""
    # must sort because hash order isn't necessarily defined.
    from pandas.core.series import Series

    if isinstance(values, Series):
        constructor = values._constructor
        values = values.values
    else:
        values = np.asanyarray(values)
        constructor = Series

    dtype = values.dtype
    if is_integer_dtype(values):
        values = _ensure_int64(values)
        result = constructor(sorted(htable.mode_int64(values)), dtype=dtype)

    elif issubclass(values.dtype.type, (np.datetime64, np.timedelta64)):
        dtype = values.dtype
        values = values.view(np.int64)
        result = constructor(sorted(htable.mode_int64(values)), dtype=dtype)

    elif is_categorical_dtype(values):
        result = constructor(values.mode())
    else:
        mask = isnull(values)
        values = _ensure_object(values)
        res = htable.mode_object(values, mask)
        try:
            res = sorted(res)
        except TypeError as e:
            warn("Unable to sort modes: %s" % e)
        result = constructor(res, dtype=dtype)

    return result


def rank(values, axis=0, method='average', na_option='keep',
         ascending=True, pct=False):
    """

    """
    if values.ndim == 1:
        f, values = _get_data_algo(values, _rank1d_functions)
        ranks = f(values, ties_method=method, ascending=ascending,
                  na_option=na_option, pct=pct)
    elif values.ndim == 2:
        f, values = _get_data_algo(values, _rank2d_functions)
        ranks = f(values, axis=axis, ties_method=method,
                  ascending=ascending, na_option=na_option, pct=pct)

    return ranks

_rank1d_functions = {
    'float64': algos.rank_1d_float64,
    'int64': algos.rank_1d_int64,
    'generic': algos.rank_1d_generic
}

_rank2d_functions = {
    'float64': algos.rank_2d_float64,
    'int64': algos.rank_2d_int64,
    'generic': algos.rank_2d_generic
}


def quantile(x, q, interpolation_method='fraction'):
    """
    Compute sample quantile or quantiles of the input array. For example, q=0.5
    computes the median.

    The `interpolation_method` parameter supports three values, namely
    `fraction` (default), `lower` and `higher`. Interpolation is done only,
    if the desired quantile lies between two data points `i` and `j`. For
    `fraction`, the result is an interpolated value between `i` and `j`;
    for `lower`, the result is `i`, for `higher` the result is `j`.

    Parameters
    ----------
    x : ndarray
        Values from which to extract score.
    q : scalar or array
        Percentile at which to extract score.
    interpolation_method : {'fraction', 'lower', 'higher'}, optional
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

        - fraction: `i + (j - i)*fraction`, where `fraction` is the
                    fractional part of the index surrounded by `i` and `j`.
        -lower: `i`.
        - higher: `j`.

    Returns
    -------
    score : float
        Score at percentile.

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5

    """
    x = np.asarray(x)
    mask = isnull(x)

    x = x[~mask]

    values = np.sort(x)

    def _get_score(at):
        if len(values) == 0:
            return np.nan

        idx = at * (len(values) - 1)
        if idx % 1 == 0:
            score = values[int(idx)]
        else:
            if interpolation_method == 'fraction':
                score = _interpolate(values[int(idx)], values[int(idx) + 1],
                                     idx % 1)
            elif interpolation_method == 'lower':
                score = values[np.floor(idx)]
            elif interpolation_method == 'higher':
                score = values[np.ceil(idx)]
            else:
                raise ValueError("interpolation_method can only be 'fraction' "
                                 ", 'lower' or 'higher'")

        return score

    if is_scalar(q):
        return _get_score(q)
    else:
        q = np.asarray(q, np.float64)
        return algos.arrmap_float64(q, _get_score)


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a) * fraction


def nsmallest(arr, n, keep='first'):
    """
    Find the indices of the n smallest values of a numpy array.

    Note: Fails silently with NaN.
    """
    if keep == 'last':
        arr = arr[::-1]

    narr = len(arr)
    n = min(n, narr)

    sdtype = str(arr.dtype)
    arr = arr.view(_dtype_map.get(sdtype, sdtype))

    kth_val = algos.kth_smallest(arr.copy(), n - 1)
    return _finalize_nsmallest(arr, kth_val, n, keep, narr)


def nlargest(arr, n, keep='first'):
    """
    Find the indices of the n largest values of a numpy array.

    Note: Fails silently with NaN.
    """
    sdtype = str(arr.dtype)
    arr = arr.view(_dtype_map.get(sdtype, sdtype))
    return nsmallest(-arr, n, keep=keep)


def select_n_slow(dropped, n, keep, method):
    reverse_it = (keep == 'last' or method == 'nlargest')
    ascending = method == 'nsmallest'
    slc = np.s_[::-1] if reverse_it else np.s_[:]
    return dropped[slc].sort_values(ascending=ascending).head(n)


_select_methods = {'nsmallest': nsmallest, 'nlargest': nlargest}


def select_n_series(series, n, keep, method):
    """Implement n largest/smallest for pandas Series

    Parameters
    ----------
    series : pandas.Series object
    n : int
    keep : {'first', 'last'}, default 'first'
    method : str, {'nlargest', 'nsmallest'}

    Returns
    -------
    nordered : Series
    """
    dtype = series.dtype
    if not issubclass(dtype.type, (np.integer, np.floating, np.datetime64,
                                   np.timedelta64)):
        raise TypeError("Cannot use method %r with dtype %s" % (method, dtype))

    if keep not in ('first', 'last'):
        raise ValueError('keep must be either "first", "last"')

    if n <= 0:
        return series[[]]

    dropped = series.dropna()

    if n >= len(series):
        return select_n_slow(dropped, n, keep, method)

    inds = _select_methods[method](dropped.values, n, keep)
    return dropped.iloc[inds]


def select_n_frame(frame, columns, n, method, keep):
    """Implement n largest/smallest for pandas DataFrame

    Parameters
    ----------
    frame : pandas.DataFrame object
    columns : list or str
    n : int
    keep : {'first', 'last'}, default 'first'
    method : str, {'nlargest', 'nsmallest'}

    Returns
    -------
    nordered : DataFrame
    """
    from pandas.core.series import Series
    if not is_list_like(columns):
        columns = [columns]
    columns = list(columns)
    ser = getattr(frame[columns[0]], method)(n, keep=keep)
    if isinstance(ser, Series):
        ser = ser.to_frame()
    return ser.merge(frame, on=columns[0], left_index=True)[frame.columns]


def _finalize_nsmallest(arr, kth_val, n, keep, narr):
    ns, = np.nonzero(arr <= kth_val)
    inds = ns[arr[ns].argsort(kind='mergesort')][:n]
    if keep == 'last':
        # reverse indices
        return narr - 1 - inds
    else:
        return inds

_dtype_map = {'datetime64[ns]': 'int64', 'timedelta64[ns]': 'int64'}


# ------- #
# helpers #
# ------- #

def _hashtable_algo(f, dtype, return_dtype=None):
    """
    f(HashTable, type_caster) -> result
    """
    if is_float_dtype(dtype):
        return f(htable.Float64HashTable, _ensure_float64)
    elif is_integer_dtype(dtype):
        return f(htable.Int64HashTable, _ensure_int64)
    elif is_datetime64_dtype(dtype):
        return_dtype = return_dtype or 'M8[ns]'
        return f(htable.Int64HashTable, _ensure_int64).view(return_dtype)
    elif is_timedelta64_dtype(dtype):
        return_dtype = return_dtype or 'm8[ns]'
        return f(htable.Int64HashTable, _ensure_int64).view(return_dtype)
    else:
        return f(htable.PyObjectHashTable, _ensure_object)

_hashtables = {
    'float64': (htable.Float64HashTable, htable.Float64Vector),
    'int64': (htable.Int64HashTable, htable.Int64Vector),
    'generic': (htable.PyObjectHashTable, htable.ObjectVector)
}


def _get_data_algo(values, func_map):
    if is_float_dtype(values):
        f = func_map['float64']
        values = _ensure_float64(values)

    elif needs_i8_conversion(values):
        f = func_map['int64']
        values = values.view('i8')

    elif is_integer_dtype(values):
        f = func_map['int64']
        values = _ensure_int64(values)
    else:
        f = func_map['generic']
        values = _ensure_object(values)
    return f, values


# ---- #
# take #
# ---- #


def _view_wrapper(f, arr_dtype=None, out_dtype=None, fill_wrap=None):
    def wrapper(arr, indexer, out, fill_value=np.nan):
        if arr_dtype is not None:
            arr = arr.view(arr_dtype)
        if out_dtype is not None:
            out = out.view(out_dtype)
        if fill_wrap is not None:
            fill_value = fill_wrap(fill_value)
        f(arr, indexer, out, fill_value=fill_value)

    return wrapper


def _convert_wrapper(f, conv_dtype):
    def wrapper(arr, indexer, out, fill_value=np.nan):
        arr = arr.astype(conv_dtype)
        f(arr, indexer, out, fill_value=fill_value)

    return wrapper


def _take_2d_multi_generic(arr, indexer, out, fill_value, mask_info):
    # this is not ideal, performance-wise, but it's better than raising
    # an exception (best to optimize in Cython to avoid getting here)
    row_idx, col_idx = indexer
    if mask_info is not None:
        (row_mask, col_mask), (row_needs, col_needs) = mask_info
    else:
        row_mask = row_idx == -1
        col_mask = col_idx == -1
        row_needs = row_mask.any()
        col_needs = col_mask.any()
    if fill_value is not None:
        if row_needs:
            out[row_mask, :] = fill_value
        if col_needs:
            out[:, col_mask] = fill_value
    for i in range(len(row_idx)):
        u_ = row_idx[i]
        for j in range(len(col_idx)):
            v = col_idx[j]
            out[i, j] = arr[u_, v]


def _take_nd_generic(arr, indexer, out, axis, fill_value, mask_info):
    if mask_info is not None:
        mask, needs_masking = mask_info
    else:
        mask = indexer == -1
        needs_masking = mask.any()
    if arr.dtype != out.dtype:
        arr = arr.astype(out.dtype)
    if arr.shape[axis] > 0:
        arr.take(_ensure_platform_int(indexer), axis=axis, out=out)
    if needs_masking:
        outindexer = [slice(None)] * arr.ndim
        outindexer[axis] = mask
        out[tuple(outindexer)] = fill_value


_take_1d_dict = {
    ('int8', 'int8'): algos.take_1d_int8_int8,
    ('int8', 'int32'): algos.take_1d_int8_int32,
    ('int8', 'int64'): algos.take_1d_int8_int64,
    ('int8', 'float64'): algos.take_1d_int8_float64,
    ('int16', 'int16'): algos.take_1d_int16_int16,
    ('int16', 'int32'): algos.take_1d_int16_int32,
    ('int16', 'int64'): algos.take_1d_int16_int64,
    ('int16', 'float64'): algos.take_1d_int16_float64,
    ('int32', 'int32'): algos.take_1d_int32_int32,
    ('int32', 'int64'): algos.take_1d_int32_int64,
    ('int32', 'float64'): algos.take_1d_int32_float64,
    ('int64', 'int64'): algos.take_1d_int64_int64,
    ('int64', 'float64'): algos.take_1d_int64_float64,
    ('float32', 'float32'): algos.take_1d_float32_float32,
    ('float32', 'float64'): algos.take_1d_float32_float64,
    ('float64', 'float64'): algos.take_1d_float64_float64,
    ('object', 'object'): algos.take_1d_object_object,
    ('bool', 'bool'): _view_wrapper(algos.take_1d_bool_bool, np.uint8,
                                    np.uint8),
    ('bool', 'object'): _view_wrapper(algos.take_1d_bool_object, np.uint8,
                                      None),
    ('datetime64[ns]', 'datetime64[ns]'): _view_wrapper(
        algos.take_1d_int64_int64, np.int64, np.int64, np.int64)
}

_take_2d_axis0_dict = {
    ('int8', 'int8'): algos.take_2d_axis0_int8_int8,
    ('int8', 'int32'): algos.take_2d_axis0_int8_int32,
    ('int8', 'int64'): algos.take_2d_axis0_int8_int64,
    ('int8', 'float64'): algos.take_2d_axis0_int8_float64,
    ('int16', 'int16'): algos.take_2d_axis0_int16_int16,
    ('int16', 'int32'): algos.take_2d_axis0_int16_int32,
    ('int16', 'int64'): algos.take_2d_axis0_int16_int64,
    ('int16', 'float64'): algos.take_2d_axis0_int16_float64,
    ('int32', 'int32'): algos.take_2d_axis0_int32_int32,
    ('int32', 'int64'): algos.take_2d_axis0_int32_int64,
    ('int32', 'float64'): algos.take_2d_axis0_int32_float64,
    ('int64', 'int64'): algos.take_2d_axis0_int64_int64,
    ('int64', 'float64'): algos.take_2d_axis0_int64_float64,
    ('float32', 'float32'): algos.take_2d_axis0_float32_float32,
    ('float32', 'float64'): algos.take_2d_axis0_float32_float64,
    ('float64', 'float64'): algos.take_2d_axis0_float64_float64,
    ('object', 'object'): algos.take_2d_axis0_object_object,
    ('bool', 'bool'): _view_wrapper(algos.take_2d_axis0_bool_bool, np.uint8,
                                    np.uint8),
    ('bool', 'object'): _view_wrapper(algos.take_2d_axis0_bool_object,
                                      np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_axis0_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}

_take_2d_axis1_dict = {
    ('int8', 'int8'): algos.take_2d_axis1_int8_int8,
    ('int8', 'int32'): algos.take_2d_axis1_int8_int32,
    ('int8', 'int64'): algos.take_2d_axis1_int8_int64,
    ('int8', 'float64'): algos.take_2d_axis1_int8_float64,
    ('int16', 'int16'): algos.take_2d_axis1_int16_int16,
    ('int16', 'int32'): algos.take_2d_axis1_int16_int32,
    ('int16', 'int64'): algos.take_2d_axis1_int16_int64,
    ('int16', 'float64'): algos.take_2d_axis1_int16_float64,
    ('int32', 'int32'): algos.take_2d_axis1_int32_int32,
    ('int32', 'int64'): algos.take_2d_axis1_int32_int64,
    ('int32', 'float64'): algos.take_2d_axis1_int32_float64,
    ('int64', 'int64'): algos.take_2d_axis1_int64_int64,
    ('int64', 'float64'): algos.take_2d_axis1_int64_float64,
    ('float32', 'float32'): algos.take_2d_axis1_float32_float32,
    ('float32', 'float64'): algos.take_2d_axis1_float32_float64,
    ('float64', 'float64'): algos.take_2d_axis1_float64_float64,
    ('object', 'object'): algos.take_2d_axis1_object_object,
    ('bool', 'bool'): _view_wrapper(algos.take_2d_axis1_bool_bool, np.uint8,
                                    np.uint8),
    ('bool', 'object'): _view_wrapper(algos.take_2d_axis1_bool_object,
                                      np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_axis1_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}

_take_2d_multi_dict = {
    ('int8', 'int8'): algos.take_2d_multi_int8_int8,
    ('int8', 'int32'): algos.take_2d_multi_int8_int32,
    ('int8', 'int64'): algos.take_2d_multi_int8_int64,
    ('int8', 'float64'): algos.take_2d_multi_int8_float64,
    ('int16', 'int16'): algos.take_2d_multi_int16_int16,
    ('int16', 'int32'): algos.take_2d_multi_int16_int32,
    ('int16', 'int64'): algos.take_2d_multi_int16_int64,
    ('int16', 'float64'): algos.take_2d_multi_int16_float64,
    ('int32', 'int32'): algos.take_2d_multi_int32_int32,
    ('int32', 'int64'): algos.take_2d_multi_int32_int64,
    ('int32', 'float64'): algos.take_2d_multi_int32_float64,
    ('int64', 'int64'): algos.take_2d_multi_int64_int64,
    ('int64', 'float64'): algos.take_2d_multi_int64_float64,
    ('float32', 'float32'): algos.take_2d_multi_float32_float32,
    ('float32', 'float64'): algos.take_2d_multi_float32_float64,
    ('float64', 'float64'): algos.take_2d_multi_float64_float64,
    ('object', 'object'): algos.take_2d_multi_object_object,
    ('bool', 'bool'): _view_wrapper(algos.take_2d_multi_bool_bool, np.uint8,
                                    np.uint8),
    ('bool', 'object'): _view_wrapper(algos.take_2d_multi_bool_object,
                                      np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_multi_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}


def _get_take_nd_function(ndim, arr_dtype, out_dtype, axis=0, mask_info=None):
    if ndim <= 2:
        tup = (arr_dtype.name, out_dtype.name)
        if ndim == 1:
            func = _take_1d_dict.get(tup, None)
        elif ndim == 2:
            if axis == 0:
                func = _take_2d_axis0_dict.get(tup, None)
            else:
                func = _take_2d_axis1_dict.get(tup, None)
        if func is not None:
            return func

        tup = (out_dtype.name, out_dtype.name)
        if ndim == 1:
            func = _take_1d_dict.get(tup, None)
        elif ndim == 2:
            if axis == 0:
                func = _take_2d_axis0_dict.get(tup, None)
            else:
                func = _take_2d_axis1_dict.get(tup, None)
        if func is not None:
            func = _convert_wrapper(func, out_dtype)
            return func

    def func(arr, indexer, out, fill_value=np.nan):
        indexer = _ensure_int64(indexer)
        _take_nd_generic(arr, indexer, out, axis=axis, fill_value=fill_value,
                         mask_info=mask_info)

    return func


def take_nd(arr, indexer, axis=0, out=None, fill_value=np.nan, mask_info=None,
            allow_fill=True):
    """
    Specialized Cython take which sets NaN values in one pass

    Parameters
    ----------
    arr : ndarray
        Input array
    indexer : ndarray
        1-D array of indices to take, subarrays corresponding to -1 value
        indicies are filed with fill_value
    axis : int, default 0
        Axis to take from
    out : ndarray or None, default None
        Optional output array, must be appropriate type to hold input and
        fill_value together, if indexer has any -1 value entries; call
        _maybe_promote to determine this type for any fill_value
    fill_value : any, default np.nan
        Fill value to replace -1 values with
    mask_info : tuple of (ndarray, boolean)
        If provided, value should correspond to:
            (indexer != -1, (indexer != -1).any())
        If not provided, it will be computed internally if necessary
    allow_fill : boolean, default True
        If False, indexer is assumed to contain no -1 values so no filling
        will be done.  This short-circuits computation of a mask.  Result is
        undefined if allow_fill == False and -1 is present in indexer.
    """

    # dispatch to internal type takes
    if is_categorical(arr):
        return arr.take_nd(indexer, fill_value=fill_value,
                           allow_fill=allow_fill)
    elif is_datetimetz(arr):
        return arr.take(indexer, fill_value=fill_value, allow_fill=allow_fill)

    if indexer is None:
        indexer = np.arange(arr.shape[axis], dtype=np.int64)
        dtype, fill_value = arr.dtype, arr.dtype.type()
    else:
        indexer = _ensure_int64(indexer)
        if not allow_fill:
            dtype, fill_value = arr.dtype, arr.dtype.type()
            mask_info = None, False
        else:
            # check for promotion based on types only (do this first because
            # it's faster than computing a mask)
            dtype, fill_value = _maybe_promote(arr.dtype, fill_value)
            if dtype != arr.dtype and (out is None or out.dtype != dtype):
                # check if promotion is actually required based on indexer
                if mask_info is not None:
                    mask, needs_masking = mask_info
                else:
                    mask = indexer == -1
                    needs_masking = mask.any()
                    mask_info = mask, needs_masking
                if needs_masking:
                    if out is not None and out.dtype != dtype:
                        raise TypeError('Incompatible type for fill_value')
                else:
                    # if not, then depromote, set fill_value to dummy
                    # (it won't be used but we don't want the cython code
                    # to crash when trying to cast it to dtype)
                    dtype, fill_value = arr.dtype, arr.dtype.type()

    flip_order = False
    if arr.ndim == 2:
        if arr.flags.f_contiguous:
            flip_order = True

    if flip_order:
        arr = arr.T
        axis = arr.ndim - axis - 1
        if out is not None:
            out = out.T

    # at this point, it's guaranteed that dtype can hold both the arr values
    # and the fill_value
    if out is None:
        out_shape = list(arr.shape)
        out_shape[axis] = len(indexer)
        out_shape = tuple(out_shape)
        if arr.flags.f_contiguous and axis == arr.ndim - 1:
            # minor tweak that can make an order-of-magnitude difference
            # for dataframes initialized directly from 2-d ndarrays
            # (s.t. df.values is c-contiguous and df._data.blocks[0] is its
            # f-contiguous transpose)
            out = np.empty(out_shape, dtype=dtype, order='F')
        else:
            out = np.empty(out_shape, dtype=dtype)

    func = _get_take_nd_function(arr.ndim, arr.dtype, out.dtype, axis=axis,
                                 mask_info=mask_info)
    indexer = _ensure_int64(indexer)
    func(arr, indexer, out, fill_value)

    if flip_order:
        out = out.T
    return out


take_1d = take_nd


def take_2d_multi(arr, indexer, out=None, fill_value=np.nan, mask_info=None,
                  allow_fill=True):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    if indexer is None or (indexer[0] is None and indexer[1] is None):
        row_idx = np.arange(arr.shape[0], dtype=np.int64)
        col_idx = np.arange(arr.shape[1], dtype=np.int64)
        indexer = row_idx, col_idx
        dtype, fill_value = arr.dtype, arr.dtype.type()
    else:
        row_idx, col_idx = indexer
        if row_idx is None:
            row_idx = np.arange(arr.shape[0], dtype=np.int64)
        else:
            row_idx = _ensure_int64(row_idx)
        if col_idx is None:
            col_idx = np.arange(arr.shape[1], dtype=np.int64)
        else:
            col_idx = _ensure_int64(col_idx)
        indexer = row_idx, col_idx
        if not allow_fill:
            dtype, fill_value = arr.dtype, arr.dtype.type()
            mask_info = None, False
        else:
            # check for promotion based on types only (do this first because
            # it's faster than computing a mask)
            dtype, fill_value = _maybe_promote(arr.dtype, fill_value)
            if dtype != arr.dtype and (out is None or out.dtype != dtype):
                # check if promotion is actually required based on indexer
                if mask_info is not None:
                    (row_mask, col_mask), (row_needs, col_needs) = mask_info
                else:
                    row_mask = row_idx == -1
                    col_mask = col_idx == -1
                    row_needs = row_mask.any()
                    col_needs = col_mask.any()
                    mask_info = (row_mask, col_mask), (row_needs, col_needs)
                if row_needs or col_needs:
                    if out is not None and out.dtype != dtype:
                        raise TypeError('Incompatible type for fill_value')
                else:
                    # if not, then depromote, set fill_value to dummy
                    # (it won't be used but we don't want the cython code
                    # to crash when trying to cast it to dtype)
                    dtype, fill_value = arr.dtype, arr.dtype.type()

    # at this point, it's guaranteed that dtype can hold both the arr values
    # and the fill_value
    if out is None:
        out_shape = len(row_idx), len(col_idx)
        out = np.empty(out_shape, dtype=dtype)

    func = _take_2d_multi_dict.get((arr.dtype.name, out.dtype.name), None)
    if func is None and arr.dtype != out.dtype:
        func = _take_2d_multi_dict.get((out.dtype.name, out.dtype.name), None)
        if func is not None:
            func = _convert_wrapper(func, out.dtype)
    if func is None:

        def func(arr, indexer, out, fill_value=np.nan):
            _take_2d_multi_generic(arr, indexer, out, fill_value=fill_value,
                                   mask_info=mask_info)

    func(arr, indexer, out=out, fill_value=fill_value)
    return out


# ---- #
# diff #
# ---- #

_diff_special = {
    'float64': algos.diff_2d_float64,
    'float32': algos.diff_2d_float32,
    'int64': algos.diff_2d_int64,
    'int32': algos.diff_2d_int32,
    'int16': algos.diff_2d_int16,
    'int8': algos.diff_2d_int8,
}


def diff(arr, n, axis=0):
    """ difference of n between self,
        analagoust to s-s.shift(n) """

    n = int(n)
    na = np.nan
    dtype = arr.dtype
    is_timedelta = False
    if needs_i8_conversion(arr):
        dtype = np.float64
        arr = arr.view('i8')
        na = tslib.iNaT
        is_timedelta = True
    elif issubclass(dtype.type, np.integer):
        dtype = np.float64
    elif issubclass(dtype.type, np.bool_):
        dtype = np.object_

    dtype = np.dtype(dtype)
    out_arr = np.empty(arr.shape, dtype=dtype)

    na_indexer = [slice(None)] * arr.ndim
    na_indexer[axis] = slice(None, n) if n >= 0 else slice(n, None)
    out_arr[tuple(na_indexer)] = na

    if arr.ndim == 2 and arr.dtype.name in _diff_special:
        f = _diff_special[arr.dtype.name]
        f(arr, out_arr, n, axis)
    else:
        res_indexer = [slice(None)] * arr.ndim
        res_indexer[axis] = slice(n, None) if n >= 0 else slice(None, n)
        res_indexer = tuple(res_indexer)

        lag_indexer = [slice(None)] * arr.ndim
        lag_indexer[axis] = slice(None, -n) if n > 0 else slice(-n, None)
        lag_indexer = tuple(lag_indexer)

        # need to make sure that we account for na for datelike/timedelta
        # we don't actually want to subtract these i8 numbers
        if is_timedelta:
            res = arr[res_indexer]
            lag = arr[lag_indexer]

            mask = (arr[res_indexer] == na) | (arr[lag_indexer] == na)
            if mask.any():
                res = res.copy()
                res[mask] = 0
                lag = lag.copy()
                lag[mask] = 0

            result = res - lag
            result[mask] = na
            out_arr[res_indexer] = result
        else:
            out_arr[res_indexer] = arr[res_indexer] - arr[lag_indexer]

    if is_timedelta:
        from pandas import TimedeltaIndex
        out_arr = TimedeltaIndex(out_arr.ravel().astype('int64')).asi8.reshape(
            out_arr.shape).astype('timedelta64[ns]')

    return out_arr
