import sys

import numpy as np

from pandas.core.common import isnull, notnull
import pandas.core.common as com
import pandas.core.config as cf
import pandas.lib as lib
import pandas.algos as algos
import pandas.hashtable as _hash
import pandas.tslib as tslib

try:
    import bottleneck as bn
    _USE_BOTTLENECK = True
except ImportError:  # pragma: no cover
    _USE_BOTTLENECK = False


def _bottleneck_switch(bn_name, alt, zero_value=None, **kwargs):
    try:
        bn_func = getattr(bn, bn_name)
    except (AttributeError, NameError):  # pragma: no cover
        bn_func = None

    def f(values, axis=None, skipna=True, **kwds):
        if len(kwargs) > 0:
            for k, v in kwargs.iteritems():
                if k not in kwds:
                    kwds[k] = v
        try:
            if zero_value is not None and values.size == 0:
                if values.ndim == 1:
                    return 0
                else:
                    result_shape = values.shape[:
                                                axis] + values.shape[axis + 1:]
                    result = np.empty(result_shape)
                    result.fill(0)
                    return result

            if _USE_BOTTLENECK and skipna and _bn_ok_dtype(values.dtype):
                result = bn_func(values, axis=axis, **kwds)
                # prefer to treat inf/-inf as NA
                if _has_infs(result):
                    result = alt(values, axis=axis, skipna=skipna, **kwds)
            else:
                result = alt(values, axis=axis, skipna=skipna, **kwds)
        except Exception:
            result = alt(values, axis=axis, skipna=skipna, **kwds)

        return result

    return f


def _bn_ok_dtype(dt):
    # Bottleneck chokes on datetime64
    return dt != np.object_ and not issubclass(dt.type, (np.datetime64,np.timedelta64))


def _has_infs(result):
    if isinstance(result, np.ndarray):
        if result.dtype == 'f8':
            return lib.has_infs_f8(result)
        elif result.dtype == 'f4':
            return lib.has_infs_f4(result)
        else:  # pragma: no cover
            raise TypeError('Only suppose float32/64 here')
    else:
        return np.isinf(result) or np.isneginf(result)

def _get_fill_value(dtype, fill_value=None, fill_value_typ=None):
    """ return the correct fill value for the dtype of the values """
    if fill_value is not None:
        return fill_value
    if _na_ok_dtype(dtype):
        if fill_value_typ is None:
            return np.nan
        else:
            if fill_value_typ == '+inf':
                return np.inf
            else:
                return -np.inf
    else:
        if fill_value_typ is None:
            return tslib.iNaT
        else:
            if fill_value_typ == '+inf':
                # need the max int here
                return np.iinfo(np.int64).max
            else:
                return tslib.iNaT

def _get_values(values, skipna, fill_value=None, fill_value_typ=None, isfinite=False, copy=True):
    """ utility to get the values view, mask, dtype
        if necessary copy and mask using the specified fill_value
        copy = True will force the copy """
    if isfinite:
        mask = _isfinite(values)
    else:
        mask = isnull(values)

    dtype    = values.dtype
    dtype_ok = _na_ok_dtype(dtype)

    # get our fill value (in case we need to provide an alternative dtype for it)
    fill_value = _get_fill_value(dtype, fill_value=fill_value, fill_value_typ=fill_value_typ)

    if skipna:
        if copy:
            values = values.copy()
        if dtype_ok:
            np.putmask(values, mask, fill_value)

        # promote if needed
        else:
            values, changed = com._maybe_upcast_putmask(values, mask, fill_value)

    elif copy:
        values = values.copy()

    values = _view_if_needed(values)
    return values, mask, dtype

def _isfinite(values):
    if issubclass(values.dtype.type, (np.timedelta64,np.datetime64)):
        return isnull(values)
    return -np.isfinite(values)

def _na_ok_dtype(dtype):
    return not issubclass(dtype.type, (np.integer, np.datetime64, np.timedelta64))

def _view_if_needed(values):
    if issubclass(values.dtype.type, (np.datetime64,np.timedelta64)):
        return values.view(np.int64)
    return values

def _wrap_results(result,dtype):
    """ wrap our results if needed """

    if issubclass(dtype.type, np.datetime64):
        if not isinstance(result, np.ndarray):
            result = lib.Timestamp(result)
        else:
            result = result.view(dtype)
    elif issubclass(dtype.type, np.timedelta64):
        if not isinstance(result, np.ndarray):

            # this is a scalar timedelta result!
            # we have series convert then take the element (scalar)
            # as series will do the right thing in py3 (and deal with numpy 1.6.2
            # bug in that it results dtype of timedelta64[us]
            from pandas import Series
            result = Series([result],dtype='timedelta64[ns]')
        else:
            result = result.view(dtype)

    return result

def nanany(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, False, copy=skipna)
    return values.any(axis)

def nanall(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, True, copy=skipna)
    return values.all(axis)

def _nansum(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, 0)
    the_sum = values.sum(axis)
    the_sum = _maybe_null_out(the_sum, axis, mask)
    return the_sum

def _nanmean(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, 0)
    the_sum = _ensure_numeric(values.sum(axis))
    count = _get_counts(mask, axis)

    if axis is not None:
        the_mean = the_sum / count
        ct_mask = count == 0
        if ct_mask.any():
            the_mean[ct_mask] = np.nan
    else:
        the_mean = the_sum / count if count > 0 else np.nan
    return the_mean


def _nanmedian(values, axis=None, skipna=True):
    def get_median(x):
        mask = notnull(x)
        if not skipna and not mask.all():
            return np.nan
        return algos.median(x[mask])

    if values.dtype != np.float64:
        values = values.astype('f8')

    if values.ndim > 1:
        return np.apply_along_axis(get_median, axis, values)
    else:
        return get_median(values)


def _nanvar(values, axis=None, skipna=True, ddof=1):
    if not isinstance(values.dtype.type, np.floating):
        values = values.astype('f8')

    mask = isnull(values)

    if axis is not None:
        count = (values.shape[axis] - mask.sum(axis)).astype(float)
    else:
        count = float(values.size - mask.sum())

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    X = _ensure_numeric(values.sum(axis))
    XX = _ensure_numeric((values ** 2).sum(axis))
    return np.fabs((XX - X ** 2 / count) / (count - ddof))


def _nanmin(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, fill_value_typ = '+inf')

    # numpy 1.6.1 workaround in Python 3.x
    if (values.dtype == np.object_
            and sys.version_info[0] >= 3):  # pragma: no cover
        import __builtin__
        if values.ndim > 1:
            apply_ax = axis if axis is not None else 0
            result = np.apply_along_axis(__builtin__.min, apply_ax, values)
        else:
            result = __builtin__.min(values)
    else:
        if ((axis is not None and values.shape[axis] == 0)
                or values.size == 0):
            result = com.ensure_float(values.sum(axis))
            result.fill(np.nan)
        else:
            result = values.min(axis)

    result = _wrap_results(result,dtype)
    return _maybe_null_out(result, axis, mask)


def _nanmax(values, axis=None, skipna=True):
    values, mask, dtype = _get_values(values, skipna, fill_value_typ ='-inf')

    # numpy 1.6.1 workaround in Python 3.x
    if (values.dtype == np.object_
            and sys.version_info[0] >= 3):  # pragma: no cover
        import __builtin__

        if values.ndim > 1:
            apply_ax = axis if axis is not None else 0
            result = np.apply_along_axis(__builtin__.max, apply_ax, values)
        else:
            result = __builtin__.max(values)
    else:
        if ((axis is not None and values.shape[axis] == 0)
                or values.size == 0):
            result = com.ensure_float(values.sum(axis))
            result.fill(np.nan)
        else:
            result = values.max(axis)

    result = _wrap_results(result,dtype)
    return _maybe_null_out(result, axis, mask)


def nanargmax(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    values, mask, dtype = _get_values(values, skipna, fill_value_typ = '-inf', isfinite=True)
    result = values.argmax(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result


def nanargmin(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    values, mask, dtype = _get_values(values, skipna, fill_value_typ = '+inf', isfinite=True)
    result = values.argmin(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result

nansum = _bottleneck_switch('nansum', _nansum, zero_value=0)
nanmean = _bottleneck_switch('nanmean', _nanmean)
nanmedian = _bottleneck_switch('nanmedian', _nanmedian)
nanvar = _bottleneck_switch('nanvar', _nanvar, ddof=1)
nanmin = _bottleneck_switch('nanmin', _nanmin)
nanmax = _bottleneck_switch('nanmax', _nanmax)


def nanskew(values, axis=None, skipna=True):
    if not isinstance(values.dtype.type, np.floating):
        values = values.astype('f8')

    mask = isnull(values)
    count = _get_counts(mask, axis)

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    A = values.sum(axis) / count
    B = (values ** 2).sum(axis) / count - A ** 2
    C = (values ** 3).sum(axis) / count - A ** 3 - 3 * A * B

    # floating point error
    B = _zero_out_fperr(B)
    C = _zero_out_fperr(C)

    result = ((np.sqrt((count ** 2 - count)) * C) /
              ((count - 2) * np.sqrt(B) ** 3))

    if isinstance(result, np.ndarray):
        result = np.where(B == 0, 0, result)
        result[count < 3] = np.nan
        return result
    else:
        result = 0 if B == 0 else result
        if count < 3:
            return np.nan
        return result


def nankurt(values, axis=None, skipna=True):
    if not isinstance(values.dtype.type, np.floating):
        values = values.astype('f8')

    mask = isnull(values)
    count = _get_counts(mask, axis)

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    A = values.sum(axis) / count
    B = (values ** 2).sum(axis) / count - A ** 2
    C = (values ** 3).sum(axis) / count - A ** 3 - 3 * A * B
    D = (values ** 4).sum(axis) / count - A ** 4 - 6 * B * A * A - 4 * C * A

    B = _zero_out_fperr(B)
    C = _zero_out_fperr(C)
    D = _zero_out_fperr(D)

    result = (((count * count - 1.) * D / (B * B) - 3 * ((count - 1.) ** 2)) /
              ((count - 2.) * (count - 3.)))
    if isinstance(result, np.ndarray):
        result = np.where(B == 0, 0, result)
        result[count < 4] = np.nan
        return result
    else:
        result = 0 if B == 0 else result
        if count < 4:
            return np.nan
        return result


def nanprod(values, axis=None, skipna=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        values[mask] = 1
    result = values.prod(axis)
    return _maybe_null_out(result, axis, mask)


def _maybe_arg_null_out(result, axis, mask, skipna):
    # helper function for nanargmin/nanargmax
    if axis is None:
        if skipna:
            if mask.all():
                result = -1
        else:
            if mask.any():
                result = -1
    else:
        if skipna:
            na_mask = mask.all(axis)
        else:
            na_mask = mask.any(axis)
        if na_mask.any():
            result[na_mask] = -1
    return result


def _get_counts(mask, axis):
    if axis is not None:
        count = (mask.shape[axis] - mask.sum(axis)).astype(float)
    else:
        count = float(mask.size - mask.sum())

    return count


def _maybe_null_out(result, axis, mask):
    if axis is not None:
        null_mask = (mask.shape[axis] - mask.sum(axis)) == 0
        if null_mask.any():
            result = result.astype('f8')
            result[null_mask] = np.nan
    else:
        null_mask = mask.size - mask.sum()
        if null_mask == 0:
            result = np.nan

    return result


def _zero_out_fperr(arg):
    if isinstance(arg, np.ndarray):
        return np.where(np.abs(arg) < 1e-14, 0, arg)
    else:
        return 0 if np.abs(arg) < 1e-14 else arg


def nancorr(a, b, method='pearson', min_periods=None):
    """
    a, b: ndarrays
    """
    if len(a) != len(b):
        raise AssertionError('Operands to nancorr must have same size')

    if min_periods is None:
        min_periods = 1

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) < min_periods:
        return np.nan

    f = get_corr_func(method)
    return f(a, b)


def get_corr_func(method):
    if method in ['kendall', 'spearman']:
        from scipy.stats import kendalltau, spearmanr

    def _pearson(a, b):
        return np.corrcoef(a, b)[0, 1]

    def _kendall(a, b):
        rs = kendalltau(a, b)
        if isinstance(rs, tuple):
            return rs[0]
        return rs

    def _spearman(a, b):
        return spearmanr(a, b)[0]

    _cor_methods = {
        'pearson': _pearson,
        'kendall': _kendall,
        'spearman': _spearman
    }
    return _cor_methods[method]


def nancov(a, b, min_periods=None):
    if len(a) != len(b):
        raise AssertionError('Operands to nancov must have same size')

    if min_periods is None:
        min_periods = 1

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) < min_periods:
        return np.nan

    return np.cov(a, b)[0, 1]


def _ensure_numeric(x):
    if isinstance(x, np.ndarray):
        if x.dtype == np.object_:
            x = x.astype(np.float64)
    elif not (com.is_float(x) or com.is_integer(x)):
        try:
            x = float(x)
        except Exception:
            raise TypeError('Could not convert %s to numeric' % str(x))

    return x

# NA-friendly array comparisons

import operator


def make_nancomp(op):
    def f(x, y):
        xmask = isnull(x)
        ymask = isnull(y)
        mask = xmask | ymask

        result = op(x, y)

        if mask.any():
            if result.dtype == np.bool_:
                result = result.astype('O')
            np.putmask(result, mask, np.nan)

        return result
    return f

nangt = make_nancomp(operator.gt)
nange = make_nancomp(operator.ge)
nanlt = make_nancomp(operator.lt)
nanle = make_nancomp(operator.le)
naneq = make_nancomp(operator.eq)
nanne = make_nancomp(operator.ne)


def unique1d(values):
    """
    Hash table-based unique
    """
    if np.issubdtype(values.dtype, np.floating):
        table = _hash.Float64HashTable(len(values))
        uniques = np.array(table.unique(com._ensure_float64(values)),
                           dtype=np.float64)
    elif np.issubdtype(values.dtype, np.datetime64):
        table = _hash.Int64HashTable(len(values))
        uniques = table.unique(com._ensure_int64(values))
        uniques = uniques.view('M8[ns]')
    elif np.issubdtype(values.dtype, np.integer):
        table = _hash.Int64HashTable(len(values))
        uniques = table.unique(com._ensure_int64(values))
    else:
        table = _hash.PyObjectHashTable(len(values))
        uniques = table.unique(com._ensure_object(values))
    return uniques
