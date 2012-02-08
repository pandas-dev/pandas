import sys

import numpy as np

from pandas.core.common import isnull, notnull
import pandas.core.common as com
import pandas._tseries as lib

try:
    import bottleneck as bn
    _USE_BOTTLENECK = True
except ImportError:  # pragma: no cover
    _USE_BOTTLENECK = False

def _bottleneck_switch(bn_name, alt, **kwargs):
    try:
        bn_func = getattr(bn, bn_name)
    except (AttributeError, NameError):  # pragma: no cover
        bn_func = None
    def f(values, axis=None, skipna=True):
        try:
            if _USE_BOTTLENECK and skipna and values.dtype != np.object_:
                result = bn_func(values, axis=axis, **kwargs)
                # prefer to treat inf/-inf as NA
                if _has_infs(result):
                    result = alt(values, axis=axis, skipna=skipna, **kwargs)
            else:
                result = alt(values, axis=axis, skipna=skipna, **kwargs)
        except Exception:
            result = alt(values, axis=axis, skipna=skipna, **kwargs)

        return result

    return f

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

def _nansum(values, axis=None, skipna=True):
    mask = isnull(values)

    if skipna and not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, 0)

    the_sum = values.sum(axis)
    the_sum = _maybe_null_out(the_sum, axis, mask)

    return the_sum

def _nanmean(values, axis=None, skipna=True):
    mask = isnull(values)

    if skipna and not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, 0)

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
        return lib.median(x[mask])

    if values.dtype != np.float64:
        values = values.astype('f8')

    if values.ndim > 1:
        return np.apply_along_axis(get_median, axis, values)
    else:
        return get_median(values)

def _nanvar(values, axis=None, skipna=True, ddof=1):
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
    return (XX - X ** 2 / count) / (count - ddof)

def _nanmin(values, axis=None, skipna=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, np.inf)
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
        result = values.min(axis)

    return _maybe_null_out(result, axis, mask)

def _nanmax(values, axis=None, skipna=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, -np.inf)
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
        result = values.max(axis)
    return _maybe_null_out(result, axis, mask)

def nanargmax(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    mask = -np.isfinite(values)
    if not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, -np.inf)
    result = values.argmax(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result

def nanargmin(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    mask = -np.isfinite(values)
    if not issubclass(values.dtype.type, np.integer):
        values = values.copy()
        np.putmask(values, mask, np.inf)
    result = values.argmin(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result

nansum = _bottleneck_switch('nansum', _nansum)
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

def nancorr(a, b, method='pearson'):
    """
    a, b: ndarrays
    """
    assert(len(a) == len(b))

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) == 0:
        return np.nan

    f = get_corr_func(method)
    return f(a, b)

def get_corr_func(method):
    if method in ['kendall', 'spearman']:
        from scipy.stats import kendalltau, spearmanr

    def _pearson(a, b):
        return np.corrcoef(a, b)[0, 1]
    def _kendall(a, b):
        return kendalltau(a, b)[0]
    def _spearman(a, b):
        return spearmanr(a, b)[0]

    _cor_methods = {
        'pearson' : _pearson,
        'kendall' : _kendall,
        'spearman' : _spearman
    }
    return _cor_methods[method]

def nancov(a, b):
    assert(len(a) == len(b))

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

    if len(a) == 0:
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
    if issubclass(values.dtype.type, np.floating):
        if values.dtype != np.float64:
            values = values.astype(np.float64)
        table = lib.Float64HashTable(len(values))
        uniques = np.array(table.unique(values), dtype=np.float64)
    elif issubclass(values.dtype.type, np.integer):
        if values.dtype != np.int64:
            values = values.astype(np.int64)
        table = lib.Int64HashTable(len(values))
        uniques = np.array(table.unique(values), dtype=np.int64)
    else:
        if not values.dtype == np.object_:
            values = values.astype(np.object_)
        table = lib.PyObjectHashTable(len(values))
        uniques = lib.list_to_object_array(table.unique(values))
        uniques = lib.maybe_convert_objects(uniques)
    return uniques

