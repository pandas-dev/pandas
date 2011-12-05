import numpy as np

from pandas.core.common import isnull, notnull
import pandas._tseries as lib

try:
    import bottleneck as bn
    _USE_BOTTLENECK = True
except ImportError:  # pragma: no cover
    _USE_BOTTLENECK = False

def nansum(values, axis=None, skipna=True, copy=True):
    if values.dtype == np.object_:
        the_sum = values.sum(axis)
    else:
        mask = isnull(values)

        if skipna and not issubclass(values.dtype.type, np.integer):
            if copy:
                values = values.copy()
            np.putmask(values, mask, 0)

        the_sum = values.sum(axis)
        the_sum = _maybe_null_out(the_sum, axis, mask)

    return the_sum

def nanmean(values, axis=None, skipna=True, copy=True):
    if values.dtype == np.object_:
        the_mean = values.sum(axis) / float(values.shape[axis])
    else:
        mask = isnull(values)

        if skipna and not issubclass(values.dtype.type, np.integer):
            if copy:
                values = values.copy()
            np.putmask(values, mask, 0)

        the_sum = values.sum(axis)
        count = _get_counts(mask, axis)

        if axis is not None:
            the_mean = the_sum / count
            ct_mask = count == 0
            if ct_mask.any():
                the_mean[ct_mask] = np.nan
        else:
            the_mean = the_sum / count if count > 0 else np.nan

    return the_mean

def nanmedian(values, axis=None, skipna=True, copy=True):
    def get_median(x):
        mask = notnull(x)
        if not skipna and not mask.all():
            return np.nan
        return lib.median(x[mask])

    if values.dtype != np.float64:
        values = values.astype('f8')

    if axis == 0:
        values = values.T

    if values.ndim > 1:
        return np.asarray([get_median(arr) for arr in values])
    else:
        return get_median(values)

def nanvar(values, axis=None, skipna=True, copy=True, ddof=1):
    mask = isnull(values)

    if axis is not None:
        count = (values.shape[axis] - mask.sum(axis)).astype(float)
    else:
        count = float(values.size - mask.sum())

    if skipna:
        if copy:
            values = values.copy()
        np.putmask(values, mask, 0)

    X = values.sum(axis)
    XX = (values ** 2).sum(axis)
    return (XX - X ** 2 / count) / (count - ddof)

def nanskew(values, axis=None, skipna=True, copy=True):
    if not isinstance(values.dtype.type, np.floating):
        values = values.astype('f8')

    mask = isnull(values)
    count = _get_counts(mask, axis)

    if skipna:
        if copy:
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
        return np.where(B == 0, 0, result)
    else:
        return 0 if B == 0 else result

def nanmin(values, axis=None, skipna=True, copy=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        if copy:
            values = values.copy()
        np.putmask(values, mask, np.inf)
    # numpy 1.6.1 workaround in Python 3.x
    if values.dtype == np.object_:  # pragma: no cover
        import __builtin__
        result = np.apply_along_axis(__builtin__.min, axis, values)
    else:
        result = values.min(axis)

    return _maybe_null_out(result, axis, mask)

def nanmax(values, axis=None, skipna=True, copy=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        if copy:
            values = values.copy()
        np.putmask(values, mask, -np.inf)
    # numpy 1.6.1 workaround in Python 3.x
    if values.dtype == np.object_:  # pragma: no cover
        import __builtin__
        result = np.apply_along_axis(__builtin__.max, axis, values)
    else:
        result = values.max(axis)
    return _maybe_null_out(result, axis, mask)

def nanprod(values, axis=None, skipna=True, copy=True):
    mask = isnull(values)
    if skipna and not issubclass(values.dtype.type, np.integer):
        if copy:
            values = values.copy()
        values[mask] = 1
    result = values.prod(axis)
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
    if len(a) == 0:
        return np.nan

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]

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
    if len(a) == 0:
        return np.nan

    valid = notnull(a) & notnull(b)
    if not valid.all():
        a = a[valid]
        b = b[valid]
    return np.cov(a, b)[0, 1]
