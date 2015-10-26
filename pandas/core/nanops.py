import itertools
import functools
import numpy as np

try:
    import bottleneck as bn
    _USE_BOTTLENECK = True
except ImportError:  # pragma: no cover
    _USE_BOTTLENECK = False

import pandas.hashtable as _hash
from pandas import compat, lib, algos, tslib
from pandas.compat import builtins
from pandas.core.common import (isnull, notnull, _values_from_object,
                                _maybe_upcast_putmask,
                                ensure_float, _ensure_float64,
                                _ensure_int64, _ensure_object,
                                is_float, is_integer, is_complex,
                                is_float_dtype,
                                is_complex_dtype, is_integer_dtype,
                                is_bool_dtype, is_object_dtype,
                                is_datetime64_dtype, is_timedelta64_dtype,
                                is_datetime_or_timedelta_dtype, _get_dtype,
                                is_int_or_datetime_dtype, is_any_int_dtype,
                                _int64_max)


class disallow(object):

    def __init__(self, *dtypes):
        super(disallow, self).__init__()
        self.dtypes = tuple(np.dtype(dtype).type for dtype in dtypes)

    def check(self, obj):
        return hasattr(obj, 'dtype') and issubclass(obj.dtype.type,
                                                    self.dtypes)

    def __call__(self, f):
        @functools.wraps(f)
        def _f(*args, **kwargs):
            obj_iter = itertools.chain(args, compat.itervalues(kwargs))
            if any(self.check(obj) for obj in obj_iter):
                raise TypeError('reduction operation {0!r} not allowed for '
                                'this dtype'.format(f.__name__.replace('nan',
                                                                       '')))
            try:
                return f(*args, **kwargs)
            except ValueError as e:
                # we want to transform an object array
                # ValueError message to the more typical TypeError
                # e.g. this is normally a disallowed function on
                # object arrays that contain strings
                if is_object_dtype(args[0]):
                    raise TypeError(e)
                raise
        return _f


class bottleneck_switch(object):

    def __init__(self, zero_value=None, **kwargs):
        self.zero_value = zero_value
        self.kwargs = kwargs

    def __call__(self, alt):
        bn_name = alt.__name__

        try:
            bn_func = getattr(bn, bn_name)
        except (AttributeError, NameError):  # pragma: no cover
            bn_func = None

        @functools.wraps(alt)
        def f(values, axis=None, skipna=True, **kwds):
            if len(self.kwargs) > 0:
                for k, v in compat.iteritems(self.kwargs):
                    if k not in kwds:
                        kwds[k] = v
            try:
                if self.zero_value is not None and values.size == 0:
                    if values.ndim == 1:

                        # wrap the 0's if needed
                        if is_timedelta64_dtype(values):
                            return lib.Timedelta(0)
                        return 0
                    else:
                        result_shape = (values.shape[:axis] +
                                        values.shape[axis + 1:])
                        result = np.empty(result_shape)
                        result.fill(0)
                        return result

                if _USE_BOTTLENECK and skipna and _bn_ok_dtype(values.dtype,
                                                               bn_name):
                    result = bn_func(values, axis=axis, **kwds)

                    # prefer to treat inf/-inf as NA, but must compute the func
                    # twice :(
                    if _has_infs(result):
                        result = alt(values, axis=axis, skipna=skipna, **kwds)
                else:
                    result = alt(values, axis=axis, skipna=skipna, **kwds)
            except Exception:
                try:
                    result = alt(values, axis=axis, skipna=skipna, **kwds)
                except ValueError as e:
                    # we want to transform an object array
                    # ValueError message to the more typical TypeError
                    # e.g. this is normally a disallowed function on
                    # object arrays that contain strings

                    if is_object_dtype(values):
                        raise TypeError(e)
                    raise

            return result

        return f


def _bn_ok_dtype(dt, name):
    # Bottleneck chokes on datetime64
    if (not is_object_dtype(dt) and
            not is_datetime_or_timedelta_dtype(dt)):

        # bottleneck does not properly upcast during the sum
        # so can overflow
        if name == 'nansum':
            if dt.itemsize < 8:
                return False

        return True
    return False


def _has_infs(result):
    if isinstance(result, np.ndarray):
        if result.dtype == 'f8':
            return lib.has_infs_f8(result.ravel())
        elif result.dtype == 'f4':
            return lib.has_infs_f4(result.ravel())
    try:
        return np.isinf(result).any()
    except (TypeError, NotImplementedError) as e:
        # if it doesn't support infs, then it can't have infs
        return False


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
                return _int64_max
            else:
                return tslib.iNaT


def _get_values(values, skipna, fill_value=None, fill_value_typ=None,
                isfinite=False, copy=True):
    """ utility to get the values view, mask, dtype
        if necessary copy and mask using the specified fill_value
        copy = True will force the copy """
    values = _values_from_object(values)
    if isfinite:
        mask = _isfinite(values)
    else:
        mask = isnull(values)

    dtype = values.dtype
    dtype_ok = _na_ok_dtype(dtype)

    # get our fill value (in case we need to provide an alternative
    # dtype for it)
    fill_value = _get_fill_value(dtype, fill_value=fill_value,
                                 fill_value_typ=fill_value_typ)

    if skipna:
        if copy:
            values = values.copy()
        if dtype_ok:
            np.putmask(values, mask, fill_value)

        # promote if needed
        else:
            values, changed = _maybe_upcast_putmask(values, mask, fill_value)

    elif copy:
        values = values.copy()

    values = _view_if_needed(values)

    # return a platform independent precision dtype
    dtype_max = dtype
    if is_integer_dtype(dtype) or is_bool_dtype(dtype):
        dtype_max = np.int64
    elif is_float_dtype(dtype):
        dtype_max = np.float64

    return values, mask, dtype, dtype_max


def _isfinite(values):
    if is_datetime_or_timedelta_dtype(values):
        return isnull(values)
    if (is_complex_dtype(values) or is_float_dtype(values) or
            is_integer_dtype(values) or is_bool_dtype(values)):
        return ~np.isfinite(values)
    return ~np.isfinite(values.astype('float64'))


def _na_ok_dtype(dtype):
    return not is_int_or_datetime_dtype(dtype)


def _view_if_needed(values):
    if is_datetime_or_timedelta_dtype(values):
        return values.view(np.int64)
    return values


def _wrap_results(result, dtype):
    """ wrap our results if needed """

    if is_datetime64_dtype(dtype):
        if not isinstance(result, np.ndarray):
            result = lib.Timestamp(result)
        else:
            result = result.view(dtype)
    elif is_timedelta64_dtype(dtype):
        if not isinstance(result, np.ndarray):

            # raise if we have a timedelta64[ns] which is too large
            if np.fabs(result) > _int64_max:
                raise ValueError("overflow in timedelta operation")

            result = lib.Timedelta(result, unit='ns')
        else:
            result = result.astype('i8').view(dtype)

    return result


def nanany(values, axis=None, skipna=True):
    values, mask, dtype, _ = _get_values(values, skipna, False, copy=skipna)
    return values.any(axis)


def nanall(values, axis=None, skipna=True):
    values, mask, dtype, _ = _get_values(values, skipna, True, copy=skipna)
    return values.all(axis)


@disallow('M8')
@bottleneck_switch(zero_value=0)
def nansum(values, axis=None, skipna=True):
    values, mask, dtype, dtype_max = _get_values(values, skipna, 0)
    dtype_sum = dtype_max
    if is_float_dtype(dtype):
        dtype_sum = dtype
    elif is_timedelta64_dtype(dtype):
        dtype_sum = np.float64
    the_sum = values.sum(axis, dtype=dtype_sum)
    the_sum = _maybe_null_out(the_sum, axis, mask)

    return _wrap_results(the_sum, dtype)


@disallow('M8')
@bottleneck_switch()
def nanmean(values, axis=None, skipna=True):
    values, mask, dtype, dtype_max = _get_values(values, skipna, 0)

    dtype_sum = dtype_max
    dtype_count = np.float64
    if is_integer_dtype(dtype) or is_timedelta64_dtype(dtype):
        dtype_sum = np.float64
    elif is_float_dtype(dtype):
        dtype_sum = dtype
        dtype_count = dtype
    count = _get_counts(mask, axis, dtype=dtype_count)
    the_sum = _ensure_numeric(values.sum(axis, dtype=dtype_sum))

    if axis is not None and getattr(the_sum, 'ndim', False):
        the_mean = the_sum / count
        ct_mask = count == 0
        if ct_mask.any():
            the_mean[ct_mask] = np.nan
    else:
        the_mean = the_sum / count if count > 0 else np.nan

    return _wrap_results(the_mean, dtype)


@disallow('M8')
@bottleneck_switch()
def nanmedian(values, axis=None, skipna=True):

    values, mask, dtype, dtype_max = _get_values(values, skipna)

    def get_median(x):
        mask = notnull(x)
        if not skipna and not mask.all():
            return np.nan
        return algos.median(_values_from_object(x[mask]))

    if not is_float_dtype(values):
        values = values.astype('f8')
        values[mask] = np.nan

    if axis is None:
        values = values.ravel()

    notempty = values.size

    # an array from a frame
    if values.ndim > 1:
        # there's a non-empty array to apply over otherwise numpy raises
        if notempty:
            return _wrap_results(np.apply_along_axis(get_median, axis, values), dtype)

        # must return the correct shape, but median is not defined for the
        # empty set so return nans of shape "everything but the passed axis"
        # since "axis" is where the reduction would occur if we had a nonempty
        # array
        shp = np.array(values.shape)
        dims = np.arange(values.ndim)
        ret = np.empty(shp[dims != axis])
        ret.fill(np.nan)
        return _wrap_results(ret, dtype)

    # otherwise return a scalar value
    return _wrap_results(get_median(values) if notempty else np.nan, dtype)


def _get_counts_nanvar(mask, axis, ddof, dtype=float):
    dtype = _get_dtype(dtype)
    count = _get_counts(mask, axis, dtype=dtype)
    d = count - dtype.type(ddof)

    # always return NaN, never inf
    if np.isscalar(count):
        if count <= ddof:
            count = np.nan
            d = np.nan
    else:
        mask2 = count <= ddof
        if mask2.any():
            np.putmask(d, mask2, np.nan)
            np.putmask(count, mask2, np.nan)
    return count, d


@disallow('M8')
@bottleneck_switch(ddof=1)
def nanstd(values, axis=None, skipna=True, ddof=1):
    result = np.sqrt(nanvar(values, axis=axis, skipna=skipna, ddof=ddof))
    return _wrap_results(result, values.dtype)


@disallow('M8')
@bottleneck_switch(ddof=1)
def nanvar(values, axis=None, skipna=True, ddof=1):

    dtype = values.dtype
    mask = isnull(values)
    if is_any_int_dtype(values):
        values = values.astype('f8')
        values[mask] = np.nan

    if is_float_dtype(values):
        count, d = _get_counts_nanvar(mask, axis, ddof, values.dtype)
    else:
        count, d = _get_counts_nanvar(mask, axis, ddof)

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    # xref GH10242
    # Compute variance via two-pass algorithm, which is stable against
    # cancellation errors and relatively accurate for small numbers of
    # observations.
    #
    # See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    avg = _ensure_numeric(values.sum(axis=axis, dtype=np.float64)) / count
    if axis is not None:
        avg = np.expand_dims(avg, axis)
    sqr = _ensure_numeric((avg - values) ** 2)
    np.putmask(sqr, mask, 0)
    result = sqr.sum(axis=axis, dtype=np.float64) / d

    # Return variance as np.float64 (the datatype used in the accumulator),
    # unless we were dealing with a float array, in which case use the same
    # precision as the original values array.
    if is_float_dtype(dtype):
        result = result.astype(dtype)
    return _wrap_results(result, values.dtype)


@disallow('M8', 'm8')
def nansem(values, axis=None, skipna=True, ddof=1):
    var = nanvar(values, axis, skipna, ddof=ddof)

    mask = isnull(values)
    if not is_float_dtype(values.dtype):
        values = values.astype('f8')
    count, _ = _get_counts_nanvar(mask, axis, ddof, values.dtype)
    var = nanvar(values, axis, skipna, ddof=ddof)

    return np.sqrt(var) / np.sqrt(count)


def _nanminmax(meth, fill_value_typ):
    @bottleneck_switch()
    def reduction(values, axis=None, skipna=True):
        values, mask, dtype, dtype_max = _get_values(
            values,
            skipna,
            fill_value_typ=fill_value_typ,
        )

        if ((axis is not None and values.shape[axis] == 0)
                or values.size == 0):
            try:
                result = getattr(values, meth)(axis, dtype=dtype_max)
                result.fill(np.nan)
            except:
                result = np.nan
        else:
            result = getattr(values, meth)(axis)

        result = _wrap_results(result, dtype)
        return _maybe_null_out(result, axis, mask)

    reduction.__name__ = 'nan' + meth
    return reduction


nanmin = _nanminmax('min', fill_value_typ='+inf')
nanmax = _nanminmax('max', fill_value_typ='-inf')


def nanargmax(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    values, mask, dtype, _ = _get_values(values, skipna, fill_value_typ='-inf',
                                         isfinite=True)
    result = values.argmax(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result


def nanargmin(values, axis=None, skipna=True):
    """
    Returns -1 in the NA case
    """
    values, mask, dtype, _ = _get_values(values, skipna, fill_value_typ='+inf',
                                         isfinite=True)
    result = values.argmin(axis)
    result = _maybe_arg_null_out(result, axis, mask, skipna)
    return result


@disallow('M8','m8')
def nanskew(values, axis=None, skipna=True):

    mask = isnull(values)
    if not is_float_dtype(values.dtype):
        values = values.astype('f8')
        count = _get_counts(mask, axis)
    else:
        count = _get_counts(mask, axis, dtype=values.dtype)

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    typ = values.dtype.type
    A = values.sum(axis) / count
    B = (values ** 2).sum(axis) / count - A ** typ(2)
    C = (values ** 3).sum(axis) / count - A ** typ(3) - typ(3) * A * B

    # floating point error
    B = _zero_out_fperr(B)
    C = _zero_out_fperr(C)

    result = ((np.sqrt(count * count - count) * C) /
              ((count - typ(2)) * np.sqrt(B) ** typ(3)))

    if isinstance(result, np.ndarray):
        result = np.where(B == 0, 0, result)
        result[count < 3] = np.nan
        return result
    else:
        result = 0 if B == 0 else result
        if count < 3:
            return np.nan
        return result


@disallow('M8','m8')
def nankurt(values, axis=None, skipna=True):

    mask = isnull(values)
    if not is_float_dtype(values.dtype):
        values = values.astype('f8')
        count = _get_counts(mask, axis)
    else:
        count = _get_counts(mask, axis, dtype=values.dtype)

    if skipna:
        values = values.copy()
        np.putmask(values, mask, 0)

    typ = values.dtype.type
    A = values.sum(axis) / count
    B = (values ** 2).sum(axis) / count - A ** typ(2)
    C = (values ** 3).sum(axis) / count - A ** typ(3) - typ(3) * A * B
    D = (values ** 4).sum(axis) / count - A ** typ(4) - typ(6) * B * A * A - typ(4) * C * A

    B = _zero_out_fperr(B)
    D = _zero_out_fperr(D)

    if not isinstance(B, np.ndarray):
        # if B is a scalar, check these corner cases first before doing division
        if count < 4:
            return np.nan
        if B == 0:
            return 0

    result = (((count * count - typ(1)) * D / (B * B) - typ(3) * ((count - typ(1)) ** typ(2))) /
              ((count - typ(2)) * (count - typ(3))))

    if isinstance(result, np.ndarray):
        result = np.where(B == 0, 0, result)
        result[count < 4] = np.nan

    return result


@disallow('M8','m8')
def nanprod(values, axis=None, skipna=True):
    mask = isnull(values)
    if skipna and not is_any_int_dtype(values):
        values = values.copy()
        values[mask] = 1
    result = values.prod(axis)
    return _maybe_null_out(result, axis, mask)


def _maybe_arg_null_out(result, axis, mask, skipna):
    # helper function for nanargmin/nanargmax
    if axis is None or not getattr(result, 'ndim', False):
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


def _get_counts(mask, axis, dtype=float):
    dtype = _get_dtype(dtype)
    if axis is None:
        return dtype.type(mask.size - mask.sum())

    count = mask.shape[axis] - mask.sum(axis)
    if np.isscalar(count):
        return dtype.type(count)
    try:
        return count.astype(dtype)
    except AttributeError:
        return np.array(count, dtype=dtype)


def _maybe_null_out(result, axis, mask):
    if axis is not None and getattr(result, 'ndim', False):
        null_mask = (mask.shape[axis] - mask.sum(axis)) == 0
        if np.any(null_mask):
            if np.iscomplexobj(result):
                result = result.astype('c16')
            else:
                result = result.astype('f8')
            result[null_mask] = np.nan
    elif result is not tslib.NaT:
        null_mask = mask.size - mask.sum()
        if null_mask == 0:
            result = np.nan

    return result


def _zero_out_fperr(arg):
    if isinstance(arg, np.ndarray):
        return np.where(np.abs(arg) < 1e-14, 0, arg)
    else:
        return arg.dtype.type(0) if np.abs(arg) < 1e-14 else arg


@disallow('M8','m8')
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


@disallow('M8','m8')
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
        if is_integer_dtype(x) or is_bool_dtype(x):
            x = x.astype(np.float64)
        elif is_object_dtype(x):
            try:
                x = x.astype(np.complex128)
            except:
                x = x.astype(np.float64)
            else:
                if not np.any(x.imag):
                    x = x.real
    elif not (is_float(x) or is_integer(x) or is_complex(x)):
        try:
            x = float(x)
        except Exception:
            try:
                x = complex(x)
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
            if is_bool_dtype(result):
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
        uniques = np.array(table.unique(_ensure_float64(values)),
                           dtype=np.float64)
    elif np.issubdtype(values.dtype, np.datetime64):
        table = _hash.Int64HashTable(len(values))
        uniques = table.unique(_ensure_int64(values))
        uniques = uniques.view('M8[ns]')
    elif np.issubdtype(values.dtype, np.timedelta64):
        table = _hash.Int64HashTable(len(values))
        uniques = table.unique(_ensure_int64(values))
        uniques = uniques.view('m8[ns]')
    elif np.issubdtype(values.dtype, np.integer):
        table = _hash.Int64HashTable(len(values))
        uniques = table.unique(_ensure_int64(values))
    else:
        table = _hash.PyObjectHashTable(len(values))
        uniques = table.unique(_ensure_object(values))
    return uniques
