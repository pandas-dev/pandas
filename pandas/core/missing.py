"""
Routines for filling missing data
"""

from functools import partial

import numpy as np

import pandas as pd
import pandas.core.common as com
import pandas.algos as algos
import pandas.lib as lib
from pandas.compat import range


def _clean_fill_method(method, allow_nearest=False):
    if method is None:
        return None
    method = method.lower()
    if method == 'ffill':
        method = 'pad'
    if method == 'bfill':
        method = 'backfill'

    valid_methods = ['pad', 'backfill']
    expecting = 'pad (ffill) or backfill (bfill)'
    if allow_nearest:
        valid_methods.append('nearest')
        expecting = 'pad (ffill), backfill (bfill) or nearest'
    if method not in valid_methods:
        msg = ('Invalid fill method. Expecting %s. Got %s'
               % (expecting, method))
        raise ValueError(msg)
    return method


def _clean_interp_method(method, **kwargs):
    order = kwargs.get('order')
    valid = ['linear', 'time', 'index', 'values', 'nearest', 'zero', 'slinear',
             'quadratic', 'cubic', 'barycentric', 'polynomial',
             'krogh', 'piecewise_polynomial',
             'pchip', 'spline']
    if method in ('spline', 'polynomial') and order is None:
        raise ValueError("You must specify the order of the spline or "
                         "polynomial.")
    if method not in valid:
        raise ValueError("method must be one of {0}."
                         "Got '{1}' instead.".format(valid, method))
    return method


def interpolate_1d(xvalues, yvalues, method='linear', limit=None,
                   limit_direction='forward',
                   fill_value=None, bounds_error=False, order=None, **kwargs):
    """
    Logic for the 1-d interpolation.  The result should be 1-d, inputs
    xvalues and yvalues will each be 1-d arrays of the same length.

    Bounds_error is currently hardcoded to False since non-scipy ones don't
    take it as an argumnet.
    """
    # Treat the original, non-scipy methods first.

    invalid = com.isnull(yvalues)
    valid = ~invalid

    if not valid.any():
        # have to call np.asarray(xvalues) since xvalues could be an Index
        # which cant be mutated
        result = np.empty_like(np.asarray(xvalues), dtype=np.float64)
        result.fill(np.nan)
        return result

    if valid.all():
        return yvalues

    if method == 'time':
        if not getattr(xvalues, 'is_all_dates', None):
        # if not issubclass(xvalues.dtype.type, np.datetime64):
            raise ValueError('time-weighted interpolation only works '
                             'on Series or DataFrames with a '
                             'DatetimeIndex')
        method = 'values'

    def _interp_limit(invalid, fw_limit, bw_limit):
        "Get idx of values that won't be filled b/c they exceed the limits."
        for x in np.where(invalid)[0]:
            if invalid[max(0, x - fw_limit):x + bw_limit + 1].all():
                yield x

    valid_limit_directions = ['forward', 'backward', 'both']
    limit_direction = limit_direction.lower()
    if limit_direction not in valid_limit_directions:
        msg = 'Invalid limit_direction: expecting one of %r, got %r.' % (
            valid_limit_directions, limit_direction)
        raise ValueError(msg)

    from pandas import Series
    ys = Series(yvalues)
    start_nans = set(range(ys.first_valid_index()))
    end_nans = set(range(1 + ys.last_valid_index(), len(valid)))

    # This is a list of the indexes in the series whose yvalue is currently NaN,
    # but whose interpolated yvalue will be overwritten with NaN after computing
    # the interpolation. For each index in this list, one of these conditions is
    # true of the corresponding NaN in the yvalues:
    #
    # a) It is one of a chain of NaNs at the beginning of the series, and either
    #    limit is not specified or limit_direction is 'forward'.
    # b) It is one of a chain of NaNs at the end of the series, and limit is
    #    specified and limit_direction is 'backward' or 'both'.
    # c) Limit is nonzero and it is further than limit from the nearest non-NaN
    #    value (with respect to the limit_direction setting).
    #
    # The default behavior is to fill forward with no limit, ignoring NaNs at
    # the beginning (see issues #9218 and #10420)
    violate_limit = sorted(start_nans)

    if limit:
        if limit_direction == 'forward':
            violate_limit = sorted(start_nans | set(_interp_limit(invalid, limit, 0)))
        if limit_direction == 'backward':
            violate_limit = sorted(end_nans | set(_interp_limit(invalid, 0, limit)))
        if limit_direction == 'both':
            violate_limit = sorted(_interp_limit(invalid, limit, limit))

    xvalues = getattr(xvalues, 'values', xvalues)
    yvalues = getattr(yvalues, 'values', yvalues)
    result = yvalues.copy()

    if method in ['linear', 'time', 'index', 'values']:
        if method in ('values', 'index'):
            inds = np.asarray(xvalues)
            # hack for DatetimeIndex, #1646
            if issubclass(inds.dtype.type, np.datetime64):
                inds = inds.view(np.int64)
            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)
        else:
            inds = xvalues
        result[invalid] = np.interp(inds[invalid], inds[valid], yvalues[valid])
        result[violate_limit] = np.nan
        return result

    sp_methods = ['nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                  'barycentric', 'krogh', 'spline', 'polynomial',
                  'piecewise_polynomial', 'pchip']
    if method in sp_methods:
        inds = np.asarray(xvalues)
        # hack for DatetimeIndex, #1646
        if issubclass(inds.dtype.type, np.datetime64):
            inds = inds.view(np.int64)
        result[invalid] = _interpolate_scipy_wrapper(
            inds[valid], yvalues[valid], inds[invalid], method=method,
            fill_value=fill_value,
            bounds_error=bounds_error, order=order, **kwargs)
        result[violate_limit] = np.nan
        return result


def _interpolate_scipy_wrapper(x, y, new_x, method, fill_value=None,
                               bounds_error=False, order=None, **kwargs):
    """
    passed off to scipy.interpolate.interp1d. method is scipy's kind.
    Returns an array interpolated at new_x.  Add any new methods to
    the list in _clean_interp_method
    """
    try:
        from scipy import interpolate
        from pandas import DatetimeIndex
    except ImportError:
        raise ImportError('{0} interpolation requires Scipy'.format(method))

    new_x = np.asarray(new_x)

    # ignores some kwargs that could be passed along.
    alt_methods = {
        'barycentric': interpolate.barycentric_interpolate,
        'krogh': interpolate.krogh_interpolate,
        'piecewise_polynomial': interpolate.piecewise_polynomial_interpolate,
    }

    if getattr(x, 'is_all_dates', False):
        # GH 5975, scipy.interp1d can't hande datetime64s
        x, new_x = x._values.astype('i8'), new_x.astype('i8')

    try:
        alt_methods['pchip'] = interpolate.pchip_interpolate
    except AttributeError:
        if method == 'pchip':
            raise ImportError("Your version of scipy does not support "
                              "PCHIP interpolation.")

    interp1d_methods = ['nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                        'polynomial']
    if method in interp1d_methods:
        if method == 'polynomial':
            method = order
        terp = interpolate.interp1d(x, y, kind=method, fill_value=fill_value,
                                    bounds_error=bounds_error)
        new_y = terp(new_x)
    elif method == 'spline':
        # GH #10633
        if not order:
            raise ValueError("order needs to be specified and greater than 0")
        terp = interpolate.UnivariateSpline(x, y, k=order, **kwargs)
        new_y = terp(new_x)
    else:
        # GH 7295: need to be able to write for some reason
        # in some circumstances: check all three
        if not x.flags.writeable:
            x = x.copy()
        if not y.flags.writeable:
            y = y.copy()
        if not new_x.flags.writeable:
            new_x = new_x.copy()
        method = alt_methods[method]
        new_y = method(x, y, new_x, **kwargs)
    return new_y


def interpolate_2d(values, method='pad', axis=0, limit=None, fill_value=None, dtype=None):
    """ perform an actual interpolation of values, values will be make 2-d if
    needed fills inplace, returns the result
    """

    transf = (lambda x: x) if axis == 0 else (lambda x: x.T)

    # reshape a 1 dim if needed
    ndim = values.ndim
    if values.ndim == 1:
        if axis != 0:  # pragma: no cover
            raise AssertionError("cannot interpolate on a ndim == 1 with "
                                 "axis != 0")
        values = values.reshape(tuple((1,) + values.shape))

    if fill_value is None:
        mask = None
    else:  # todo create faster fill func without masking
        mask = com.mask_missing(transf(values), fill_value)

    method = _clean_fill_method(method)
    if method == 'pad':
        values = transf(pad_2d(transf(values), limit=limit, mask=mask, dtype=dtype))
    else:
        values = transf(backfill_2d(transf(values), limit=limit, mask=mask, dtype=dtype))

    # reshape back
    if ndim == 1:
        values = values[0]

    return values


def _interp_wrapper(f, wrap_dtype, na_override=None):
    def wrapper(arr, mask, limit=None):
        view = arr.view(wrap_dtype)
        f(view, mask, limit=limit)
    return wrapper


_pad_1d_datetime = _interp_wrapper(algos.pad_inplace_int64, np.int64)
_pad_2d_datetime = _interp_wrapper(algos.pad_2d_inplace_int64, np.int64)
_backfill_1d_datetime = _interp_wrapper(algos.backfill_inplace_int64,
                                        np.int64)
_backfill_2d_datetime = _interp_wrapper(algos.backfill_2d_inplace_int64,
                                        np.int64)


def pad_1d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if com.is_float_dtype(values):
        _method = getattr(algos, 'pad_inplace_%s' % dtype.name, None)
    elif dtype in com._DATELIKE_DTYPES or com.is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif com.is_integer_dtype(values):
        values = com._ensure_float64(values)
        _method = algos.pad_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_1d [%s]' % dtype.name)

    if mask is None:
        mask = com.isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)
    return values


def backfill_1d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if com.is_float_dtype(values):
        _method = getattr(algos, 'backfill_inplace_%s' % dtype.name, None)
    elif dtype in com._DATELIKE_DTYPES or com.is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif com.is_integer_dtype(values):
        values = com._ensure_float64(values)
        _method = algos.backfill_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_1d [%s]' % dtype.name)

    if mask is None:
        mask = com.isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)
    return values


def pad_2d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if com.is_float_dtype(values):
        _method = getattr(algos, 'pad_2d_inplace_%s' % dtype.name, None)
    elif dtype in com._DATELIKE_DTYPES or com.is_datetime64_dtype(values):
        _method = _pad_2d_datetime
    elif com.is_integer_dtype(values):
        values = com._ensure_float64(values)
        _method = algos.pad_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_2d [%s]' % dtype.name)

    if mask is None:
        mask = com.isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def backfill_2d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if com.is_float_dtype(values):
        _method = getattr(algos, 'backfill_2d_inplace_%s' % dtype.name, None)
    elif dtype in com._DATELIKE_DTYPES or com.is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
    elif com.is_integer_dtype(values):
        values = com._ensure_float64(values)
        _method = algos.backfill_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_2d [%s]' % dtype.name)

    if mask is None:
        mask = com.isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


_fill_methods = {'pad': pad_1d, 'backfill': backfill_1d}


def _get_fill_func(method):
    method = _clean_fill_method(method)
    return _fill_methods[method]


def _clean_reindex_fill_method(method):
    return _clean_fill_method(method, allow_nearest=True)
