"""
Misc tools for implementing data structures
"""

import re
import collections
import numbers
from datetime import datetime, timedelta
from functools import partial

import numpy as np
import pandas as pd
import pandas.algos as algos
import pandas.lib as lib
import pandas.tslib as tslib
from pandas import compat
from pandas.compat import (long, zip, map, string_types,
                           iteritems)
from pandas.types import api as gt
from pandas.types.api import *  # noqa
from pandas.core.config import get_option


class PandasError(Exception):
    pass


class PerformanceWarning(Warning):
    pass


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


class AmbiguousIndexError(PandasError, KeyError):
    pass


class AbstractMethodError(NotImplementedError):
    """Raise this error instead of NotImplementedError for abstract methods
    while keeping compatibility with Python 2 and Python 3.
    """

    def __init__(self, class_instance):
        self.class_instance = class_instance

    def __str__(self):
        return ("This method must be defined in the concrete class of %s" %
                self.class_instance.__class__.__name__)


_POSSIBLY_CAST_DTYPES = set([np.dtype(t).name
                             for t in ['O', 'int8', 'uint8', 'int16', 'uint16',
                                       'int32', 'uint32', 'int64', 'uint64']])

_NS_DTYPE = np.dtype('M8[ns]')
_TD_DTYPE = np.dtype('m8[ns]')
_INT64_DTYPE = np.dtype(np.int64)
_DATELIKE_DTYPES = set([np.dtype(t)
                        for t in ['M8[ns]', '<M8[ns]', '>M8[ns]',
                                  'm8[ns]', '<m8[ns]', '>m8[ns]']])
_int8_max = np.iinfo(np.int8).max
_int16_max = np.iinfo(np.int16).max
_int32_max = np.iinfo(np.int32).max
_int64_max = np.iinfo(np.int64).max


def isnull(obj):
    """Detect missing values (NaN in numeric arrays, None/NaN in object arrays)

    Parameters
    ----------
    arr : ndarray or object value
        Object to check for null-ness

    Returns
    -------
    isnulled : array-like of bool or bool
        Array or bool indicating whether an object is null or if an array is
        given which of the element is null.

    See also
    --------
    pandas.notnull: boolean inverse of pandas.isnull
    """
    return _isnull(obj)


def _isnull_new(obj):
    if lib.isscalar(obj):
        return lib.checknull(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, pd.MultiIndex):
        raise NotImplementedError("isnull is not defined for MultiIndex")
    elif isinstance(obj, (gt.ABCSeries, np.ndarray, pd.Index)):
        return _isnull_ndarraylike(obj)
    elif isinstance(obj, gt.ABCGeneric):
        return obj._constructor(obj._data.isnull(func=isnull))
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike(np.asarray(obj))
    else:
        return obj is None


def _isnull_old(obj):
    """Detect missing values. Treat None, NaN, INF, -INF as null.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    """
    if lib.isscalar(obj):
        return lib.checknull_old(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, pd.MultiIndex):
        raise NotImplementedError("isnull is not defined for MultiIndex")
    elif isinstance(obj, (gt.ABCSeries, np.ndarray, pd.Index)):
        return _isnull_ndarraylike_old(obj)
    elif isinstance(obj, gt.ABCGeneric):
        return obj._constructor(obj._data.isnull(func=_isnull_old))
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike_old(np.asarray(obj))
    else:
        return obj is None


_isnull = _isnull_new


def _use_inf_as_null(key):
    """Option change callback for null/inf behaviour
    Choose which replacement for numpy.isnan / -numpy.isfinite is used.

    Parameters
    ----------
    flag: bool
        True means treat None, NaN, INF, -INF as null (old way),
        False means None and NaN are null, but INF, -INF are not null
        (new way).

    Notes
    -----
    This approach to setting global module values is discussed and
    approved here:

    * http://stackoverflow.com/questions/4859217/
      programmatically-creating-variables-in-python/4859312#4859312
    """
    flag = get_option(key)
    if flag:
        globals()['_isnull'] = _isnull_old
    else:
        globals()['_isnull'] = _isnull_new


def _isnull_ndarraylike(obj):

    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if is_string_dtype(dtype):
        if is_categorical_dtype(values):
            from pandas import Categorical
            if not isinstance(values, Categorical):
                values = values.values
            result = values.isnull()
        else:

            # Working around NumPy ticket 1542
            shape = values.shape

            if is_string_like_dtype(dtype):
                result = np.zeros(values.shape, dtype=bool)
            else:
                result = np.empty(shape, dtype=bool)
                vec = lib.isnullobj(values.ravel())
                result[...] = vec.reshape(shape)

    elif is_datetimelike(obj):
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        result = np.isnan(values)

    # box
    if isinstance(obj, gt.ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result


def _isnull_ndarraylike_old(obj):
    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if is_string_dtype(dtype):
        # Working around NumPy ticket 1542
        shape = values.shape

        if is_string_like_dtype(dtype):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj_old(values.ravel())
            result[:] = vec.reshape(shape)

    elif dtype in _DATELIKE_DTYPES:
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        result = ~np.isfinite(values)

    # box
    if isinstance(obj, gt.ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result


def notnull(obj):
    """Replacement for numpy.isfinite / -numpy.isnan which is suitable for use
    on object arrays.

    Parameters
    ----------
    arr : ndarray or object value
        Object to check for *not*-null-ness

    Returns
    -------
    isnulled : array-like of bool or bool
        Array or bool indicating whether an object is *not* null or if an array
        is given which of the element is *not* null.

    See also
    --------
    pandas.isnull : boolean inverse of pandas.notnull
    """
    res = isnull(obj)
    if lib.isscalar(res):
        return not res
    return ~res


def is_null_datelike_scalar(other):
    """ test whether the object is a null datelike, e.g. Nat
    but guard against passing a non-scalar """
    if other is pd.NaT or other is None:
        return True
    elif lib.isscalar(other):

        # a timedelta
        if hasattr(other, 'dtype'):
            return other.view('i8') == tslib.iNaT
        elif is_integer(other) and other == tslib.iNaT:
            return True
        return isnull(other)
    return False


def array_equivalent(left, right, strict_nan=False):
    """
    True if two arrays, left and right, have equal non-NaN elements, and NaNs
    in corresponding locations.  False otherwise. It is assumed that left and
    right are NumPy arrays of the same dtype. The behavior of this function
    (particularly with respect to NaNs) is not defined if the dtypes are
    different.

    Parameters
    ----------
    left, right : ndarrays
    strict_nan : bool, default False
        If True, consider NaN and None to be different.

    Returns
    -------
    b : bool
        Returns True if the arrays are equivalent.

    Examples
    --------
    >>> array_equivalent(
    ...     np.array([1, 2, np.nan]),
    ...     np.array([1, 2, np.nan]))
    True
    >>> array_equivalent(
    ...     np.array([1, np.nan, 2]),
    ...     np.array([1, 2, np.nan]))
    False
    """

    left, right = np.asarray(left), np.asarray(right)

    # shape compat
    if left.shape != right.shape:
        return False

    # Object arrays can contain None, NaN and NaT.
    if is_object_dtype(left) or is_object_dtype(right):

        if not strict_nan:
            # pd.isnull considers NaN and None to be equivalent.
            return lib.array_equivalent_object(
                _ensure_object(left.ravel()), _ensure_object(right.ravel()))

        for left_value, right_value in zip(left, right):
            if left_value is tslib.NaT and right_value is not tslib.NaT:
                return False

            elif isinstance(left_value, float) and np.isnan(left_value):
                if (not isinstance(right_value, float) or
                        not np.isnan(right_value)):
                    return False
            else:
                if left_value != right_value:
                    return False
        return True

    # NaNs can occur in float and complex arrays.
    if is_float_dtype(left) or is_complex_dtype(left):
        return ((left == right) | (np.isnan(left) & np.isnan(right))).all()

    # numpy will will not allow this type of datetimelike vs integer comparison
    elif is_datetimelike_v_numeric(left, right):
        return False

    # M8/m8
    elif needs_i8_conversion(left) and needs_i8_conversion(right):
        if not is_dtype_equal(left.dtype, right.dtype):
            return False

        left = left.view('i8')
        right = right.view('i8')

    # NaNs cannot occur otherwise.
    return np.array_equal(left, right)


def _iterable_not_string(x):
    return (isinstance(x, collections.Iterable) and
            not isinstance(x, compat.string_types))


def flatten(l):
    """Flatten an arbitrarily nested sequence.

    Parameters
    ----------
    l : sequence
        The non string sequence to flatten

    Notes
    -----
    This doesn't consider strings sequences.

    Returns
    -------
    flattened : generator
    """
    for el in l:
        if _iterable_not_string(el):
            for s in flatten(el):
                yield s
        else:
            yield el


def _coerce_indexer_dtype(indexer, categories):
    """ coerce the indexer input array to the smallest dtype possible """
    l = len(categories)
    if l < _int8_max:
        return _ensure_int8(indexer)
    elif l < _int16_max:
        return _ensure_int16(indexer)
    elif l < _int32_max:
        return _ensure_int32(indexer)
    return _ensure_int64(indexer)


def _coerce_to_dtypes(result, dtypes):
    """ given a dtypes and a result set, coerce the result elements to the
    dtypes
    """
    if len(result) != len(dtypes):
        raise AssertionError("_coerce_to_dtypes requires equal len arrays")

    from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type

    def conv(r, dtype):
        try:
            if isnull(r):
                pass
            elif dtype == _NS_DTYPE:
                r = lib.Timestamp(r)
            elif dtype == _TD_DTYPE:
                r = _coerce_scalar_to_timedelta_type(r)
            elif dtype == np.bool_:
                # messy. non 0/1 integers do not get converted.
                if is_integer(r) and r not in [0, 1]:
                    return int(r)
                r = bool(r)
            elif dtype.kind == 'f':
                r = float(r)
            elif dtype.kind == 'i':
                r = int(r)
        except:
            pass

        return r

    return [conv(r, dtype) for r, dtype in zip(result, dtypes)]


def _infer_fill_value(val):
    """
    infer the fill value for the nan/NaT from the provided
    scalar/ndarray/list-like if we are a NaT, return the correct dtyped
    element to provide proper block construction
    """

    if not is_list_like(val):
        val = [val]
    val = np.array(val, copy=False)
    if is_datetimelike(val):
        return np.array('NaT', dtype=val.dtype)
    elif is_object_dtype(val.dtype):
        dtype = lib.infer_dtype(_ensure_object(val))
        if dtype in ['datetime', 'datetime64']:
            return np.array('NaT', dtype=_NS_DTYPE)
        elif dtype in ['timedelta', 'timedelta64']:
            return np.array('NaT', dtype=_TD_DTYPE)
    return np.nan


def _infer_dtype_from_scalar(val):
    """ interpret the dtype from a scalar """

    dtype = np.object_

    # a 1-element ndarray
    if isinstance(val, np.ndarray):
        if val.ndim != 0:
            raise ValueError(
                "invalid ndarray passed to _infer_dtype_from_scalar")

        dtype = val.dtype
        val = val.item()

    elif isinstance(val, compat.string_types):

        # If we create an empty array using a string to infer
        # the dtype, NumPy will only allocate one character per entry
        # so this is kind of bad. Alternately we could use np.repeat
        # instead of np.empty (but then you still don't want things
        # coming out as np.str_!

        dtype = np.object_

    elif isinstance(val, (np.datetime64,
                          datetime)) and getattr(val, 'tzinfo', None) is None:
        val = lib.Timestamp(val).value
        dtype = np.dtype('M8[ns]')

    elif isinstance(val, (np.timedelta64, timedelta)):
        val = lib.Timedelta(val).value
        dtype = np.dtype('m8[ns]')

    elif is_bool(val):
        dtype = np.bool_

    elif is_integer(val):
        if isinstance(val, np.integer):
            dtype = type(val)
        else:
            dtype = np.int64

    elif is_float(val):
        if isinstance(val, np.floating):
            dtype = type(val)
        else:
            dtype = np.float64

    elif is_complex(val):
        dtype = np.complex_

    return dtype, val


def _is_na_compat(arr, fill_value=np.nan):
    """
    Parameters
    ----------
    arr: a numpy array
    fill_value: fill value, default to np.nan

    Returns
    -------
    True if we can fill using this fill_value
    """
    dtype = arr.dtype
    if isnull(fill_value):
        return not (is_bool_dtype(dtype) or
                    is_integer_dtype(dtype))
    return True


def _maybe_fill(arr, fill_value=np.nan):
    """
    if we have a compatiable fill_value and arr dtype, then fill
    """
    if _is_na_compat(arr, fill_value):
        arr.fill(fill_value)
    return arr


def _maybe_promote(dtype, fill_value=np.nan):

    # if we passed an array here, determine the fill value by dtype
    if isinstance(fill_value, np.ndarray):
        if issubclass(fill_value.dtype.type, (np.datetime64, np.timedelta64)):
            fill_value = tslib.iNaT
        else:

            # we need to change to object type as our
            # fill_value is of object type
            if fill_value.dtype == np.object_:
                dtype = np.dtype(np.object_)
            fill_value = np.nan

    # returns tuple of (dtype, fill_value)
    if issubclass(dtype.type, (np.datetime64, np.timedelta64)):
        # for now: refuse to upcast datetime64
        # (this is because datetime64 will not implicitly upconvert
        #  to object correctly as of numpy 1.6.1)
        if isnull(fill_value):
            fill_value = tslib.iNaT
        else:
            if issubclass(dtype.type, np.datetime64):
                try:
                    fill_value = lib.Timestamp(fill_value).value
                except:
                    # the proper thing to do here would probably be to upcast
                    # to object (but numpy 1.6.1 doesn't do this properly)
                    fill_value = tslib.iNaT
            elif issubclass(dtype.type, np.timedelta64):
                try:
                    fill_value = lib.Timedelta(fill_value).value
                except:
                    # as for datetimes, cannot upcast to object
                    fill_value = tslib.iNaT
            else:
                fill_value = tslib.iNaT
    elif is_datetimetz(dtype):
        if isnull(fill_value):
            fill_value = tslib.iNaT
    elif is_float(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, np.integer):
            dtype = np.float64
    elif is_bool(fill_value):
        if not issubclass(dtype.type, np.bool_):
            dtype = np.object_
    elif is_integer(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, np.integer):
            # upcast to prevent overflow
            arr = np.asarray(fill_value)
            if arr != arr.astype(dtype):
                dtype = arr.dtype
    elif is_complex(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, (np.integer, np.floating)):
            dtype = np.complex128
    elif fill_value is None:
        if is_float_dtype(dtype) or is_complex_dtype(dtype):
            fill_value = np.nan
        elif is_integer_dtype(dtype):
            dtype = np.float64
            fill_value = np.nan
        elif is_datetime_or_timedelta_dtype(dtype):
            fill_value = tslib.iNaT
        else:
            dtype = np.object_
    else:
        dtype = np.object_

    # in case we have a string that looked like a number
    if is_categorical_dtype(dtype):
        pass
    elif is_datetimetz(dtype):
        pass
    elif issubclass(np.dtype(dtype).type, compat.string_types):
        dtype = np.object_

    return dtype, fill_value


def _maybe_upcast_putmask(result, mask, other):
    """
    A safe version of putmask that potentially upcasts the result

    Parameters
    ----------
    result : ndarray
        The destination array. This will be mutated in-place if no upcasting is
        necessary.
    mask : boolean ndarray
    other : ndarray or scalar
        The source array or value

    Returns
    -------
    result : ndarray
    changed : boolean
        Set to true if the result array was upcasted
    """

    if mask.any():
        # Two conversions for date-like dtypes that can't be done automatically
        # in np.place:
        #   NaN -> NaT
        #   integer or integer array -> date-like array
        if result.dtype in _DATELIKE_DTYPES:
            if lib.isscalar(other):
                if isnull(other):
                    other = result.dtype.type('nat')
                elif is_integer(other):
                    other = np.array(other, dtype=result.dtype)
            elif is_integer_dtype(other):
                other = np.array(other, dtype=result.dtype)

        def changeit():

            # try to directly set by expanding our array to full
            # length of the boolean
            try:
                om = other[mask]
                om_at = om.astype(result.dtype)
                if (om == om_at).all():
                    new_result = result.values.copy()
                    new_result[mask] = om_at
                    result[:] = new_result
                    return result, False
            except:
                pass

            # we are forced to change the dtype of the result as the input
            # isn't compatible
            r, _ = _maybe_upcast(result, fill_value=other, copy=True)
            np.place(r, mask, other)

            return r, True

        # we want to decide whether place will work
        # if we have nans in the False portion of our mask then we need to
        # upcast (possibly), otherwise we DON't want to upcast (e.g. if we
        # have values, say integers, in the success portion then it's ok to not
        # upcast)
        new_dtype, _ = _maybe_promote(result.dtype, other)
        if new_dtype != result.dtype:

            # we have a scalar or len 0 ndarray
            # and its nan and we are changing some values
            if (lib.isscalar(other) or
                    (isinstance(other, np.ndarray) and other.ndim < 1)):
                if isnull(other):
                    return changeit()

            # we have an ndarray and the masking has nans in it
            else:

                if isnull(other[mask]).any():
                    return changeit()

        try:
            np.place(result, mask, other)
        except:
            return changeit()

    return result, False


def _maybe_upcast(values, fill_value=np.nan, dtype=None, copy=False):
    """ provide explict type promotion and coercion

    Parameters
    ----------
    values : the ndarray that we want to maybe upcast
    fill_value : what we want to fill with
    dtype : if None, then use the dtype of the values, else coerce to this type
    copy : if True always make a copy even if no upcast is required
    """

    if is_extension_type(values):
        if copy:
            values = values.copy()
    else:
        if dtype is None:
            dtype = values.dtype
        new_dtype, fill_value = _maybe_promote(dtype, fill_value)
        if new_dtype != values.dtype:
            values = values.astype(new_dtype)
        elif copy:
            values = values.copy()

    return values, fill_value


def _possibly_cast_item(obj, item, dtype):
    chunk = obj[item]

    if chunk.values.dtype != dtype:
        if dtype in (np.object_, np.bool_):
            obj[item] = chunk.astype(np.object_)
        elif not issubclass(dtype, (np.integer, np.bool_)):  # pragma: no cover
            raise ValueError("Unexpected dtype encountered: %s" % dtype)


def _possibly_downcast_to_dtype(result, dtype):
    """ try to cast to the specified dtype (e.g. convert back to bool/int
    or could be an astype of float64->float32
    """

    if lib.isscalar(result):
        return result

    def trans(x):
        return x

    if isinstance(dtype, compat.string_types):
        if dtype == 'infer':
            inferred_type = lib.infer_dtype(_ensure_object(result.ravel()))
            if inferred_type == 'boolean':
                dtype = 'bool'
            elif inferred_type == 'integer':
                dtype = 'int64'
            elif inferred_type == 'datetime64':
                dtype = 'datetime64[ns]'
            elif inferred_type == 'timedelta64':
                dtype = 'timedelta64[ns]'

            # try to upcast here
            elif inferred_type == 'floating':
                dtype = 'int64'
                if issubclass(result.dtype.type, np.number):

                    def trans(x):  # noqa
                        return x.round()
            else:
                dtype = 'object'

    if isinstance(dtype, compat.string_types):
        dtype = np.dtype(dtype)

    try:

        # don't allow upcasts here (except if empty)
        if dtype.kind == result.dtype.kind:
            if (result.dtype.itemsize <= dtype.itemsize and
                    np.prod(result.shape)):
                return result

        if issubclass(dtype.type, np.floating):
            return result.astype(dtype)
        elif dtype == np.bool_ or issubclass(dtype.type, np.integer):

            # if we don't have any elements, just astype it
            if not np.prod(result.shape):
                return trans(result).astype(dtype)

            # do a test on the first element, if it fails then we are done
            r = result.ravel()
            arr = np.array([r[0]])

            # if we have any nulls, then we are done
            if isnull(arr).any() or not np.allclose(arr,
                                                    trans(arr).astype(dtype)):
                return result

            # a comparable, e.g. a Decimal may slip in here
            elif not isinstance(r[0], (np.integer, np.floating, np.bool, int,
                                       float, bool)):
                return result

            if (issubclass(result.dtype.type, (np.object_, np.number)) and
                    notnull(result).all()):
                new_result = trans(result).astype(dtype)
                try:
                    if np.allclose(new_result, result):
                        return new_result
                except:

                    # comparison of an object dtype with a number type could
                    # hit here
                    if (new_result == result).all():
                        return new_result

        # a datetimelike
        elif dtype.kind in ['M', 'm'] and result.dtype.kind in ['i']:
            try:
                result = result.astype(dtype)
            except:
                if dtype.tz:
                    # convert to datetime and change timezone
                    result = pd.to_datetime(result).tz_localize(dtype.tz)

    except:
        pass

    return result


def _maybe_convert_string_to_object(values):
    """

    Convert string-like and string-like array to convert object dtype.
    This is to avoid numpy to handle the array as str dtype.
    """
    if isinstance(values, string_types):
        values = np.array([values], dtype=object)
    elif (isinstance(values, np.ndarray) and
          issubclass(values.dtype.type, (np.string_, np.unicode_))):
        values = values.astype(object)
    return values


def _maybe_convert_scalar(values):
    """
    Convert a python scalar to the appropriate numpy dtype if possible
    This avoids numpy directly converting according to platform preferences
    """
    if lib.isscalar(values):
        dtype, values = _infer_dtype_from_scalar(values)
        try:
            values = dtype(values)
        except TypeError:
            pass
    return values


def _lcd_dtypes(a_dtype, b_dtype):
    """ return the lcd dtype to hold these types """

    if is_datetime64_dtype(a_dtype) or is_datetime64_dtype(b_dtype):
        return _NS_DTYPE
    elif is_timedelta64_dtype(a_dtype) or is_timedelta64_dtype(b_dtype):
        return _TD_DTYPE
    elif is_complex_dtype(a_dtype):
        if is_complex_dtype(b_dtype):
            return a_dtype
        return np.float64
    elif is_integer_dtype(a_dtype):
        if is_integer_dtype(b_dtype):
            if a_dtype.itemsize == b_dtype.itemsize:
                return a_dtype
            return np.int64
        return np.float64
    elif is_float_dtype(a_dtype):
        if is_float_dtype(b_dtype):
            if a_dtype.itemsize == b_dtype.itemsize:
                return a_dtype
            else:
                return np.float64
        elif is_integer(b_dtype):
            return np.float64
    return np.object


def _consensus_name_attr(objs):
    name = objs[0].name
    for obj in objs[1:]:
        if obj.name != name:
            return None
    return name

# ----------------------------------------------------------------------
# Lots of little utilities


def _validate_date_like_dtype(dtype):
    try:
        typ = np.datetime_data(dtype)[0]
    except ValueError as e:
        raise TypeError('%s' % e)
    if typ != 'generic' and typ != 'ns':
        raise ValueError('%r is too specific of a frequency, try passing %r' %
                         (dtype.name, dtype.type.__name__))


def _invalidate_string_dtypes(dtype_set):
    """Change string like dtypes to object for
    ``DataFrame.select_dtypes()``.
    """
    non_string_dtypes = dtype_set - _string_dtypes
    if non_string_dtypes != dtype_set:
        raise TypeError("string dtypes are not allowed, use 'object' instead")


def _get_dtype_from_object(dtype):
    """Get a numpy dtype.type-style object. This handles the datetime64[ns]
    and datetime64[ns, TZ] compat

    Notes
    -----
    If nothing can be found, returns ``object``.
    """
    # type object from a dtype
    if isinstance(dtype, type) and issubclass(dtype, np.generic):
        return dtype
    elif is_categorical(dtype):
        return gt.CategoricalDtype().type
    elif is_datetimetz(dtype):
        return gt.DatetimeTZDtype(dtype).type
    elif isinstance(dtype, np.dtype):  # dtype object
        try:
            _validate_date_like_dtype(dtype)
        except TypeError:
            # should still pass if we don't have a datelike
            pass
        return dtype.type
    elif isinstance(dtype, compat.string_types):
        if dtype == 'datetime' or dtype == 'timedelta':
            dtype += '64'

        try:
            return _get_dtype_from_object(getattr(np, dtype))
        except (AttributeError, TypeError):
            # handles cases like _get_dtype(int)
            # i.e., python objects that are valid dtypes (unlike user-defined
            # types, in general)
            # TypeError handles the float16 typecode of 'e'
            # further handle internal types
            pass

    return _get_dtype_from_object(np.dtype(dtype))


def _get_info_slice(obj, indexer):
    """Slice the info axis of `obj` with `indexer`."""
    if not hasattr(obj, '_info_axis_number'):
        raise TypeError('object of type %r has no info axis' %
                        type(obj).__name__)
    slices = [slice(None)] * obj.ndim
    slices[obj._info_axis_number] = indexer
    return tuple(slices)


def _maybe_box(indexer, values, obj, key):

    # if we have multiples coming back, box em
    if isinstance(values, np.ndarray):
        return obj[indexer.get_loc(key)]

    # return the value
    return values


def _maybe_box_datetimelike(value):
    # turn a datetime like into a Timestamp/timedelta as needed

    if isinstance(value, (np.datetime64, datetime)):
        value = tslib.Timestamp(value)
    elif isinstance(value, (np.timedelta64, timedelta)):
        value = tslib.Timedelta(value)

    return value


_values_from_object = lib.values_from_object


def _possibly_castable(arr):
    # return False to force a non-fastpath

    # check datetime64[ns]/timedelta64[ns] are valid
    # otherwise try to coerce
    kind = arr.dtype.kind
    if kind == 'M' or kind == 'm':
        return arr.dtype in _DATELIKE_DTYPES

    return arr.dtype.name not in _POSSIBLY_CAST_DTYPES


def _possibly_convert_platform(values):
    """ try to do platform conversion, allow ndarray or list here """

    if isinstance(values, (list, tuple)):
        values = lib.list_to_object_array(values)
    if getattr(values, 'dtype', None) == np.object_:
        if hasattr(values, '_values'):
            values = values._values
        values = lib.maybe_convert_objects(values)

    return values


def _possibly_cast_to_datetime(value, dtype, errors='raise'):
    """ try to cast the array/value to a datetimelike dtype, converting float
    nan to iNaT
    """
    from pandas.tseries.timedeltas import to_timedelta
    from pandas.tseries.tools import to_datetime

    if dtype is not None:
        if isinstance(dtype, compat.string_types):
            dtype = np.dtype(dtype)

        is_datetime64 = is_datetime64_dtype(dtype)
        is_datetime64tz = is_datetime64tz_dtype(dtype)
        is_timedelta64 = is_timedelta64_dtype(dtype)

        if is_datetime64 or is_datetime64tz or is_timedelta64:

            # force the dtype if needed
            if is_datetime64 and not is_dtype_equal(dtype, _NS_DTYPE):
                if dtype.name == 'datetime64[ns]':
                    dtype = _NS_DTYPE
                else:
                    raise TypeError("cannot convert datetimelike to "
                                    "dtype [%s]" % dtype)
            elif is_datetime64tz:

                # our NaT doesn't support tz's
                # this will coerce to DatetimeIndex with
                # a matching dtype below
                if lib.isscalar(value) and isnull(value):
                    value = [value]

            elif is_timedelta64 and not is_dtype_equal(dtype, _TD_DTYPE):
                if dtype.name == 'timedelta64[ns]':
                    dtype = _TD_DTYPE
                else:
                    raise TypeError("cannot convert timedeltalike to "
                                    "dtype [%s]" % dtype)

            if lib.isscalar(value):
                if value == tslib.iNaT or isnull(value):
                    value = tslib.iNaT
            else:
                value = np.array(value, copy=False)

                # have a scalar array-like (e.g. NaT)
                if value.ndim == 0:
                    value = tslib.iNaT

                # we have an array of datetime or timedeltas & nulls
                elif np.prod(value.shape) or not is_dtype_equal(value.dtype,
                                                                dtype):
                    try:
                        if is_datetime64:
                            value = to_datetime(value, errors=errors)._values
                        elif is_datetime64tz:
                            # input has to be UTC at this point, so just
                            # localize
                            value = to_datetime(
                                value,
                                errors=errors).tz_localize(dtype.tz)
                        elif is_timedelta64:
                            value = to_timedelta(value, errors=errors)._values
                    except (AttributeError, ValueError, TypeError):
                        pass

        # coerce datetimelike to object
        elif is_datetime64_dtype(value) and not is_datetime64_dtype(dtype):
            if is_object_dtype(dtype):
                ints = np.asarray(value).view('i8')
                return tslib.ints_to_pydatetime(ints)

            # we have a non-castable dtype that was passed
            raise TypeError('Cannot cast datetime64 to %s' % dtype)

    else:

        is_array = isinstance(value, np.ndarray)

        # catch a datetime/timedelta that is not of ns variety
        # and no coercion specified
        if is_array and value.dtype.kind in ['M', 'm']:
            dtype = value.dtype

            if dtype.kind == 'M' and dtype != _NS_DTYPE:
                value = value.astype(_NS_DTYPE)

            elif dtype.kind == 'm' and dtype != _TD_DTYPE:
                value = to_timedelta(value)

        # only do this if we have an array and the dtype of the array is not
        # setup already we are not an integer/object, so don't bother with this
        # conversion
        elif not (is_array and not (issubclass(value.dtype.type, np.integer) or
                                    value.dtype == np.object_)):
            value = _possibly_infer_to_datetimelike(value)

    return value


def _possibly_infer_to_datetimelike(value, convert_dates=False):
    """
    we might have a array (or single object) that is datetime like,
    and no dtype is passed don't change the value unless we find a
    datetime/timedelta set

    this is pretty strict in that a datetime/timedelta is REQUIRED
    in addition to possible nulls/string likes

    ONLY strings are NOT datetimelike

    Parameters
    ----------
    value : np.array / Series / Index / list-like
    convert_dates : boolean, default False
       if True try really hard to convert dates (such as datetime.date), other
       leave inferred dtype 'date' alone

    """

    if isinstance(value, (gt.ABCDatetimeIndex, gt.ABCPeriodIndex)):
        return value
    elif isinstance(value, gt.ABCSeries):
        if isinstance(value._values, gt.ABCDatetimeIndex):
            return value._values

    v = value
    if not is_list_like(v):
        v = [v]
    v = np.array(v, copy=False)
    shape = v.shape
    if not v.ndim == 1:
        v = v.ravel()

    if len(v):

        def _try_datetime(v):
            # safe coerce to datetime64
            try:
                v = tslib.array_to_datetime(v, errors='raise')
            except ValueError:

                # we might have a sequence of the same-datetimes with tz's
                # if so coerce to a DatetimeIndex; if they are not the same,
                # then these stay as object dtype
                try:
                    from pandas import to_datetime
                    return to_datetime(v)
                except:
                    pass

            except:
                pass

            return v.reshape(shape)

        def _try_timedelta(v):
            # safe coerce to timedelta64

            # will try first with a string & object conversion
            from pandas.tseries.timedeltas import to_timedelta
            try:
                return to_timedelta(v)._values.reshape(shape)
            except:
                return v

        # do a quick inference for perf
        sample = v[:min(3, len(v))]
        inferred_type = lib.infer_dtype(sample)

        if (inferred_type in ['datetime', 'datetime64'] or
                (convert_dates and inferred_type in ['date'])):
            value = _try_datetime(v)
        elif inferred_type in ['timedelta', 'timedelta64']:
            value = _try_timedelta(v)

        # It's possible to have nulls intermixed within the datetime or
        # timedelta.  These will in general have an inferred_type of 'mixed',
        # so have to try both datetime and timedelta.

        # try timedelta first to avoid spurious datetime conversions
        # e.g. '00:00:01' is a timedelta but technically is also a datetime
        elif inferred_type in ['mixed']:

            if lib.is_possible_datetimelike_array(_ensure_object(v)):
                value = _try_timedelta(v)
                if lib.infer_dtype(value) in ['mixed']:
                    value = _try_datetime(v)

    return value


def is_bool_indexer(key):
    if isinstance(key, (gt.ABCSeries, np.ndarray)):
        if key.dtype == np.object_:
            key = np.asarray(_values_from_object(key))

            if not lib.is_bool_array(key):
                if isnull(key).any():
                    raise ValueError('cannot index with vector containing '
                                     'NA / NaN values')
                return False
            return True
        elif key.dtype == np.bool_:
            return True
    elif isinstance(key, list):
        try:
            arr = np.asarray(key)
            return arr.dtype == np.bool_ and len(arr) == len(key)
        except TypeError:  # pragma: no cover
            return False

    return False


def _default_index(n):
    from pandas.core.index import RangeIndex
    return RangeIndex(0, n, name=None)


def ensure_float(arr):
    if issubclass(arr.dtype.type, (np.integer, np.bool_)):
        arr = arr.astype(float)
    return arr


def _mut_exclusive(**kwargs):
    item1, item2 = kwargs.items()
    label1, val1 = item1
    label2, val2 = item2
    if val1 is not None and val2 is not None:
        raise TypeError('mutually exclusive arguments: %r and %r' %
                        (label1, label2))
    elif val1 is not None:
        return val1
    else:
        return val2


def _not_none(*args):
    return (arg for arg in args if arg is not None)


def _any_none(*args):
    for arg in args:
        if arg is None:
            return True
    return False


def _all_not_none(*args):
    for arg in args:
        if arg is None:
            return False
    return True


def _try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed


def _count_not_none(*args):
    return sum(x is not None for x in args)


def iterpairs(seq):
    """
    Parameters
    ----------
    seq : sequence

    Returns
    -------
    iterator returning overlapping pairs of elements

    Examples
    --------
    >>> list(iterpairs([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """
    # input may not be sliceable
    seq_it = iter(seq)
    seq_it_next = iter(seq)
    next(seq_it_next)

    return zip(seq_it, seq_it_next)


def split_ranges(mask):
    """ Generates tuples of ranges which cover all True value in mask

    >>> list(split_ranges([1,0,0,1,0]))
    [(0, 1), (3, 4)]
    """
    ranges = [(0, len(mask))]

    for pos, val in enumerate(mask):
        if not val:  # this pos should be ommited, split off the prefix range
            r = ranges.pop()
            if pos > r[0]:  # yield non-zero range
                yield (r[0], pos)
            if pos + 1 < len(mask):  # save the rest for processing
                ranges.append((pos + 1, len(mask)))
    if ranges:
        yield ranges[-1]


def _long_prod(vals):
    result = long(1)
    for x in vals:
        result *= x
    return result


class groupby(dict):
    """
    A simple groupby different from the one in itertools.

    Does not require the sequence elements to be sorted by keys,
    however it is slower.
    """

    def __init__(self, seq, key=lambda x: x):
        for value in seq:
            k = key(value)
            self.setdefault(k, []).append(value)

    try:
        __iter__ = dict.iteritems
    except AttributeError:  # pragma: no cover
        # Python 3
        def __iter__(self):
            return iter(dict.items(self))


def map_indices_py(arr):
    """
    Returns a dictionary with (element, index) pairs for each element in the
    given array/list
    """
    return dict([(x, i) for i, x in enumerate(arr)])


def union(*seqs):
    result = set([])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result |= seq
    return type(seqs[0])(list(result))


def difference(a, b):
    return type(a)(list(set(a) - set(b)))


def intersection(*seqs):
    result = set(seqs[0])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result &= seq
    return type(seqs[0])(list(result))


def _asarray_tuplesafe(values, dtype=None):
    from pandas.core.index import Index

    if not (isinstance(values, (list, tuple)) or hasattr(values, '__array__')):
        values = list(values)
    elif isinstance(values, Index):
        return values.values

    if isinstance(values, list) and dtype in [np.object_, object]:
        return lib.list_to_object_array(values)

    result = np.asarray(values, dtype=dtype)

    if issubclass(result.dtype.type, compat.string_types):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        if isinstance(values, list):
            return lib.list_to_object_array(values)
        else:
            # Making a 1D array that safely contains tuples is a bit tricky
            # in numpy, leading to the following
            try:
                result = np.empty(len(values), dtype=object)
                result[:] = values
            except ValueError:
                # we have a list-of-list
                result[:] = [tuple(x) for x in values]

    return result


def _index_labels_to_array(labels):
    if isinstance(labels, (compat.string_types, tuple)):
        labels = [labels]

    if not isinstance(labels, (list, np.ndarray)):
        try:
            labels = list(labels)
        except TypeError:  # non-iterable
            labels = [labels]

    labels = _asarray_tuplesafe(labels)

    return labels


def _maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj

# TYPE TESTING

is_bool = lib.is_bool

is_integer = lib.is_integer

is_float = lib.is_float

is_complex = lib.is_complex


def is_string_like(obj):
    return isinstance(obj, (compat.text_type, compat.string_types))


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_number(obj):
    return isinstance(obj, (numbers.Number, np.number))


def is_period_arraylike(arr):
    """ return if we are period arraylike / PeriodIndex """
    if isinstance(arr, pd.PeriodIndex):
        return True
    elif isinstance(arr, (np.ndarray, gt.ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'period'
    return getattr(arr, 'inferred_type', None) == 'period'


def is_datetime_arraylike(arr):
    """ return if we are datetime arraylike / DatetimeIndex """
    if isinstance(arr, gt.ABCDatetimeIndex):
        return True
    elif isinstance(arr, (np.ndarray, gt.ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'datetime'
    return getattr(arr, 'inferred_type', None) == 'datetime'


def is_datetimelike(arr):
    return (arr.dtype in _DATELIKE_DTYPES or
            isinstance(arr, gt.ABCPeriodIndex) or
            is_datetimetz(arr))


def _coerce_to_dtype(dtype):
    """ coerce a string / np.dtype to a dtype """
    if is_categorical_dtype(dtype):
        dtype = gt.CategoricalDtype()
    elif is_datetime64tz_dtype(dtype):
        dtype = gt.DatetimeTZDtype(dtype)
    else:
        dtype = np.dtype(dtype)
    return dtype


def _get_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype)
    elif isinstance(arr_or_dtype, gt.CategoricalDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, gt.DatetimeTZDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, compat.string_types):
        if is_categorical_dtype(arr_or_dtype):
            return gt.CategoricalDtype.construct_from_string(arr_or_dtype)
        elif is_datetime64tz_dtype(arr_or_dtype):
            return gt.DatetimeTZDtype.construct_from_string(arr_or_dtype)

    if hasattr(arr_or_dtype, 'dtype'):
        arr_or_dtype = arr_or_dtype.dtype
    return np.dtype(arr_or_dtype)


def _get_dtype_type(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype).type
    elif isinstance(arr_or_dtype, gt.CategoricalDtype):
        return gt.CategoricalDtypeType
    elif isinstance(arr_or_dtype, gt.DatetimeTZDtype):
        return gt.DatetimeTZDtypeType
    elif isinstance(arr_or_dtype, compat.string_types):
        if is_categorical_dtype(arr_or_dtype):
            return gt.CategoricalDtypeType
        elif is_datetime64tz_dtype(arr_or_dtype):
            return gt.DatetimeTZDtypeType
        return _get_dtype_type(np.dtype(arr_or_dtype))
    try:
        return arr_or_dtype.dtype.type
    except AttributeError:
        return type(None)


def is_dtype_equal(source, target):
    """ return a boolean if the dtypes are equal """
    try:
        source = _get_dtype(source)
        target = _get_dtype(target)
        return source == target
    except (TypeError, AttributeError):

        # invalid comparison
        # object == category will hit this
        return False


def is_any_int_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.integer)


def is_integer_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_int64_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.int64)


def is_int_or_datetime_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) or
            issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_datetime64_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except TypeError:
        return False
    return issubclass(tipo, np.datetime64)


def is_datetime64tz_dtype(arr_or_dtype):
    return gt.DatetimeTZDtype.is_dtype(arr_or_dtype)


def is_datetime64_any_dtype(arr_or_dtype):
    return (is_datetime64_dtype(arr_or_dtype) or
            is_datetime64tz_dtype(arr_or_dtype))


def is_datetime64_ns_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype(arr_or_dtype)
    except TypeError:
        return False
    return tipo == _NS_DTYPE


def is_timedelta64_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.timedelta64)


def is_timedelta64_ns_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return tipo == _TD_DTYPE


def is_datetime_or_timedelta_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, (np.datetime64, np.timedelta64))


def is_numeric_v_string_like(a, b):
    """
    numpy doesn't like to compare numeric arrays vs scalar string-likes

    return a boolean result if this is the case for a,b or b,a

    """
    is_a_array = isinstance(a, np.ndarray)
    is_b_array = isinstance(b, np.ndarray)

    is_a_numeric_array = is_a_array and is_numeric_dtype(a)
    is_b_numeric_array = is_b_array and is_numeric_dtype(b)
    is_a_string_array = is_a_array and is_string_like_dtype(a)
    is_b_string_array = is_b_array and is_string_like_dtype(b)

    is_a_scalar_string_like = not is_a_array and is_string_like(a)
    is_b_scalar_string_like = not is_b_array and is_string_like(b)

    return ((is_a_numeric_array and is_b_scalar_string_like) or
            (is_b_numeric_array and is_a_scalar_string_like) or
            (is_a_numeric_array and is_b_string_array) or
            (is_b_numeric_array and is_a_string_array))


def is_datetimelike_v_numeric(a, b):
    # return if we have an i8 convertible and numeric comparison
    if not hasattr(a, 'dtype'):
        a = np.asarray(a)
    if not hasattr(b, 'dtype'):
        b = np.asarray(b)

    def is_numeric(x):
        return is_integer_dtype(x) or is_float_dtype(x)

    is_datetimelike = needs_i8_conversion
    return ((is_datetimelike(a) and is_numeric(b)) or
            (is_datetimelike(b) and is_numeric(a)))


def is_datetimelike_v_object(a, b):
    # return if we have an i8 convertible and object comparsion
    if not hasattr(a, 'dtype'):
        a = np.asarray(a)
    if not hasattr(b, 'dtype'):
        b = np.asarray(b)

    def f(x):
        return is_object_dtype(x)

    def is_object(x):
        return is_integer_dtype(x) or is_float_dtype(x)

    is_datetimelike = needs_i8_conversion
    return ((is_datetimelike(a) and is_object(b)) or
            (is_datetimelike(b) and is_object(a)))


def needs_i8_conversion(arr_or_dtype):
    return (is_datetime_or_timedelta_dtype(arr_or_dtype) or
            is_datetime64tz_dtype(arr_or_dtype))


def is_numeric_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, (np.number, np.bool_)) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_string_dtype(arr_or_dtype):
    dtype = _get_dtype(arr_or_dtype)
    return dtype.kind in ('O', 'S', 'U')


def is_string_like_dtype(arr_or_dtype):
    # exclude object as its a mixed dtype
    dtype = _get_dtype(arr_or_dtype)
    return dtype.kind in ('S', 'U')


def is_float_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.floating)


def is_floating_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return isinstance(tipo, np.floating)


def is_bool_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except ValueError:
        # this isn't even a dtype
        return False
    return issubclass(tipo, np.bool_)


def is_sparse(array):
    """ return if we are a sparse array """
    return isinstance(array, (gt.ABCSparseArray, gt.ABCSparseSeries))


def is_datetimetz(array):
    """ return if we are a datetime with tz array """
    return ((isinstance(array, gt.ABCDatetimeIndex) and
             getattr(array, 'tz', None) is not None) or
            is_datetime64tz_dtype(array))


def is_extension_type(value):
    """
    if we are a klass that is preserved by the internals
    these are internal klasses that we represent (and don't use a np.array)
    """
    if is_categorical(value):
        return True
    elif is_sparse(value):
        return True
    elif is_datetimetz(value):
        return True
    return False


def is_categorical(array):
    """ return if we are a categorical possibility """
    return isinstance(array, gt.ABCCategorical) or is_categorical_dtype(array)


def is_categorical_dtype(arr_or_dtype):
    return gt.CategoricalDtype.is_dtype(arr_or_dtype)


def is_complex_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.complexfloating)


def is_object_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.object_)


def is_re(obj):
    return isinstance(obj, re._pattern_type)


def is_re_compilable(obj):
    try:
        re.compile(obj)
    except TypeError:
        return False
    else:
        return True


def is_list_like(arg):
    return (hasattr(arg, '__iter__') and
            not isinstance(arg, compat.string_and_binary_types))


def is_dict_like(arg):
    return hasattr(arg, '__getitem__') and hasattr(arg, 'keys')


def is_named_tuple(arg):
    return isinstance(arg, tuple) and hasattr(arg, '_fields')


def is_null_slice(obj):
    """ we have a null slice """
    return (isinstance(obj, slice) and obj.start is None and
            obj.stop is None and obj.step is None)


def is_full_slice(obj, l):
    """ we have a full length slice """
    return (isinstance(obj, slice) and obj.start == 0 and obj.stop == l and
            obj.step is None)


def is_hashable(arg):
    """Return True if hash(arg) will succeed, False otherwise.

    Some types will pass a test against collections.Hashable but fail when they
    are actually hashed with hash().

    Distinguish between these and other types by trying the call to hash() and
    seeing if they raise TypeError.

    Examples
    --------
    >>> a = ([],)
    >>> isinstance(a, collections.Hashable)
    True
    >>> is_hashable(a)
    False
    """
    # unfortunately, we can't use isinstance(arg, collections.Hashable), which
    # can be faster than calling hash, because numpy scalars on Python 3 fail
    # this test

    # reconsider this decision once this numpy bug is fixed:
    # https://github.com/numpy/numpy/issues/5562

    try:
        hash(arg)
    except TypeError:
        return False
    else:
        return True


def is_sequence(x):
    try:
        iter(x)
        len(x)  # it has a length
        return not isinstance(x, compat.string_and_binary_types)
    except (TypeError, AttributeError):
        return False


def _get_callable_name(obj):
    # typical case has name
    if hasattr(obj, '__name__'):
        return getattr(obj, '__name__')
    # some objects don't; could recurse
    if isinstance(obj, partial):
        return _get_callable_name(obj.func)
    # fall back to class name
    if hasattr(obj, '__call__'):
        return obj.__class__.__name__
    # everything failed (probably because the argument
    # wasn't actually callable); we return None
    # instead of the empty string in this case to allow
    # distinguishing between no name and a name of ''
    return None


def _apply_if_callable(maybe_callable, obj, **kwargs):
    """
    Evaluate possibly callable input using obj and kwargs if it is callable,
    otherwise return as it is
    """
    if callable(maybe_callable):
        return maybe_callable(obj, **kwargs)
    return maybe_callable


_string_dtypes = frozenset(map(_get_dtype_from_object, (compat.binary_type,
                                                        compat.text_type)))

_ensure_float64 = algos.ensure_float64
_ensure_float32 = algos.ensure_float32
_ensure_int64 = algos.ensure_int64
_ensure_int32 = algos.ensure_int32
_ensure_int16 = algos.ensure_int16
_ensure_int8 = algos.ensure_int8
_ensure_platform_int = algos.ensure_platform_int
_ensure_object = algos.ensure_object


def _astype_nansafe(arr, dtype, copy=True):
    """ return a view if copy is False, but
        need to be very careful as the result shape could change! """
    if not isinstance(dtype, np.dtype):
        dtype = _coerce_to_dtype(dtype)

    if issubclass(dtype.type, compat.text_type):
        # in Py3 that's str, in Py2 that's unicode
        return lib.astype_unicode(arr.ravel()).reshape(arr.shape)
    elif issubclass(dtype.type, compat.string_types):
        return lib.astype_str(arr.ravel()).reshape(arr.shape)
    elif is_datetime64_dtype(arr):
        if dtype == object:
            return tslib.ints_to_pydatetime(arr.view(np.int64))
        elif dtype == np.int64:
            return arr.view(dtype)
        elif dtype != _NS_DTYPE:
            raise TypeError("cannot astype a datetimelike from [%s] to [%s]" %
                            (arr.dtype, dtype))
        return arr.astype(_NS_DTYPE)
    elif is_timedelta64_dtype(arr):
        if dtype == np.int64:
            return arr.view(dtype)
        elif dtype == object:
            return tslib.ints_to_pytimedelta(arr.view(np.int64))

        # in py3, timedelta64[ns] are int64
        elif ((compat.PY3 and dtype not in [_INT64_DTYPE, _TD_DTYPE]) or
              (not compat.PY3 and dtype != _TD_DTYPE)):

            # allow frequency conversions
            if dtype.kind == 'm':
                mask = isnull(arr)
                result = arr.astype(dtype).astype(np.float64)
                result[mask] = np.nan
                return result

            raise TypeError("cannot astype a timedelta from [%s] to [%s]" %
                            (arr.dtype, dtype))

        return arr.astype(_TD_DTYPE)
    elif (np.issubdtype(arr.dtype, np.floating) and
          np.issubdtype(dtype, np.integer)):

        if np.isnan(arr).any():
            raise ValueError('Cannot convert NA to integer')
    elif arr.dtype == np.object_ and np.issubdtype(dtype.type, np.integer):
        # work around NumPy brokenness, #1987
        return lib.astype_intsafe(arr.ravel(), dtype).reshape(arr.shape)

    if copy:
        return arr.astype(dtype)
    return arr.view(dtype)


def _all_none(*args):
    for arg in args:
        if arg is not None:
            return False
    return True


def _where_compat(mask, arr1, arr2):
    if arr1.dtype == _NS_DTYPE and arr2.dtype == _NS_DTYPE:
        new_vals = np.where(mask, arr1.view('i8'), arr2.view('i8'))
        return new_vals.view(_NS_DTYPE)

    import pandas.tslib as tslib
    if arr1.dtype == _NS_DTYPE:
        arr1 = tslib.ints_to_pydatetime(arr1.view('i8'))
    if arr2.dtype == _NS_DTYPE:
        arr2 = tslib.ints_to_pydatetime(arr2.view('i8'))

    return np.where(mask, arr1, arr2)


def _dict_compat(d):
    """
    Helper function to convert datetimelike-keyed dicts to Timestamp-keyed dict

    Parameters
    ----------
    d: dict like object

    Returns
    -------
    dict

    """
    return dict((_maybe_box_datetimelike(key), value)
                for key, value in iteritems(d))


def sentinel_factory():
    class Sentinel(object):
        pass

    return Sentinel()


def in_interactive_session():
    """ check if we're running in an interactive shell

    returns True if running under python/ipython interactive shell
    """

    def check_main():
        import __main__ as main
        return (not hasattr(main, '__file__') or
                get_option('mode.sim_interactive'))

    try:
        return __IPYTHON__ or check_main()  # noqa
    except:
        return check_main()


def in_qtconsole():
    """
    check if we're inside an IPython qtconsole

    DEPRECATED: This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'qtconsole' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipnb():
    """
    check if we're inside an IPython Notebook

    DEPRECATED: This is no longer used in pandas, and won't work in IPython 3
    and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'notebook' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipython_frontend():
    """
    check if we're inside an an IPython zmq frontend
    """
    try:
        ip = get_ipython()  # noqa
        return 'zmq' in str(type(ip)).lower()
    except:
        pass

    return False


def _maybe_match_name(a, b):
    a_has = hasattr(a, 'name')
    b_has = hasattr(b, 'name')
    if a_has and b_has:
        if a.name == b.name:
            return a.name
        else:
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None


def _random_state(state=None):
    """
    Helper function for processing random_state arguments.

    Parameters
    ----------
    state : int, np.random.RandomState, None.
        If receives an int, passes to np.random.RandomState() as seed.
        If receives an np.random.RandomState object, just returns object.
        If receives `None`, returns an np.random.RandomState object.
        If receives anything else, raises an informative ValueError.
        Default None.

    Returns
    -------
    np.random.RandomState
    """

    if is_integer(state):
        return np.random.RandomState(state)
    elif isinstance(state, np.random.RandomState):
        return state
    elif state is None:
        return np.random.RandomState()
    else:
        raise ValueError("random_state must be an integer, a numpy "
                         "RandomState, or None")
