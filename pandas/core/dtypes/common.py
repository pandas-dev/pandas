""" common type operations """

import numpy as np
from pandas.compat import (string_types, text_type, binary_type,
                           PY3, PY36)
from pandas._libs import algos, lib
from .dtypes import (CategoricalDtype, CategoricalDtypeType,
                     DatetimeTZDtype, DatetimeTZDtypeType,
                     PeriodDtype, PeriodDtypeType,
                     IntervalDtype, IntervalDtypeType,
                     ExtensionDtype)
from .generic import (ABCCategorical, ABCPeriodIndex,
                      ABCDatetimeIndex, ABCSeries,
                      ABCSparseArray, ABCSparseSeries)
from .inference import is_string_like
from .inference import *  # noqa


_POSSIBLY_CAST_DTYPES = set([np.dtype(t).name
                             for t in ['O', 'int8', 'uint8', 'int16', 'uint16',
                                       'int32', 'uint32', 'int64', 'uint64']])

_NS_DTYPE = np.dtype('M8[ns]')
_TD_DTYPE = np.dtype('m8[ns]')
_INT64_DTYPE = np.dtype(np.int64)

# oh the troubles to reduce import time
_is_scipy_sparse = None

_ensure_float64 = algos.ensure_float64
_ensure_float32 = algos.ensure_float32


def _ensure_float(arr):
    """
    Ensure that an array object has a float dtype if possible.

    Parameters
    ----------
    arr : ndarray, Series
        The array whose data type we want to enforce as float.

    Returns
    -------
    float_arr : The original array cast to the float dtype if
                possible. Otherwise, the original array is returned.
    """

    if issubclass(arr.dtype.type, (np.integer, np.bool_)):
        arr = arr.astype(float)
    return arr


_ensure_uint64 = algos.ensure_uint64
_ensure_int64 = algos.ensure_int64
_ensure_int32 = algos.ensure_int32
_ensure_int16 = algos.ensure_int16
_ensure_int8 = algos.ensure_int8
_ensure_platform_int = algos.ensure_platform_int
_ensure_object = algos.ensure_object


def _ensure_categorical(arr):
    """
    Ensure that an array-like object is a Categorical (if not already).

    Parameters
    ----------
    arr : array-like
        The array that we want to convert into a Categorical.

    Returns
    -------
    cat_arr : The original array cast as a Categorical. If it already
              is a Categorical, we return as is.
    """

    if not is_categorical(arr):
        from pandas import Categorical
        arr = Categorical(arr)
    return arr


def is_object_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.object_)


def is_sparse(array):
    """ return if we are a sparse array """
    return isinstance(array, (ABCSparseArray, ABCSparseSeries))


def is_scipy_sparse(array):
    """ return if we are a scipy.sparse.spmatrix """
    global _is_scipy_sparse
    if _is_scipy_sparse is None:
        try:
            from scipy.sparse import issparse as _is_scipy_sparse
        except ImportError:
            _is_scipy_sparse = lambda _: False
    return _is_scipy_sparse(array)


def is_categorical(array):
    """ return if we are a categorical possibility """
    return isinstance(array, ABCCategorical) or is_categorical_dtype(array)


def is_datetimetz(array):
    """ return if we are a datetime with tz array """
    return ((isinstance(array, ABCDatetimeIndex) and
             getattr(array, 'tz', None) is not None) or
            is_datetime64tz_dtype(array))


def is_period(array):
    """ return if we are a period array """
    return isinstance(array, ABCPeriodIndex) or is_period_arraylike(array)


def is_datetime64_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except TypeError:
        return False
    return issubclass(tipo, np.datetime64)


def is_datetime64tz_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    return DatetimeTZDtype.is_dtype(arr_or_dtype)


def is_timedelta64_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.timedelta64)


def is_period_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    return PeriodDtype.is_dtype(arr_or_dtype)


def is_interval_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    return IntervalDtype.is_dtype(arr_or_dtype)


def is_categorical_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    return CategoricalDtype.is_dtype(arr_or_dtype)


def is_string_dtype(arr_or_dtype):
    """
    Check whether the provided array or dtype is of the string dtype.

    Parameters
    ----------
    arr_or_dtype : ndarray, dtype, type
        The array or dtype to check.

    Returns
    -------
    boolean : Whether or not the array or dtype is of the string dtype.

    Examples
    --------
    >>> is_string_dtype(str)
    True
    >>> is_string_dtype(object)
    True
    >>> is_string_dtype(int)
    False
    >>>
    >>> is_string_dtype(np.array(['a', 'b']))
    True
    >>> is_string_dtype(np.array([1, 2]))
    False
    """

    # TODO: gh-15585: consider making the checks stricter.

    if arr_or_dtype is None:
        return False
    try:
        dtype = _get_dtype(arr_or_dtype)
        return dtype.kind in ('O', 'S', 'U') and not is_period_dtype(dtype)
    except TypeError:
        return False


def is_period_arraylike(arr):
    """ return if we are period arraylike / PeriodIndex """
    if isinstance(arr, ABCPeriodIndex):
        return True
    elif isinstance(arr, (np.ndarray, ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'period'
    return getattr(arr, 'inferred_type', None) == 'period'


def is_datetime_arraylike(arr):
    """ return if we are datetime arraylike / DatetimeIndex """
    if isinstance(arr, ABCDatetimeIndex):
        return True
    elif isinstance(arr, (np.ndarray, ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'datetime'
    return getattr(arr, 'inferred_type', None) == 'datetime'


def is_datetimelike(arr):
    return (is_datetime64_dtype(arr) or is_datetime64tz_dtype(arr) or
            is_timedelta64_dtype(arr) or
            isinstance(arr, ABCPeriodIndex) or
            is_datetimetz(arr))


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
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.integer)


def is_integer_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_signed_integer_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.signedinteger) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_unsigned_integer_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.unsignedinteger) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_int64_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.int64)


def is_int_or_datetime_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) or
            issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_datetime64_any_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    return (is_datetime64_dtype(arr_or_dtype) or
            is_datetime64tz_dtype(arr_or_dtype))


def is_datetime64_ns_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    try:
        tipo = _get_dtype(arr_or_dtype)
    except TypeError:
        if is_datetime64tz_dtype(arr_or_dtype):
            tipo = _get_dtype(arr_or_dtype.dtype)
        else:
            return False
    return tipo == _NS_DTYPE or getattr(tipo, 'base', None) == _NS_DTYPE


def is_timedelta64_ns_dtype(arr_or_dtype):
    """
    Check whether the provided array or dtype is of the timedelta64[ns] dtype.

    This is a very specific dtype, so generic ones like `np.timedelta64`
    will return False if passed into this function.

    Parameters
    ----------
    arr_or_dtype : ndarray, dtype, type
        The array or dtype to check.

    Returns
    -------
    boolean : Whether or not the array or dtype
              is of the timedelta64[ns] dtype.

    Examples
    --------
    >>> is_timedelta64_ns_dtype(np.dtype('m8[ns]')
    True
    >>> is_timedelta64_ns_dtype(np.dtype('m8[ps]')  # Wrong frequency
    False
    >>>
    >>> is_timedelta64_ns_dtype(np.array([1, 2], dtype='m8[ns]'))
    True
    >>> is_timedelta64_ns_dtype(np.array([1, 2], dtype=np.timedelta64))
    False
    """

    if arr_or_dtype is None:
        return False
    try:
        tipo = _get_dtype(arr_or_dtype)
        return tipo == _TD_DTYPE
    except TypeError:
        return False


def is_datetime_or_timedelta_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, (np.datetime64, np.timedelta64))


def _is_unorderable_exception(e):
    """
    Check if the exception raised is an unorderable exception.

    The error message differs for 3 <= PY <= 3.5 and PY >= 3.6, so
    we need to condition based on Python version.

    Parameters
    ----------
    e : Exception or sub-class
        The exception object to check.

    Returns
    -------
    boolean : Whether or not the exception raised is an unorderable exception.
    """

    if PY36:
        return "'>' not supported between instances of" in str(e)

    elif PY3:
        return 'unorderable' in str(e)
    return False


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
    if arr_or_dtype is None:
        return False
    return (is_datetime_or_timedelta_dtype(arr_or_dtype) or
            is_datetime64tz_dtype(arr_or_dtype) or
            is_period_dtype(arr_or_dtype))


def is_numeric_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, (np.number, np.bool_)) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_string_like_dtype(arr_or_dtype):
    """
    Check whether the provided array or dtype is of a string-like dtype.

    Unlike `is_string_dtype`, the object dtype is excluded because it
    is a mixed dtype.

    Parameters
    ----------
    arr_or_dtype : ndarray, dtype, type
        The array or dtype to check.

    Returns
    -------
    boolean : Whether or not the array or dtype is of the string dtype.

    Examples
    --------
    >>> is_string_like_dtype(str)
    True
    >>> is_string_like_dtype(object)
    False
    >>>
    >>> is_string_like_dtype(np.array(['a', 'b']))
    True
    >>> is_string_like_dtype(np.array([1, 2]))
    False
    """

    if arr_or_dtype is None:
        return False
    try:
        dtype = _get_dtype(arr_or_dtype)
        return dtype.kind in ('S', 'U')
    except TypeError:
        return False


def is_float_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.floating)


def is_floating_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return isinstance(tipo, np.floating)


def is_bool_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except ValueError:
        # this isn't even a dtype
        return False
    return issubclass(tipo, np.bool_)


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


def is_complex_dtype(arr_or_dtype):
    if arr_or_dtype is None:
        return False
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.complexfloating)


def _coerce_to_dtype(dtype):
    """
    Coerce a string or np.dtype to a pandas or numpy
    dtype if possible.

    If we cannot convert to a pandas dtype initially,
    we convert to a numpy dtype.

    Parameters
    ----------
    dtype : The dtype that we want to coerce.

    Returns
    -------
    pd_or_np_dtype : The coerced dtype.
    """

    if is_categorical_dtype(dtype):
        dtype = CategoricalDtype()
    elif is_datetime64tz_dtype(dtype):
        dtype = DatetimeTZDtype(dtype)
    elif is_period_dtype(dtype):
        dtype = PeriodDtype(dtype)
    elif is_interval_dtype(dtype):
        dtype = IntervalDtype(dtype)
    else:
        dtype = np.dtype(dtype)
    return dtype


def _get_dtype(arr_or_dtype):
    """
    Get the dtype instance associated with an array
    or dtype object.

    Parameters
    ----------
    arr_or_dtype : ndarray, Series, dtype, type
        The array-like or dtype object whose dtype we want to extract.

    Returns
    -------
    obj_dtype : The extract dtype instance from the
                passed in array or dtype object.

    Raises
    ------
    TypeError : The passed in object is None.
    """

    if arr_or_dtype is None:
        raise TypeError("Cannot deduce dtype from null object")
    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype)
    elif isinstance(arr_or_dtype, CategoricalDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, DatetimeTZDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, PeriodDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, IntervalDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, string_types):
        if is_categorical_dtype(arr_or_dtype):
            return CategoricalDtype.construct_from_string(arr_or_dtype)
        elif is_datetime64tz_dtype(arr_or_dtype):
            return DatetimeTZDtype.construct_from_string(arr_or_dtype)
        elif is_period_dtype(arr_or_dtype):
            return PeriodDtype.construct_from_string(arr_or_dtype)
        elif is_interval_dtype(arr_or_dtype):
            return IntervalDtype.construct_from_string(arr_or_dtype)

    if hasattr(arr_or_dtype, 'dtype'):
        arr_or_dtype = arr_or_dtype.dtype
    return np.dtype(arr_or_dtype)


def _get_dtype_type(arr_or_dtype):
    """
    Get the type (NOT dtype) instance associated with
    an array or dtype object.

    Parameters
    ----------
    arr_or_dtype : ndarray, Series, dtype, type
        The array-like or dtype object whose type we want to extract.

    Returns
    -------
    obj_type : The extract type instance from the
               passed in array or dtype object.
    """

    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype).type
    elif isinstance(arr_or_dtype, CategoricalDtype):
        return CategoricalDtypeType
    elif isinstance(arr_or_dtype, DatetimeTZDtype):
        return DatetimeTZDtypeType
    elif isinstance(arr_or_dtype, IntervalDtype):
        return IntervalDtypeType
    elif isinstance(arr_or_dtype, PeriodDtype):
        return PeriodDtypeType
    elif isinstance(arr_or_dtype, string_types):
        if is_categorical_dtype(arr_or_dtype):
            return CategoricalDtypeType
        elif is_datetime64tz_dtype(arr_or_dtype):
            return DatetimeTZDtypeType
        elif is_period_dtype(arr_or_dtype):
            return PeriodDtypeType
        elif is_interval_dtype(arr_or_dtype):
            return IntervalDtypeType
        return _get_dtype_type(np.dtype(arr_or_dtype))
    try:
        return arr_or_dtype.dtype.type
    except AttributeError:
        return type(None)


def _get_dtype_from_object(dtype):
    """
    Get a numpy dtype.type-style object for a dtype object.

    This methods also includes handling of the datetime64[ns] and
    datetime64[ns, TZ] objects.

    If no dtype can be found, we return ``object``.

    Parameters
    ----------
    dtype : dtype, type
        The dtype object whose numpy dtype.type-style
        object we want to extract.

    Returns
    -------
    dtype_object : The extracted numpy dtype.type-style object.
    """

    if isinstance(dtype, type) and issubclass(dtype, np.generic):
        # Type object from a dtype
        return dtype
    elif is_categorical(dtype):
        return CategoricalDtype().type
    elif is_datetimetz(dtype):
        return DatetimeTZDtype(dtype).type
    elif isinstance(dtype, np.dtype):  # dtype object
        try:
            _validate_date_like_dtype(dtype)
        except TypeError:
            # Should still pass if we don't have a date-like
            pass
        return dtype.type
    elif isinstance(dtype, string_types):
        if dtype in ['datetimetz', 'datetime64tz']:
            return DatetimeTZDtype.type
        elif dtype in ['period']:
            raise NotImplementedError

        if dtype == 'datetime' or dtype == 'timedelta':
            dtype += '64'

        try:
            return _get_dtype_from_object(getattr(np, dtype))
        except (AttributeError, TypeError):
            # Handles cases like _get_dtype(int) i.e.,
            # Python objects that are valid dtypes
            # (unlike user-defined types, in general)
            #
            # TypeError handles the float16 type code of 'e'
            # further handle internal types
            pass

    return _get_dtype_from_object(np.dtype(dtype))


def _validate_date_like_dtype(dtype):
    """
    Check whether the dtype is a date-like dtype. Raises an error if invalid.

    Parameters
    ----------
    dtype : dtype, type
        The dtype to check.

    Raises
    ------
    TypeError : The dtype could not be casted to a date-like dtype.
    ValueError : The dtype is an illegal date-like dtype (e.g. the
                 the frequency provided is too specific)
    """

    try:
        typ = np.datetime_data(dtype)[0]
    except ValueError as e:
        raise TypeError('%s' % e)
    if typ != 'generic' and typ != 'ns':
        raise ValueError('%r is too specific of a frequency, try passing %r' %
                         (dtype.name, dtype.type.__name__))


_string_dtypes = frozenset(map(_get_dtype_from_object, (binary_type,
                                                        text_type)))


def pandas_dtype(dtype):
    """
    Converts input into a pandas only dtype object or a numpy dtype object.

    Parameters
    ----------
    dtype : object to be converted

    Returns
    -------
    np.dtype or a pandas dtype
    """
    if isinstance(dtype, DatetimeTZDtype):
        return dtype
    elif isinstance(dtype, PeriodDtype):
        return dtype
    elif isinstance(dtype, CategoricalDtype):
        return dtype
    elif isinstance(dtype, IntervalDtype):
        return dtype
    elif isinstance(dtype, string_types):
        try:
            return DatetimeTZDtype.construct_from_string(dtype)
        except TypeError:
            pass

        if dtype.startswith('period[') or dtype.startswith('Period['):
            # do not parse string like U as period[U]
            try:
                return PeriodDtype.construct_from_string(dtype)
            except TypeError:
                pass

        elif dtype.startswith('interval[') or dtype.startswith('Interval['):
            try:
                return IntervalDtype.construct_from_string(dtype)
            except TypeError:
                pass

        try:
            return CategoricalDtype.construct_from_string(dtype)
        except TypeError:
            pass
    elif isinstance(dtype, ExtensionDtype):
        return dtype

    try:
        npdtype = np.dtype(dtype)
    except (TypeError, ValueError):
        raise

    # Any invalid dtype (such as pd.Timestamp) should raise an error.
    # np.dtype(invalid_type).kind = 0 for such objects. However, this will
    # also catch some valid dtypes such as object, np.object_ and 'object'
    # which we safeguard against by catching them earlier and returning
    # np.dtype(valid_dtype) before this condition is evaluated.
    if dtype in [object, np.object_, 'object', 'O']:
        return npdtype
    elif npdtype.kind == 'O':
        raise TypeError('dtype {0} not understood'.format(dtype))

    return npdtype
