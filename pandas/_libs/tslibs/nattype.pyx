# -*- coding: utf-8 -*-
# cython: profile=False
import warnings

from cpython cimport (
    PyFloat_Check, PyComplex_Check,
    PyObject_RichCompare,
    Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE)

from cpython.datetime cimport (datetime,
                               PyDateTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

import numpy as np
cimport numpy as np
from numpy cimport int64_t
np.import_array()

from util cimport (get_nat,
                   is_integer_object, is_float_object,
                   is_datetime64_object, is_timedelta64_object)

# ----------------------------------------------------------------------
# Constants
_nat_strings = set(['NaT', 'nat', 'NAT', 'nan', 'NaN', 'NAN'])

cdef int64_t NPY_NAT = get_nat()

cdef bint _nat_scalar_rules[6]
_nat_scalar_rules[Py_EQ] = False
_nat_scalar_rules[Py_NE] = True
_nat_scalar_rules[Py_LT] = False
_nat_scalar_rules[Py_LE] = False
_nat_scalar_rules[Py_GT] = False
_nat_scalar_rules[Py_GE] = False

# ----------------------------------------------------------------------


def _make_nan_func(func_name, cls):
    def f(*args, **kwargs):
        return np.nan
    f.__name__ = func_name
    f.__doc__ = getattr(cls, func_name).__doc__
    return f


def _make_nat_func(func_name, cls):
    def f(*args, **kwargs):
        return NaT

    f.__name__ = func_name
    if isinstance(cls, str):
        # passed the literal docstring directly
        f.__doc__ = cls
    else:
        f.__doc__ = getattr(cls, func_name).__doc__
    return f


def _make_error_func(func_name, cls):
    def f(*args, **kwargs):
        raise ValueError("NaTType does not support " + func_name)

    f.__name__ = func_name
    if isinstance(cls, str):
        # passed the literal docstring directly
        f.__doc__ = cls
    elif cls is not None:
        f.__doc__ = getattr(cls, func_name).__doc__
    return f


cdef _nat_divide_op(self, other):
    if PyDelta_Check(other) or is_timedelta64_object(other) or other is NaT:
        return np.nan
    if is_integer_object(other) or is_float_object(other):
        return NaT
    return NotImplemented


cdef _nat_rdivide_op(self, other):
    if PyDelta_Check(other):
        return np.nan
    return NotImplemented


def __nat_unpickle(*args):
    # return constant defined in the module
    return NaT

# ----------------------------------------------------------------------


cdef class _NaT(datetime):
    cdef readonly:
        int64_t value
        object freq

    def __hash__(_NaT self):
        # py3k needs this defined here
        return hash(self.value)

    def __richcmp__(_NaT self, object other, int op):
        cdef int ndim = getattr(other, 'ndim', -1)

        if ndim == -1:
            return _nat_scalar_rules[op]

        if ndim == 0:
            if is_datetime64_object(other):
                return _nat_scalar_rules[op]
            else:
                raise TypeError('Cannot compare type %r with type %r' %
                                (type(self).__name__, type(other).__name__))
        # Note: instead of passing "other, self, _reverse_ops[op]", we observe
        # that `_nat_scalar_rules` is invariant under `_reverse_ops`,
        # rendering it unnecessary.
        return PyObject_RichCompare(other, self, op)

    def __add__(self, other):
        if PyDateTime_Check(other):
            return NaT

        elif hasattr(other, 'delta'):
            # Timedelta, offsets.Tick, offsets.Week
            return NaT
        elif getattr(other, '_typ', None) in ['dateoffset', 'series',
                                              'period', 'datetimeindex',
                                              'timedeltaindex']:
            # Duplicate logic in _Timestamp.__add__ to avoid needing
            # to subclass; allows us to @final(_Timestamp.__add__)
            return NotImplemented
        return NaT

    def __sub__(self, other):
        # Duplicate some logic from _Timestamp.__sub__ to avoid needing
        # to subclass; allows us to @final(_Timestamp.__sub__)
        if PyDateTime_Check(other):
            return  NaT
        elif PyDelta_Check(other):
            return NaT

        elif getattr(other, '_typ', None) == 'datetimeindex':
            # a Timestamp-DatetimeIndex -> yields a negative TimedeltaIndex
            return -other.__sub__(self)

        elif getattr(other, '_typ', None) == 'timedeltaindex':
            # a Timestamp-TimedeltaIndex -> yields a negative TimedeltaIndex
            return (-other).__add__(self)

        elif hasattr(other, 'delta'):
            # offsets.Tick, offsets.Week
            neg_other = -other
            return self + neg_other

        elif getattr(other, '_typ', None) in ['period',
                                              'periodindex', 'dateoffset']:
            return NotImplemented

        return NaT

    def __pos__(self):
        return NaT

    def __neg__(self):
        return NaT

    def __div__(self, other):
        return _nat_divide_op(self, other)

    def __truediv__(self, other):
        return _nat_divide_op(self, other)

    def __floordiv__(self, other):
        return _nat_divide_op(self, other)

    def __mul__(self, other):
        if is_integer_object(other) or is_float_object(other):
            return NaT
        return NotImplemented

    @property
    def asm8(self):
        return np.datetime64(NPY_NAT, 'ns')

    def to_datetime64(self):
        """ Returns a numpy.datetime64 object with 'ns' precision """
        return np.datetime64('NaT')


class NaTType(_NaT):
    """(N)ot-(A)-(T)ime, the time equivalent of NaN"""

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        base.value = NPY_NAT
        base.freq = None

        return base

    def __repr__(self):
        return 'NaT'

    def __str__(self):
        return 'NaT'

    def isoformat(self, sep='T'):
        # This allows Timestamp(ts.isoformat()) to always correctly roundtrip.
        return 'NaT'

    def __hash__(self):
        return NPY_NAT

    def __int__(self):
        return NPY_NAT

    def __long__(self):
        return NPY_NAT

    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # http://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    def __reduce__(self):
        return (__nat_unpickle, (None, ))

    def total_seconds(self):
        """
        Total duration of timedelta in seconds (to ns precision)
        """
        # GH 10939
        return np.nan

    @property
    def is_leap_year(self):
        return False

    @property
    def is_month_start(self):
        return False

    @property
    def is_quarter_start(self):
        return False

    @property
    def is_year_start(self):
        return False

    @property
    def is_month_end(self):
        return False

    @property
    def is_quarter_end(self):
        return False

    @property
    def is_year_end(self):
        return False

    def __rdiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rtruediv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rfloordiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rmul__(self, other):
        if is_integer_object(other) or is_float_object(other):
            return NaT
        return NotImplemented

    # ----------------------------------------------------------------------
    # inject the Timestamp field properties
    # these by definition return np.nan

    year = property(fget=lambda self: np.nan)
    quarter = property(fget=lambda self: np.nan)
    month = property(fget=lambda self: np.nan)
    day = property(fget=lambda self: np.nan)
    hour = property(fget=lambda self: np.nan)
    minute = property(fget=lambda self: np.nan)
    second = property(fget=lambda self: np.nan)
    millisecond = property(fget=lambda self: np.nan)
    microsecond = property(fget=lambda self: np.nan)
    nanosecond = property(fget=lambda self: np.nan)

    week = property(fget=lambda self: np.nan)
    dayofyear = property(fget=lambda self: np.nan)
    weekofyear = property(fget=lambda self: np.nan)
    days_in_month = property(fget=lambda self: np.nan)
    daysinmonth = property(fget=lambda self: np.nan)
    dayofweek = property(fget=lambda self: np.nan)
    weekday_name = property(fget=lambda self: np.nan)

    # inject Timedelta properties
    days = property(fget=lambda self: np.nan)
    seconds = property(fget=lambda self: np.nan)
    microseconds = property(fget=lambda self: np.nan)
    nanoseconds = property(fget=lambda self: np.nan)

    # inject pd.Period properties
    qyear = property(fget=lambda self: np.nan)

    # ----------------------------------------------------------------------
    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    # These are the ones that can get their docstrings from datetime.

    # nan methods
    weekday = _make_nan_func('weekday', datetime)
    isoweekday = _make_nan_func('isoweekday', datetime)

    # _nat_methods
    date = _make_nat_func('date', datetime)

    utctimetuple = _make_error_func('utctimetuple', datetime)
    timetz = _make_error_func('timetz', datetime)
    timetuple = _make_error_func('timetuple', datetime)
    strptime = _make_error_func('strptime', datetime)
    strftime = _make_error_func('strftime', datetime)
    isocalendar = _make_error_func('isocalendar', datetime)
    dst = _make_error_func('dst', datetime)
    ctime = _make_error_func('ctime', datetime)
    time = _make_error_func('time', datetime)
    toordinal = _make_error_func('toordinal', datetime)
    tzname = _make_error_func('tzname', datetime)
    utcoffset = _make_error_func('utcoffset', datetime)

    # Timestamp has empty docstring for some methods.
    utcfromtimestamp = _make_error_func('utcfromtimestamp', None) 
    fromtimestamp = _make_error_func('fromtimestamp', None)
    combine = _make_error_func('combine', None)
    utcnow = _make_error_func('utcnow', None)

    # ----------------------------------------------------------------------
    # The remaining methods are created with empty docstrings that will
    # be patched with the `Timestamp` versions once those are imported.

    timestamp = _make_error_func('timestamp', '')

    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    astimezone = _make_error_func('astimezone', '')
    fromordinal = _make_error_func('fromordinal', '')

    # _nat_methods
    to_pydatetime = _make_nat_func('to_pydatetime', '')

    now = _make_nat_func('now', '')
    today = _make_nat_func('today', '')
    round = _make_nat_func('round', '')
    floor = _make_nat_func('floor', '')
    ceil = _make_nat_func('ceil', '')

    tz_convert = _make_nat_func('tz_convert', '')
    tz_localize = _make_nat_func('tz_localize', '')
    replace = _make_nat_func('replace', '')

    def to_datetime(self):
        """
        DEPRECATED: use :meth:`to_pydatetime` instead.

        Convert a Timestamp object to a native Python datetime object.
        """
        warnings.warn("to_datetime is deprecated. Use self.to_pydatetime()",
                      FutureWarning, stacklevel=2)
        return self.to_pydatetime(warn=False)


NaT = NaTType()


# ----------------------------------------------------------------------

cdef inline bint _checknull_with_nat(object val):
    """ utility to check if a value is a nat or not """
    return val is None or (
        PyFloat_Check(val) and val != val) or val is NaT
