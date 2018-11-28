# -*- coding: utf-8 -*-
from __future__ import division

from datetime import timedelta
import operator
import warnings

import numpy as np

from pandas._libs import tslibs
from pandas._libs.tslibs import NaT, Timedelta, Timestamp, iNaT
from pandas._libs.tslibs.fields import get_timedelta_field
from pandas._libs.tslibs.timedeltas import (
    array_to_timedelta64, parse_timedelta_unit)
import pandas.compat as compat
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    _TD_DTYPE, ensure_int64, is_datetime64_dtype, is_float_dtype,
    is_integer_dtype, is_list_like, is_object_dtype, is_string_dtype,
    is_timedelta64_dtype)
from pandas.core.dtypes.generic import (
    ABCDataFrame, ABCIndexClass, ABCSeries, ABCTimedeltaIndex)
from pandas.core.dtypes.missing import isna

from pandas.core import ops
from pandas.core.algorithms import checked_add_with_arr
import pandas.core.common as com

from pandas.tseries.frequencies import to_offset
from pandas.tseries.offsets import Tick

from . import datetimelike as dtl


def _to_m8(key):
    """
    Timedelta-like => dt64
    """
    if not isinstance(key, Timedelta):
        # this also converts strings
        key = Timedelta(key)

    # return an type that can be compared
    return np.int64(key.value).view(_TD_DTYPE)


def _is_convertible_to_td(key):
    return isinstance(key, (Tick, timedelta,
                            np.timedelta64, compat.string_types))


def _field_accessor(name, alias, docstring=None):
    def f(self):
        values = self.asi8
        result = get_timedelta_field(values, alias)
        if self.hasnans:
            result = self._maybe_mask_results(result, fill_value=None,
                                              convert='float64')

        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _td_array_cmp(cls, op):
    """
    Wrap comparison operations to convert timedelta-like to timedelta64
    """
    opname = '__{name}__'.format(name=op.__name__)
    nat_result = True if opname == '__ne__' else False

    def wrapper(self, other):
        msg = "cannot compare a {cls} with type {typ}"
        meth = getattr(dtl.DatetimeLikeArrayMixin, opname)
        if _is_convertible_to_td(other) or other is NaT:
            try:
                other = _to_m8(other)
            except ValueError:
                # failed to parse as timedelta
                raise TypeError(msg.format(cls=type(self).__name__,
                                           typ=type(other).__name__))
            result = meth(self, other)
            if isna(other):
                result.fill(nat_result)

        elif not is_list_like(other):
            raise TypeError(msg.format(cls=type(self).__name__,
                                       typ=type(other).__name__))
        else:
            other = type(self)(other)._data
            result = meth(self, other)
            result = com.values_from_object(result)

            o_mask = np.array(isna(other))
            if o_mask.any():
                result[o_mask] = nat_result

        if self.hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


def _wrap_tdi_op(op):
    """
    Instead of re-implementing multiplication/division etc operations
    in the Array class, for now we dispatch to the TimedeltaIndex
    implementations.
    """
    # TODO: implement directly here and wrap in TimedeltaIndex, instead of
    #  the other way around
    def method(self, other):
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        from pandas import TimedeltaIndex
        obj = TimedeltaIndex(self)
        result = op(obj, other)
        if is_timedelta64_dtype(result):
            return type(self)(result)
        return np.array(result)

    method.__name__ = '__{name}__'.format(name=op.__name__)
    return method


class TimedeltaArrayMixin(dtl.DatetimeLikeArrayMixin):
    _typ = "timedeltaarray"
    __array_priority__ = 1000

    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit='ns')

    @property
    def dtype(self):
        return _TD_DTYPE

    # ----------------------------------------------------------------
    # Constructors
    _attributes = ["freq"]

    @classmethod
    def _simple_new(cls, values, freq=None, dtype=_TD_DTYPE):
        # `dtype` is passed by _shallow_copy in corner cases, should always
        #  be timedelta64[ns] if present
        assert dtype == _TD_DTYPE
        assert isinstance(values, np.ndarray), type(values)

        if values.dtype == 'i8':
            values = values.view('m8[ns]')

        assert values.dtype == 'm8[ns]'

        result = object.__new__(cls)
        result._data = values
        result._freq = freq
        return result

    def __new__(cls, values, freq=None, dtype=_TD_DTYPE):

        freq, freq_infer = dtl.maybe_infer_freq(freq)

        values = np.array(values, copy=False)
        if values.dtype == np.object_:
            values = array_to_timedelta64(values)

        result = cls._simple_new(values, freq=freq)
        if freq_infer:
            result.freq = to_offset(result.inferred_freq)

        return result

    @classmethod
    def _generate_range(cls, start, end, periods, freq, closed=None):

        periods = dtl.validate_periods(periods)
        if freq is None and any(x is None for x in [periods, start, end]):
            raise ValueError('Must provide freq argument if no data is '
                             'supplied')

        if com.count_not_none(start, end, periods, freq) != 3:
            raise ValueError('Of the four parameters: start, end, periods, '
                             'and freq, exactly three must be specified')

        if start is not None:
            start = Timedelta(start)

        if end is not None:
            end = Timedelta(end)

        if start is None and end is None:
            if closed is not None:
                raise ValueError("Closed has to be None if not both of start"
                                 "and end are defined")

        left_closed, right_closed = dtl.validate_endpoints(closed)

        if freq is not None:
            index = _generate_regular_range(start, end, periods, freq)
        else:
            index = np.linspace(start.value, end.value, periods).astype('i8')

        if not left_closed:
            index = index[1:]
        if not right_closed:
            index = index[:-1]

        return cls._simple_new(index, freq=freq)

    # ----------------------------------------------------------------
    # Array-Like / EA-Interface Methods

    @Appender(dtl.DatetimeLikeArrayMixin._validate_fill_value.__doc__)
    def _validate_fill_value(self, fill_value):
        if isna(fill_value):
            fill_value = iNaT
        elif isinstance(fill_value, (timedelta, np.timedelta64, Tick)):
            fill_value = Timedelta(fill_value).value
        else:
            raise ValueError("'fill_value' should be a Timedelta. "
                             "Got '{got}'.".format(got=fill_value))
        return fill_value

    # ----------------------------------------------------------------
    # Arithmetic Methods

    _create_comparison_method = classmethod(_td_array_cmp)

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        raise TypeError("cannot add the type {typ} to a {cls}"
                        .format(typ=type(other).__name__,
                                cls=type(self).__name__))

    def _add_delta(self, delta):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self, yielding a new TimedeltaArray.

        Parameters
        ----------
        other : {timedelta, np.timedelta64, Tick,
                 TimedeltaIndex, ndarray[timedelta64]}

        Returns
        -------
        result : TimedeltaArray
        """
        new_values = dtl.DatetimeLikeArrayMixin._add_delta(self, delta)
        return type(self)(new_values, freq='infer')

    def _add_datetime_arraylike(self, other):
        """
        Add DatetimeArray/Index or ndarray[datetime64] to TimedeltaArray.
        """
        if isinstance(other, np.ndarray):
            # At this point we have already checked that dtype is datetime64
            from pandas.core.arrays import DatetimeArrayMixin
            other = DatetimeArrayMixin(other)

        # defer to implementation in DatetimeArray
        return other + self

    def _add_datetimelike_scalar(self, other):
        # adding a timedeltaindex to a datetimelike
        from pandas.core.arrays import DatetimeArrayMixin

        assert other is not NaT
        other = Timestamp(other)
        if other is NaT:
            # In this case we specifically interpret NaT as a datetime, not
            # the timedelta interpretation we would get by returning self + NaT
            result = self.asi8.view('m8[ms]') + NaT.to_datetime64()
            return DatetimeArrayMixin(result)

        i8 = self.asi8
        result = checked_add_with_arr(i8, other.value,
                                      arr_mask=self._isnan)
        result = self._maybe_mask_results(result)
        return DatetimeArrayMixin(result, tz=other.tz)

    def _addsub_offset_array(self, other, op):
        # Add or subtract Array-like of DateOffset objects
        try:
            # TimedeltaIndex can only operate with a subset of DateOffset
            # subclasses.  Incompatible classes will raise AttributeError,
            # which we re-raise as TypeError
            return dtl.DatetimeLikeArrayMixin._addsub_offset_array(self, other,
                                                                   op)
        except AttributeError:
            raise TypeError("Cannot add/subtract non-tick DateOffset to {cls}"
                            .format(cls=type(self).__name__))

    def _evaluate_with_timedelta_like(self, other, op):
        if isinstance(other, ABCSeries):
            # GH#19042
            return NotImplemented

        opstr = '__{opname}__'.format(opname=op.__name__).replace('__r', '__')
        # allow division by a timedelta
        if opstr in ['__div__', '__truediv__', '__floordiv__']:
            if _is_convertible_to_td(other):
                other = Timedelta(other)
                if isna(other):
                    raise NotImplementedError(
                        "division by pd.NaT not implemented")

                i8 = self.asi8
                left, right = i8, other.value

                if opstr in ['__floordiv__']:
                    result = op(left, right)
                else:
                    result = op(left, np.float64(right))
                result = self._maybe_mask_results(result, fill_value=None,
                                                  convert='float64')
                return result

        return NotImplemented

    __mul__ = _wrap_tdi_op(operator.mul)
    __rmul__ = __mul__
    __truediv__ = _wrap_tdi_op(operator.truediv)
    __floordiv__ = _wrap_tdi_op(operator.floordiv)
    __rfloordiv__ = _wrap_tdi_op(ops.rfloordiv)

    if compat.PY2:
        __div__ = __truediv__

    # Note: TimedeltaIndex overrides this in call to cls._add_numeric_methods
    def __neg__(self):
        if self.freq is not None:
            return type(self)(-self._data, freq=-self.freq)
        return type(self)(-self._data)

    def __abs__(self):
        # Note: freq is not preserved
        return type(self)(np.abs(self._data))

    # ----------------------------------------------------------------
    # Conversion Methods - Vectorized analogues of Timedelta methods

    def total_seconds(self):
        """
        Return total duration of each element expressed in seconds.

        This method is available directly on TimedeltaArray, TimedeltaIndex
        and on Series containing timedelta values under the ``.dt`` namespace.

        Returns
        -------
        seconds : [ndarray, Float64Index, Series]
            When the calling object is a TimedeltaArray, the return type
            is ndarray.  When the calling object is a TimedeltaIndex,
            the return type is a Float64Index. When the calling object
            is a Series, the return type is Series of type `float64` whose
            index is the same as the original.

        See Also
        --------
        datetime.timedelta.total_seconds : Standard library version
            of this method.
        TimedeltaIndex.components : Return a DataFrame with components of
            each Timedelta.

        Examples
        --------
        **Series**

        >>> s = pd.Series(pd.to_timedelta(np.arange(5), unit='d'))
        >>> s
        0   0 days
        1   1 days
        2   2 days
        3   3 days
        4   4 days
        dtype: timedelta64[ns]

        >>> s.dt.total_seconds()
        0         0.0
        1     86400.0
        2    172800.0
        3    259200.0
        4    345600.0
        dtype: float64

        **TimedeltaIndex**

        >>> idx = pd.to_timedelta(np.arange(5), unit='d')
        >>> idx
        TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                       dtype='timedelta64[ns]', freq=None)

        >>> idx.total_seconds()
        Float64Index([0.0, 86400.0, 172800.0, 259200.00000000003, 345600.0],
                     dtype='float64')
        """
        return self._maybe_mask_results(1e-9 * self.asi8, fill_value=None)

    def to_pytimedelta(self):
        """
        Return Timedelta Array/Index as object ndarray of datetime.timedelta
        objects.

        Returns
        -------
        datetimes : ndarray
        """
        return tslibs.ints_to_pytimedelta(self.asi8)

    days = _field_accessor("days", "days",
                           " Number of days for each element. ")
    seconds = _field_accessor("seconds", "seconds",
                              " Number of seconds (>= 0 and less than 1 day) "
                              "for each element. ")
    microseconds = _field_accessor("microseconds", "microseconds",
                                   "\nNumber of microseconds (>= 0 and less "
                                   "than 1 second) for each\nelement. ")
    nanoseconds = _field_accessor("nanoseconds", "nanoseconds",
                                  "\nNumber of nanoseconds (>= 0 and less "
                                  "than 1 microsecond) for each\nelement.\n")

    @property
    def components(self):
        """
        Return a dataframe of the components (days, hours, minutes,
        seconds, milliseconds, microseconds, nanoseconds) of the Timedeltas.

        Returns
        -------
        a DataFrame
        """
        from pandas import DataFrame

        columns = ['days', 'hours', 'minutes', 'seconds',
                   'milliseconds', 'microseconds', 'nanoseconds']
        hasnans = self.hasnans
        if hasnans:
            def f(x):
                if isna(x):
                    return [np.nan] * len(columns)
                return x.components
        else:
            def f(x):
                return x.components

        result = DataFrame([f(x) for x in self], columns=columns)
        if not hasnans:
            result = result.astype('int64')
        return result


TimedeltaArrayMixin._add_comparison_ops()
TimedeltaArrayMixin._add_datetimelike_methods()


# ---------------------------------------------------------------------
# Constructor Helpers

def sequence_to_td64ns(data, copy=False, unit="ns", errors="raise"):
    """
    Parameters
    ----------
    array : list-like
    copy : bool, default False
    unit : str, default "ns"
        The timedelta unit to treat integers as multiples of.
    errors : {"raise", "coerce", "ignore"}, default "raise"
        How to handle elements that cannot be converted to timedelta64[ns].
        See ``pandas.to_timedelta`` for details.

    Returns
    -------
    converted : numpy.ndarray
        The sequence converted to a numpy array with dtype ``timedelta64[ns]``.
    inferred_freq : Tick or None
        The inferred frequency of the sequence.

    Raises
    ------
    ValueError : Data cannot be converted to timedelta64[ns].

    Notes
    -----
    Unlike `pandas.to_timedelta`, if setting ``errors=ignore`` will not cause
    errors to be ignored; they are caught and subsequently ignored at a
    higher level.
    """
    inferred_freq = None
    unit = parse_timedelta_unit(unit)

    # Unwrap whatever we have into a np.ndarray
    if not hasattr(data, 'dtype'):
        # e.g. list, tuple
        if np.ndim(data) == 0:
            # i.e. generator
            data = list(data)
        data = np.array(data, copy=False)
    elif isinstance(data, ABCSeries):
        data = data._values
    elif isinstance(data, (ABCTimedeltaIndex, TimedeltaArrayMixin)):
        inferred_freq = data.freq
        data = data._data

    # Convert whatever we have into timedelta64[ns] dtype
    if is_object_dtype(data) or is_string_dtype(data):
        # no need to make a copy, need to convert if string-dtyped
        data = objects_to_td64ns(data, unit=unit, errors=errors)
        copy = False

    elif is_integer_dtype(data):
        # treat as multiples of the given unit
        data, copy_made = ints_to_td64ns(data, unit=unit)
        copy = copy and not copy_made

    elif is_float_dtype(data):
        # treat as multiples of the given unit.  If after converting to nanos,
        #  there are fractional components left, these are truncated
        #  (i.e. NOT rounded)
        mask = np.isnan(data)
        coeff = np.timedelta64(1, unit) / np.timedelta64(1, 'ns')
        data = (coeff * data).astype(np.int64).view('timedelta64[ns]')
        data[mask] = iNaT
        copy = False

    elif is_timedelta64_dtype(data):
        if data.dtype != _TD_DTYPE:
            # non-nano unit
            # TODO: watch out for overflows
            data = data.astype(_TD_DTYPE)
            copy = False

    elif is_datetime64_dtype(data):
        # GH#23539
        warnings.warn("Passing datetime64-dtype data to TimedeltaIndex is "
                      "deprecated, will raise a TypeError in a future "
                      "version",
                      FutureWarning, stacklevel=3)
        data = ensure_int64(data).view(_TD_DTYPE)

    else:
        raise TypeError("dtype {dtype} cannot be converted to timedelta64[ns]"
                        .format(dtype=data.dtype))

    data = np.array(data, copy=copy)
    assert data.dtype == 'm8[ns]', data
    return data, inferred_freq


def ints_to_td64ns(data, unit="ns"):
    """
    Convert an ndarray with integer-dtype to timedelta64[ns] dtype, treating
    the integers as multiples of the given timedelta unit.

    Parameters
    ----------
    data : numpy.ndarray with integer-dtype
    unit : str, default "ns"
        The timedelta unit to treat integers as multiples of.

    Returns
    -------
    numpy.ndarray : timedelta64[ns] array converted from data
    bool : whether a copy was made
    """
    copy_made = False
    unit = unit if unit is not None else "ns"

    if data.dtype != np.int64:
        # converting to int64 makes a copy, so we can avoid
        # re-copying later
        data = data.astype(np.int64)
        copy_made = True

    if unit != "ns":
        dtype_str = "timedelta64[{unit}]".format(unit=unit)
        data = data.view(dtype_str)

        # TODO: watch out for overflows when converting from lower-resolution
        data = data.astype("timedelta64[ns]")
        # the astype conversion makes a copy, so we can avoid re-copying later
        copy_made = True

    else:
        data = data.view("timedelta64[ns]")

    return data, copy_made


def objects_to_td64ns(data, unit="ns", errors="raise"):
    """
    Convert a object-dtyped or string-dtyped array into an
    timedelta64[ns]-dtyped array.

    Parameters
    ----------
    data : ndarray or Index
    unit : str, default "ns"
        The timedelta unit to treat integers as multiples of.
    errors : {"raise", "coerce", "ignore"}, default "raise"
        How to handle elements that cannot be converted to timedelta64[ns].
        See ``pandas.to_timedelta`` for details.

    Returns
    -------
    numpy.ndarray : timedelta64[ns] array converted from data

    Raises
    ------
    ValueError : Data cannot be converted to timedelta64[ns].

    Notes
    -----
    Unlike `pandas.to_timedelta`, if setting `errors=ignore` will not cause
    errors to be ignored; they are caught and subsequently ignored at a
    higher level.
    """
    # coerce Index to np.ndarray, converting string-dtype if necessary
    values = np.array(data, dtype=np.object_, copy=False)

    result = array_to_timedelta64(values,
                                  unit=unit, errors=errors)
    return result.view('timedelta64[ns]')


def _generate_regular_range(start, end, periods, offset):
    stride = offset.nanos
    if periods is None:
        b = Timedelta(start).value
        e = Timedelta(end).value
        e += stride - e % stride
    elif start is not None:
        b = Timedelta(start).value
        e = b + periods * stride
    elif end is not None:
        e = Timedelta(end).value + stride
        b = e - periods * stride
    else:
        raise ValueError("at least 'start' or 'end' should be specified "
                         "if a 'period' is given.")

    data = np.arange(b, e, stride, dtype=np.int64)
    return data
