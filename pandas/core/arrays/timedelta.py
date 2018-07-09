# -*- coding: utf-8 -*-
from datetime import timedelta

import numpy as np

from pandas._libs import tslibs
from pandas._libs.tslibs import Timedelta, NaT
from pandas._libs.tslibs.fields import get_timedelta_field
from pandas._libs.tslibs.timedeltas import array_to_timedelta64

from pandas import compat

from pandas.core.dtypes.common import (
    _TD_DTYPE, _ensure_int64, is_timedelta64_dtype)
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.tseries.offsets import Tick, DateOffset
from pandas.tseries.frequencies import to_offset

from .datetimelike import DatetimeLikeArrayMixin


def _is_convertible_to_td(key):
    return isinstance(key, (Tick, timedelta,
                            np.timedelta64, compat.string_types))


def _field_accessor(name, alias, docstring=None):
    def f(self):
        values = self.asi8
        result = get_timedelta_field(values, alias)
        if self.hasnans:
            result = self._maybe_mask_results(result, convert='float64')

        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


class TimedeltaArrayMixin(DatetimeLikeArrayMixin):
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
    def _simple_new(cls, values, freq=None, **kwargs):
        values = np.array(values, copy=False)
        if values.dtype == np.object_:
            values = array_to_timedelta64(values)
        if values.dtype != _TD_DTYPE:
            if is_timedelta64_dtype(values):
                # non-nano unit
                values = values.astype(_TD_DTYPE)
            else:
                values = _ensure_int64(values).view(_TD_DTYPE)

        result = object.__new__(cls)
        result._data = values
        result._freq = freq
        return result

    def __new__(cls, values, freq=None):
        if (freq is not None and not isinstance(freq, DateOffset) and
                freq != 'infer'):
            freq = to_offset(freq)

        result = cls._simple_new(values, freq=freq)
        if freq == 'infer':
            inferred = result.inferred_freq
            if inferred:
                result._freq = to_offset(inferred)

        return result

    # ----------------------------------------------------------------
    # Arithmetic Methods

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        raise TypeError("cannot add the type {typ} to a {cls}"
                        .format(typ=type(other).__name__,
                                cls=type(self).__name__))

    def _sub_datelike(self, other):
        assert other is not NaT
        raise TypeError("cannot subtract a datelike from a {cls}"
                        .format(cls=type(self).__name__))

    def _add_delta(self, delta):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self.

        Parameters
        ----------
        delta : timedelta, np.timedelta64, Tick, TimedeltaArray, TimedeltaIndex

        Returns
        -------
        result : same type as self

        Notes
        -----
        The result's name is set outside of _add_delta by the calling
        method (__add__ or __sub__)
        """
        if isinstance(delta, (Tick, timedelta, np.timedelta64)):
            new_values = self._add_delta_td(delta)
        elif isinstance(delta, TimedeltaArrayMixin):
            new_values = self._add_delta_tdi(delta)
        elif is_timedelta64_dtype(delta):
            # ndarray[timedelta64] --> wrap in TimedeltaArray/Index
            delta = type(self)(delta)
            new_values = self._add_delta_tdi(delta)
        else:
            raise TypeError("cannot add the type {0} to a TimedeltaIndex"
                            .format(type(delta)))

        return type(self)(new_values, freq='infer')

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
                result = self._maybe_mask_results(result, convert='float64')
                return result

        return NotImplemented

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
        """
        return self._maybe_mask_results(1e-9 * self.asi8)

    def to_pytimedelta(self):
        """
        Return Timedelta Array/Index as object ndarray of datetime.timedelta
        objects

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
