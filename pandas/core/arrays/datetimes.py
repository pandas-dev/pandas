# -*- coding: utf-8 -*-
from datetime import datetime, timedelta, time
import warnings

import numpy as np
from pytz import utc

from pandas._libs import lib, tslib
from pandas._libs.tslib import Timestamp, NaT, iNaT
from pandas._libs.tslibs import (
    normalize_date,
    conversion, fields, timezones,
    resolution as libresolution)

from pandas.util._decorators import cache_readonly
from pandas.errors import PerformanceWarning
from pandas import compat

from pandas.core.dtypes.common import (
    _NS_DTYPE,
    is_object_dtype,
    is_datetime64tz_dtype,
    is_datetime64_dtype,
    is_timedelta64_dtype,
    ensure_int64)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import isna
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries

import pandas.core.common as com
from pandas.core.algorithms import checked_add_with_arr
from pandas.core import ops

from pandas.tseries.frequencies import to_offset, get_period_alias
from pandas.tseries.offsets import Tick, generate_range

from pandas.core.arrays import datetimelike as dtl


_midnight = time(0, 0)


def _to_m8(key, tz=None):
    """
    Timestamp-like => dt64
    """
    if not isinstance(key, Timestamp):
        # this also converts strings
        key = Timestamp(key, tz=tz)

    return np.int64(conversion.pydt_to_i8(key)).view(_NS_DTYPE)


def _field_accessor(name, field, docstring=None):
    def f(self):
        values = self.asi8
        if self.tz is not None:
            if self.tz is not utc:
                values = self._local_timestamps()

        if field in self._bool_ops:
            if field.endswith(('start', 'end')):
                freq = self.freq
                month_kw = 12
                if freq:
                    kwds = freq.kwds
                    month_kw = kwds.get('startingMonth', kwds.get('month', 12))

                result = fields.get_start_end_field(values, field,
                                                    self.freqstr, month_kw)
            else:
                result = fields.get_date_field(values, field)

            # these return a boolean by-definition
            return result

        if field in self._object_ops:
            result = fields.get_date_name_field(values, field)
            result = self._maybe_mask_results(result)

        else:
            result = fields.get_date_field(values, field)
            result = self._maybe_mask_results(result, convert='float64')

        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _dt_array_cmp(cls, op):
    """
    Wrap comparison operations to convert datetime-like to datetime64
    """
    opname = '__{name}__'.format(name=op.__name__)
    nat_result = True if opname == '__ne__' else False

    def wrapper(self, other):
        meth = getattr(dtl.DatetimeLikeArrayMixin, opname)

        if isinstance(other, (datetime, np.datetime64, compat.string_types)):
            if isinstance(other, (datetime, np.datetime64)):
                # GH#18435 strings get a pass from tzawareness compat
                self._assert_tzawareness_compat(other)

            try:
                other = _to_m8(other, tz=self.tz)
            except ValueError:
                # string that cannot be parsed to Timestamp
                return ops.invalid_comparison(self, other, op)

            result = meth(self, other)
            if isna(other):
                result.fill(nat_result)
        elif lib.is_scalar(other):
            return ops.invalid_comparison(self, other, op)
        else:
            if isinstance(other, list):
                # FIXME: This can break for object-dtype with mixed types
                other = type(self)(other)
            elif not isinstance(other, (np.ndarray, ABCIndexClass, ABCSeries)):
                # Following Timestamp convention, __eq__ is all-False
                # and __ne__ is all True, others raise TypeError.
                return ops.invalid_comparison(self, other, op)

            if is_object_dtype(other):
                result = op(self.astype('O'), np.array(other))
            elif not (is_datetime64_dtype(other) or
                      is_datetime64tz_dtype(other)):
                # e.g. is_timedelta64_dtype(other)
                return ops.invalid_comparison(self, other, op)
            else:
                self._assert_tzawareness_compat(other)
                result = meth(self, np.asarray(other))

            result = com.values_from_object(result)

            # Make sure to pass an array to result[...]; indexing with
            # Series breaks with older version of numpy
            o_mask = np.array(isna(other))
            if o_mask.any():
                result[o_mask] = nat_result

        if self.hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


class DatetimeArrayMixin(dtl.DatetimeLikeArrayMixin):
    """
    Assumes that subclass __new__/__init__ defines:
        tz
        _freq
        _data
    """
    _bool_ops = ['is_month_start', 'is_month_end',
                 'is_quarter_start', 'is_quarter_end', 'is_year_start',
                 'is_year_end', 'is_leap_year']
    _object_ops = ['weekday_name', 'freq', 'tz']

    # dummy attribute so that datetime.__eq__(DatetimeArray) defers
    # by returning NotImplemented
    timetuple = None

    # -----------------------------------------------------------------
    # Constructors

    _attributes = ["freq", "tz"]

    @classmethod
    def _simple_new(cls, values, freq=None, tz=None, **kwargs):
        """
        we require the we have a dtype compat for the values
        if we are passed a non-dtype compat, then coerce using the constructor
        """

        if getattr(values, 'dtype', None) is None:
            # empty, but with dtype compat
            if values is None:
                values = np.empty(0, dtype=_NS_DTYPE)
                return cls(values, freq=freq, tz=tz, **kwargs)
            values = np.array(values, copy=False)

        if not is_datetime64_dtype(values):
            values = ensure_int64(values).view(_NS_DTYPE)

        result = object.__new__(cls)
        result._data = values
        result._freq = freq
        tz = timezones.maybe_get_tz(tz)
        result._tz = timezones.tz_standardize(tz)
        return result

    def __new__(cls, values, freq=None, tz=None, dtype=None):
        if tz is None and hasattr(values, 'tz'):
            # e.g. DatetimeIndex
            tz = values.tz

        if freq is None and hasattr(values, "freq"):
            # i.e. DatetimeArray, DatetimeIndex
            freq = values.freq

        freq, freq_infer = dtl.maybe_infer_freq(freq)

        # if dtype has an embedded tz, capture it
        tz = dtl.validate_tz_from_dtype(dtype, tz)

        result = cls._simple_new(values, freq=freq, tz=tz)
        if freq_infer:
            inferred = result.inferred_freq
            if inferred:
                result.freq = to_offset(inferred)

        # NB: Among other things not yet ported from the DatetimeIndex
        # constructor, this does not call _deepcopy_if_needed
        return result

    @classmethod
    def _generate_range(cls, start, end, periods, freq, tz=None,
                        normalize=False, ambiguous='raise', closed=None):

        periods = dtl.validate_periods(periods)
        if freq is None and any(x is None for x in [periods, start, end]):
            raise ValueError('Must provide freq argument if no data is '
                             'supplied')

        if com.count_not_none(start, end, periods, freq) != 3:
            raise ValueError('Of the four parameters: start, end, periods, '
                             'and freq, exactly three must be specified')
        freq = to_offset(freq)

        if start is not None:
            start = Timestamp(start)

        if end is not None:
            end = Timestamp(end)

        if start is None and end is None:
            if closed is not None:
                raise ValueError("Closed has to be None if not both of start"
                                 "and end are defined")

        left_closed, right_closed = dtl.validate_endpoints(closed)

        start, end, _normalized = _maybe_normalize_endpoints(start, end,
                                                             normalize)

        tz, _ = _infer_tz_from_endpoints(start, end, tz)

        if tz is not None:
            # Localize the start and end arguments
            start = _maybe_localize_point(
                start, getattr(start, 'tz', None), start, freq, tz
            )
            end = _maybe_localize_point(
                end, getattr(end, 'tz', None), end, freq, tz
            )
        if start and end:
            # Make sure start and end have the same tz
            start = _maybe_localize_point(
                start, start.tz, end.tz, freq, tz
            )
            end = _maybe_localize_point(
                end, end.tz, start.tz, freq, tz
            )
        if freq is not None:
            # TODO: consider re-implementing _cached_range; GH#17914
            index = _generate_regular_range(cls, start, end, periods, freq)

            if tz is not None and getattr(index, 'tz', None) is None:
                arr = conversion.tz_localize_to_utc(
                    ensure_int64(index.values),
                    tz, ambiguous=ambiguous)

                index = cls(arr)

                # index is localized datetime64 array -> have to convert
                # start/end as well to compare
                if start is not None:
                    start = start.tz_localize(tz).asm8
                if end is not None:
                    end = end.tz_localize(tz).asm8
        else:
            # Create a linearly spaced date_range in local time
            arr = np.linspace(start.value, end.value, periods)
            index = cls._simple_new(
                arr.astype('M8[ns]', copy=False), freq=None, tz=tz
            )

        if not left_closed and len(index) and index[0] == start:
            index = index[1:]
        if not right_closed and len(index) and index[-1] == end:
            index = index[:-1]

        return cls._simple_new(index.values, freq=freq, tz=tz)

    # -----------------------------------------------------------------
    # Descriptive Properties

    @property
    def _box_func(self):
        return lambda x: Timestamp(x, freq=self.freq, tz=self.tz)

    @cache_readonly
    def dtype(self):
        if self.tz is None:
            return _NS_DTYPE
        return DatetimeTZDtype('ns', self.tz)

    @property
    def tz(self):
        # GH 18595
        return self._tz

    @tz.setter
    def tz(self, value):
        # GH 3746: Prevent localizing or converting the index by setting tz
        raise AttributeError("Cannot directly set timezone. Use tz_localize() "
                             "or tz_convert() as appropriate")

    @property
    def tzinfo(self):
        """
        Alias for tz attribute
        """
        return self.tz

    @property  # NB: override with cache_readonly in immutable subclasses
    def _timezone(self):
        """ Comparable timezone both for pytz / dateutil"""
        return timezones.get_timezone(self.tzinfo)

    @property
    def offset(self):
        """get/set the frequency of the instance"""
        msg = ('{cls}.offset has been deprecated and will be removed '
               'in a future version; use {cls}.freq instead.'
               .format(cls=type(self).__name__))
        warnings.warn(msg, FutureWarning, stacklevel=2)
        return self.freq

    @offset.setter
    def offset(self, value):
        """get/set the frequency of the instance"""
        msg = ('{cls}.offset has been deprecated and will be removed '
               'in a future version; use {cls}.freq instead.'
               .format(cls=type(self).__name__))
        warnings.warn(msg, FutureWarning, stacklevel=2)
        self.freq = value

    @property  # NB: override with cache_readonly in immutable subclasses
    def is_normalized(self):
        """
        Returns True if all of the dates are at midnight ("no time")
        """
        return conversion.is_date_array_normalized(self.asi8, self.tz)

    @property  # NB: override with cache_readonly in immutable subclasses
    def _resolution(self):
        return libresolution.resolution(self.asi8, self.tz)

    # ----------------------------------------------------------------
    # Array-like Methods

    def __iter__(self):
        """
        Return an iterator over the boxed values

        Yields
        -------
        tstamp : Timestamp
        """

        # convert in chunks of 10k for efficiency
        data = self.asi8
        length = len(self)
        chunksize = 10000
        chunks = int(length / chunksize) + 1
        for i in range(chunks):
            start_i = i * chunksize
            end_i = min((i + 1) * chunksize, length)
            converted = tslib.ints_to_pydatetime(data[start_i:end_i],
                                                 tz=self.tz, freq=self.freq,
                                                 box="timestamp")
            for v in converted:
                yield v

    # -----------------------------------------------------------------
    # Comparison Methods

    _create_comparison_method = classmethod(_dt_array_cmp)

    def _has_same_tz(self, other):
        zzone = self._timezone

        # vzone sholdn't be None if value is non-datetime like
        if isinstance(other, np.datetime64):
            # convert to Timestamp as np.datetime64 doesn't have tz attr
            other = Timestamp(other)
        vzone = timezones.get_timezone(getattr(other, 'tzinfo', '__no_tz__'))
        return zzone == vzone

    def _assert_tzawareness_compat(self, other):
        # adapted from _Timestamp._assert_tzawareness_compat
        other_tz = getattr(other, 'tzinfo', None)
        if is_datetime64tz_dtype(other):
            # Get tzinfo from Series dtype
            other_tz = other.dtype.tz
        if other is NaT:
            # pd.NaT quacks both aware and naive
            pass
        elif self.tz is None:
            if other_tz is not None:
                raise TypeError('Cannot compare tz-naive and tz-aware '
                                'datetime-like objects.')
        elif other_tz is None:
            raise TypeError('Cannot compare tz-naive and tz-aware '
                            'datetime-like objects')

    # -----------------------------------------------------------------
    # Arithmetic Methods

    def _sub_datelike_dti(self, other):
        """subtraction of DatetimeArray/Index or ndarray[datetime64]"""
        if len(self) != len(other):
            raise ValueError("cannot add indices of unequal length")

        if not isinstance(other, DatetimeArrayMixin):
            # i.e. np.ndarray; we assume it is datetime64-dtype
            other = type(self)(other)

        if not self._has_same_tz(other):
            # require tz compat
            raise TypeError("{cls} subtraction must have the same "
                            "timezones or no timezones"
                            .format(cls=type(self).__name__))

        self_i8 = self.asi8
        other_i8 = other.asi8
        new_values = checked_add_with_arr(self_i8, -other_i8,
                                          arr_mask=self._isnan)
        if self.hasnans or other.hasnans:
            mask = (self._isnan) | (other._isnan)
            new_values[mask] = iNaT
        return new_values.view('timedelta64[ns]')

    def _add_offset(self, offset):
        assert not isinstance(offset, Tick)
        try:
            if self.tz is not None:
                values = self.tz_localize(None)
            else:
                values = self
            result = offset.apply_index(values)
            if self.tz is not None:
                result = result.tz_localize(self.tz)

        except NotImplementedError:
            warnings.warn("Non-vectorized DateOffset being applied to Series "
                          "or DatetimeIndex", PerformanceWarning)
            result = self.astype('O') + offset

        return type(self)(result, freq='infer')

    def _sub_datelike(self, other):
        # subtract a datetime from myself, yielding a ndarray[timedelta64[ns]]
        if isinstance(other, (DatetimeArrayMixin, np.ndarray)):
            return self._sub_datelike_dti(other)

        assert isinstance(other, (datetime, np.datetime64))
        assert other is not NaT
        other = Timestamp(other)
        if other is NaT:
            return self - NaT
        elif not self._has_same_tz(other):
            # require tz compat
            raise TypeError("Timestamp subtraction must have the same "
                            "timezones or no timezones")
        i8 = self.asi8
        result = checked_add_with_arr(i8, -other.value,
                                      arr_mask=self._isnan)
        result = self._maybe_mask_results(result,
                                          fill_value=iNaT)
        return result.view('timedelta64[ns]')

    def _add_delta(self, delta):
        """
        Add a timedelta-like, DateOffset, or TimedeltaIndex-like object
        to self.

        Parameters
        ----------
        delta : {timedelta, np.timedelta64, DateOffset,
                 TimedeltaIndex, ndarray[timedelta64]}

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
        elif is_timedelta64_dtype(delta):
            new_values = self._add_delta_tdi(delta)

        return type(self)(new_values, tz=self.tz, freq='infer')

    # -----------------------------------------------------------------
    # Timezone Conversion and Localization Methods

    def _local_timestamps(self):
        """
        Convert to an i8 (unix-like nanosecond timestamp) representation
        while keeping the local timezone and not using UTC.
        This is used to calculate time-of-day information as if the timestamps
        were timezone-naive.
        """
        return conversion.tz_convert(self.asi8, utc, self.tz)

    def tz_convert(self, tz):
        """
        Convert tz-aware Datetime Array/Index from one time zone to another.

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time. Corresponding timestamps would be converted
            to this time zone of the Datetime Array/Index. A `tz` of None will
            convert to UTC and remove the timezone information.

        Returns
        -------
        normalized : same type as self

        Raises
        ------
        TypeError
            If Datetime Array/Index is tz-naive.

        See Also
        --------
        DatetimeIndex.tz : A timezone that has a variable offset from UTC
        DatetimeIndex.tz_localize : Localize tz-naive DatetimeIndex to a
            given time zone, or remove timezone from a tz-aware DatetimeIndex.

        Examples
        --------
        With the `tz` parameter, we can change the DatetimeIndex
        to other time zones:

        >>> dti = pd.DatetimeIndex(start='2014-08-01 09:00',
        ...                        freq='H', periods=3, tz='Europe/Berlin')

        >>> dti
        DatetimeIndex(['2014-08-01 09:00:00+02:00',
                       '2014-08-01 10:00:00+02:00',
                       '2014-08-01 11:00:00+02:00'],
                      dtype='datetime64[ns, Europe/Berlin]', freq='H')

        >>> dti.tz_convert('US/Central')
        DatetimeIndex(['2014-08-01 02:00:00-05:00',
                       '2014-08-01 03:00:00-05:00',
                       '2014-08-01 04:00:00-05:00'],
                      dtype='datetime64[ns, US/Central]', freq='H')

        With the ``tz=None``, we can remove the timezone (after converting
        to UTC if necessary):

        >>> dti = pd.DatetimeIndex(start='2014-08-01 09:00',freq='H',
        ...                        periods=3, tz='Europe/Berlin')

        >>> dti
        DatetimeIndex(['2014-08-01 09:00:00+02:00',
                       '2014-08-01 10:00:00+02:00',
                       '2014-08-01 11:00:00+02:00'],
                        dtype='datetime64[ns, Europe/Berlin]', freq='H')

        >>> dti.tz_convert(None)
        DatetimeIndex(['2014-08-01 07:00:00',
                       '2014-08-01 08:00:00',
                       '2014-08-01 09:00:00'],
                        dtype='datetime64[ns]', freq='H')
        """
        tz = timezones.maybe_get_tz(tz)

        if self.tz is None:
            # tz naive, use tz_localize
            raise TypeError('Cannot convert tz-naive timestamps, use '
                            'tz_localize to localize')

        # No conversion since timestamps are all UTC to begin with
        return self._shallow_copy(tz=tz)

    def tz_localize(self, tz, ambiguous='raise', errors='raise'):
        """
        Localize tz-naive Datetime Array/Index to tz-aware
        Datetime Array/Index.

        This method takes a time zone (tz) naive Datetime Array/Index object
        and makes this time zone aware. It does not move the time to another
        time zone.
        Time zone localization helps to switch from time zone aware to time
        zone unaware objects.

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone to convert timestamps to. Passing ``None`` will
            remove the time zone information preserving local time.
        ambiguous : str {'infer', 'NaT', 'raise'} or bool array,
            default 'raise'

            - 'infer' will attempt to infer fall dst-transition hours based on
              order
            - bool-ndarray where True signifies a DST time, False signifies a
              non-DST time (note that this flag is only applicable for
              ambiguous times)
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an AmbiguousTimeError if there are ambiguous
              times

        errors : {'raise', 'coerce'}, default 'raise'

            - 'raise' will raise a NonExistentTimeError if a timestamp is not
              valid in the specified time zone (e.g. due to a transition from
              or to DST time)
            - 'coerce' will return NaT if the timestamp can not be converted
              to the specified time zone

            .. versionadded:: 0.19.0

        Returns
        -------
        result : same type as self
            Array/Index converted to the specified time zone.

        Raises
        ------
        TypeError
            If the Datetime Array/Index is tz-aware and tz is not None.

        See Also
        --------
        DatetimeIndex.tz_convert : Convert tz-aware DatetimeIndex from
            one time zone to another.

        Examples
        --------
        >>> tz_naive = pd.date_range('2018-03-01 09:00', periods=3)
        >>> tz_naive
        DatetimeIndex(['2018-03-01 09:00:00', '2018-03-02 09:00:00',
                       '2018-03-03 09:00:00'],
                      dtype='datetime64[ns]', freq='D')

        Localize DatetimeIndex in US/Eastern time zone:

        >>> tz_aware = tz_naive.tz_localize(tz='US/Eastern')
        >>> tz_aware
        DatetimeIndex(['2018-03-01 09:00:00-05:00',
                       '2018-03-02 09:00:00-05:00',
                       '2018-03-03 09:00:00-05:00'],
                      dtype='datetime64[ns, US/Eastern]', freq='D')

        With the ``tz=None``, we can remove the time zone information
        while keeping the local time (not converted to UTC):

        >>> tz_aware.tz_localize(None)
        DatetimeIndex(['2018-03-01 09:00:00', '2018-03-02 09:00:00',
                       '2018-03-03 09:00:00'],
                      dtype='datetime64[ns]', freq='D')
        """
        if self.tz is not None:
            if tz is None:
                new_dates = conversion.tz_convert(self.asi8, 'UTC', self.tz)
            else:
                raise TypeError("Already tz-aware, use tz_convert to convert.")
        else:
            tz = timezones.maybe_get_tz(tz)
            # Convert to UTC

            new_dates = conversion.tz_localize_to_utc(self.asi8, tz,
                                                      ambiguous=ambiguous,
                                                      errors=errors)
        new_dates = new_dates.view(_NS_DTYPE)
        return self._shallow_copy(new_dates, tz=tz)

    # ----------------------------------------------------------------
    # Conversion Methods - Vectorized analogues of Timestamp methods

    def to_pydatetime(self):
        """
        Return Datetime Array/Index as object ndarray of datetime.datetime
        objects

        Returns
        -------
        datetimes : ndarray
        """
        return tslib.ints_to_pydatetime(self.asi8, tz=self.tz)

    def normalize(self):
        """
        Convert times to midnight.

        The time component of the date-time is converted to midnight i.e.
        00:00:00. This is useful in cases, when the time does not matter.
        Length is unaltered. The timezones are unaffected.

        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on Datetime Array/Index.

        Returns
        -------
        DatetimeArray, DatetimeIndex or Series
            The same type as the original data. Series will have the same
            name and index. DatetimeIndex will have the same name.

        See Also
        --------
        floor : Floor the datetimes to the specified freq.
        ceil : Ceil the datetimes to the specified freq.
        round : Round the datetimes to the specified freq.

        Examples
        --------
        >>> idx = pd.DatetimeIndex(start='2014-08-01 10:00', freq='H',
        ...                        periods=3, tz='Asia/Calcutta')
        >>> idx
        DatetimeIndex(['2014-08-01 10:00:00+05:30',
                       '2014-08-01 11:00:00+05:30',
                       '2014-08-01 12:00:00+05:30'],
                        dtype='datetime64[ns, Asia/Calcutta]', freq='H')
        >>> idx.normalize()
        DatetimeIndex(['2014-08-01 00:00:00+05:30',
                       '2014-08-01 00:00:00+05:30',
                       '2014-08-01 00:00:00+05:30'],
                       dtype='datetime64[ns, Asia/Calcutta]', freq=None)
        """
        new_values = conversion.normalize_i8_timestamps(self.asi8, self.tz)
        return type(self)(new_values, freq='infer').tz_localize(self.tz)

    def to_period(self, freq=None):
        """
        Cast to PeriodArray/Index at a particular frequency.

        Converts DatetimeArray/Index to PeriodArray/Index.

        Parameters
        ----------
        freq : string or Offset, optional
            One of pandas' :ref:`offset strings <timeseries.offset_aliases>`
            or an Offset object. Will be inferred by default.

        Returns
        -------
        PeriodArray/Index

        Raises
        ------
        ValueError
            When converting a DatetimeArray/Index with non-regular values,
            so that a frequency cannot be inferred.

        Examples
        --------
        >>> df = pd.DataFrame({"y": [1,2,3]},
        ...                   index=pd.to_datetime(["2000-03-31 00:00:00",
        ...                                         "2000-05-31 00:00:00",
        ...                                         "2000-08-31 00:00:00"]))
        >>> df.index.to_period("M")
        PeriodIndex(['2000-03', '2000-05', '2000-08'],
                    dtype='period[M]', freq='M')

        Infer the daily frequency

        >>> idx = pd.date_range("2017-01-01", periods=2)
        >>> idx.to_period()
        PeriodIndex(['2017-01-01', '2017-01-02'],
                    dtype='period[D]', freq='D')

        See also
        --------
        pandas.PeriodIndex: Immutable ndarray holding ordinal values
        pandas.DatetimeIndex.to_pydatetime: Return DatetimeIndex as object
        """
        from pandas.core.arrays.period import PeriodArrayMixin

        if self.tz is not None:
            warnings.warn("Converting to PeriodArray/Index representation "
                          "will drop timezone information.", UserWarning)

        if freq is None:
            freq = self.freqstr or self.inferred_freq

            if freq is None:
                raise ValueError("You must pass a freq argument as "
                                 "current index has none.")

            freq = get_period_alias(freq)

        return PeriodArrayMixin(self.values, freq=freq)

    # -----------------------------------------------------------------
    # Properties - Vectorized Timestamp Properties/Methods

    def month_name(self, locale=None):
        """
        Return the month names of the DateTimeIndex with specified locale.

        .. versionadded:: 0.23.0

        Parameters
        ----------
        locale : str, optional
            Locale determining the language in which to return the month name.
            Default is English locale.

        Returns
        -------
        Index
            Index of month names.

        Examples
        --------
        >>> idx = pd.DatetimeIndex(start='2018-01', freq='M', periods=3)
        >>> idx
        DatetimeIndex(['2018-01-31', '2018-02-28', '2018-03-31'],
                      dtype='datetime64[ns]', freq='M')
        >>> idx.month_name()
        Index(['January', 'February', 'March'], dtype='object')
        """
        if self.tz is not None and self.tz is not utc:
            values = self._local_timestamps()
        else:
            values = self.asi8

        result = fields.get_date_name_field(values, 'month_name',
                                            locale=locale)
        result = self._maybe_mask_results(result)
        return result

    def day_name(self, locale=None):
        """
        Return the day names of the DateTimeIndex with specified locale.

        .. versionadded:: 0.23.0

        Parameters
        ----------
        locale : str, optional
            Locale determining the language in which to return the day name.
            Default is English locale.

        Returns
        -------
        Index
            Index of day names.

        Examples
        --------
        >>> idx = pd.DatetimeIndex(start='2018-01-01', freq='D', periods=3)
        >>> idx
        DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03'],
                      dtype='datetime64[ns]', freq='D')
        >>> idx.day_name()
        Index(['Monday', 'Tuesday', 'Wednesday'], dtype='object')
        """
        if self.tz is not None and self.tz is not utc:
            values = self._local_timestamps()
        else:
            values = self.asi8

        result = fields.get_date_name_field(values, 'day_name',
                                            locale=locale)
        result = self._maybe_mask_results(result)
        return result

    @property
    def time(self):
        """
        Returns numpy array of datetime.time. The time part of the Timestamps.
        """
        # If the Timestamps have a timezone that is not UTC,
        # convert them into their i8 representation while
        # keeping their timezone and not using UTC
        if self.tz is not None and self.tz is not utc:
            timestamps = self._local_timestamps()
        else:
            timestamps = self.asi8

        return tslib.ints_to_pydatetime(timestamps, box="time")

    @property
    def timetz(self):
        """
        Returns numpy array of datetime.time also containing timezone
        information. The time part of the Timestamps.
        """
        return tslib.ints_to_pydatetime(self.asi8, self.tz, box="time")

    @property
    def date(self):
        """
        Returns numpy array of python datetime.date objects (namely, the date
        part of Timestamps without timezone information).
        """
        # If the Timestamps have a timezone that is not UTC,
        # convert them into their i8 representation while
        # keeping their timezone and not using UTC
        if self.tz is not None and self.tz is not utc:
            timestamps = self._local_timestamps()
        else:
            timestamps = self.asi8

        return tslib.ints_to_pydatetime(timestamps, box="date")

    year = _field_accessor('year', 'Y', "The year of the datetime")
    month = _field_accessor('month', 'M',
                            "The month as January=1, December=12")
    day = _field_accessor('day', 'D', "The days of the datetime")
    hour = _field_accessor('hour', 'h', "The hours of the datetime")
    minute = _field_accessor('minute', 'm', "The minutes of the datetime")
    second = _field_accessor('second', 's', "The seconds of the datetime")
    microsecond = _field_accessor('microsecond', 'us',
                                  "The microseconds of the datetime")
    nanosecond = _field_accessor('nanosecond', 'ns',
                                 "The nanoseconds of the datetime")
    weekofyear = _field_accessor('weekofyear', 'woy',
                                 "The week ordinal of the year")
    week = weekofyear
    _dayofweek_doc = """
    The day of the week with Monday=0, Sunday=6.

    Return the day of the week. It is assumed the week starts on
    Monday, which is denoted by 0 and ends on Sunday which is denoted
    by 6. This method is available on both Series with datetime
    values (using the `dt` accessor) or DatetimeIndex.

    See Also
    --------
    Series.dt.dayofweek : Alias.
    Series.dt.weekday : Alias.
    Series.dt.day_name : Returns the name of the day of the week.

    Returns
    -------
    Series or Index
        Containing integers indicating the day number.

    Examples
    --------
    >>> s = pd.date_range('2016-12-31', '2017-01-08', freq='D').to_series()
    >>> s.dt.dayofweek
    2016-12-31    5
    2017-01-01    6
    2017-01-02    0
    2017-01-03    1
    2017-01-04    2
    2017-01-05    3
    2017-01-06    4
    2017-01-07    5
    2017-01-08    6
    Freq: D, dtype: int64
    """
    dayofweek = _field_accessor('dayofweek', 'dow', _dayofweek_doc)
    weekday = dayofweek

    weekday_name = _field_accessor(
        'weekday_name',
        'weekday_name',
        "The name of day in a week (ex: Friday)\n\n.. deprecated:: 0.23.0")

    dayofyear = _field_accessor('dayofyear', 'doy',
                                "The ordinal day of the year")
    quarter = _field_accessor('quarter', 'q', "The quarter of the date")
    days_in_month = _field_accessor(
        'days_in_month',
        'dim',
        "The number of days in the month")
    daysinmonth = days_in_month
    is_month_start = _field_accessor(
        'is_month_start',
        'is_month_start',
        "Logical indicating if first day of month (defined by frequency)")
    is_month_end = _field_accessor(
        'is_month_end',
        'is_month_end',
        """
        Indicator for whether the date is the last day of the month.

        Returns
        -------
        Series or array
            For Series, returns a Series with boolean values. For
            DatetimeIndex, returns a boolean array.

        See Also
        --------
        is_month_start : Indicator for whether the date is the first day
            of the month.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> dates = pd.Series(pd.date_range("2018-02-27", periods=3))
        >>> dates
        0   2018-02-27
        1   2018-02-28
        2   2018-03-01
        dtype: datetime64[ns]
        >>> dates.dt.is_month_end
        0    False
        1    True
        2    False
        dtype: bool

        >>> idx = pd.date_range("2018-02-27", periods=3)
        >>> idx.is_month_end
        array([False,  True, False], dtype=bool)
        """)
    is_quarter_start = _field_accessor(
        'is_quarter_start',
        'is_quarter_start',
        """
        Indicator for whether the date is the first day of a quarter.

        Returns
        -------
        is_quarter_start : Series or DatetimeIndex
            The same type as the original data with boolean values. Series will
            have the same name and index. DatetimeIndex will have the same
            name.

        See Also
        --------
        quarter : Return the quarter of the date.
        is_quarter_end : Similar property for indicating the quarter start.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> df = pd.DataFrame({'dates': pd.date_range("2017-03-30",
        ...                   periods=4)})
        >>> df.assign(quarter=df.dates.dt.quarter,
        ...           is_quarter_start=df.dates.dt.is_quarter_start)
               dates  quarter  is_quarter_start
        0 2017-03-30        1             False
        1 2017-03-31        1             False
        2 2017-04-01        2              True
        3 2017-04-02        2             False

        >>> idx = pd.date_range('2017-03-30', periods=4)
        >>> idx
        DatetimeIndex(['2017-03-30', '2017-03-31', '2017-04-01', '2017-04-02'],
                      dtype='datetime64[ns]', freq='D')

        >>> idx.is_quarter_start
        array([False, False,  True, False])
        """)
    is_quarter_end = _field_accessor(
        'is_quarter_end',
        'is_quarter_end',
        """
        Indicator for whether the date is the last day of a quarter.

        Returns
        -------
        is_quarter_end : Series or DatetimeIndex
            The same type as the original data with boolean values. Series will
            have the same name and index. DatetimeIndex will have the same
            name.

        See Also
        --------
        quarter : Return the quarter of the date.
        is_quarter_start : Similar property indicating the quarter start.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> df = pd.DataFrame({'dates': pd.date_range("2017-03-30",
        ...                    periods=4)})
        >>> df.assign(quarter=df.dates.dt.quarter,
        ...           is_quarter_end=df.dates.dt.is_quarter_end)
               dates  quarter    is_quarter_end
        0 2017-03-30        1             False
        1 2017-03-31        1              True
        2 2017-04-01        2             False
        3 2017-04-02        2             False

        >>> idx = pd.date_range('2017-03-30', periods=4)
        >>> idx
        DatetimeIndex(['2017-03-30', '2017-03-31', '2017-04-01', '2017-04-02'],
                      dtype='datetime64[ns]', freq='D')

        >>> idx.is_quarter_end
        array([False,  True, False, False])
        """)
    is_year_start = _field_accessor(
        'is_year_start',
        'is_year_start',
        """
        Indicate whether the date is the first day of a year.

        Returns
        -------
        Series or DatetimeIndex
            The same type as the original data with boolean values. Series will
            have the same name and index. DatetimeIndex will have the same
            name.

        See Also
        --------
        is_year_end : Similar property indicating the last day of the year.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> dates = pd.Series(pd.date_range("2017-12-30", periods=3))
        >>> dates
        0   2017-12-30
        1   2017-12-31
        2   2018-01-01
        dtype: datetime64[ns]

        >>> dates.dt.is_year_start
        0    False
        1    False
        2    True
        dtype: bool

        >>> idx = pd.date_range("2017-12-30", periods=3)
        >>> idx
        DatetimeIndex(['2017-12-30', '2017-12-31', '2018-01-01'],
                      dtype='datetime64[ns]', freq='D')

        >>> idx.is_year_start
        array([False, False,  True])
        """)
    is_year_end = _field_accessor(
        'is_year_end',
        'is_year_end',
        """
        Indicate whether the date is the last day of the year.

        Returns
        -------
        Series or DatetimeIndex
            The same type as the original data with boolean values. Series will
            have the same name and index. DatetimeIndex will have the same
            name.

        See Also
        --------
        is_year_start : Similar property indicating the start of the year.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> dates = pd.Series(pd.date_range("2017-12-30", periods=3))
        >>> dates
        0   2017-12-30
        1   2017-12-31
        2   2018-01-01
        dtype: datetime64[ns]

        >>> dates.dt.is_year_end
        0    False
        1     True
        2    False
        dtype: bool

        >>> idx = pd.date_range("2017-12-30", periods=3)
        >>> idx
        DatetimeIndex(['2017-12-30', '2017-12-31', '2018-01-01'],
                      dtype='datetime64[ns]', freq='D')

        >>> idx.is_year_end
        array([False,  True, False])
        """)
    is_leap_year = _field_accessor(
        'is_leap_year',
        'is_leap_year',
        """
        Boolean indicator if the date belongs to a leap year.

        A leap year is a year, which has 366 days (instead of 365) including
        29th of February as an intercalary day.
        Leap years are years which are multiples of four with the exception
        of years divisible by 100 but not by 400.

        Returns
        -------
        Series or ndarray
             Booleans indicating if dates belong to a leap year.

        Examples
        --------
        This method is available on Series with datetime values under
        the ``.dt`` accessor, and directly on DatetimeIndex.

        >>> idx = pd.date_range("2012-01-01", "2015-01-01", freq="Y")
        >>> idx
        DatetimeIndex(['2012-12-31', '2013-12-31', '2014-12-31'],
                      dtype='datetime64[ns]', freq='A-DEC')
        >>> idx.is_leap_year
        array([ True, False, False], dtype=bool)

        >>> dates = pd.Series(idx)
        >>> dates_series
        0   2012-12-31
        1   2013-12-31
        2   2014-12-31
        dtype: datetime64[ns]
        >>> dates_series.dt.is_leap_year
        0     True
        1    False
        2    False
        dtype: bool
        """)

    def to_julian_date(self):
        """
        Convert Datetime Array to float64 ndarray of Julian Dates.
        0 Julian date is noon January 1, 4713 BC.
        http://en.wikipedia.org/wiki/Julian_day
        """

        # http://mysite.verizon.net/aesir_research/date/jdalg2.htm
        year = np.asarray(self.year)
        month = np.asarray(self.month)
        day = np.asarray(self.day)
        testarr = month < 3
        year[testarr] -= 1
        month[testarr] += 12
        return (day +
                np.fix((153 * month - 457) / 5) +
                365 * year +
                np.floor(year / 4) -
                np.floor(year / 100) +
                np.floor(year / 400) +
                1721118.5 +
                (self.hour +
                 self.minute / 60.0 +
                 self.second / 3600.0 +
                 self.microsecond / 3600.0 / 1e+6 +
                 self.nanosecond / 3600.0 / 1e+9
                 ) / 24.0)


DatetimeArrayMixin._add_comparison_ops()
DatetimeArrayMixin._add_datetimelike_methods()


def _generate_regular_range(cls, start, end, periods, freq):
    if isinstance(freq, Tick):
        stride = freq.nanos
        if periods is None:
            b = Timestamp(start).value
            # cannot just use e = Timestamp(end) + 1 because arange breaks when
            # stride is too large, see GH10887
            e = (b + (Timestamp(end).value - b) // stride * stride +
                 stride // 2 + 1)
            # end.tz == start.tz by this point due to _generate implementation
            tz = start.tz
        elif start is not None:
            b = Timestamp(start).value
            e = b + np.int64(periods) * stride
            tz = start.tz
        elif end is not None:
            e = Timestamp(end).value + stride
            b = e - np.int64(periods) * stride
            tz = end.tz
        else:
            raise ValueError("at least 'start' or 'end' should be specified "
                             "if a 'period' is given.")

        data = np.arange(b, e, stride, dtype=np.int64)
        data = cls._simple_new(data.view(_NS_DTYPE), None, tz=tz)
    else:
        tz = None
        # start and end should have the same timezone by this point
        if isinstance(start, Timestamp):
            tz = start.tz
        elif isinstance(end, Timestamp):
            tz = end.tz

        xdr = generate_range(start=start, end=end,
                             periods=periods, offset=freq)

        values = np.array([x.value for x in xdr])
        data = cls._simple_new(values, freq=freq, tz=tz)

    return data


def _infer_tz_from_endpoints(start, end, tz):
    """
    If a timezone is not explicitly given via `tz`, see if one can
    be inferred from the `start` and `end` endpoints.  If more than one
    of these inputs provides a timezone, require that they all agree.

    Parameters
    ----------
    start : Timestamp
    end : Timestamp
    tz : tzinfo or None

    Returns
    -------
    tz : tzinfo or None
    inferred_tz : tzinfo or None

    Raises
    ------
    TypeError : if start and end timezones do not agree
    """
    try:
        inferred_tz = timezones.infer_tzinfo(start, end)
    except Exception:
        raise TypeError('Start and end cannot both be tz-aware with '
                        'different timezones')

    inferred_tz = timezones.maybe_get_tz(inferred_tz)
    tz = timezones.maybe_get_tz(tz)

    if tz is not None and inferred_tz is not None:
        if not timezones.tz_compare(inferred_tz, tz):
            raise AssertionError("Inferred time zone not equal to passed "
                                 "time zone")

    elif inferred_tz is not None:
        tz = inferred_tz

    return tz, inferred_tz


def _maybe_normalize_endpoints(start, end, normalize):
    _normalized = True

    if start is not None:
        if normalize:
            start = normalize_date(start)
            _normalized = True
        else:
            _normalized = _normalized and start.time() == _midnight

    if end is not None:
        if normalize:
            end = normalize_date(end)
            _normalized = True
        else:
            _normalized = _normalized and end.time() == _midnight

    return start, end, _normalized


def _maybe_localize_point(ts, is_none, is_not_none, freq, tz):
    """
    Localize a start or end Timestamp to the timezone of the corresponding
    start or end Timestamp

    Parameters
    ----------
    ts : start or end Timestamp to potentially localize
    is_none : argument that should be None
    is_not_none : argument that should not be None
    freq : Tick, DateOffset, or None
    tz : str, timezone object or None

    Returns
    -------
    ts : Timestamp
    """
    # Make sure start and end are timezone localized if:
    # 1) freq = a Timedelta-like frequency (Tick)
    # 2) freq = None i.e. generating a linspaced range
    if isinstance(freq, Tick) or freq is None:
        localize_args = {'tz': tz, 'ambiguous': False}
    else:
        localize_args = {'tz': None}
    if is_none is None and is_not_none is not None:
        ts = ts.tz_localize(**localize_args)
    return ts
