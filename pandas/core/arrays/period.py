# -*- coding: utf-8 -*-
from datetime import timedelta
import operator

import numpy as np

from pandas import compat
from pandas.compat.numpy import function as nv
from pandas._libs import lib
from pandas._libs.tslib import NaT, iNaT
from pandas._libs.tslibs.period import (
    Period, IncompatibleFrequency, DIFFERENT_FREQ_INDEX,
    get_period_field_arr, period_asfreq_arr,
)
from pandas._libs.tslibs import period as libperiod
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds, Timedelta
from pandas._libs.tslibs.fields import isleapyear_arr
from pandas.util._decorators import cache_readonly
from pandas.core.dtypes.common import (
    is_integer_dtype, is_float_dtype, is_period_dtype,
    is_float, is_integer, pandas_dtype, is_scalar,
    is_datetime64_dtype,
    is_categorical_dtype,
    is_timedelta64_dtype,
    is_object_dtype,
    is_string_dtype,
    is_datetime_or_timedelta_dtype,
    is_dtype_equal,
    ensure_object,
    _TD_DTYPE,
)

from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.dtypes.generic import (
    ABCSeries, ABCIndexClass,
)

import pandas.core.common as com

from pandas.tseries import frequencies
from pandas.tseries.frequencies import get_freq_code as _gfc
from pandas.tseries.offsets import Tick, DateOffset

from pandas.core.arrays import ExtensionArray
from pandas.core.arrays import datetimelike as dtl


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = frequencies.get_freq_code(self.freq)
        result = get_period_field_arr(alias, self._ndarray_values, base)
        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _period_array_cmp(cls, op):
    """
    Wrap comparison operations to convert Period-like to PeriodDtype
    """
    opname = '__{name}__'.format(name=op.__name__)
    nat_result = True if opname == '__ne__' else False

    def wrapper(self, other):
        op = getattr(self._ndarray_values, opname)
        if isinstance(other, (ABCSeries, ABCIndexClass)):
            other = other.values

        if isinstance(other, Period):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            result = op(other.ordinal)
        elif isinstance(other, cls):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            result = op(other._ndarray_values)

            mask = self._isnan | other._isnan
            if mask.any():
                result[mask] = nat_result

            return result
        elif other is NaT:
            result = np.empty(len(self._ndarray_values), dtype=bool)
            result.fill(nat_result)
        elif isinstance(other, (list, np.ndarray)):
            # XXX: is this correct? Why not convert the
            # sequence to a PeriodArray?
            return NotImplemented
        else:
            other = Period(other, freq=self.freq)
            result = op(other.ordinal)

        if self.hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


class PeriodArray(dtl.DatetimeLikeArrayMixin, ExtensionArray):
    """
    Pandas ExtensionArray for storing Period data.

    Users should use the :func:`period_array` function to create
    new instances of PeriodArray.

    Notes
    -----
    There are two components to a PeriodArray

    - ordinals : integer ndarray
    - freq : pd.tseries.offsets.Tick

    The values are physically stored as a 1-D ndarray of integers. These are
    called "ordinals" and represent some kind of offset from a base.

    The `freq` indicates the span covered by each element of the array.
    All elements in the PeriodArray have the same `freq`.

    See Also
    --------
    period_array : Create a new PeriodArray
    pandas.PeriodIndex : Immutable Index for period data
    """
    _attributes = ["freq"]
    _typ = "periodarray"  # ABCPeriodArray

    # Names others delegate to us
    _other_ops = []
    _bool_ops = ['is_leap_year']
    _object_ops = ['start_time', 'end_time', 'freq']
    _field_ops = ['year', 'month', 'day', 'hour', 'minute', 'second',
                  'weekofyear', 'weekday', 'week', 'dayofweek',
                  'dayofyear', 'quarter', 'qyear',
                  'days_in_month', 'daysinmonth']
    _datetimelike_ops = _field_ops + _object_ops + _bool_ops
    _datetimelike_methods = ['strftime', 'to_timestamp', 'asfreq']

    # --------------------------------------------------------------------
    # Constructors
    def __init__(self, values, freq=None):
        # type: (Union[PeriodArray, np.ndarray], Union[str, Tick]) -> None
        if isinstance(values, type(self)):
            values, freq = values._data, values.freq

        values = np.array(values, dtype='int64', copy=False)
        self._data = values
        if freq is None:
            raise ValueError('freq is not specified and cannot be inferred')
        freq = Period._maybe_convert_freq(freq)
        self._dtype = PeriodDtype(freq)

    @classmethod
    def _complex_new(cls, data=None, ordinal=None, freq=None, start=None,
                     end=None, periods=None, tz=None, dtype=None, copy=False,
                     **fields):
        from pandas import PeriodIndex, DatetimeIndex, Int64Index

        # copy-pase from PeriodIndex.__new__ with slight adjustments.
        #
        # - removed all uses of name
        # - refactored to smaller, more dedicated constructors.

        # TODO: move fields validation to range init
        valid_field_set = {'year', 'month', 'day', 'quarter',
                           'hour', 'minute', 'second'}

        if not set(fields).issubset(valid_field_set):
            raise TypeError('__new__() got an unexpected keyword argument {}'.
                            format(list(set(fields) - valid_field_set)[0]))

        if periods is not None:
            if is_float(periods):
                periods = int(periods)
            elif not is_integer(periods):
                msg = 'periods must be a number, got {periods}'
                raise TypeError(msg.format(periods=periods))

            periods = dtl.validate_periods(periods)

        if dtype is not None:
            dtype = pandas_dtype(dtype)
            if not is_period_dtype(dtype):
                raise ValueError('dtype must be PeriodDtype')
            if freq is None:
                freq = dtype.freq
            elif freq != dtype.freq:
                msg = 'specified freq and dtype are different'
                raise IncompatibleFrequency(msg)

        # coerce freq to freq object, otherwise it can be coerced elementwise
        # which is slow
        if freq:
            freq = Period._maybe_convert_freq(freq)

        if data is None:
            if ordinal is not None:
                data = np.asarray(ordinal, dtype=np.int64)
            else:
                data, freq = cls._generate_range(start, end, periods,
                                                 freq, fields)
            return cls._from_ordinals(data, freq=freq)

        if isinstance(data, (cls, PeriodIndex)):
            if freq is None or freq == data.freq:  # no freq change
                freq = data.freq
                data = data._ndarray_values
            else:
                base1, _ = _gfc(data.freq)
                base2, _ = _gfc(freq)
                data = libperiod.period_asfreq_arr(data._ndarray_values,
                                                   base1, base2, 1)
            if copy:
                data = data.copy()
            return cls._simple_new(data, freq=freq)

        # not array / index
        if not isinstance(data, (np.ndarray, PeriodIndex,
                                 DatetimeIndex, Int64Index)):
            if is_scalar(data):
                raise TypeError('{0}(...) must be called with a '
                                'collection of some '
                                'kind, {1} was passed'.format(cls.__name__,
                                                              repr(data)))

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            data = np.asarray(data)

        # datetime other than period
        if is_datetime64_dtype(data.dtype):
            data = dt64arr_to_periodarr(data, freq, tz)
            return cls._from_ordinals(data, freq=freq)

        # check not floats
        if lib.infer_dtype(data) == 'floating' and len(data) > 0:
            raise TypeError("PeriodIndex does not allow "
                            "floating point in construction")

        # anything else, likely an array of strings or periods
        data = ensure_object(data)
        return cls._from_periods(data, freq=freq)

    @classmethod
    def _simple_new(cls, values, freq=None, **kwargs):
        """
        Values can be any type that can be coerced to Periods.
        Ordinals in an ndarray are fastpath-ed to `_from_ordinals`
        """
        if not is_integer_dtype(values):
            values = np.array(values, copy=False)
            if len(values) > 0 and is_float_dtype(values):
                raise TypeError("{cls} can't take floats"
                                .format(cls=cls.__name__))
            return cls(values, freq=freq, **kwargs)

        return cls(values, freq=freq)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        # type: (Sequence[Optional[Period]], Dtype, bool) -> PeriodArray
        if dtype:
            freq = dtype.freq
        else:
            freq = None
        scalars = np.asarray(scalars, dtype=object)
        return cls._from_periods(scalars, freq=freq)

    def _values_for_factorize(self):
        return self.values, iNaT

    @classmethod
    def _from_factorized(cls, values, original):
        # type: (Sequence[Optional[Period]], PeriodArray) -> PeriodArray
        return cls._simple_new(values, freq=original.freq)

    @classmethod
    def _from_ordinals(cls, values, freq=None):
        # type: (ndarray[int], Optional[Tick]) -> PeriodArray
        """
        Values should be int ordinals
        `__new__` & `_simple_new` coerce to ordinals and call this method
        """
        return cls(values, freq=freq)

    @classmethod
    def _from_periods(cls, periods, freq=None):
        # type: (np.ndarray[Optional[Period]], Optional[Tick]) -> PeriodArray
        periods = np.asarray(periods, dtype=object)
        freq = freq or libperiod.extract_freq(periods)
        ordinals = libperiod.extract_ordinals(periods, freq)
        return cls._from_ordinals(ordinals, freq=freq)

    @classmethod
    def _from_datetime64(cls, data, freq, tz=None):
        """Construct a PeriodArray from a datetime64 array

        Parameters
        ----------
        data : ndarray[datetime64[ns], datetime64[ns, tz]]
        freq : str or Tick
        tz : tzinfo, optional

        Returns
        -------
        PeriodArray[freq]
        """
        data = dt64arr_to_periodarr(data, freq, tz)
        return cls._simple_new(data, freq=freq)

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        if freq is not None:
            freq = Period._maybe_convert_freq(freq)

        field_count = len(fields)
        if com.count_not_none(start, end) > 0:
            if field_count > 0:
                raise ValueError('Can either instantiate from fields '
                                 'or endpoints, but not both')
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
        elif field_count > 0:
            subarr, freq = _range_from_fields(freq=freq, **fields)
        else:
            raise ValueError('Not enough parameters to construct '
                             'Period range')

        return subarr, freq

    @classmethod
    def _concat_same_type(cls, to_concat):
        freq = {x.freq for x in to_concat}
        assert len(freq) == 1
        freq = list(freq)[0]
        values = np.concatenate([x._data for x in to_concat])
        return cls._from_ordinals(values, freq=freq)

    @property
    def asi8(self):
        return self._ndarray_values.view('i8')

    # --------------------------------------------------------------------
    # Data / Attributes
    @property
    def nbytes(self):
        # TODO(DatetimeArray): remove
        return self._data.nbytes

    @cache_readonly
    def dtype(self):
        return self._dtype

    @property
    def _ndarray_values(self):
        # Ordinals
        return self._data

    @property
    def freq(self):
        """Return the frequency object for this PeriodArray."""
        return self.dtype.freq

    @property
    def flags(self):
        """Deprecated"""
        # Just here to support Index.flags deprecation.
        # could also override PeriodIndex.flags if we don't want a
        # version with PeriodArray.flags
        return self.values.flags

    @property
    def base(self):
        return self.values.base

    @property
    def data(self):
        return self.astype(object).data

    # --------------------------------------------------------------------
    # Vectorized analogues of Period properties

    year = _field_accessor('year', 0, "The year of the period")
    month = _field_accessor('month', 3, "The month as January=1, December=12")
    day = _field_accessor('day', 4, "The days of the period")
    hour = _field_accessor('hour', 5, "The hour of the period")
    minute = _field_accessor('minute', 6, "The minute of the period")
    second = _field_accessor('second', 7, "The second of the period")
    weekofyear = _field_accessor('week', 8, "The week ordinal of the year")
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 10,
                                "The day of the week with Monday=0, Sunday=6")
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor('dayofyear', 9,
                                              "The ordinal day of the year")
    quarter = _field_accessor('quarter', 2, "The quarter of the date")
    qyear = _field_accessor('qyear', 1)
    days_in_month = _field_accessor('days_in_month', 11,
                                    "The number of days in the month")
    daysinmonth = days_in_month

    @property
    def is_leap_year(self):
        """ Logical indicating if the date belongs to a leap year """
        return isleapyear_arr(np.asarray(self.year))

    @property
    def start_time(self):
        return self.to_timestamp(how='start')

    @property
    def end_time(self):
        return self.to_timestamp(how='end')

    def __repr__(self):
        return '<{}>\n{}\nLength: {}, dtype: {}'.format(
            self.__class__.__name__,
            [str(s) for s in self],
            len(self),
            self.dtype
        )

    def __len__(self):
        return len(self._data)

    def __setitem__(self, key, value):
        from pandas.core.dtypes.missing import isna

        if isinstance(value, (compat.Sequence, type(self))):
            if len(key) != len(value) and not com.is_bool_indexer(key):
                msg = ("shape mismatch: value array of length '{}' does not "
                       "match indexing result of length '{}'.")
                raise ValueError(msg.format(len(key), len(value)))
            if len(key) == 0:
                return

            value = type(self)._complex_new(value)
            if self.freqstr != value.freqstr:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, value.freqstr)
                raise IncompatibleFrequency(msg)

            value = value.asi8
        elif isinstance(value, Period):

            if self.freqstr != value.freqstr:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, value.freqstr)
                raise IncompatibleFrequency(msg)

            value = value.ordinal
        elif isna(value):
            # Previously we allowed setting np.nan on a Series[object]
            # do we still want to allow that, or should we require None / NaT?
            value = iNaT
        else:
            msg = ("'value' should be a 'Period', 'NaT', or array of those. "
                   "Got '{}' instead.".format(type(value).__name__))
            raise TypeError(msg)
        self._data[key] = value

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take
        from pandas import isna

        if allow_fill:
            if isna(fill_value):
                fill_value = iNaT
            elif isinstance(fill_value, Period):
                fill_value = fill_value.ordinal
            else:
                msg = "'fill_value' should be a Period. Got '{}'."
                raise ValueError(msg.format(fill_value))

        new_values = take(self._data,
                          indices,
                          allow_fill=allow_fill,
                          fill_value=fill_value)

        return self._from_ordinals(new_values, self.freq)

    def isna(self):
        return self._data == iNaT

    def fillna(self, value=None, method=None, limit=None):
        # TODO(#20300)
        # To avoid converting to object, we re-implement here with the changes
        # 1. Passing `_ndarray_values` to func instead of self.astype(object)
        # 2. Re-boxing with `_from_ordinals`
        # #20300 should let us do this kind of logic on ExtensionArray.fillna
        # and we can use it.
        from pandas.api.types import is_array_like
        from pandas.util._validators import validate_fillna_kwargs
        from pandas.core.missing import pad_1d, backfill_1d

        if isinstance(value, ABCSeries):
            value = value.values

        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError("Length of 'value' does not match. Got ({}) "
                                 " expected {}".format(len(value), len(self)))
            value = value[mask]

        if mask.any():
            if method is not None:
                func = pad_1d if method == 'pad' else backfill_1d
                new_values = func(self._ndarray_values, limit=limit,
                                  mask=mask)
                new_values = self._from_ordinals(new_values, freq=self.freq)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def copy(self, deep=False):
        return self._from_ordinals(self._data.copy(), freq=self.freq)

    def value_counts(self, dropna=False):
        from pandas.core.algorithms import value_counts
        from pandas.core.indexes.period import PeriodIndex

        if dropna:
            values = self[~self.isna()]._data
        else:
            values = self._data

        result = value_counts(values, sort=False)
        index = PeriodIndex._from_ordinals(result.index,
                                           name=result.index.name,
                                           freq=self.freq)
        return type(result)(result.values,
                            index=index,
                            name=result.name)

    def shift(self, periods=1):
        """
        Shift values by desired number.

        Newly introduced missing values are filled with
        ``self.dtype.na_value``.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        periods : int, default 1
            The number of periods to shift. Negative values are allowed
            for shifting backwards.

        Returns
        -------
        shifted : PeriodArray
        """
        # TODO(DatetimeArray): remove
        # The semantics for Index.shift differ from EA.shift
        # then just call super.
        return ExtensionArray.shift(self, periods)

    def _time_shift(self, n, freq=None):
        """
        Shift each value by `periods`.

        Note this is different from ExtensionArray.shift, which
        shifts the *position* of each element, padding the end with
        missing values.

        Parameters
        ----------
        periods : int
            Number of periods to shift by.
        freq : pandas.DateOffset, pandas.Timedelta, or string
            Frequency increment to shift by.
        """
        values = self.values + n * self.freq.n
        if self.hasnans:
            values[self._isnan] = iNaT
        return self._simple_new(values, freq=self.freq)

    @property
    def _box_func(self):
        # Used in DatelikeArray.__iter__
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    def asfreq(self, freq=None, how='E'):
        """
        Convert the Period Array/Index to the specified frequency `freq`.

        Parameters
        ----------
        freq : str
            a frequency
        how : str {'E', 'S'}
            'E', 'END', or 'FINISH' for end,
            'S', 'START', or 'BEGIN' for start.
            Whether the elements should be aligned to the end
            or start within pa period. January 31st ('END') vs.
            January 1st ('START') for example.

        Returns
        -------
        new : Period Array/Index with the new frequency

        Examples
        --------
        >>> pidx = pd.period_range('2010-01-01', '2015-01-01', freq='A')
        >>> pidx
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010, ..., 2015]
        Length: 6, Freq: A-DEC

        >>> pidx.asfreq('M')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-12, ..., 2015-12]
        Length: 6, Freq: M

        >>> pidx.asfreq('M', how='S')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-01, ..., 2015-01]
        Length: 6, Freq: M
        """
        how = libperiod._validate_end_alias(how)

        freq = Period._maybe_convert_freq(freq)

        base1, mult1 = frequencies.get_freq_code(self.freq)
        base2, mult2 = frequencies.get_freq_code(freq)

        asi8 = self.asi8
        # mult1 can't be negative or 0
        end = how == 'E'
        if end:
            ordinal = asi8 + mult1 - 1
        else:
            ordinal = asi8

        new_data = period_asfreq_arr(ordinal, base1, base2, end)

        if self.hasnans:
            new_data[self._isnan] = iNaT

        return self._shallow_copy(new_data, freq=freq)

    def to_timestamp(self, freq=None, how='start'):
        """
        Cast to DatetimeArray/Index

        Parameters
        ----------
        freq : string or DateOffset, optional
            Target frequency. The default is 'D' for week or longer,
            'S' otherwise
        how : {'s', 'e', 'start', 'end'}

        Returns
        -------
        DatetimeArray/Index
        """
        from pandas.core.arrays.datetimes import DatetimeArrayMixin

        how = libperiod._validate_end_alias(how)

        end = how == 'E'
        if end:
            if freq == 'B':
                # roll forward to ensure we land on B date
                adjust = Timedelta(1, 'D') - Timedelta(1, 'ns')
                return self.to_timestamp(how='start') + adjust
            else:
                adjust = Timedelta(1, 'ns')
                return (self + 1).to_timestamp(how='start') - adjust

        if freq is None:
            base, mult = frequencies.get_freq_code(self.freq)
            freq = frequencies.get_to_timestamp_base(base)
        else:
            freq = Period._maybe_convert_freq(freq)

        base, mult = frequencies.get_freq_code(freq)
        new_data = self.asfreq(freq, how=how)

        new_data = libperiod.periodarr_to_dt64arr(new_data._ndarray_values,
                                                  base)
        return DatetimeArrayMixin(new_data, freq='infer')

    # ------------------------------------------------------------------
    # Arithmetic Methods
    _create_comparison_method = classmethod(_period_array_cmp)

    def _sub_datelike(self, other):
        assert other is not NaT
        return NotImplemented

    def _sub_period(self, other):
        # If the operation is well-defined, we return an object-Index
        # of DateOffsets.  Null entries are filled with pd.NaT
        if self.freq != other.freq:
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)

        asi8 = self.asi8
        new_data = asi8 - other.ordinal
        new_data = np.array([self.freq * x for x in new_data])

        if self.hasnans:
            new_data[self._isnan] = NaT

        return new_data

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        base = frequencies.get_base_alias(other.rule_code)
        if base != self.freq.rule_code:
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        return self._time_shift(other.n)

    def _add_delta_td(self, other):
        assert isinstance(self.freq, Tick)  # checked by calling function
        assert isinstance(other, (timedelta, np.timedelta64, Tick))

        delta = self._check_timedeltalike_freq_compat(other)

        # Note: when calling parent class's _add_delta_td, it will call
        #  delta_to_nanoseconds(delta).  Because delta here is an integer,
        #  delta_to_nanoseconds will return it unchanged.
        return type(self)._add_delta_td(self, delta)

    def _add_delta_tdi(self, other):
        assert isinstance(self.freq, Tick)  # checked by calling function

        delta = self._check_timedeltalike_freq_compat(other)
        return self._addsub_int_array(delta, operator.add)

    def _add_delta(self, other):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self.

        Parameters
        ----------
        other : {timedelta, np.timedelta64, Tick,
                 TimedeltaIndex, ndarray[timedelta64]}

        Returns
        -------
        result : same type as self
        """
        if not isinstance(self.freq, Tick):
            # We cannot add timedelta-like to non-tick PeriodArray
            raise IncompatibleFrequency("Input has different freq from "
                                        "{cls}(freq={freqstr})"
                                        .format(cls=type(self).__name__,
                                                freqstr=self.freqstr))

        # TODO: standardize across datetimelike subclasses whether to return
        #  i8 view or _shallow_copy
        if isinstance(other, (Tick, timedelta, np.timedelta64)):
            new_values = self._add_delta_td(other)
            return self._shallow_copy(new_values)
        elif is_timedelta64_dtype(other):
            # ndarray[timedelta64] or TimedeltaArray/index
            new_values = self._add_delta_tdi(other)
            return self._shallow_copy(new_values)
        else:  # pragma: no cover
            raise TypeError(type(other).__name__)

    def _maybe_convert_timedelta(self, other):
        """
        Convert timedelta-like input to an integer multiple of self.freq

        Parameters
        ----------
        other : timedelta, np.timedelta64, DateOffset, int, np.ndarray

        Returns
        -------
        converted : int, np.ndarray[int64]

        Raises
        ------
        IncompatibleFrequency : if the input cannot be written as a multiple
            of self.freq.  Note IncompatibleFrequency subclasses ValueError.
        """
        if isinstance(
                other, (timedelta, np.timedelta64, Tick, np.ndarray)):
            offset = frequencies.to_offset(self.freq.rule_code)
            if isinstance(offset, Tick):
                # _check_timedeltalike_freq_compat will raise if incompatible
                delta = self._check_timedeltalike_freq_compat(other)
                return delta
        elif isinstance(other, DateOffset):
            freqstr = other.rule_code
            base = frequencies.get_base_alias(freqstr)
            if base == self.freq.rule_code:
                return other.n
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        elif lib.is_integer(other):
            # integer is passed to .shift via
            # _add_datetimelike_methods basically
            # but ufunc may pass integer to _add_delta
            return other

        # raise when input doesn't have freq
        msg = "Input has different freq from {cls}(freq={freqstr})"
        raise IncompatibleFrequency(msg.format(cls=type(self).__name__,
                                               freqstr=self.freqstr))

    # ------------------------------------------------------------------
    # Formatting
    def _format_native_types(self, na_rep=u'NaT', date_format=None,
                             **kwargs):
        # TODO(DatetimeArray): remove
        values = self.astype(object)

        if date_format:
            formatter = lambda dt: dt.strftime(date_format)
        else:
            formatter = lambda dt: u'%s' % dt

        if self.hasnans:
            mask = self._isnan
            values[mask] = na_rep
            imask = ~mask
            values[imask] = np.array([formatter(dt) for dt
                                      in values[imask]])
        else:
            values = np.array([formatter(dt) for dt in values])
        return values

    def _check_timedeltalike_freq_compat(self, other):
        """
        Arithmetic operations with timedelta-like scalars or array `other`
        are only valid if `other` is an integer multiple of `self.freq`.
        If the operation is valid, find that integer multiple.  Otherwise,
        raise because the operation is invalid.

        Parameters
        ----------
        other : timedelta, np.timedelta64, Tick,
                ndarray[timedelta64], TimedeltaArray, TimedeltaIndex

        Returns
        -------
        multiple : int or ndarray[int64]

        Raises
        ------
        IncompatibleFrequency
        """
        assert isinstance(self.freq, Tick)  # checked by calling function
        own_offset = frequencies.to_offset(self.freq.rule_code)
        base_nanos = delta_to_nanoseconds(own_offset)

        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            nanos = delta_to_nanoseconds(other)

        elif isinstance(other, np.ndarray):
            # numpy timedelta64 array; all entries must be compatible
            assert other.dtype.kind == 'm'
            if other.dtype != _TD_DTYPE:
                # i.e. non-nano unit
                # TODO: disallow unit-less timedelta64
                other = other.astype(_TD_DTYPE)
            nanos = other.view('i8')
        else:
            # TimedeltaArray/Index
            nanos = other.asi8

        if np.all(nanos % base_nanos == 0):
            # nanos being added is an integer multiple of the
            #  base-frequency to self.freq
            delta = nanos // base_nanos
            # delta is the integer (or integer-array) number of periods
            # by which will be added to self.
            return delta

        raise IncompatibleFrequency("Input has different freq from "
                                    "{cls}(freq={freqstr})"
                                    .format(cls=type(self).__name__,
                                            freqstr=self.freqstr))

    def repeat(self, repeats, *args, **kwargs):
        """
        Repeat elements of a Categorical.

        See also
        --------
        numpy.ndarray.repeat
        """
        # TODO: Share with Categorical.repeat?
        # need to use ndarray_values in Categorical
        # and some kind of _constructor (from_ordinals, from_codes).
        nv.validate_repeat(args, kwargs)
        values = self._ndarray_values.repeat(repeats)
        return self._from_ordinals(values, self.freq)

    # Delegation...
    def strftime(self, date_format):
        return self._format_native_types(date_format=date_format)

    def astype(self, dtype, copy=True):
        # TODO: Figure out something better here...
        # We have DatetimeLikeArrayMixin ->
        #     super(...), which ends up being... DatetimeIndexOpsMixin?
        # this is complicated.
        # need a pandas_astype(arr, dtype).
        from pandas import Categorical

        dtype = pandas_dtype(dtype)

        if is_object_dtype(dtype):
            return np.asarray(self, dtype=object)
        elif is_string_dtype(dtype) and not is_categorical_dtype(dtype):
            return self._format_native_types()
        elif is_integer_dtype(dtype):
            return self.values.astype("i8", copy=copy)
        elif (is_datetime_or_timedelta_dtype(dtype) and
              not is_dtype_equal(self.dtype, dtype)) or is_float_dtype(dtype):
            # disallow conversion between datetime/timedelta,
            # and conversions for any datetimelike to float
            msg = 'Cannot cast {name} to dtype {dtype}'
            raise TypeError(msg.format(name=type(self).__name__, dtype=dtype))
        elif is_categorical_dtype(dtype):
            return Categorical(self, dtype=dtype)
        elif is_period_dtype(dtype):
            return self.asfreq(dtype.freq)
        else:
            return np.asarray(self, dtype=dtype)

    def item(self):
        if len(self) == 1:
            return Period._from_ordinal(self.values[0], self.freq)
        else:
            raise ValueError('can only convert an array of size 1 to a '
                             'Python scalar')


PeriodArray._add_comparison_ops()
PeriodArray._add_datetimelike_methods()


# -------------------------------------------------------------------
# Constructor Helpers

def period_array(data, freq=None):
    # type: (Sequence[Optional[Period]], Optional[Tick]) -> PeriodArray
    """
    Construct a new PeriodArray from a sequence of Period scalars.

    Parameters
    ----------
    data : Sequence of Period objects
        A sequence of Period objects. These are required to all have
        the same ``freq.`` Missing values can be indicated by ``None``
        or ``pandas.NaT``.
    freq : str, Tick, or Offset
        The frequency of every element of the array. This can be specified
        to avoid inferring the `freq` from `data`.

    Returns
    -------
    PeriodArray

    See Also
    --------
    PeriodArray
    pandas.PeriodIndex

    Examples
    --------
    >>> period_array([pd.Period('2017', freq='A'),
    ...               pd.Period('2018', freq='A')])
    <PeriodArray>
    ['2017', '2018']
    Length: 2, dtype: period[A-DEC]

    >>> period_array([pd.Period('2017', freq='A'),
    ...               pd.Period('2018', freq='A'),
    ...               pd.NaT])
    <PeriodArray>
    ['2017', '2018', 'NaT']
    Length: 3, dtype: period[A-DEC]
    """
    return PeriodArray._from_periods(data, freq=freq)


def dt64arr_to_periodarr(data, freq, tz=None):
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    freq = Period._maybe_convert_freq(freq)
    base, mult = frequencies.get_freq_code(freq)
    return libperiod.dt64arr_to_periodarr(data.view('i8'), base, tz)


def _get_ordinal_range(start, end, periods, freq, mult=1):
    if com.count_not_none(start, end, periods) != 2:
        raise ValueError('Of the three parameters: start, end, and periods, '
                         'exactly two must be specified')

    if freq is not None:
        _, mult = frequencies.get_freq_code(freq)

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError('start and end must have same freq')
    if (start is NaT or end is NaT):
        raise ValueError('start and end must not be NaT')

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:  # pragma: no cover
            raise ValueError('Could not infer freq from start/end')

    if periods is not None:
        periods = periods * mult
        if start is None:
            data = np.arange(end.ordinal - periods + mult,
                             end.ordinal + 1, mult,
                             dtype=np.int64)
        else:
            data = np.arange(start.ordinal, start.ordinal + periods, mult,
                             dtype=np.int64)
    else:
        data = np.arange(start.ordinal, end.ordinal + 1, mult, dtype=np.int64)

    return data, freq


def _range_from_fields(year=None, month=None, quarter=None, day=None,
                       hour=None, minute=None, second=None, freq=None):
    if hour is None:
        hour = 0
    if minute is None:
        minute = 0
    if second is None:
        second = 0
    if day is None:
        day = 1

    ordinals = []

    if quarter is not None:
        if freq is None:
            freq = 'Q'
            base = frequencies.FreqGroup.FR_QTR
        else:
            base, mult = frequencies.get_freq_code(freq)
            if base != frequencies.FreqGroup.FR_QTR:
                raise AssertionError("base must equal FR_QTR")

        year, quarter = _make_field_arrays(year, quarter)
        for y, q in compat.zip(year, quarter):
            y, m = libperiod.quarter_to_myear(y, q, freq)
            val = libperiod.period_ordinal(y, m, 1, 1, 1, 1, 0, 0, base)
            ordinals.append(val)
    else:
        base, mult = frequencies.get_freq_code(freq)
        arrays = _make_field_arrays(year, month, day, hour, minute, second)
        for y, mth, d, h, mn, s in compat.zip(*arrays):
            ordinals.append(libperiod.period_ordinal(
                y, mth, d, h, mn, s, 0, 0, base))

    return np.array(ordinals, dtype=np.int64), freq


def _make_field_arrays(*fields):
    length = None
    for x in fields:
        if isinstance(x, (list, np.ndarray, ABCSeries)):
            if length is not None and len(x) != length:
                raise ValueError('Mismatched Period array lengths')
            elif length is None:
                length = len(x)

    arrays = [np.asarray(x) if isinstance(x, (np.ndarray, list, ABCSeries))
              else np.repeat(x, length) for x in fields]

    return arrays
