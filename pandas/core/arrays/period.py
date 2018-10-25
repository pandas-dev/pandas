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
from pandas.util._validators import validate_fillna_kwargs
import pandas.core.algorithms as algos
from pandas.core.dtypes.common import (
    is_integer_dtype, is_float_dtype, is_period_dtype,
    pandas_dtype,
    is_datetime64_dtype,
    is_categorical_dtype,
    is_timedelta64_dtype,
    is_list_like,
    is_array_like,
    is_object_dtype,
    is_string_dtype,
    is_datetime_or_timedelta_dtype,
    is_dtype_equal,
    ensure_object,
    _TD_DTYPE,
)


from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.dtypes.generic import (
    ABCSeries, ABCIndexClass, ABCPeriodIndex
)
from pandas.core.dtypes.missing import isna
from pandas.core.missing import pad_1d, backfill_1d

import pandas.core.common as com

from pandas.tseries import frequencies
from pandas.tseries.offsets import Tick, DateOffset

from pandas.core.arrays import ExtensionArray
from pandas.core.arrays import datetimelike as dtl


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = frequencies.get_freq_code(self.freq)
        result = get_period_field_arr(alias, self.asi8, base)
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
        op = getattr(self.asi8, opname)
        # We want to eventually defer to the Series or PeriodIndex (which will
        # return here with an unboxed PeriodArray). But before we do that,
        # we do a bit of validation on type (Period) and freq, so that our
        # error messages are sensible
        not_implemented = isinstance(other, (ABCSeries, ABCIndexClass))
        if not_implemented:
            other = other._values

        if isinstance(other, Period):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            result = op(other.ordinal)
        elif isinstance(other, cls):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            if not_implemented:
                return NotImplemented
            result = op(other.asi8)

            mask = self._isnan | other._isnan
            if mask.any():
                result[mask] = nat_result

            return result
        elif other is NaT:
            result = np.empty(len(self.asi8), dtype=bool)
            result.fill(nat_result)
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

    Users should use :func:`period_array` to create new instances.

    Parameters
    ----------
    values : Union[PeriodArray, Series[period], ndarary[int], PeriodIndex]
        The data to store. These should be arrays that can be directly
        converted to ordinals without inference or copy (PeriodArray,
        ndarray[int64]), or a box around such an array (Series[period],
        PeriodIndex).
    freq : str or DateOffset
        The `freq` to use for the array. Mostly applicable when `values`
        is an ndarray of integers, when `freq` is required. When `values`
        is a PeriodArray (or box around), it's checked that ``values.freq``
        matches `freq`.
    copy : bool, default False
        Whether to copy the ordinals before storing.

    Notes
    -----
    There are two components to a PeriodArray

    - ordinals : integer ndarray
    - freq : pd.tseries.offsets.Offset

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
    def __init__(self, values, freq=None, copy=False):
        if freq is not None:
            freq = Period._maybe_convert_freq(freq)

        if isinstance(values, ABCSeries):
            values = values._values
            if not isinstance(values, type(self)):
                raise TypeError("Incorrect dtype")

        elif isinstance(values, ABCPeriodIndex):
            values = values._values

        if isinstance(values, type(self)):
            if freq is not None and freq != values.freq:
                msg = DIFFERENT_FREQ_INDEX.format(values.freq.freqstr,
                                                  freq.freqstr)
                raise IncompatibleFrequency(msg)
            values, freq = values._data, values.freq

        values = np.array(values, dtype='int64', copy=copy)
        self._data = values
        if freq is None:
            raise ValueError('freq is not specified and cannot be inferred')
        self._dtype = PeriodDtype(freq)

    @classmethod
    def _simple_new(cls, values, freq=None, **kwargs):
        # TODO(DatetimeArray): remove once all constructors are aligned.
        # alias from PeriodArray.__init__
        return cls(values, freq=freq, **kwargs)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        # type: (Sequence[Optional[Period]], PeriodDtype, bool) -> PeriodArray
        if dtype:
            freq = dtype.freq
        else:
            freq = None
        periods = np.asarray(scalars, dtype=object)
        if copy:
            periods = periods.copy()

        freq = freq or libperiod.extract_freq(periods)
        ordinals = libperiod.extract_ordinals(periods, freq)
        return cls(ordinals, freq=freq)

    def _values_for_factorize(self):
        return self.asi8, iNaT

    @classmethod
    def _from_factorized(cls, values, original):
        # type: (Sequence[Optional[Period]], PeriodArray) -> PeriodArray
        return cls(values, freq=original.freq)

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
        data, freq = dt64arr_to_periodarr(data, freq, tz)
        return cls(data, freq=freq)

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        periods = dtl.validate_periods(periods)

        if freq is not None:
            freq = Period._maybe_convert_freq(freq)

        field_count = len(fields)
        if start is not None or end is not None:
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
        return cls(values, freq=freq)

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
    def asi8(self):
        return self._data

    @property
    def freq(self):
        """Return the frequency object for this PeriodArray."""
        return self.dtype.freq

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

    def __setitem__(
            self,
            key,   # type: Union[int, Sequence[int], Sequence[bool]]
            value  # type: Union[NaTType, Period, Sequence[Period]]
    ):
        # type: (...) -> None
        # n.b. the type on `value` is a bit too restrictive.
        # we also accept a sequence of stuff coercible to a PeriodArray
        # by period_array, which includes things like ndarray[object],
        # ndarray[datetime64ns]. I think ndarray[int] / ndarray[str] won't
        # work, since the freq can't be inferred.
        if is_list_like(value):
            if len(key) != len(value) and not com.is_bool_indexer(key):
                msg = ("shape mismatch: value array of length '{}' does not "
                       "match indexing result of length '{}'.")
                raise ValueError(msg.format(len(key), len(value)))
            if len(key) == 0:
                return

            value = period_array(value)

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
            value = iNaT
        else:
            msg = ("'value' should be a 'Period', 'NaT', or array of those. "
                   "Got '{}' instead.".format(type(value).__name__))
            raise TypeError(msg)
        self._data[key] = value

    def take(self, indices, allow_fill=False, fill_value=None):
        if allow_fill:
            if isna(fill_value):
                fill_value = iNaT
            elif isinstance(fill_value, Period):
                if self.freq != fill_value.freq:
                    msg = DIFFERENT_FREQ_INDEX.format(
                        self.freq.freqstr,
                        fill_value.freqstr
                    )
                    raise IncompatibleFrequency(msg)

                fill_value = fill_value.ordinal
            else:
                msg = "'fill_value' should be a Period. Got '{}'."
                raise ValueError(msg.format(fill_value))

        new_values = algos.take(self._data,
                                indices,
                                allow_fill=allow_fill,
                                fill_value=fill_value)

        return type(self)(new_values, self.freq)

    def isna(self):
        return self._data == iNaT

    def fillna(self, value=None, method=None, limit=None):
        # TODO(#20300)
        # To avoid converting to object, we re-implement here with the changes
        # 1. Passing `_data` to func instead of self.astype(object)
        # 2. Re-boxing output of 1.
        # #20300 should let us do this kind of logic on ExtensionArray.fillna
        # and we can use it.

        if isinstance(value, ABCSeries):
            value = value._values

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
                new_values = func(self._data, limit=limit,
                                  mask=mask)
                new_values = type(self)(new_values, freq=self.freq)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def copy(self, deep=False):
        return type(self)(self._data.copy(), freq=self.freq)

    def value_counts(self, dropna=False):
        from pandas import Series, PeriodIndex

        if dropna:
            values = self[~self.isna()]._data
        else:
            values = self._data

        cls = type(self)

        result = algos.value_counts(values, sort=False)
        index = PeriodIndex(cls(result.index, freq=self.freq),
                            name=result.index.name)
        return Series(result.values, index=index, name=result.name)

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
        values = self._data + n * self.freq.n
        if self.hasnans:
            values[self._isnan] = iNaT
        return type(self)(values, freq=self.freq)

    @property
    def _box_func(self):
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

        return type(self)(new_data, freq=freq)

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
        from pandas.core.arrays import DatetimeArrayMixin

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

        new_data = libperiod.periodarr_to_dt64arr(new_data.asi8, base)
        return DatetimeArrayMixin(new_data, freq='infer')

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
        """ actually format my specific types """
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
        # TODO(DatetimeArray): remove
        nv.validate_repeat(args, kwargs)
        values = self._data.repeat(repeats)
        return type(self)(values, self.freq)

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
            values = self._data

            if values.dtype != dtype:
                # int32 vs. int64
                values = values.astype(dtype)

            elif copy:
                values = values.copy()

            return values
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

    @property
    def flags(self):
        # TODO: remove
        # We need this since reduction.SeriesBinGrouper uses values.flags
        # Ideally, we wouldn't be passing objects down there in the first
        # place.
        return self._data.flags

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

    def _addsub_int_array(
            self,
            other,  # type: Union[Index, ExtensionArray, np.ndarray[int]]
            op,     # type: Callable[Any, Any]
    ):
        # type: (...) -> PeriodArray
        assert op in [operator.add, operator.sub]
        # easy case for PeriodIndex
        if op is operator.sub:
            other = -other
        res_values = algos.checked_add_with_arr(self.asi8, other,
                                                arr_mask=self._isnan)
        res_values = res_values.view('i8')
        res_values[self._isnan] = iNaT
        return type(self)(res_values, freq=self.freq)

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
        ordinals = super(PeriodArray, self)._add_delta_td(delta)
        return type(self)(ordinals, self.freq)

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
            return self._add_delta_td(other)
        elif is_timedelta64_dtype(other):
            # ndarray[timedelta64] or TimedeltaArray/index
            return self._add_delta_tdi(other)
        else:  # pragma: no cover
            raise TypeError(type(other).__name__)


PeriodArray._add_comparison_ops()
PeriodArray._add_datetimelike_methods()


# -------------------------------------------------------------------
# Constructor Helpers

def period_array(data, freq=None, copy=False):
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
    copy : bool, default False
        Whether to ensure a copy of the data is made.

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

    Integers that look like years are handled

    >>> period_array([2000, 2001, 2002], freq='D')
    ['2000-01-01', '2001-01-01', '2002-01-01']
    Length: 3, dtype: period[D]

    Datetime-like strings may also be passed

    >>> period_array(['2000-Q1', '2000-Q2', '2000-Q3', '2000-Q4'], freq='Q')
    <PeriodArray>
    ['2000Q1', '2000Q2', '2000Q3', '2000Q4']
    Length: 4, dtype: period[Q-DEC]
    """
    if is_datetime64_dtype(data):
        return PeriodArray._from_datetime64(data, freq)
    if isinstance(data, (ABCPeriodIndex, ABCSeries, PeriodArray)):
        return PeriodArray(data, freq)

    # other iterable of some kind
    if not isinstance(data, (np.ndarray, list, tuple)):
        data = list(data)

    data = np.asarray(data)

    if freq:
        dtype = PeriodDtype(freq)
    else:
        dtype = None

    if is_float_dtype(data) and len(data) > 0:
        raise TypeError("PeriodIndex does not allow "
                        "floating point in construction")

    data = ensure_object(data)

    return PeriodArray._from_sequence(data, dtype=dtype)


def dt64arr_to_periodarr(data, freq, tz=None):
    """
    Convert an datetime-like array to values Period ordinals.

    Parameters
    ----------
    data : Union[Series[datetime64[ns]], DatetimeIndex, ndarray[datetime64ns]]
    freq : Optional[Union[str, Tick]]
        Must match the `freq` on the `data` if `data` is a DatetimeIndex
        or Series.
    tz : Optional[tzinfo]

    Returns
    -------
    ordinals : ndarray[int]
    freq : Tick
        The frequencey extracted from the Series or DatetimeIndex if that's
        used.

    """
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    if freq is not None:
        freq = Period._maybe_convert_freq(freq)

    if isinstance(data, ABCIndexClass):
        if freq is None:
            freq = data.freq
        elif freq != data.freq:
            msg = DIFFERENT_FREQ_INDEX.format(freq.freqstr, data.freq.freqstr)
            raise IncompatibleFrequency(msg)
        data = data._values

    elif isinstance(data, ABCSeries):
        if freq is None:
            freq = data.dt.freq
        elif freq != data.dt.freq:
            msg = DIFFERENT_FREQ_INDEX.format(freq.freqstr,
                                              data.dt.freq.freqstr)
            raise IncompatibleFrequency(msg)
        data = data._values

    base, mult = frequencies.get_freq_code(freq)
    return libperiod.dt64arr_to_periodarr(data.view('i8'), base, tz), freq


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
