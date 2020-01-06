from datetime import timedelta
import operator
from typing import Any, Callable, List, Optional, Sequence, Union

import numpy as np

from pandas._libs.tslibs import (
    NaT,
    NaTType,
    frequencies as libfrequencies,
    iNaT,
    period as libperiod,
)
from pandas._libs.tslibs.fields import isleapyear_arr
from pandas._libs.tslibs.period import (
    DIFFERENT_FREQ,
    IncompatibleFrequency,
    Period,
    get_period_field_arr,
    period_asfreq_arr,
)
from pandas._libs.tslibs.timedeltas import Timedelta, delta_to_nanoseconds
import pandas.compat as compat
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    _TD_DTYPE,
    ensure_object,
    is_datetime64_dtype,
    is_float_dtype,
    is_list_like,
    is_object_dtype,
    is_period_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.dtypes.generic import (
    ABCIndexClass,
    ABCPeriodArray,
    ABCPeriodIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna, notna

from pandas.core import ops
import pandas.core.algorithms as algos
from pandas.core.arrays import datetimelike as dtl
import pandas.core.common as com
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.ops.invalid import invalid_comparison

from pandas.tseries import frequencies
from pandas.tseries.offsets import DateOffset, Tick, _delta_to_tick


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = libfrequencies.get_freq_code(self.freq)
        result = get_period_field_arr(alias, self.asi8, base)
        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _period_array_cmp(cls, op):
    """
    Wrap comparison operations to convert Period-like to PeriodDtype
    """
    opname = f"__{op.__name__}__"
    nat_result = opname == "__ne__"

    @unpack_zerodim_and_defer(opname)
    def wrapper(self, other):

        if isinstance(other, str):
            try:
                other = self._scalar_from_string(other)
            except ValueError:
                # string that can't be parsed as Period
                return invalid_comparison(self, other, op)

        if isinstance(other, self._recognized_scalars) or other is NaT:
            other = self._scalar_type(other)
            self._check_compatible_with(other)

            other_i8 = self._unbox_scalar(other)

            result = op(self.view("i8"), other_i8)
            if isna(other):
                result.fill(nat_result)

        elif not is_list_like(other):
            return invalid_comparison(self, other, op)

        elif len(other) != len(self):
            raise ValueError("Lengths must match")

        else:
            if isinstance(other, list):
                # TODO: could use pd.Index to do inference?
                other = np.array(other)

            if not isinstance(other, (np.ndarray, cls)):
                return invalid_comparison(self, other, op)

            if is_object_dtype(other):
                with np.errstate(all="ignore"):
                    result = ops.comp_method_OBJECT_ARRAY(
                        op, self.astype(object), other
                    )
                o_mask = isna(other)

            elif not cls._is_recognized_dtype(other.dtype):
                # e.g. is_timedelta64_dtype(other)
                return invalid_comparison(self, other, op)

            else:
                assert isinstance(other, cls), type(other)

                self._check_compatible_with(other)

                result = op(self.view("i8"), other.view("i8"))
                o_mask = other._isnan

            if o_mask.any():
                result[o_mask] = nat_result

        if self._hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


class PeriodArray(dtl.DatetimeLikeArrayMixin, dtl.DatelikeOps):
    """
    Pandas ExtensionArray for storing Period data.

    Users should use :func:`period_array` to create new instances.

    Parameters
    ----------
    values : Union[PeriodArray, Series[period], ndarray[int], PeriodIndex]
        The data to store. These should be arrays that can be directly
        converted to ordinals without inference or copy (PeriodArray,
        ndarray[int64]), or a box around such an array (Series[period],
        PeriodIndex).
    freq : str or DateOffset
        The `freq` to use for the array. Mostly applicable when `values`
        is an ndarray of integers, when `freq` is required. When `values`
        is a PeriodArray (or box around), it's checked that ``values.freq``
        matches `freq`.
    dtype : PeriodDtype, optional
        A PeriodDtype instance from which to extract a `freq`. If both
        `freq` and `dtype` are specified, then the frequencies must match.
    copy : bool, default False
        Whether to copy the ordinals before storing.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    period_array : Create a new PeriodArray.
    PeriodIndex : Immutable Index for period data.

    Notes
    -----
    There are two components to a PeriodArray

    - ordinals : integer ndarray
    - freq : pd.tseries.offsets.Offset

    The values are physically stored as a 1-D ndarray of integers. These are
    called "ordinals" and represent some kind of offset from a base.

    The `freq` indicates the span covered by each element of the array.
    All elements in the PeriodArray have the same `freq`.
    """

    # array priority higher than numpy scalars
    __array_priority__ = 1000
    _typ = "periodarray"  # ABCPeriodArray
    _scalar_type = Period
    _recognized_scalars = (Period,)
    _is_recognized_dtype = is_period_dtype

    # Names others delegate to us
    _other_ops: List[str] = []
    _bool_ops = ["is_leap_year"]
    _object_ops = ["start_time", "end_time", "freq"]
    _field_ops = [
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "weekofyear",
        "weekday",
        "week",
        "dayofweek",
        "dayofyear",
        "quarter",
        "qyear",
        "days_in_month",
        "daysinmonth",
    ]
    _datetimelike_ops = _field_ops + _object_ops + _bool_ops
    _datetimelike_methods = ["strftime", "to_timestamp", "asfreq"]

    # --------------------------------------------------------------------
    # Constructors

    def __init__(self, values, freq=None, dtype=None, copy=False):
        freq = validate_dtype_freq(dtype, freq)

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
                raise raise_on_incompatible(values, freq)
            values, freq = values._data, values.freq

        values = np.array(values, dtype="int64", copy=copy)
        self._data = values
        if freq is None:
            raise ValueError("freq is not specified and cannot be inferred")
        self._dtype = PeriodDtype(freq)

    @classmethod
    def _simple_new(cls, values, freq=None, **kwargs):
        # alias for PeriodArray.__init__
        return cls(values, freq=freq, **kwargs)

    @classmethod
    def _from_sequence(
        cls,
        scalars: Sequence[Optional[Period]],
        dtype: Optional[PeriodDtype] = None,
        copy: bool = False,
    ) -> ABCPeriodArray:
        if dtype:
            freq = dtype.freq
        else:
            freq = None

        if isinstance(scalars, cls):
            validate_dtype_freq(scalars.dtype, freq)
            if copy:
                scalars = scalars.copy()
            return scalars

        periods = np.asarray(scalars, dtype=object)
        if copy:
            periods = periods.copy()

        freq = freq or libperiod.extract_freq(periods)
        ordinals = libperiod.extract_ordinals(periods, freq)
        return cls(ordinals, freq=freq)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        return cls._from_sequence(strings, dtype, copy)

    @classmethod
    def _from_datetime64(cls, data, freq, tz=None):
        """
        Construct a PeriodArray from a datetime64 array

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
                raise ValueError(
                    "Can either instantiate from fields or endpoints, but not both"
                )
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
        elif field_count > 0:
            subarr, freq = _range_from_fields(freq=freq, **fields)
        else:
            raise ValueError("Not enough parameters to construct Period range")

        return subarr, freq

    # -----------------------------------------------------------------
    # DatetimeLike Interface

    def _unbox_scalar(self, value: Union[Period, NaTType]) -> int:
        if value is NaT:
            return value.value
        elif isinstance(value, self._scalar_type):
            if not isna(value):
                self._check_compatible_with(value)
            return value.ordinal
        else:
            raise ValueError(f"'value' should be a Period. Got '{value}' instead.")

    def _scalar_from_string(self, value: str) -> Period:
        return Period(value, freq=self.freq)

    def _check_compatible_with(self, other, setitem: bool = False):
        if other is NaT:
            return
        if self.freqstr != other.freqstr:
            raise raise_on_incompatible(self, other)

    # --------------------------------------------------------------------
    # Data / Attributes

    @cache_readonly
    def dtype(self):
        return self._dtype

    # error: Read-only property cannot override read-write property  [misc]
    @property  # type: ignore
    def freq(self):
        """
        Return the frequency object for this PeriodArray.
        """
        return self.dtype.freq

    def __array__(self, dtype=None):
        # overriding DatetimelikeArray
        return np.array(list(self), dtype=object)

    # --------------------------------------------------------------------
    # Vectorized analogues of Period properties

    year = _field_accessor(
        "year",
        0,
        """
        The year of the period.
        """,
    )
    month = _field_accessor(
        "month",
        3,
        """
        The month as January=1, December=12.
        """,
    )
    day = _field_accessor(
        "day",
        4,
        """
        The days of the period.
        """,
    )
    hour = _field_accessor(
        "hour",
        5,
        """
        The hour of the period.
        """,
    )
    minute = _field_accessor(
        "minute",
        6,
        """
        The minute of the period.
        """,
    )
    second = _field_accessor(
        "second",
        7,
        """
        The second of the period.
        """,
    )
    weekofyear = _field_accessor(
        "week",
        8,
        """
        The week ordinal of the year.
        """,
    )
    week = weekofyear
    dayofweek = _field_accessor(
        "dayofweek",
        10,
        """
        The day of the week with Monday=0, Sunday=6.
        """,
    )
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor(
        "dayofyear",
        9,
        """
        The ordinal day of the year.
        """,
    )
    quarter = _field_accessor(
        "quarter",
        2,
        """
        The quarter of the date.
        """,
    )
    qyear = _field_accessor("qyear", 1)
    days_in_month = _field_accessor(
        "days_in_month",
        11,
        """
        The number of days in the month.
        """,
    )
    daysinmonth = days_in_month

    @property
    def is_leap_year(self):
        """
        Logical indicating if the date belongs to a leap year.
        """
        return isleapyear_arr(np.asarray(self.year))

    @property
    def start_time(self):
        return self.to_timestamp(how="start")

    @property
    def end_time(self):
        return self.to_timestamp(how="end")

    def to_timestamp(self, freq=None, how="start"):
        """
        Cast to DatetimeArray/Index.

        Parameters
        ----------
        freq : str or DateOffset, optional
            Target frequency. The default is 'D' for week or longer,
            'S' otherwise.
        how : {'s', 'e', 'start', 'end'}
            Whether to use the start or end of the time period being converted.

        Returns
        -------
        DatetimeArray/Index
        """
        from pandas.core.arrays import DatetimeArray

        how = libperiod._validate_end_alias(how)

        end = how == "E"
        if end:
            if freq == "B":
                # roll forward to ensure we land on B date
                adjust = Timedelta(1, "D") - Timedelta(1, "ns")
                return self.to_timestamp(how="start") + adjust
            else:
                adjust = Timedelta(1, "ns")
                return (self + self.freq).to_timestamp(how="start") - adjust

        if freq is None:
            base, mult = libfrequencies.get_freq_code(self.freq)
            freq = libfrequencies.get_to_timestamp_base(base)
        else:
            freq = Period._maybe_convert_freq(freq)

        base, mult = libfrequencies.get_freq_code(freq)
        new_data = self.asfreq(freq, how=how)

        new_data = libperiod.periodarr_to_dt64arr(new_data.asi8, base)
        return DatetimeArray._from_sequence(new_data, freq="infer")

    # --------------------------------------------------------------------
    # Array-like / EA-Interface Methods

    def _values_for_argsort(self):
        return self._data

    # --------------------------------------------------------------------

    def _time_shift(self, periods, freq=None):
        """
        Shift each value by `periods`.

        Note this is different from ExtensionArray.shift, which
        shifts the *position* of each element, padding the end with
        missing values.

        Parameters
        ----------
        periods : int
            Number of periods to shift by.
        freq : pandas.DateOffset, pandas.Timedelta, or str
            Frequency increment to shift by.
        """
        if freq is not None:
            raise TypeError(
                "`freq` argument is not supported for "
                f"{type(self).__name__}._time_shift"
            )
        values = self.asi8 + periods * self.freq.n
        if self._hasnans:
            values[self._isnan] = iNaT
        return type(self)(values, freq=self.freq)

    @property
    def _box_func(self):
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    def asfreq(self, freq=None, how="E"):
        """
        Convert the Period Array/Index to the specified frequency `freq`.

        Parameters
        ----------
        freq : str
            A frequency.
        how : str {'E', 'S'}
            Whether the elements should be aligned to the end
            or start within pa period.

            * 'E', 'END', or 'FINISH' for end,
            * 'S', 'START', or 'BEGIN' for start.

            January 31st ('END') vs. January 1st ('START') for example.

        Returns
        -------
        Period Array/Index
            Constructed with the new frequency.

        Examples
        --------
        >>> pidx = pd.period_range('2010-01-01', '2015-01-01', freq='A')
        >>> pidx
        PeriodIndex(['2010', '2011', '2012', '2013', '2014', '2015'],
        dtype='period[A-DEC]', freq='A-DEC')

        >>> pidx.asfreq('M')
        PeriodIndex(['2010-12', '2011-12', '2012-12', '2013-12', '2014-12',
        '2015-12'], dtype='period[M]', freq='M')

        >>> pidx.asfreq('M', how='S')
        PeriodIndex(['2010-01', '2011-01', '2012-01', '2013-01', '2014-01',
        '2015-01'], dtype='period[M]', freq='M')
        """
        how = libperiod._validate_end_alias(how)

        freq = Period._maybe_convert_freq(freq)

        base1, mult1 = libfrequencies.get_freq_code(self.freq)
        base2, mult2 = libfrequencies.get_freq_code(freq)

        asi8 = self.asi8
        # mult1 can't be negative or 0
        end = how == "E"
        if end:
            ordinal = asi8 + mult1 - 1
        else:
            ordinal = asi8

        new_data = period_asfreq_arr(ordinal, base1, base2, end)

        if self._hasnans:
            new_data[self._isnan] = iNaT

        return type(self)(new_data, freq=freq)

    # ------------------------------------------------------------------
    # Rendering Methods

    def _formatter(self, boxed=False):
        if boxed:
            return str
        return "'{}'".format

    def _format_native_types(self, na_rep="NaT", date_format=None, **kwargs):
        """
        actually format my specific types
        """
        values = self.astype(object)

        if date_format:
            formatter = lambda dt: dt.strftime(date_format)
        else:
            formatter = lambda dt: str(dt)

        if self._hasnans:
            mask = self._isnan
            values[mask] = na_rep
            imask = ~mask
            values[imask] = np.array([formatter(dt) for dt in values[imask]])
        else:
            values = np.array([formatter(dt) for dt in values])
        return values

    # ------------------------------------------------------------------

    def astype(self, dtype, copy=True):
        # We handle Period[T] -> Period[U]
        # Our parent handles everything else.
        dtype = pandas_dtype(dtype)

        if is_period_dtype(dtype):
            return self.asfreq(dtype.freq)
        return super().astype(dtype, copy=copy)

    # ------------------------------------------------------------------
    # Arithmetic Methods
    _create_comparison_method = classmethod(_period_array_cmp)

    def _sub_datelike(self, other):
        assert other is not NaT
        return NotImplemented

    def _sub_period(self, other):
        # If the operation is well-defined, we return an object-Index
        # of DateOffsets.  Null entries are filled with pd.NaT
        self._check_compatible_with(other)
        asi8 = self.asi8
        new_data = asi8 - other.ordinal
        new_data = np.array([self.freq * x for x in new_data])

        if self._hasnans:
            new_data[self._isnan] = NaT

        return new_data

    def _addsub_int_array(
        self, other: np.ndarray, op: Callable[[Any, Any], Any],
    ) -> "PeriodArray":
        """
        Add or subtract array of integers; equivalent to applying
        `_time_shift` pointwise.

        Parameters
        ----------
        other : np.ndarray[integer-dtype]
        op : {operator.add, operator.sub}

        Returns
        -------
        result : PeriodArray
        """

        assert op in [operator.add, operator.sub]
        if op is operator.sub:
            other = -other
        res_values = algos.checked_add_with_arr(self.asi8, other, arr_mask=self._isnan)
        res_values = res_values.view("i8")
        res_values[self._isnan] = iNaT
        return type(self)(res_values, freq=self.freq)

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        base = libfrequencies.get_base_alias(other.rule_code)
        if base != self.freq.rule_code:
            raise raise_on_incompatible(self, other)

        # Note: when calling parent class's _add_timedeltalike_scalar,
        #  it will call delta_to_nanoseconds(delta).  Because delta here
        #  is an integer, delta_to_nanoseconds will return it unchanged.
        result = super()._add_timedeltalike_scalar(other.n)
        return type(self)(result, freq=self.freq)

    def _add_timedeltalike_scalar(self, other):
        """
        Parameters
        ----------
        other : timedelta, Tick, np.timedelta64

        Returns
        -------
        result : ndarray[int64]
        """
        assert isinstance(self.freq, Tick)  # checked by calling function
        assert isinstance(other, (timedelta, np.timedelta64, Tick))

        if notna(other):
            # special handling for np.timedelta64("NaT"), avoid calling
            #  _check_timedeltalike_freq_compat as that would raise TypeError
            other = self._check_timedeltalike_freq_compat(other)

        # Note: when calling parent class's _add_timedeltalike_scalar,
        #  it will call delta_to_nanoseconds(delta).  Because delta here
        #  is an integer, delta_to_nanoseconds will return it unchanged.
        ordinals = super()._add_timedeltalike_scalar(other)
        return ordinals

    def _add_delta_tdi(self, other):
        """
        Parameters
        ----------
        other : TimedeltaArray or ndarray[timedelta64]

        Returns
        -------
        result : ndarray[int64]
        """
        assert isinstance(self.freq, Tick)  # checked by calling function

        if not np.all(isna(other)):
            delta = self._check_timedeltalike_freq_compat(other)
        else:
            # all-NaT TimedeltaIndex is equivalent to a single scalar td64 NaT
            return self + np.timedelta64("NaT")

        return self._addsub_int_array(delta, operator.add).asi8

    def _add_delta(self, other):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self, yielding a new PeriodArray

        Parameters
        ----------
        other : {timedelta, np.timedelta64, Tick,
                 TimedeltaIndex, ndarray[timedelta64]}

        Returns
        -------
        result : PeriodArray
        """
        if not isinstance(self.freq, Tick):
            # We cannot add timedelta-like to non-tick PeriodArray
            raise raise_on_incompatible(self, other)

        new_ordinals = super()._add_delta(other)
        return type(self)(new_ordinals, freq=self.freq)

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
            assert other.dtype.kind == "m"
            if other.dtype != _TD_DTYPE:
                # i.e. non-nano unit
                # TODO: disallow unit-less timedelta64
                other = other.astype(_TD_DTYPE)
            nanos = other.view("i8")
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

        raise raise_on_incompatible(self, other)


PeriodArray._add_comparison_ops()


def raise_on_incompatible(left, right):
    """
    Helper function to render a consistent error message when raising
    IncompatibleFrequency.

    Parameters
    ----------
    left : PeriodArray
    right : None, DateOffset, Period, ndarray, or timedelta-like

    Returns
    ------
    IncompatibleFrequency
        Exception to be raised by the caller.
    """
    # GH#24283 error message format depends on whether right is scalar
    if isinstance(right, np.ndarray) or right is None:
        other_freq = None
    elif isinstance(right, (ABCPeriodIndex, PeriodArray, Period, DateOffset)):
        other_freq = right.freqstr
    else:
        other_freq = _delta_to_tick(Timedelta(right)).freqstr

    msg = DIFFERENT_FREQ.format(
        cls=type(left).__name__, own_freq=left.freqstr, other_freq=other_freq
    )
    return IncompatibleFrequency(msg)


# -------------------------------------------------------------------
# Constructor Helpers


def period_array(
    data: Sequence[Optional[Period]],
    freq: Optional[Union[str, Tick]] = None,
    copy: bool = False,
) -> PeriodArray:
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

    dtype: Optional[PeriodDtype]
    if freq:
        dtype = PeriodDtype(freq)
    else:
        dtype = None

    if is_float_dtype(data) and len(data) > 0:
        raise TypeError("PeriodIndex does not allow floating point in construction")

    data = ensure_object(data)

    return PeriodArray._from_sequence(data, dtype=dtype)


def validate_dtype_freq(dtype, freq):
    """
    If both a dtype and a freq are available, ensure they match.  If only
    dtype is available, extract the implied freq.

    Parameters
    ----------
    dtype : dtype
    freq : DateOffset or None

    Returns
    -------
    freq : DateOffset

    Raises
    ------
    ValueError : non-period dtype
    IncompatibleFrequency : mismatch between dtype and freq
    """
    if freq is not None:
        freq = frequencies.to_offset(freq)

    if dtype is not None:
        dtype = pandas_dtype(dtype)
        if not is_period_dtype(dtype):
            raise ValueError("dtype must be PeriodDtype")
        if freq is None:
            freq = dtype.freq
        elif freq != dtype.freq:
            raise IncompatibleFrequency("specified freq and dtype are different")
    return freq


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
        The frequency extracted from the Series or DatetimeIndex if that's
        used.

    """
    if data.dtype != np.dtype("M8[ns]"):
        raise ValueError(f"Wrong dtype: {data.dtype}")

    if freq is None:
        if isinstance(data, ABCIndexClass):
            data, freq = data._values, data.freq
        elif isinstance(data, ABCSeries):
            data, freq = data._values, data.dt.freq

    freq = Period._maybe_convert_freq(freq)

    if isinstance(data, (ABCIndexClass, ABCSeries)):
        data = data._values

    base, mult = libfrequencies.get_freq_code(freq)
    return libperiod.dt64arr_to_periodarr(data.view("i8"), base, tz), freq


def _get_ordinal_range(start, end, periods, freq, mult=1):
    if com.count_not_none(start, end, periods) != 2:
        raise ValueError(
            "Of the three parameters: start, end, and periods, "
            "exactly two must be specified"
        )

    if freq is not None:
        _, mult = libfrequencies.get_freq_code(freq)

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError("start and end must have same freq")
    if start is NaT or end is NaT:
        raise ValueError("start and end must not be NaT")

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:  # pragma: no cover
            raise ValueError("Could not infer freq from start/end")

    if periods is not None:
        periods = periods * mult
        if start is None:
            data = np.arange(
                end.ordinal - periods + mult, end.ordinal + 1, mult, dtype=np.int64
            )
        else:
            data = np.arange(
                start.ordinal, start.ordinal + periods, mult, dtype=np.int64
            )
    else:
        data = np.arange(start.ordinal, end.ordinal + 1, mult, dtype=np.int64)

    return data, freq


def _range_from_fields(
    year=None,
    month=None,
    quarter=None,
    day=None,
    hour=None,
    minute=None,
    second=None,
    freq=None,
):
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
            freq = "Q"
            base = libfrequencies.FreqGroup.FR_QTR
        else:
            base, mult = libfrequencies.get_freq_code(freq)
            if base != libfrequencies.FreqGroup.FR_QTR:
                raise AssertionError("base must equal FR_QTR")

        year, quarter = _make_field_arrays(year, quarter)
        for y, q in zip(year, quarter):
            y, m = libperiod.quarter_to_myear(y, q, freq)
            val = libperiod.period_ordinal(y, m, 1, 1, 1, 1, 0, 0, base)
            ordinals.append(val)
    else:
        base, mult = libfrequencies.get_freq_code(freq)
        arrays = _make_field_arrays(year, month, day, hour, minute, second)
        for y, mth, d, h, mn, s in zip(*arrays):
            ordinals.append(libperiod.period_ordinal(y, mth, d, h, mn, s, 0, 0, base))

    return np.array(ordinals, dtype=np.int64), freq


def _make_field_arrays(*fields):
    length = None
    for x in fields:
        if isinstance(x, (list, np.ndarray, ABCSeries)):
            if length is not None and len(x) != length:
                raise ValueError("Mismatched Period array lengths")
            elif length is None:
                length = len(x)

    arrays = [
        np.asarray(x)
        if isinstance(x, (np.ndarray, list, ABCSeries))
        else np.repeat(x, length)
        for x in fields
    ]

    return arrays
