from datetime import timedelta
from typing import List

import numpy as np

from pandas._libs import lib, tslibs
from pandas._libs.tslibs import NaT, Timedelta, Timestamp, iNaT
from pandas._libs.tslibs.fields import get_timedelta_field
from pandas._libs.tslibs.timedeltas import (
    array_to_timedelta64,
    parse_timedelta_unit,
    precision_from_unit,
)
import pandas.compat as compat
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    _NS_DTYPE,
    _TD_DTYPE,
    is_dtype_equal,
    is_float_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
    is_string_dtype,
    is_timedelta64_dtype,
    is_timedelta64_ns_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndexClass,
    ABCSeries,
    ABCTimedeltaIndex,
)
from pandas.core.dtypes.missing import isna

from pandas.core import nanops, ops
from pandas.core.algorithms import checked_add_with_arr
import pandas.core.common as com
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.ops.invalid import invalid_comparison

from pandas.tseries.frequencies import to_offset
from pandas.tseries.offsets import Tick

from . import datetimelike as dtl

_BAD_DTYPE = "dtype {dtype} cannot be converted to timedelta64[ns]"


def _is_convertible_to_td(key):
    return isinstance(key, (Tick, timedelta, np.timedelta64, str))


def _field_accessor(name, alias, docstring=None):
    def f(self):
        values = self.asi8
        result = get_timedelta_field(values, alias)
        if self._hasnans:
            result = self._maybe_mask_results(
                result, fill_value=None, convert="float64"
            )

        return result

    f.__name__ = name
    f.__doc__ = f"\n{docstring}\n"
    return property(f)


def _td_array_cmp(cls, op):
    """
    Wrap comparison operations to convert timedelta-like to timedelta64
    """
    opname = f"__{op.__name__}__"
    nat_result = opname == "__ne__"

    @unpack_zerodim_and_defer(opname)
    def wrapper(self, other):

        if isinstance(other, str):
            try:
                other = self._scalar_from_string(other)
            except ValueError:
                # failed to parse as timedelta
                return invalid_comparison(self, other, op)

        if _is_convertible_to_td(other) or other is NaT:
            other = Timedelta(other)

            result = op(self.view("i8"), other.value)
            if isna(other):
                result.fill(nat_result)

        elif not is_list_like(other):
            return invalid_comparison(self, other, op)

        elif len(other) != len(self):
            raise ValueError("Lengths must match")

        else:
            if isinstance(other, list):
                other = np.array(other)

            if not isinstance(other, (np.ndarray, cls)):
                return invalid_comparison(self, other, op)

            if is_object_dtype(other):
                with np.errstate(all="ignore"):
                    result = ops.comp_method_OBJECT_ARRAY(
                        op, self.astype(object), other
                    )
                o_mask = isna(other)

            elif not is_timedelta64_dtype(other):
                # e.g. other is datetimearray
                return invalid_comparison(self, other, op)

            else:
                other = type(self)._from_sequence(other)

                result = op(self.view("i8"), other.view("i8"))
                o_mask = other._isnan

            if o_mask.any():
                result[o_mask] = nat_result

        if self._hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


class TimedeltaArray(dtl.DatetimeLikeArrayMixin, dtl.TimelikeOps):
    """
    Pandas ExtensionArray for timedelta data.

    .. versionadded:: 0.24.0

    .. warning::

       TimedeltaArray is currently experimental, and its API may change
       without warning. In particular, :attr:`TimedeltaArray.dtype` is
       expected to change to be an instance of an ``ExtensionDtype``
       subclass.

    Parameters
    ----------
    values : array-like
        The timedelta data.

    dtype : numpy.dtype
        Currently, only ``numpy.dtype("timedelta64[ns]")`` is accepted.
    freq : Offset, optional
    copy : bool, default False
        Whether to copy the underlying array of data.

    Attributes
    ----------
    None

    Methods
    -------
    None
    """

    _typ = "timedeltaarray"
    _scalar_type = Timedelta
    __array_priority__ = 1000
    # define my properties & methods for delegation
    _other_ops: List[str] = []
    _bool_ops: List[str] = []
    _object_ops = ["freq"]
    _field_ops = ["days", "seconds", "microseconds", "nanoseconds"]
    _datetimelike_ops = _field_ops + _object_ops + _bool_ops
    _datetimelike_methods = [
        "to_pytimedelta",
        "total_seconds",
        "round",
        "floor",
        "ceil",
    ]

    # Note: ndim must be defined to ensure NaT.__richcmp(TimedeltaArray)
    #  operates pointwise.

    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit="ns")

    @property
    def dtype(self):
        """
        The dtype for the TimedeltaArray.

        .. warning::

           A future version of pandas will change dtype to be an instance
           of a :class:`pandas.api.extensions.ExtensionDtype` subclass,
           not a ``numpy.dtype``.

        Returns
        -------
        numpy.dtype
        """
        return _TD_DTYPE

    # ----------------------------------------------------------------
    # Constructors

    def __init__(self, values, dtype=_TD_DTYPE, freq=None, copy=False):
        if isinstance(values, (ABCSeries, ABCIndexClass)):
            values = values._values

        inferred_freq = getattr(values, "_freq", None)

        if isinstance(values, type(self)):
            if freq is None:
                freq = values.freq
            elif freq and values.freq:
                freq = to_offset(freq)
                freq, _ = dtl.validate_inferred_freq(freq, values.freq, False)
            values = values._data

        if not isinstance(values, np.ndarray):
            msg = (
                f"Unexpected type '{type(values).__name__}'. 'values' must be a"
                " TimedeltaArray ndarray, or Series or Index containing one of those."
            )
            raise ValueError(msg)
        if values.ndim not in [1, 2]:
            raise ValueError("Only 1-dimensional input arrays are supported.")

        if values.dtype == "i8":
            # for compat with datetime/timedelta/period shared methods,
            #  we can sometimes get here with int64 values.  These represent
            #  nanosecond UTC (or tz-naive) unix timestamps
            values = values.view(_TD_DTYPE)

        _validate_td64_dtype(values.dtype)
        dtype = _validate_td64_dtype(dtype)

        if freq == "infer":
            msg = (
                "Frequency inference not allowed in TimedeltaArray.__init__. "
                "Use 'pd.array()' instead."
            )
            raise ValueError(msg)

        if copy:
            values = values.copy()
        if freq:
            freq = to_offset(freq)

        self._data = values
        self._dtype = dtype
        self._freq = freq

        if inferred_freq is None and freq is not None:
            type(self)._validate_frequency(self, freq)

    @classmethod
    def _simple_new(cls, values, freq=None, dtype=_TD_DTYPE):
        assert dtype == _TD_DTYPE, dtype
        assert isinstance(values, np.ndarray), type(values)

        result = object.__new__(cls)
        result._data = values.view(_TD_DTYPE)
        result._freq = to_offset(freq)
        result._dtype = _TD_DTYPE
        return result

    @classmethod
    def _from_sequence(cls, data, dtype=_TD_DTYPE, copy=False, freq=None, unit=None):
        if dtype:
            _validate_td64_dtype(dtype)
        freq, freq_infer = dtl.maybe_infer_freq(freq)

        data, inferred_freq = sequence_to_td64ns(data, copy=copy, unit=unit)
        freq, freq_infer = dtl.validate_inferred_freq(freq, inferred_freq, freq_infer)

        result = cls._simple_new(data, freq=freq)

        if inferred_freq is None and freq is not None:
            # this condition precludes `freq_infer`
            cls._validate_frequency(result, freq)

        elif freq_infer:
            # Set _freq directly to bypass duplicative _validate_frequency
            # check.
            result._freq = to_offset(result.inferred_freq)

        return result

    @classmethod
    def _generate_range(cls, start, end, periods, freq, closed=None):

        periods = dtl.validate_periods(periods)
        if freq is None and any(x is None for x in [periods, start, end]):
            raise ValueError("Must provide freq argument if no data is supplied")

        if com.count_not_none(start, end, periods, freq) != 3:
            raise ValueError(
                "Of the four parameters: start, end, periods, "
                "and freq, exactly three must be specified"
            )

        if start is not None:
            start = Timedelta(start)

        if end is not None:
            end = Timedelta(end)

        if start is None and end is None:
            if closed is not None:
                raise ValueError(
                    "Closed has to be None if not both of startand end are defined"
                )

        left_closed, right_closed = dtl.validate_endpoints(closed)

        if freq is not None:
            index = _generate_regular_range(start, end, periods, freq)
        else:
            index = np.linspace(start.value, end.value, periods).astype("i8")

        if not left_closed:
            index = index[1:]
        if not right_closed:
            index = index[:-1]

        return cls._simple_new(index, freq=freq)

    # ----------------------------------------------------------------
    # DatetimeLike Interface

    def _unbox_scalar(self, value):
        if not isinstance(value, self._scalar_type) and value is not NaT:
            raise ValueError("'value' should be a Timedelta.")
        self._check_compatible_with(value)
        return value.value

    def _scalar_from_string(self, value):
        return Timedelta(value)

    def _check_compatible_with(self, other):
        # we don't have anything to validate.
        pass

    def _maybe_clear_freq(self):
        self._freq = None

    # ----------------------------------------------------------------
    # Array-Like / EA-Interface Methods

    @Appender(dtl.DatetimeLikeArrayMixin._validate_fill_value.__doc__)
    def _validate_fill_value(self, fill_value):
        if isna(fill_value):
            fill_value = iNaT
        elif isinstance(fill_value, (timedelta, np.timedelta64, Tick)):
            fill_value = Timedelta(fill_value).value
        else:
            raise ValueError(f"'fill_value' should be a Timedelta. Got '{fill_value}'.")
        return fill_value

    def astype(self, dtype, copy=True):
        # We handle
        #   --> timedelta64[ns]
        #   --> timedelta64
        # DatetimeLikeArrayMixin super call handles other cases
        dtype = pandas_dtype(dtype)

        if is_timedelta64_dtype(dtype) and not is_timedelta64_ns_dtype(dtype):
            # by pandas convention, converting to non-nano timedelta64
            #  returns an int64-dtyped array with ints representing multiples
            #  of the desired timedelta unit.  This is essentially division
            if self._hasnans:
                # avoid double-copying
                result = self._data.astype(dtype, copy=False)
                values = self._maybe_mask_results(
                    result, fill_value=None, convert="float64"
                )
                return values
            result = self._data.astype(dtype, copy=copy)
            return result.astype("i8")
        elif is_timedelta64_ns_dtype(dtype):
            if copy:
                return self.copy()
            return self
        return dtl.DatetimeLikeArrayMixin.astype(self, dtype, copy=copy)

    # ----------------------------------------------------------------
    # Reductions

    def sum(
        self,
        axis=None,
        dtype=None,
        out=None,
        keepdims: bool = False,
        initial=None,
        skipna: bool = True,
        min_count: int = 0,
    ):
        nv.validate_sum(
            (), dict(dtype=dtype, out=out, keepdims=keepdims, initial=initial)
        )
        if not len(self):
            return NaT
        if not skipna and self._hasnans:
            return NaT

        result = nanops.nansum(
            self._data, axis=axis, skipna=skipna, min_count=min_count
        )
        return Timedelta(result)

    def std(
        self,
        axis=None,
        dtype=None,
        out=None,
        ddof: int = 1,
        keepdims: bool = False,
        skipna: bool = True,
    ):
        nv.validate_stat_ddof_func(
            (), dict(dtype=dtype, out=out, keepdims=keepdims), fname="std"
        )
        if not len(self):
            return NaT
        if not skipna and self._hasnans:
            return NaT

        result = nanops.nanstd(self._data, axis=axis, skipna=skipna, ddof=ddof)
        return Timedelta(result)

    def median(
        self,
        axis=None,
        out=None,
        overwrite_input: bool = False,
        keepdims: bool = False,
        skipna: bool = True,
    ):
        nv.validate_median(
            (), dict(out=out, overwrite_input=overwrite_input, keepdims=keepdims)
        )
        return nanops.nanmedian(self._data, axis=axis, skipna=skipna)

    # ----------------------------------------------------------------
    # Rendering Methods

    def _formatter(self, boxed=False):
        from pandas.io.formats.format import _get_format_timedelta64

        return _get_format_timedelta64(self, box=True)

    def _format_native_types(self, na_rep="NaT", date_format=None, **kwargs):
        from pandas.io.formats.format import _get_format_timedelta64

        formatter = _get_format_timedelta64(self._data, na_rep)
        return np.array([formatter(x) for x in self._data])

    # ----------------------------------------------------------------
    # Arithmetic Methods

    _create_comparison_method = classmethod(_td_array_cmp)

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        raise TypeError(
            f"cannot add the type {type(other).__name__} to a {type(self).__name__}"
        )

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
        new_values = super()._add_delta(delta)
        return type(self)._from_sequence(new_values, freq="infer")

    def _add_datetime_arraylike(self, other):
        """
        Add DatetimeArray/Index or ndarray[datetime64] to TimedeltaArray.
        """
        if isinstance(other, np.ndarray):
            # At this point we have already checked that dtype is datetime64
            from pandas.core.arrays import DatetimeArray

            other = DatetimeArray(other)

        # defer to implementation in DatetimeArray
        return other + self

    def _add_datetimelike_scalar(self, other):
        # adding a timedeltaindex to a datetimelike
        from pandas.core.arrays import DatetimeArray

        assert other is not NaT
        other = Timestamp(other)
        if other is NaT:
            # In this case we specifically interpret NaT as a datetime, not
            # the timedelta interpretation we would get by returning self + NaT
            result = self.asi8.view("m8[ms]") + NaT.to_datetime64()
            return DatetimeArray(result)

        i8 = self.asi8
        result = checked_add_with_arr(i8, other.value, arr_mask=self._isnan)
        result = self._maybe_mask_results(result)
        dtype = DatetimeTZDtype(tz=other.tz) if other.tz else _NS_DTYPE
        return DatetimeArray(result, dtype=dtype, freq=self.freq)

    def _addsub_object_array(self, other, op):
        # Add or subtract Array-like of objects
        try:
            # TimedeltaIndex can only operate with a subset of DateOffset
            # subclasses.  Incompatible classes will raise AttributeError,
            # which we re-raise as TypeError
            return super()._addsub_object_array(other, op)
        except AttributeError:
            raise TypeError(
                f"Cannot add/subtract non-tick DateOffset to {type(self).__name__}"
            )

    def __mul__(self, other):
        other = lib.item_from_zerodim(other)

        if isinstance(other, (ABCDataFrame, ABCSeries, ABCIndexClass)):
            return NotImplemented

        if is_scalar(other):
            # numpy will accept float and int, raise TypeError for others
            result = self._data * other
            freq = None
            if self.freq is not None and not isna(other):
                freq = self.freq * other
            return type(self)(result, freq=freq)

        if not hasattr(other, "dtype"):
            # list, tuple
            other = np.array(other)
        if len(other) != len(self) and not is_timedelta64_dtype(other):
            # Exclude timedelta64 here so we correctly raise TypeError
            #  for that instead of ValueError
            raise ValueError("Cannot multiply with unequal lengths")

        if is_object_dtype(other.dtype):
            # this multiplication will succeed only if all elements of other
            #  are int or float scalars, so we will end up with
            #  timedelta64[ns]-dtyped result
            result = [self[n] * other[n] for n in range(len(self))]
            result = np.array(result)
            return type(self)(result)

        # numpy will accept float or int dtype, raise TypeError for others
        result = self._data * other
        return type(self)(result)

    __rmul__ = __mul__

    def __truediv__(self, other):
        # timedelta / X is well-defined for timedelta-like or numeric X
        other = lib.item_from_zerodim(other)

        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)
            if other is NaT:
                # specifically timedelta64-NaT
                result = np.empty(self.shape, dtype=np.float64)
                result.fill(np.nan)
                return result

            # otherwise, dispatch to Timedelta implementation
            return self._data / other

        elif lib.is_scalar(other):
            # assume it is numeric
            result = self._data / other
            freq = None
            if self.freq is not None:
                # Tick division is not implemented, so operate on Timedelta
                freq = self.freq.delta / other
            return type(self)(result, freq=freq)

        if not hasattr(other, "dtype"):
            # e.g. list, tuple
            other = np.array(other)

        if len(other) != len(self):
            raise ValueError("Cannot divide vectors with unequal lengths")

        elif is_timedelta64_dtype(other.dtype):
            # let numpy handle it
            return self._data / other

        elif is_object_dtype(other.dtype):
            # Note: we do not do type inference on the result, so either
            #  an object array or numeric-dtyped (if numpy does inference)
            #  will be returned.  GH#23829
            result = [self[n] / other[n] for n in range(len(self))]
            result = np.array(result)
            return result

        else:
            result = self._data / other
            return type(self)(result)

    def __rtruediv__(self, other):
        # X / timedelta is defined only for timedelta-like X
        other = lib.item_from_zerodim(other)

        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)
            if other is NaT:
                # specifically timedelta64-NaT
                result = np.empty(self.shape, dtype=np.float64)
                result.fill(np.nan)
                return result

            # otherwise, dispatch to Timedelta implementation
            return other / self._data

        elif lib.is_scalar(other):
            raise TypeError(
                f"Cannot divide {type(other).__name__} by {type(self).__name__}"
            )

        if not hasattr(other, "dtype"):
            # e.g. list, tuple
            other = np.array(other)

        if len(other) != len(self):
            raise ValueError("Cannot divide vectors with unequal lengths")

        elif is_timedelta64_dtype(other.dtype):
            # let numpy handle it
            return other / self._data

        elif is_object_dtype(other.dtype):
            # Note: unlike in __truediv__, we do not _need_ to do type
            #  inference on the result.  It does not raise, a numeric array
            #  is returned.  GH#23829
            result = [other[n] / self[n] for n in range(len(self))]
            return np.array(result)

        else:
            raise TypeError(
                f"Cannot divide {other.dtype} data by {type(self).__name__}"
            )

    def __floordiv__(self, other):
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if is_scalar(other):
            if isinstance(other, (timedelta, np.timedelta64, Tick)):
                other = Timedelta(other)
                if other is NaT:
                    # treat this specifically as timedelta-NaT
                    result = np.empty(self.shape, dtype=np.float64)
                    result.fill(np.nan)
                    return result

                # dispatch to Timedelta implementation
                result = other.__rfloordiv__(self._data)
                return result

            # at this point we should only have numeric scalars; anything
            #  else will raise
            result = self.asi8 // other
            result[self._isnan] = iNaT
            freq = None
            if self.freq is not None:
                # Note: freq gets division, not floor-division
                freq = self.freq / other
            return type(self)(result.view("m8[ns]"), freq=freq)

        if not hasattr(other, "dtype"):
            # list, tuple
            other = np.array(other)
        if len(other) != len(self):
            raise ValueError("Cannot divide with unequal lengths")

        elif is_timedelta64_dtype(other.dtype):
            other = type(self)(other)

            # numpy timedelta64 does not natively support floordiv, so operate
            #  on the i8 values
            result = self.asi8 // other.asi8
            mask = self._isnan | other._isnan
            if mask.any():
                result = result.astype(np.int64)
                result[mask] = np.nan
            return result

        elif is_object_dtype(other.dtype):
            result = [self[n] // other[n] for n in range(len(self))]
            result = np.array(result)
            if lib.infer_dtype(result, skipna=False) == "timedelta":
                result, _ = sequence_to_td64ns(result)
                return type(self)(result)
            return result

        elif is_integer_dtype(other.dtype) or is_float_dtype(other.dtype):
            result = self._data // other
            return type(self)(result)

        else:
            dtype = getattr(other, "dtype", type(other).__name__)
            raise TypeError(f"Cannot divide {dtype} by {type(self).__name__}")

    def __rfloordiv__(self, other):
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if is_scalar(other):
            if isinstance(other, (timedelta, np.timedelta64, Tick)):
                other = Timedelta(other)
                if other is NaT:
                    # treat this specifically as timedelta-NaT
                    result = np.empty(self.shape, dtype=np.float64)
                    result.fill(np.nan)
                    return result

                # dispatch to Timedelta implementation
                result = other.__floordiv__(self._data)
                return result

            raise TypeError(
                f"Cannot divide {type(other).__name__} by {type(self).__name__}"
            )

        if not hasattr(other, "dtype"):
            # list, tuple
            other = np.array(other)
        if len(other) != len(self):
            raise ValueError("Cannot divide with unequal lengths")

        elif is_timedelta64_dtype(other.dtype):
            other = type(self)(other)

            # numpy timedelta64 does not natively support floordiv, so operate
            #  on the i8 values
            result = other.asi8 // self.asi8
            mask = self._isnan | other._isnan
            if mask.any():
                result = result.astype(np.int64)
                result[mask] = np.nan
            return result

        elif is_object_dtype(other.dtype):
            result = [other[n] // self[n] for n in range(len(self))]
            result = np.array(result)
            return result

        else:
            dtype = getattr(other, "dtype", type(other).__name__)
            raise TypeError(f"Cannot divide {dtype} by {type(self).__name__}")

    def __mod__(self, other):
        # Note: This is a naive implementation, can likely be optimized
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)
        return self - (self // other) * other

    def __rmod__(self, other):
        # Note: This is a naive implementation, can likely be optimized
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)
        return other - (other // self) * self

    def __divmod__(self, other):
        # Note: This is a naive implementation, can likely be optimized
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)

        res1 = self // other
        res2 = self - res1 * other
        return res1, res2

    def __rdivmod__(self, other):
        # Note: This is a naive implementation, can likely be optimized
        if isinstance(other, (ABCSeries, ABCDataFrame, ABCIndexClass)):
            return NotImplemented

        other = lib.item_from_zerodim(other)
        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            other = Timedelta(other)

        res1 = other // self
        res2 = other - res1 * self
        return res1, res2

    def __neg__(self):
        if self.freq is not None:
            return type(self)(-self._data, freq=-self.freq)
        return type(self)(-self._data)

    def __pos__(self):
        return type(self)(self._data, freq=self.freq)

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

    days = _field_accessor("days", "days", "Number of days for each element.")
    seconds = _field_accessor(
        "seconds",
        "seconds",
        "Number of seconds (>= 0 and less than 1 day) for each element.",
    )
    microseconds = _field_accessor(
        "microseconds",
        "microseconds",
        "Number of microseconds (>= 0 and less than 1 second) for each element.",
    )
    nanoseconds = _field_accessor(
        "nanoseconds",
        "nanoseconds",
        "Number of nanoseconds (>= 0 and less than 1 microsecond) for each element.",
    )

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

        columns = [
            "days",
            "hours",
            "minutes",
            "seconds",
            "milliseconds",
            "microseconds",
            "nanoseconds",
        ]
        hasnans = self._hasnans
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
            result = result.astype("int64")
        return result


TimedeltaArray._add_comparison_ops()


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
    if not hasattr(data, "dtype"):
        # e.g. list, tuple
        if np.ndim(data) == 0:
            # i.e. generator
            data = list(data)
        data = np.array(data, copy=False)
    elif isinstance(data, ABCSeries):
        data = data._values
    elif isinstance(data, (ABCTimedeltaIndex, TimedeltaArray)):
        inferred_freq = data.freq
        data = data._data

    # Convert whatever we have into timedelta64[ns] dtype
    if is_object_dtype(data.dtype) or is_string_dtype(data.dtype):
        # no need to make a copy, need to convert if string-dtyped
        data = objects_to_td64ns(data, unit=unit, errors=errors)
        copy = False

    elif is_integer_dtype(data.dtype):
        # treat as multiples of the given unit
        data, copy_made = ints_to_td64ns(data, unit=unit)
        copy = copy and not copy_made

    elif is_float_dtype(data.dtype):
        # cast the unit, multiply base/frace separately
        # to avoid precision issues from float -> int
        mask = np.isnan(data)
        m, p = precision_from_unit(unit)
        base = data.astype(np.int64)
        frac = data - base
        if p:
            frac = np.round(frac, p)
        data = (base * m + (frac * m).astype(np.int64)).view("timedelta64[ns]")
        data[mask] = iNaT
        copy = False

    elif is_timedelta64_dtype(data.dtype):
        if data.dtype != _TD_DTYPE:
            # non-nano unit
            # TODO: watch out for overflows
            data = data.astype(_TD_DTYPE)
            copy = False

    else:
        # This includes datetime64-dtype, see GH#23539, GH#29794
        raise TypeError(f"dtype {data.dtype} cannot be converted to timedelta64[ns]")

    data = np.array(data, copy=copy)

    assert data.dtype == "m8[ns]", data
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
        dtype_str = f"timedelta64[{unit}]"
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

    result = array_to_timedelta64(values, unit=unit, errors=errors)
    return result.view("timedelta64[ns]")


def _validate_td64_dtype(dtype):
    dtype = pandas_dtype(dtype)
    if is_dtype_equal(dtype, np.dtype("timedelta64")):
        # no precision disallowed GH#24806
        msg = (
            "Passing in 'timedelta' dtype with no precision is not allowed. "
            "Please pass in 'timedelta64[ns]' instead."
        )
        raise ValueError(msg)

    if not is_dtype_equal(dtype, _TD_DTYPE):
        raise ValueError(_BAD_DTYPE.format(dtype=dtype))

    return dtype


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
        raise ValueError(
            "at least 'start' or 'end' should be specified if a 'period' is given."
        )

    data = np.arange(b, e, stride, dtype=np.int64)
    return data
