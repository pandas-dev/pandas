from datetime import datetime, timedelta
import operator
from typing import Any, Callable, Optional, Sequence, Tuple, Type, TypeVar, Union, cast
import warnings

import numpy as np

from pandas._libs import algos, lib
from pandas._libs.tslibs import (
    BaseOffset,
    NaT,
    NaTType,
    Period,
    Resolution,
    Tick,
    Timestamp,
    delta_to_nanoseconds,
    iNaT,
    to_offset,
)
from pandas._libs.tslibs.timestamps import (
    RoundTo,
    integer_op_not_supported,
    round_nsint64,
)
from pandas._typing import DatetimeLikeScalar, DtypeObj
from pandas.compat import set_function_name
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError, NullFrequencyError, PerformanceWarning
from pandas.util._decorators import Appender, Substitution
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_datetime_or_timedelta_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_period_dtype,
    is_string_dtype,
    is_timedelta64_dtype,
    is_unsigned_integer_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.inference import is_array_like
from pandas.core.dtypes.missing import is_valid_nat_for_dtype, isna

from pandas.core import missing, nanops, ops
from pandas.core.algorithms import checked_add_with_arr, unique1d, value_counts
from pandas.core.array_algos.transforms import shift
from pandas.core.arrays._mixins import _T, NDArrayBackedExtensionArray
from pandas.core.arrays.base import ExtensionArray, ExtensionOpsMixin
import pandas.core.common as com
from pandas.core.construction import array, extract_array
from pandas.core.indexers import check_array_indexer
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.ops.invalid import invalid_comparison, make_invalid_op

from pandas.tseries import frequencies

DTScalarOrNaT = Union[DatetimeLikeScalar, NaTType]


def _datetimelike_array_cmp(cls, op):
    """
    Wrap comparison operations to convert Timestamp/Timedelta/Period-like to
    boxed scalars/arrays.
    """
    opname = f"__{op.__name__}__"
    nat_result = opname == "__ne__"

    class InvalidComparison(Exception):
        pass

    def _validate_comparison_value(self, other):
        if isinstance(other, str):
            try:
                # GH#18435 strings get a pass from tzawareness compat
                other = self._scalar_from_string(other)
            except ValueError:
                # failed to parse as Timestamp/Timedelta/Period
                raise InvalidComparison(other)

        if isinstance(other, self._recognized_scalars) or other is NaT:
            other = self._scalar_type(other)
            self._check_compatible_with(other)

        elif not is_list_like(other):
            raise InvalidComparison(other)

        elif len(other) != len(self):
            raise ValueError("Lengths must match")

        else:
            try:
                other = self._validate_listlike(other, opname, allow_object=True)
            except TypeError as err:
                raise InvalidComparison(other) from err

        return other

    @unpack_zerodim_and_defer(opname)
    def wrapper(self, other):
        if self.ndim > 1 and getattr(other, "shape", None) == self.shape:
            # TODO: handle 2D-like listlikes
            return op(self.ravel(), other.ravel()).reshape(self.shape)

        try:
            other = _validate_comparison_value(self, other)
        except InvalidComparison:
            return invalid_comparison(self, other, op)

        dtype = getattr(other, "dtype", None)
        if is_object_dtype(dtype):
            # We have to use comp_method_OBJECT_ARRAY instead of numpy
            #  comparison otherwise it would fail to raise when
            #  comparing tz-aware and tz-naive
            with np.errstate(all="ignore"):
                result = ops.comp_method_OBJECT_ARRAY(op, self.astype(object), other)
            return result

        other_i8 = self._unbox(other)
        result = op(self.asi8, other_i8)

        o_mask = isna(other)
        if self._hasnans | np.any(o_mask):
            result[self._isnan | o_mask] = nat_result

        return result

    return set_function_name(wrapper, opname, cls)


class AttributesMixin:
    _data: np.ndarray

    @classmethod
    def _simple_new(cls, values: np.ndarray, **kwargs):
        raise AbstractMethodError(cls)

    @property
    def _scalar_type(self) -> Type[DatetimeLikeScalar]:
        """
        The scalar associated with this datelike

        * PeriodArray : Period
        * DatetimeArray : Timestamp
        * TimedeltaArray : Timedelta
        """
        raise AbstractMethodError(self)

    def _scalar_from_string(self, value: str) -> DTScalarOrNaT:
        """
        Construct a scalar type from a string.

        Parameters
        ----------
        value : str

        Returns
        -------
        Period, Timestamp, or Timedelta, or NaT
            Whatever the type of ``self._scalar_type`` is.

        Notes
        -----
        This should call ``self._check_compatible_with`` before
        unboxing the result.
        """
        raise AbstractMethodError(self)

    def _unbox_scalar(self, value: DTScalarOrNaT) -> int:
        """
        Unbox the integer value of a scalar `value`.

        Parameters
        ----------
        value : Period, Timestamp, Timedelta, or NaT
            Depending on subclass.

        Returns
        -------
        int

        Examples
        --------
        >>> self._unbox_scalar(Timedelta("10s"))  # doctest: +SKIP
        10000000000
        """
        raise AbstractMethodError(self)

    def _check_compatible_with(
        self, other: DTScalarOrNaT, setitem: bool = False
    ) -> None:
        """
        Verify that `self` and `other` are compatible.

        * DatetimeArray verifies that the timezones (if any) match
        * PeriodArray verifies that the freq matches
        * Timedelta has no verification

        In each case, NaT is considered compatible.

        Parameters
        ----------
        other
        setitem : bool, default False
            For __setitem__ we may have stricter compatibility restrictions than
            for comparisons.

        Raises
        ------
        Exception
        """
        raise AbstractMethodError(self)


class DatelikeOps:
    """
    Common ops for DatetimeIndex/PeriodIndex, but not TimedeltaIndex.
    """

    @Substitution(
        URL="https://docs.python.org/3/library/datetime.html"
        "#strftime-and-strptime-behavior"
    )
    def strftime(self, date_format):
        """
        Convert to Index using specified date_format.

        Return an Index of formatted strings specified by date_format, which
        supports the same string format as the python standard library. Details
        of the string format can be found in `python string format
        doc <%(URL)s>`__.

        Parameters
        ----------
        date_format : str
            Date format string (e.g. "%%Y-%%m-%%d").

        Returns
        -------
        ndarray
            NumPy ndarray of formatted strings.

        See Also
        --------
        to_datetime : Convert the given argument to datetime.
        DatetimeIndex.normalize : Return DatetimeIndex with times to midnight.
        DatetimeIndex.round : Round the DatetimeIndex to the specified freq.
        DatetimeIndex.floor : Floor the DatetimeIndex to the specified freq.

        Examples
        --------
        >>> rng = pd.date_range(pd.Timestamp("2018-03-10 09:00"),
        ...                     periods=3, freq='s')
        >>> rng.strftime('%%B %%d, %%Y, %%r')
        Index(['March 10, 2018, 09:00:00 AM', 'March 10, 2018, 09:00:01 AM',
               'March 10, 2018, 09:00:02 AM'],
              dtype='object')
        """
        result = self._format_native_types(date_format=date_format, na_rep=np.nan)
        return result.astype(object)


class TimelikeOps:
    """
    Common ops for TimedeltaIndex/DatetimeIndex, but not PeriodIndex.
    """

    _round_doc = """
        Perform {op} operation on the data to the specified `freq`.

        Parameters
        ----------
        freq : str or Offset
            The frequency level to {op} the index to. Must be a fixed
            frequency like 'S' (second) not 'ME' (month end). See
            :ref:`frequency aliases <timeseries.offset_aliases>` for
            a list of possible `freq` values.
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            Only relevant for DatetimeIndex:

            - 'infer' will attempt to infer fall dst-transition hours based on
              order
            - bool-ndarray where True signifies a DST time, False designates
              a non-DST time (note that this flag is only applicable for
              ambiguous times)
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an AmbiguousTimeError if there are ambiguous
              times.

            .. versionadded:: 0.24.0

        nonexistent : 'shift_forward', 'shift_backward', 'NaT', timedelta, \
default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST.

            - 'shift_forward' will shift the nonexistent time forward to the
              closest existing time
            - 'shift_backward' will shift the nonexistent time backward to the
              closest existing time
            - 'NaT' will return NaT where there are nonexistent times
            - timedelta objects will shift nonexistent times by the timedelta
            - 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

            .. versionadded:: 0.24.0

        Returns
        -------
        DatetimeIndex, TimedeltaIndex, or Series
            Index of the same type for a DatetimeIndex or TimedeltaIndex,
            or a Series with the same index for a Series.

        Raises
        ------
        ValueError if the `freq` cannot be converted.

        Examples
        --------
        **DatetimeIndex**

        >>> rng = pd.date_range('1/1/2018 11:59:00', periods=3, freq='min')
        >>> rng
        DatetimeIndex(['2018-01-01 11:59:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:01:00'],
                      dtype='datetime64[ns]', freq='T')
        """

    _round_example = """>>> rng.round('H')
        DatetimeIndex(['2018-01-01 12:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.round("H")
        0   2018-01-01 12:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 12:00:00
        dtype: datetime64[ns]
        """

    _floor_example = """>>> rng.floor('H')
        DatetimeIndex(['2018-01-01 11:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.floor("H")
        0   2018-01-01 11:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 12:00:00
        dtype: datetime64[ns]
        """

    _ceil_example = """>>> rng.ceil('H')
        DatetimeIndex(['2018-01-01 12:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 13:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.ceil("H")
        0   2018-01-01 12:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 13:00:00
        dtype: datetime64[ns]
        """

    def _round(self, freq, mode, ambiguous, nonexistent):
        # round the local times
        if is_datetime64tz_dtype(self.dtype):
            # operate on naive timestamps, then convert back to aware
            naive = self.tz_localize(None)
            result = naive._round(freq, mode, ambiguous, nonexistent)
            aware = result.tz_localize(
                self.tz, ambiguous=ambiguous, nonexistent=nonexistent
            )
            return aware

        values = self.view("i8")
        result = round_nsint64(values, mode, freq)
        result = self._maybe_mask_results(result, fill_value=NaT)
        return self._simple_new(result, dtype=self.dtype)

    @Appender((_round_doc + _round_example).format(op="round"))
    def round(self, freq, ambiguous="raise", nonexistent="raise"):
        return self._round(freq, RoundTo.NEAREST_HALF_EVEN, ambiguous, nonexistent)

    @Appender((_round_doc + _floor_example).format(op="floor"))
    def floor(self, freq, ambiguous="raise", nonexistent="raise"):
        return self._round(freq, RoundTo.MINUS_INFTY, ambiguous, nonexistent)

    @Appender((_round_doc + _ceil_example).format(op="ceil"))
    def ceil(self, freq, ambiguous="raise", nonexistent="raise"):
        return self._round(freq, RoundTo.PLUS_INFTY, ambiguous, nonexistent)

    def _with_freq(self, freq):
        """
        Helper to get a view on the same data, with a new freq.

        Parameters
        ----------
        freq : DateOffset, None, or "infer"

        Returns
        -------
        Same type as self
        """
        # GH#29843
        if freq is None:
            # Always valid
            pass
        elif len(self) == 0 and isinstance(freq, BaseOffset):
            # Always valid.  In the TimedeltaArray case, we assume this
            #  is a Tick offset.
            pass
        else:
            # As an internal method, we can ensure this assertion always holds
            assert freq == "infer"
            freq = to_offset(self.inferred_freq)

        arr = self.view()
        arr._freq = freq
        return arr


DatetimeLikeArrayT = TypeVar("DatetimeLikeArrayT", bound="DatetimeLikeArrayMixin")


class DatetimeLikeArrayMixin(
    ExtensionOpsMixin, AttributesMixin, NDArrayBackedExtensionArray
):
    """
    Shared Base/Mixin class for DatetimeArray, TimedeltaArray, PeriodArray

    Assumes that __new__/__init__ defines:
        _data
        _freq

    and that the inheriting class has methods:
        _generate_range
    """

    _is_recognized_dtype: Callable[[DtypeObj], bool]
    _recognized_scalars: Tuple[Type, ...]

    # ------------------------------------------------------------------
    # NDArrayBackedExtensionArray compat

    # TODO: make this a cache_readonly; need to get around _index_data
    #  kludge in libreduction
    @property
    def _ndarray(self) -> np.ndarray:
        # NB: A bunch of Interval tests fail if we use ._data
        return self.asi8

    def _from_backing_data(self: _T, arr: np.ndarray) -> _T:
        # Note: we do not retain `freq`
        return type(self)(arr, dtype=self.dtype)  # type: ignore

    # ------------------------------------------------------------------

    @property
    def _box_func(self):
        """
        box function to get object from internal representation
        """
        raise AbstractMethodError(self)

    def _box_values(self, values):
        """
        apply box func to passed values
        """
        return lib.map_infer(values, self._box_func)

    def __iter__(self):
        return (self._box_func(v) for v in self.asi8)

    @property
    def asi8(self) -> np.ndarray:
        """
        Integer representation of the values.

        Returns
        -------
        ndarray
            An ndarray with int64 dtype.
        """
        # do not cache or you'll create a memory leak
        return self._data.view("i8")

    # ----------------------------------------------------------------
    # Rendering Methods

    def _format_native_types(self, na_rep="NaT", date_format=None):
        """
        Helper method for astype when converting to strings.

        Returns
        -------
        ndarray[str]
        """
        raise AbstractMethodError(self)

    def _formatter(self, boxed=False):
        # TODO: Remove Datetime & DatetimeTZ formatters.
        return "'{}'".format

    # ----------------------------------------------------------------
    # Array-Like / EA-Interface Methods

    def __array__(self, dtype=None) -> np.ndarray:
        # used for Timedelta/DatetimeArray, overwritten by PeriodArray
        if is_object_dtype(dtype):
            return np.array(list(self), dtype=object)
        return self._data

    def __getitem__(self, key):
        """
        This getitem defers to the underlying array, which by-definition can
        only handle list-likes, slices, and integer scalars
        """

        if lib.is_integer(key):
            # fast-path
            result = self._data[key]
            if self.ndim == 1:
                return self._box_func(result)
            return self._simple_new(result, dtype=self.dtype)

        if com.is_bool_indexer(key):
            # first convert to boolean, because check_array_indexer doesn't
            # allow object dtype
            if is_object_dtype(key):
                key = np.asarray(key, dtype=bool)

            key = check_array_indexer(self, key)
            key = lib.maybe_booleans_to_slice(key.view(np.uint8))
        elif isinstance(key, list) and len(key) == 1 and isinstance(key[0], slice):
            # see https://github.com/pandas-dev/pandas/issues/31299, need to allow
            # this for now (would otherwise raise in check_array_indexer)
            pass
        else:
            key = check_array_indexer(self, key)

        freq = self._get_getitem_freq(key)
        result = self._data[key]
        if lib.is_scalar(result):
            return self._box_func(result)
        return self._simple_new(result, dtype=self.dtype, freq=freq)

    def _get_getitem_freq(self, key):
        """
        Find the `freq` attribute to assign to the result of a __getitem__ lookup.
        """
        is_period = is_period_dtype(self.dtype)
        if is_period:
            freq = self.freq
        else:
            freq = None
            if isinstance(key, slice):
                if self.freq is not None and key.step is not None:
                    freq = key.step * self.freq
                else:
                    freq = self.freq
            elif key is Ellipsis:
                # GH#21282 indexing with Ellipsis is similar to a full slice,
                #  should preserve `freq` attribute
                freq = self.freq
        return freq

    def __setitem__(
        self,
        key: Union[int, Sequence[int], Sequence[bool], slice],
        value: Union[NaTType, Any, Sequence[Any]],
    ) -> None:
        # I'm fudging the types a bit here. "Any" above really depends
        # on type(self). For PeriodArray, it's Period (or stuff coercible
        # to a period in from_sequence). For DatetimeArray, it's Timestamp...
        # I don't know if mypy can do that, possibly with Generics.
        # https://mypy.readthedocs.io/en/latest/generics.html
        if is_list_like(value):
            is_slice = isinstance(key, slice)

            if lib.is_scalar(key):
                raise ValueError("setting an array element with a sequence.")

            if not is_slice:
                key = cast(Sequence, key)
                if len(key) != len(value) and not com.is_bool_indexer(key):
                    msg = (
                        f"shape mismatch: value array of length '{len(key)}' "
                        "does not match indexing result of length "
                        f"'{len(value)}'."
                    )
                    raise ValueError(msg)
                elif not len(key):
                    return

        value = self._validate_setitem_value(value)
        key = check_array_indexer(self, key)
        self._data[key] = value
        self._maybe_clear_freq()

    def _maybe_clear_freq(self):
        # inplace operations like __setitem__ may invalidate the freq of
        # DatetimeArray and TimedeltaArray
        pass

    def astype(self, dtype, copy=True):
        # Some notes on cases we don't have to handle here in the base class:
        #   1. PeriodArray.astype handles period -> period
        #   2. DatetimeArray.astype handles conversion between tz.
        #   3. DatetimeArray.astype handles datetime -> period
        dtype = pandas_dtype(dtype)

        if is_object_dtype(dtype):
            return self._box_values(self.asi8.ravel()).reshape(self.shape)
        elif is_string_dtype(dtype) and not is_categorical_dtype(dtype):
            if is_extension_array_dtype(dtype):
                arr_cls = dtype.construct_array_type()
                return arr_cls._from_sequence(self, dtype=dtype)
            else:
                return self._format_native_types()
        elif is_integer_dtype(dtype):
            # we deliberately ignore int32 vs. int64 here.
            # See https://github.com/pandas-dev/pandas/issues/24381 for more.
            values = self.asi8

            if is_unsigned_integer_dtype(dtype):
                # Again, we ignore int32 vs. int64
                values = values.view("uint64")

            if copy:
                values = values.copy()
            return values
        elif (
            is_datetime_or_timedelta_dtype(dtype)
            and not is_dtype_equal(self.dtype, dtype)
        ) or is_float_dtype(dtype):
            # disallow conversion between datetime/timedelta,
            # and conversions for any datetimelike to float
            msg = f"Cannot cast {type(self).__name__} to dtype {dtype}"
            raise TypeError(msg)
        elif is_categorical_dtype(dtype):
            arr_cls = dtype.construct_array_type()
            return arr_cls(self, dtype=dtype)
        else:
            return np.asarray(self, dtype=dtype)

    def view(self, dtype=None):
        if dtype is None or dtype is self.dtype:
            return type(self)(self._data, dtype=self.dtype)
        return self._data.view(dtype=dtype)

    # ------------------------------------------------------------------
    # ExtensionArray Interface

    @classmethod
    def _concat_same_type(cls, to_concat, axis: int = 0):

        # do not pass tz to set because tzlocal cannot be hashed
        dtypes = {str(x.dtype) for x in to_concat}
        if len(dtypes) != 1:
            raise ValueError("to_concat must have the same dtype (tz)", dtypes)

        obj = to_concat[0]
        dtype = obj.dtype

        i8values = [x.asi8 for x in to_concat]
        values = np.concatenate(i8values, axis=axis)

        new_freq = None
        if is_period_dtype(dtype):
            new_freq = obj.freq
        elif axis == 0:
            # GH 3232: If the concat result is evenly spaced, we can retain the
            # original frequency
            to_concat = [x for x in to_concat if len(x)]

            if obj.freq is not None and all(x.freq == obj.freq for x in to_concat):
                pairs = zip(to_concat[:-1], to_concat[1:])
                if all(pair[0][-1] + obj.freq == pair[1][0] for pair in pairs):
                    new_freq = obj.freq

        return cls._simple_new(values, dtype=dtype, freq=new_freq)

    def copy(self: DatetimeLikeArrayT) -> DatetimeLikeArrayT:
        values = self.asi8.copy()
        return type(self)._simple_new(values, dtype=self.dtype, freq=self.freq)

    def _values_for_factorize(self):
        return self.asi8, iNaT

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values, dtype=original.dtype)

    def _values_for_argsort(self):
        return self._data

    @Appender(ExtensionArray.shift.__doc__)
    def shift(self, periods=1, fill_value=None, axis=0):

        fill_value = self._validate_shift_value(fill_value)
        new_values = shift(self._data, periods, axis, fill_value)

        return type(self)._simple_new(new_values, dtype=self.dtype)

    # ------------------------------------------------------------------
    # Validation Methods
    # TODO: try to de-duplicate these, ensure identical behavior

    def _validate_fill_value(self, fill_value):
        """
        If a fill_value is passed to `take` convert it to an i8 representation,
        raising ValueError if this is not possible.

        Parameters
        ----------
        fill_value : object

        Returns
        -------
        fill_value : np.int64

        Raises
        ------
        ValueError
        """
        msg = (
            f"'fill_value' should be a {self._scalar_type}. "
            f"Got '{str(fill_value)}'."
        )
        try:
            fill_value = self._validate_scalar(fill_value, msg)
        except TypeError as err:
            raise ValueError(msg) from err
        return self._unbox(fill_value)

    def _validate_shift_value(self, fill_value):
        # TODO(2.0): once this deprecation is enforced, use _validate_fill_value
        if is_valid_nat_for_dtype(fill_value, self.dtype):
            fill_value = NaT
        elif isinstance(fill_value, self._recognized_scalars):
            fill_value = self._scalar_type(fill_value)
        else:
            # only warn if we're not going to raise
            if self._scalar_type is Period and lib.is_integer(fill_value):
                # kludge for #31971 since Period(integer) tries to cast to str
                new_fill = Period._from_ordinal(fill_value, freq=self.dtype.freq)
            else:
                new_fill = self._scalar_type(fill_value)

            # stacklevel here is chosen to be correct when called from
            #  DataFrame.shift or Series.shift
            warnings.warn(
                f"Passing {type(fill_value)} to shift is deprecated and "
                "will raise in a future version, pass "
                f"{self._scalar_type.__name__} instead.",
                FutureWarning,
                stacklevel=8,
            )
            fill_value = new_fill

        return self._unbox(fill_value)

    def _validate_scalar(
        self, value, msg: Optional[str] = None, cast_str: bool = False
    ):
        """
        Validate that the input value can be cast to our scalar_type.

        Parameters
        ----------
        value : object
        msg : str, optional.
            Message to raise in TypeError on invalid input.
            If not provided, `value` is cast to a str and used
            as the message.
        cast_str : bool, default False
            Whether to try to parse string input to scalar_type.

        Returns
        -------
        self._scalar_type or NaT
        """
        if cast_str and isinstance(value, str):
            # NB: Careful about tzawareness
            try:
                value = self._scalar_from_string(value)
            except ValueError as err:
                raise TypeError(msg) from err

        elif is_valid_nat_for_dtype(value, self.dtype):
            # GH#18295
            value = NaT

        elif isinstance(value, self._recognized_scalars):
            value = self._scalar_type(value)  # type: ignore

        else:
            if msg is None:
                msg = str(value)
            raise TypeError(msg)

        return value

    def _validate_listlike(
        self, value, opname: str, cast_str: bool = False, allow_object: bool = False
    ):
        if isinstance(value, type(self)):
            return value

        # Do type inference if necessary up front
        # e.g. we passed PeriodIndex.values and got an ndarray of Periods
        value = array(value)
        value = extract_array(value, extract_numpy=True)

        if cast_str and is_dtype_equal(value.dtype, "string"):
            # We got a StringArray
            try:
                # TODO: Could use from_sequence_of_strings if implemented
                # Note: passing dtype is necessary for PeriodArray tests
                value = type(self)._from_sequence(value, dtype=self.dtype)
            except ValueError:
                pass

        if is_categorical_dtype(value.dtype):
            # e.g. we have a Categorical holding self.dtype
            if is_dtype_equal(value.categories.dtype, self.dtype):
                # TODO: do we need equal dtype or just comparable?
                value = value._internal_get_values()

        if allow_object and is_object_dtype(value.dtype):
            pass

        elif not type(self)._is_recognized_dtype(value.dtype):
            raise TypeError(
                f"{opname} requires compatible dtype or scalar, "
                f"not {type(value).__name__}"
            )

        return value

    def _validate_searchsorted_value(self, value):
        msg = "searchsorted requires compatible dtype or scalar"
        if not is_list_like(value):
            value = self._validate_scalar(value, msg, cast_str=True)
        else:
            # TODO: cast_str?  we accept it for scalar
            value = self._validate_listlike(value, "searchsorted")

        return self._unbox(value)

    def _validate_setitem_value(self, value):
        msg = (
            f"'value' should be a '{self._scalar_type.__name__}', 'NaT', "
            f"or array of those. Got '{type(value).__name__}' instead."
        )
        if is_list_like(value):
            value = self._validate_listlike(value, "setitem", cast_str=True)
        else:
            # TODO: cast_str for consistency?
            value = self._validate_scalar(value, msg, cast_str=False)

        self._check_compatible_with(value, setitem=True)
        return self._unbox(value)

    def _validate_insert_value(self, value):
        msg = f"cannot insert {type(self).__name__} with incompatible label"
        value = self._validate_scalar(value, msg, cast_str=False)

        self._check_compatible_with(value, setitem=True)
        # TODO: if we dont have compat, should we raise or astype(object)?
        #  PeriodIndex does astype(object)
        return value

    def _validate_where_value(self, other):
        msg = f"Where requires matching dtype, not {type(other)}"
        if not is_list_like(other):
            other = self._validate_scalar(other, msg)
        else:
            other = self._validate_listlike(other, "where")
            self._check_compatible_with(other, setitem=True)

        self._check_compatible_with(other, setitem=True)
        return self._unbox(other)

    def _unbox(self, other) -> Union[np.int64, np.ndarray]:
        """
        Unbox either a scalar with _unbox_scalar or an instance of our own type.
        """
        if lib.is_scalar(other):
            other = self._unbox_scalar(other)
        else:
            # same type as self
            self._check_compatible_with(other)
            other = other.view("i8")
        return other

    # ------------------------------------------------------------------
    # Additional array methods
    #  These are not part of the EA API, but we implement them because
    #  pandas assumes they're there.

    def searchsorted(self, value, side="left", sorter=None):
        """
        Find indices where elements should be inserted to maintain order.

        Find the indices into a sorted array `self` such that, if the
        corresponding elements in `value` were inserted before the indices,
        the order of `self` would be preserved.

        Parameters
        ----------
        value : array_like
            Values to insert into `self`.
        side : {'left', 'right'}, optional
            If 'left', the index of the first suitable location found is given.
            If 'right', return the last such index.  If there is no suitable
            index, return either 0 or N (where N is the length of `self`).
        sorter : 1-D array_like, optional
            Optional array of integer indices that sort `self` into ascending
            order. They are typically the result of ``np.argsort``.

        Returns
        -------
        indices : array of ints
            Array of insertion points with the same shape as `value`.
        """
        value = self._validate_searchsorted_value(value)

        # TODO: Use datetime64 semantics for sorting, xref GH#29844
        return self.asi8.searchsorted(value, side=side, sorter=sorter)

    def value_counts(self, dropna=False):
        """
        Return a Series containing counts of unique values.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of NaT values.

        Returns
        -------
        Series
        """
        from pandas import Index, Series

        if dropna:
            values = self[~self.isna()]._data
        else:
            values = self._data

        cls = type(self)

        result = value_counts(values, sort=False, dropna=dropna)
        index = Index(
            cls(result.index.view("i8"), dtype=self.dtype), name=result.index.name
        )
        return Series(result._values, index=index, name=result.name)

    def map(self, mapper):
        # TODO(GH-23179): Add ExtensionArray.map
        # Need to figure out if we want ExtensionArray.map first.
        # If so, then we can refactor IndexOpsMixin._map_values to
        # a standalone function and call from here..
        # Else, just rewrite _map_infer_values to do the right thing.
        from pandas import Index

        return Index(self).map(mapper).array

    # ------------------------------------------------------------------
    # Null Handling

    def isna(self):
        return self._isnan

    @property  # NB: override with cache_readonly in immutable subclasses
    def _isnan(self):
        """
        return if each value is nan
        """
        return self.asi8 == iNaT

    @property  # NB: override with cache_readonly in immutable subclasses
    def _hasnans(self):
        """
        return if I have any nans; enables various perf speedups
        """
        return bool(self._isnan.any())

    def _maybe_mask_results(self, result, fill_value=iNaT, convert=None):
        """
        Parameters
        ----------
        result : a ndarray
        fill_value : object, default iNaT
        convert : str, dtype or None

        Returns
        -------
        result : ndarray with values replace by the fill_value

        mask the result if needed, convert to the provided dtype if its not
        None

        This is an internal routine.
        """
        if self._hasnans:
            if convert:
                result = result.astype(convert)
            if fill_value is None:
                fill_value = np.nan
            result[self._isnan] = fill_value
        return result

    def fillna(self, value=None, method=None, limit=None):
        # TODO(GH-20300): remove this
        # Just overriding to ensure that we avoid an astype(object).
        # Either 20300 or a `_values_for_fillna` would avoid this duplication.
        if isinstance(value, ABCSeries):
            value = value.array

        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError(
                    f"Length of 'value' does not match. Got ({len(value)}) "
                    f" expected {len(self)}"
                )
            value = value[mask]

        if mask.any():
            if method is not None:
                if method == "pad":
                    func = missing.pad_1d
                else:
                    func = missing.backfill_1d

                values = self._data
                if not is_period_dtype(self.dtype):
                    # For PeriodArray self._data is i8, which gets copied
                    #  by `func`.  Otherwise we need to make a copy manually
                    # to avoid modifying `self` in-place.
                    values = values.copy()

                new_values = func(values, limit=limit, mask=mask)
                if is_datetime64tz_dtype(self.dtype):
                    # we need to pass int64 values to the constructor to avoid
                    #  re-localizing incorrectly
                    new_values = new_values.view("i8")
                new_values = type(self)(new_values, dtype=self.dtype)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    # ------------------------------------------------------------------
    # Frequency Properties/Methods

    @property
    def freq(self):
        """
        Return the frequency object if it is set, otherwise None.
        """
        return self._freq

    @freq.setter
    def freq(self, value):
        if value is not None:
            value = to_offset(value)
            self._validate_frequency(self, value)

        self._freq = value

    @property
    def freqstr(self):
        """
        Return the frequency object as a string if its set, otherwise None.
        """
        if self.freq is None:
            return None
        return self.freq.freqstr

    @property  # NB: override with cache_readonly in immutable subclasses
    def inferred_freq(self):
        """
        Tries to return a string representing a frequency guess,
        generated by infer_freq.  Returns None if it can't autodetect the
        frequency.
        """
        if self.ndim != 1:
            return None
        try:
            return frequencies.infer_freq(self)
        except ValueError:
            return None

    @property  # NB: override with cache_readonly in immutable subclasses
    def _resolution_obj(self) -> Optional[Resolution]:
        try:
            return Resolution.get_reso_from_freq(self.freqstr)
        except KeyError:
            return None

    @property  # NB: override with cache_readonly in immutable subclasses
    def resolution(self) -> str:
        """
        Returns day, hour, minute, second, millisecond or microsecond
        """
        return self._resolution_obj.attrname  # type: ignore

    @classmethod
    def _validate_frequency(cls, index, freq, **kwargs):
        """
        Validate that a frequency is compatible with the values of a given
        Datetime Array/Index or Timedelta Array/Index

        Parameters
        ----------
        index : DatetimeIndex or TimedeltaIndex
            The index on which to determine if the given frequency is valid
        freq : DateOffset
            The frequency to validate
        """
        # TODO: this is not applicable to PeriodArray, move to correct Mixin
        inferred = index.inferred_freq
        if index.size == 0 or inferred == freq.freqstr:
            return None

        try:
            on_freq = cls._generate_range(
                start=index[0], end=None, periods=len(index), freq=freq, **kwargs
            )
            if not np.array_equal(index.asi8, on_freq.asi8):
                raise ValueError
        except ValueError as e:
            if "non-fixed" in str(e):
                # non-fixed frequencies are not meaningful for timedelta64;
                #  we retain that error message
                raise e
            # GH#11587 the main way this is reached is if the `np.array_equal`
            #  check above is False.  This can also be reached if index[0]
            #  is `NaT`, in which case the call to `cls._generate_range` will
            #  raise a ValueError, which we re-raise with a more targeted
            #  message.
            raise ValueError(
                f"Inferred frequency {inferred} from passed values "
                f"does not conform to passed frequency {freq.freqstr}"
            ) from e

    # monotonicity/uniqueness properties are called via frequencies.infer_freq,
    #  see GH#23789

    @property
    def _is_monotonic_increasing(self):
        return algos.is_monotonic(self.asi8, timelike=True)[0]

    @property
    def _is_monotonic_decreasing(self):
        return algos.is_monotonic(self.asi8, timelike=True)[1]

    @property
    def _is_unique(self):
        return len(unique1d(self.asi8)) == len(self)

    # ------------------------------------------------------------------
    # Arithmetic Methods
    _create_comparison_method = classmethod(_datetimelike_array_cmp)

    # pow is invalid for all three subclasses; TimedeltaArray will override
    #  the multiplication and division ops
    __pow__ = make_invalid_op("__pow__")
    __rpow__ = make_invalid_op("__rpow__")
    __mul__ = make_invalid_op("__mul__")
    __rmul__ = make_invalid_op("__rmul__")
    __truediv__ = make_invalid_op("__truediv__")
    __rtruediv__ = make_invalid_op("__rtruediv__")
    __floordiv__ = make_invalid_op("__floordiv__")
    __rfloordiv__ = make_invalid_op("__rfloordiv__")
    __mod__ = make_invalid_op("__mod__")
    __rmod__ = make_invalid_op("__rmod__")
    __divmod__ = make_invalid_op("__divmod__")
    __rdivmod__ = make_invalid_op("__rdivmod__")

    def _add_datetimelike_scalar(self, other):
        # Overridden by TimedeltaArray
        raise TypeError(f"cannot add {type(self).__name__} and {type(other).__name__}")

    _add_datetime_arraylike = _add_datetimelike_scalar

    def _sub_datetimelike_scalar(self, other):
        # Overridden by DatetimeArray
        assert other is not NaT
        raise TypeError(f"cannot subtract a datelike from a {type(self).__name__}")

    _sub_datetime_arraylike = _sub_datetimelike_scalar

    def _sub_period(self, other):
        # Overridden by PeriodArray
        raise TypeError(f"cannot subtract Period from a {type(self).__name__}")

    def _add_period(self, other: Period):
        # Overriden by TimedeltaArray
        raise TypeError(f"cannot add Period to a {type(self).__name__}")

    def _add_offset(self, offset):
        raise AbstractMethodError(self)

    def _add_timedeltalike_scalar(self, other):
        """
        Add a delta of a timedeltalike

        Returns
        -------
        Same type as self
        """
        if isna(other):
            # i.e np.timedelta64("NaT"), not recognized by delta_to_nanoseconds
            new_values = np.empty(self.shape, dtype="i8")
            new_values[:] = iNaT
            return type(self)(new_values, dtype=self.dtype)

        inc = delta_to_nanoseconds(other)
        new_values = checked_add_with_arr(self.asi8, inc, arr_mask=self._isnan).view(
            "i8"
        )
        new_values = self._maybe_mask_results(new_values)

        new_freq = None
        if isinstance(self.freq, Tick) or is_period_dtype(self.dtype):
            # adding a scalar preserves freq
            new_freq = self.freq

        return type(self)._simple_new(new_values, dtype=self.dtype, freq=new_freq)

    def _add_timedelta_arraylike(self, other):
        """
        Add a delta of a TimedeltaIndex

        Returns
        -------
        Same type as self
        """
        # overridden by PeriodArray

        if len(self) != len(other):
            raise ValueError("cannot add indices of unequal length")

        if isinstance(other, np.ndarray):
            # ndarray[timedelta64]; wrap in TimedeltaIndex for op
            from pandas.core.arrays import TimedeltaArray

            other = TimedeltaArray._from_sequence(other)

        self_i8 = self.asi8
        other_i8 = other.asi8
        new_values = checked_add_with_arr(
            self_i8, other_i8, arr_mask=self._isnan, b_mask=other._isnan
        )
        if self._hasnans or other._hasnans:
            mask = (self._isnan) | (other._isnan)
            new_values[mask] = iNaT

        return type(self)(new_values, dtype=self.dtype)

    def _add_nat(self):
        """
        Add pd.NaT to self
        """
        if is_period_dtype(self.dtype):
            raise TypeError(
                f"Cannot add {type(self).__name__} and {type(NaT).__name__}"
            )

        # GH#19124 pd.NaT is treated like a timedelta for both timedelta
        # and datetime dtypes
        result = np.zeros(self.shape, dtype=np.int64)
        result.fill(iNaT)
        return type(self)(result, dtype=self.dtype, freq=None)

    def _sub_nat(self):
        """
        Subtract pd.NaT from self
        """
        # GH#19124 Timedelta - datetime is not in general well-defined.
        # We make an exception for pd.NaT, which in this case quacks
        # like a timedelta.
        # For datetime64 dtypes by convention we treat NaT as a datetime, so
        # this subtraction returns a timedelta64 dtype.
        # For period dtype, timedelta64 is a close-enough return dtype.
        result = np.zeros(self.shape, dtype=np.int64)
        result.fill(iNaT)
        return result.view("timedelta64[ns]")

    def _sub_period_array(self, other):
        # Overridden by PeriodArray
        raise TypeError(
            f"cannot subtract {other.dtype}-dtype from {type(self).__name__}"
        )

    def _addsub_object_array(self, other: np.ndarray, op):
        """
        Add or subtract array-like of DateOffset objects

        Parameters
        ----------
        other : np.ndarray[object]
        op : {operator.add, operator.sub}

        Returns
        -------
        result : same class as self
        """
        assert op in [operator.add, operator.sub]
        if len(other) == 1:
            # If both 1D then broadcasting is unambiguous
            # TODO(EA2D): require self.ndim == other.ndim here
            return op(self, other[0])

        warnings.warn(
            "Adding/subtracting object-dtype array to "
            f"{type(self).__name__} not vectorized",
            PerformanceWarning,
        )

        # Caller is responsible for broadcasting if necessary
        assert self.shape == other.shape, (self.shape, other.shape)

        res_values = op(self.astype("O"), np.asarray(other))
        result = array(res_values.ravel())
        result = extract_array(result, extract_numpy=True).reshape(self.shape)
        return result

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
        if freq is not None and freq != self.freq:
            if isinstance(freq, str):
                freq = to_offset(freq)
            offset = periods * freq
            result = self + offset
            return result

        if periods == 0:
            # immutable so OK
            return self.copy()

        if self.freq is None:
            raise NullFrequencyError("Cannot shift with no freq")

        start = self[0] + periods * self.freq
        end = self[-1] + periods * self.freq

        # Note: in the DatetimeTZ case, _generate_range will infer the
        #  appropriate timezone from `start` and `end`, so tz does not need
        #  to be passed explicitly.
        return self._generate_range(start=start, end=end, periods=None, freq=self.freq)

    @unpack_zerodim_and_defer("__add__")
    def __add__(self, other):
        other_dtype = getattr(other, "dtype", None)

        # scalar others
        if other is NaT:
            result = self._add_nat()
        elif isinstance(other, (Tick, timedelta, np.timedelta64)):
            result = self._add_timedeltalike_scalar(other)
        elif isinstance(other, BaseOffset):
            # specifically _not_ a Tick
            result = self._add_offset(other)
        elif isinstance(other, (datetime, np.datetime64)):
            result = self._add_datetimelike_scalar(other)
        elif isinstance(other, Period) and is_timedelta64_dtype(self.dtype):
            result = self._add_period(other)
        elif lib.is_integer(other):
            # This check must come after the check for np.timedelta64
            # as is_integer returns True for these
            if not is_period_dtype(self.dtype):
                raise integer_op_not_supported(self)
            result = self._time_shift(other)

        # array-like others
        elif is_timedelta64_dtype(other_dtype):
            # TimedeltaIndex, ndarray[timedelta64]
            result = self._add_timedelta_arraylike(other)
        elif is_object_dtype(other_dtype):
            # e.g. Array/Index of DateOffset objects
            result = self._addsub_object_array(other, operator.add)
        elif is_datetime64_dtype(other_dtype) or is_datetime64tz_dtype(other_dtype):
            # DatetimeIndex, ndarray[datetime64]
            return self._add_datetime_arraylike(other)
        elif is_integer_dtype(other_dtype):
            if not is_period_dtype(self.dtype):
                raise integer_op_not_supported(self)
            result = self._addsub_int_array(other, operator.add)
        else:
            # Includes Categorical, other ExtensionArrays
            # For PeriodDtype, if self is a TimedeltaArray and other is a
            #  PeriodArray with  a timedelta-like (i.e. Tick) freq, this
            #  operation is valid.  Defer to the PeriodArray implementation.
            #  In remaining cases, this will end up raising TypeError.
            return NotImplemented

        if isinstance(result, np.ndarray) and is_timedelta64_dtype(result.dtype):
            from pandas.core.arrays import TimedeltaArray

            return TimedeltaArray(result)
        return result

    def __radd__(self, other):
        # alias for __add__
        return self.__add__(other)

    @unpack_zerodim_and_defer("__sub__")
    def __sub__(self, other):

        other_dtype = getattr(other, "dtype", None)

        # scalar others
        if other is NaT:
            result = self._sub_nat()
        elif isinstance(other, (Tick, timedelta, np.timedelta64)):
            result = self._add_timedeltalike_scalar(-other)
        elif isinstance(other, BaseOffset):
            # specifically _not_ a Tick
            result = self._add_offset(-other)
        elif isinstance(other, (datetime, np.datetime64)):
            result = self._sub_datetimelike_scalar(other)
        elif lib.is_integer(other):
            # This check must come after the check for np.timedelta64
            # as is_integer returns True for these
            if not is_period_dtype(self.dtype):
                raise integer_op_not_supported(self)
            result = self._time_shift(-other)

        elif isinstance(other, Period):
            result = self._sub_period(other)

        # array-like others
        elif is_timedelta64_dtype(other_dtype):
            # TimedeltaIndex, ndarray[timedelta64]
            result = self._add_timedelta_arraylike(-other)
        elif is_object_dtype(other_dtype):
            # e.g. Array/Index of DateOffset objects
            result = self._addsub_object_array(other, operator.sub)
        elif is_datetime64_dtype(other_dtype) or is_datetime64tz_dtype(other_dtype):
            # DatetimeIndex, ndarray[datetime64]
            result = self._sub_datetime_arraylike(other)
        elif is_period_dtype(other_dtype):
            # PeriodIndex
            result = self._sub_period_array(other)
        elif is_integer_dtype(other_dtype):
            if not is_period_dtype(self.dtype):
                raise integer_op_not_supported(self)
            result = self._addsub_int_array(other, operator.sub)
        else:
            # Includes ExtensionArrays, float_dtype
            return NotImplemented

        if isinstance(result, np.ndarray) and is_timedelta64_dtype(result.dtype):
            from pandas.core.arrays import TimedeltaArray

            return TimedeltaArray(result)
        return result

    def __rsub__(self, other):
        other_dtype = getattr(other, "dtype", None)

        if is_datetime64_any_dtype(other_dtype) and is_timedelta64_dtype(self.dtype):
            # ndarray[datetime64] cannot be subtracted from self, so
            # we need to wrap in DatetimeArray/Index and flip the operation
            if lib.is_scalar(other):
                # i.e. np.datetime64 object
                return Timestamp(other) - self
            if not isinstance(other, DatetimeLikeArrayMixin):
                # Avoid down-casting DatetimeIndex
                from pandas.core.arrays import DatetimeArray

                other = DatetimeArray(other)
            return other - self
        elif (
            is_datetime64_any_dtype(self.dtype)
            and hasattr(other, "dtype")
            and not is_datetime64_any_dtype(other.dtype)
        ):
            # GH#19959 datetime - datetime is well-defined as timedelta,
            # but any other type - datetime is not well-defined.
            raise TypeError(
                f"cannot subtract {type(self).__name__} from {type(other).__name__}"
            )
        elif is_period_dtype(self.dtype) and is_timedelta64_dtype(other_dtype):
            # TODO: Can we simplify/generalize these cases at all?
            raise TypeError(f"cannot subtract {type(self).__name__} from {other.dtype}")
        elif is_timedelta64_dtype(self.dtype):
            return (-self) + other

        # We get here with e.g. datetime objects
        return -(self - other)

    def __iadd__(self, other):
        result = self + other
        self[:] = result[:]

        if not is_period_dtype(self.dtype):
            # restore freq, which is invalidated by setitem
            self._freq = result._freq
        return self

    def __isub__(self, other):
        result = self - other
        self[:] = result[:]

        if not is_period_dtype(self.dtype):
            # restore freq, which is invalidated by setitem
            self._freq = result._freq
        return self

    # --------------------------------------------------------------
    # Reductions

    def _reduce(self, name: str, skipna: bool = True, **kwargs):
        op = getattr(self, name, None)
        if op:
            return op(skipna=skipna, **kwargs)
        else:
            return super()._reduce(name, skipna, **kwargs)

    def min(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the minimum value of the Array or minimum along
        an axis.

        See Also
        --------
        numpy.ndarray.min
        Index.min : Return the minimum value in an Index.
        Series.min : Return the minimum value in a Series.
        """
        nv.validate_min(args, kwargs)
        nv.validate_minmax_axis(axis)

        result = nanops.nanmin(self.asi8, skipna=skipna, mask=self.isna())
        if isna(result):
            # Period._from_ordinal does not handle np.nan gracefully
            return NaT
        return self._box_func(result)

    def max(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the maximum value of the Array or maximum along
        an axis.

        See Also
        --------
        numpy.ndarray.max
        Index.max : Return the maximum value in an Index.
        Series.max : Return the maximum value in a Series.
        """
        # TODO: skipna is broken with max.
        # See https://github.com/pandas-dev/pandas/issues/24265
        nv.validate_max(args, kwargs)
        nv.validate_minmax_axis(axis)

        mask = self.isna()
        if skipna:
            values = self[~mask].asi8
        elif mask.any():
            return NaT
        else:
            values = self.asi8

        if not len(values):
            # short-circuit for empty max / min
            return NaT

        result = nanops.nanmax(values, skipna=skipna)
        # Don't have to worry about NA `result`, since no NA went in.
        return self._box_func(result)

    def mean(self, skipna=True):
        """
        Return the mean value of the Array.

        .. versionadded:: 0.25.0

        Parameters
        ----------
        skipna : bool, default True
            Whether to ignore any NaT elements.

        Returns
        -------
        scalar
            Timestamp or Timedelta.

        See Also
        --------
        numpy.ndarray.mean : Returns the average of array elements along a given axis.
        Series.mean : Return the mean value in a Series.

        Notes
        -----
        mean is only defined for Datetime and Timedelta dtypes, not for Period.
        """
        if is_period_dtype(self.dtype):
            # See discussion in GH#24757
            raise TypeError(
                f"mean is not implemented for {type(self).__name__} since the "
                "meaning is ambiguous.  An alternative is "
                "obj.to_timestamp(how='start').mean()"
            )

        mask = self.isna()
        if skipna:
            values = self[~mask]
        elif mask.any():
            return NaT
        else:
            values = self

        if not len(values):
            # short-circuit for empty max / min
            return NaT

        result = nanops.nanmean(values.view("i8"), skipna=skipna)
        # Don't have to worry about NA `result`, since no NA went in.
        return self._box_func(result)

    # --------------------------------------------------------------

    def factorize(self, na_sentinel=-1, sort: bool = False):
        if self.freq is not None:
            # We must be unique, so can short-circuit (and retain freq)
            codes = np.arange(len(self), dtype=np.intp)
            uniques = self.copy()  # TODO: copy or view?
            if sort and self.freq.n < 0:
                codes = codes[::-1]
                uniques = uniques[::-1]
            return codes, uniques
        # FIXME: shouldn't get here; we are ignoring sort
        return super().factorize(na_sentinel=na_sentinel)


DatetimeLikeArrayMixin._add_comparison_ops()

# -------------------------------------------------------------------
# Shared Constructor Helpers


def validate_periods(periods):
    """
    If a `periods` argument is passed to the Datetime/Timedelta Array/Index
    constructor, cast it to an integer.

    Parameters
    ----------
    periods : None, float, int

    Returns
    -------
    periods : None or int

    Raises
    ------
    TypeError
        if periods is None, float, or int
    """
    if periods is not None:
        if lib.is_float(periods):
            periods = int(periods)
        elif not lib.is_integer(periods):
            raise TypeError(f"periods must be a number, got {periods}")
    return periods


def validate_endpoints(closed):
    """
    Check that the `closed` argument is among [None, "left", "right"]

    Parameters
    ----------
    closed : {None, "left", "right"}

    Returns
    -------
    left_closed : bool
    right_closed : bool

    Raises
    ------
    ValueError : if argument is not among valid values
    """
    left_closed = False
    right_closed = False

    if closed is None:
        left_closed = True
        right_closed = True
    elif closed == "left":
        left_closed = True
    elif closed == "right":
        right_closed = True
    else:
        raise ValueError("Closed has to be either 'left', 'right' or None")

    return left_closed, right_closed


def validate_inferred_freq(freq, inferred_freq, freq_infer):
    """
    If the user passes a freq and another freq is inferred from passed data,
    require that they match.

    Parameters
    ----------
    freq : DateOffset or None
    inferred_freq : DateOffset or None
    freq_infer : bool

    Returns
    -------
    freq : DateOffset or None
    freq_infer : bool

    Notes
    -----
    We assume at this point that `maybe_infer_freq` has been called, so
    `freq` is either a DateOffset object or None.
    """
    if inferred_freq is not None:
        if freq is not None and freq != inferred_freq:
            raise ValueError(
                f"Inferred frequency {inferred_freq} from passed "
                "values does not conform to passed frequency "
                f"{freq.freqstr}"
            )
        elif freq is None:
            freq = inferred_freq
        freq_infer = False

    return freq, freq_infer


def maybe_infer_freq(freq):
    """
    Comparing a DateOffset to the string "infer" raises, so we need to
    be careful about comparisons.  Make a dummy variable `freq_infer` to
    signify the case where the given freq is "infer" and set freq to None
    to avoid comparison trouble later on.

    Parameters
    ----------
    freq : {DateOffset, None, str}

    Returns
    -------
    freq : {DateOffset, None}
    freq_infer : bool
        Whether we should inherit the freq of passed data.
    """
    freq_infer = False
    if not isinstance(freq, BaseOffset):
        # if a passed freq is None, don't infer automatically
        if freq != "infer":
            freq = to_offset(freq)
        else:
            freq_infer = True
            freq = None
    return freq, freq_infer
