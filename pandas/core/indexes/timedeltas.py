"""implement the TimedeltaIndex"""

from __future__ import annotations

import operator
from typing import (
    TYPE_CHECKING,
    Self,
    cast,
)
import warnings

from pandas._libs import (
    index as libindex,
    lib,
)
from pandas._libs.tslibs import (
    Day,
    Resolution,
    Tick,
    Timedelta,
    Timestamp,
    timezones,
    to_offset,
)
from pandas._libs.tslibs.dtypes import abbrev_to_npy_unit
from pandas._libs.tslibs.timedeltas import parse_timedelta_string_reso
from pandas.util._decorators import set_module

from pandas.core.dtypes.common import (
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import ArrowDtype
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core.arrays.timedeltas import TimedeltaArray
import pandas.core.common as com
from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)
from pandas.core.indexes.datetimelike import DatetimeTimedeltaMixin
from pandas.core.roperator import (
    rmul,
    rsub,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    import numpy as np

    from pandas._libs import NaTType
    from pandas._typing import (
        AxisInt,
        DtypeObj,
        NpDtype,
        TimeAmbiguous,
        TimeNonexistent,
        TimeUnit,
        npt,
    )

    from pandas import DataFrame


def _new_TimedeltaIndex(cls, d):
    """
    This is called upon unpickling, rather than the default which doesn't
    have arguments and breaks __new__
    """
    if "data" in d and not isinstance(d["data"], TimedeltaIndex):
        data = d.pop("data")
        if isinstance(data, TimedeltaArray):
            tdarr = data
        else:
            tdarr = TimedeltaArray._simple_new(data, dtype=data.dtype)
        # Legacy pickles stored freq on the TimedeltaArray; current pickles
        # include it in ``d``. Migrate either up onto the Index.
        legacy_freq = vars(tdarr).pop("_freq", None)
        freq = d.pop("freq", legacy_freq)
        result = cls._simple_new(tdarr, **d)
        result._freq = freq
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cls.__new__(cls, **d)
    return result


@set_module("pandas")
class TimedeltaIndex(DatetimeTimedeltaMixin):
    """
    Immutable Index of timedelta64 data.

    Represented internally as int64, and scalars returned Timedelta objects.

    Parameters
    ----------
    data : array-like (1-dimensional), optional
        Optional timedelta-like data to construct index with.
    freq : str or pandas offset object, optional
        One of pandas date offset strings or corresponding objects. The string
        ``'infer'`` can be passed in order to set the frequency of the index as
        the inferred frequency upon creation.
    dtype : numpy.dtype or str, default None
        Valid ``numpy`` dtypes are ``timedelta64[ns]``, ``timedelta64[us]``,
        ``timedelta64[ms]``, and ``timedelta64[s]``.
    copy : bool, default None
        Whether to copy input data, only relevant for array, Series, and Index
        inputs (for other input, e.g. a list, a new array is created anyway).
        Defaults to True for array input and False for Index/Series.
        Set to False to avoid copying array input at your own risk (if you
        know the input data won't be modified elsewhere).
        Set to True to force copying Series/Index input up front.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    days
    seconds
    microseconds
    nanoseconds
    components
    inferred_freq

    Methods
    -------
    to_pytimedelta
    to_series
    round
    floor
    ceil
    to_frame
    mean

    See Also
    --------
    Index : The base pandas Index type.
    Timedelta : Represents a duration between two dates or times.
    DatetimeIndex : Index of datetime64 data.
    PeriodIndex : Index of Period data.
    timedelta_range : Create a fixed-frequency TimedeltaIndex.

    Notes
    -----
    To learn more about the frequency strings, please see
    :ref:`this link<timeseries.offset_aliases>`.

    Examples
    --------
    >>> pd.TimedeltaIndex(["0 days", "1 days", "2 days", "3 days", "4 days"])
    TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq=None)

    We can also let pandas infer the frequency when possible.

    >>> pd.TimedeltaIndex(np.arange(5) * 24 * 3600 * 1e9, freq="infer")
    TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')
    """

    _typ = "timedeltaindex"

    _data_cls = TimedeltaArray

    @property
    def _engine_type(self) -> type[libindex.TimedeltaEngine]:
        return libindex.TimedeltaEngine

    _data: TimedeltaArray

    # Use base class method instead of DatetimeTimedeltaMixin._get_string_slice
    _get_string_slice = Index._get_string_slice

    # -------------------------------------------------------------------
    # Methods that dispatch to TimedeltaArray and wrap the result

    def _wrap_field_result(self, result: np.ndarray) -> Index:
        """Wrap an ndarray result computed from ``self._data`` in an Index."""
        return Index(result, name=self.name, dtype=result.dtype, copy=False)

    def _wrap_td_result(self, result: TimedeltaArray) -> Self:
        """Wrap a TimedeltaArray result computed from ``self._data``."""
        return type(self)._simple_new(result, name=self.name)

    def total_seconds(self) -> Index:
        """
        Return total duration of each element expressed in seconds.

        This method is available directly on TimedeltaArray, TimedeltaIndex
        and on Series containing timedelta values under the ``.dt`` namespace.

        Returns
        -------
        Index
            An Index with a float64 dtype.

        See Also
        --------
        datetime.timedelta.total_seconds : Standard library version
            of this method.
        TimedeltaIndex.components : Return a DataFrame with components of
            each Timedelta.

        Examples
        --------
        >>> idx = pd.to_timedelta(np.arange(5), unit="D")
        >>> idx
        TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                       dtype='timedelta64[us]', freq=None)

        >>> idx.total_seconds()
        Index([0.0, 86400.0, 172800.0, 259200.0, 345600.0], dtype='float64')
        """
        return self._wrap_field_result(self._data.total_seconds())

    # error: Signature of "round" incompatible with supertype "Index"
    # TimedeltaIndex.round rounds to a freq like Index.floor/ceil, unlike
    # Index.round which rounds numeric values to a number of decimals.
    def round(  # type: ignore[override]
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ) -> Self:
        """
        Perform round operation on the data to the specified `freq`.

        This method rounds each timedelta value in the Series/Index to the
        nearest specified frequency using standard rounding rules (round half
        to even).

        Parameters
        ----------
        freq : str or Offset
            The frequency level to round the index to. Must be a fixed
            frequency like 's' (second) not 'ME' (month end). See
            :ref:`frequency aliases <timeseries.offset_aliases>` for
            a list of possible `freq` values.
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.round`;
            has no effect since TimedeltaIndex values are never tz-aware.
        nonexistent : 'shift_forward', 'shift_backward', 'NaT', timedelta, \
            default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.round`;
            has no effect since TimedeltaIndex values are never tz-aware.

        Returns
        -------
        TimedeltaIndex
            Values rounded to the given `freq`.

        Raises
        ------
        ValueError if the `freq` cannot be converted.

        See Also
        --------
        TimedeltaIndex.floor : Perform floor operation on the data to the
            specified `freq`.
        TimedeltaIndex.ceil : Perform ceil operation on the data to the
            specified `freq`.

        Examples
        --------
        >>> tdelta_idx = pd.timedelta_range("1 day", periods=3, freq="6h20min")
        >>> tdelta_idx
        TimedeltaIndex(['1 days 00:00:00', '1 days 06:20:00', '1 days 12:40:00'],
                       dtype='timedelta64[us]', freq='380min')
        >>> tdelta_idx.round("h")
        TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 13:00:00'],
                       dtype='timedelta64[us]', freq=None)
        """
        return self._wrap_td_result(self._data.round(freq, ambiguous, nonexistent))

    def floor(
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ) -> Self:
        """
        Perform floor operation on the data to the specified `freq`.

        This method rounds each timedelta value in the Series/Index down to
        the specified frequency (i.e., towards negative infinity).

        Parameters
        ----------
        freq : str or Offset
            The frequency level to floor the index to. Must be a fixed
            frequency like 's' (second) not 'ME' (month end). See
            :ref:`frequency aliases <timeseries.offset_aliases>` for
            a list of possible `freq` values.
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.floor`;
            has no effect since TimedeltaIndex values are never tz-aware.
        nonexistent : 'shift_forward', 'shift_backward', 'NaT', timedelta, \
            default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.floor`;
            has no effect since TimedeltaIndex values are never tz-aware.

        Returns
        -------
        TimedeltaIndex
            Values rounded down to the given `freq`.

        Raises
        ------
        ValueError if the `freq` cannot be converted.

        See Also
        --------
        TimedeltaIndex.round : Perform round operation on the data to the
            specified `freq`.
        TimedeltaIndex.ceil : Perform ceil operation on the data to the
            specified `freq`.

        Examples
        --------
        >>> tdelta_idx = pd.timedelta_range("1 day", periods=3, freq="6h20min")
        >>> tdelta_idx
        TimedeltaIndex(['1 days 00:00:00', '1 days 06:20:00', '1 days 12:40:00'],
                       dtype='timedelta64[us]', freq='380min')
        >>> tdelta_idx.floor("h")
        TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 12:00:00'],
                       dtype='timedelta64[us]', freq=None)
        """
        return self._wrap_td_result(self._data.floor(freq, ambiguous, nonexistent))

    def ceil(
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ) -> Self:
        """
        Perform ceil operation on the data to the specified `freq`.

        This method rounds each timedelta value in the Series/Index up to
        the specified frequency (i.e., towards positive infinity).

        Parameters
        ----------
        freq : str or Offset
            The frequency level to ceil the index to. Must be a fixed
            frequency like 's' (second) not 'ME' (month end). See
            :ref:`frequency aliases <timeseries.offset_aliases>` for
            a list of possible `freq` values.
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.ceil`;
            has no effect since TimedeltaIndex values are never tz-aware.
        nonexistent : 'shift_forward', 'shift_backward', 'NaT', timedelta, \
            default 'raise'
            Accepted for compatibility with :meth:`DatetimeIndex.ceil`;
            has no effect since TimedeltaIndex values are never tz-aware.

        Returns
        -------
        TimedeltaIndex
            Values rounded up to the given `freq`.

        Raises
        ------
        ValueError if the `freq` cannot be converted.

        See Also
        --------
        TimedeltaIndex.round : Perform round operation on the data to the
            specified `freq`.
        TimedeltaIndex.floor : Perform floor operation on the data to the
            specified `freq`.

        Examples
        --------
        >>> tdelta_idx = pd.timedelta_range("1 day", periods=3, freq="6h20min")
        >>> tdelta_idx
        TimedeltaIndex(['1 days 00:00:00', '1 days 06:20:00', '1 days 12:40:00'],
                       dtype='timedelta64[us]', freq='380min')
        >>> tdelta_idx.ceil("h")
        TimedeltaIndex(['1 days 00:00:00', '1 days 07:00:00', '1 days 13:00:00'],
                       dtype='timedelta64[us]', freq=None)
        """
        return self._wrap_td_result(self._data.ceil(freq, ambiguous, nonexistent))

    @property
    def days(self) -> Index:
        """
        Number of days for each element.

        This attribute returns the number of whole days in each timedelta value.
        It represents the days component of the duration, not the total duration
        expressed in days.

        See Also
        --------
        Series.dt.seconds : Return number of seconds for each element.
        Series.dt.microseconds : Return number of microseconds for each element.
        Series.dt.nanoseconds : Return number of nanoseconds for each element.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta(["0 days", "10 days", "20 days"])
        >>> tdelta_idx
        TimedeltaIndex(['0 days', '10 days', '20 days'],
                        dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.days
        Index([0, 10, 20], dtype='int64')
        """
        return self._wrap_field_result(self._data.days)

    @property
    def seconds(self) -> Index:
        """
        Number of seconds (>= 0 and less than 1 day) for each element.

        This attribute returns the seconds component of each timedelta value,
        which is the number of seconds remaining after subtracting whole days.
        Values range from 0 to 86399.

        See Also
        --------
        Series.dt.seconds : Return number of seconds for each element.
        Series.dt.nanoseconds : Return number of nanoseconds for each element.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="s")
        >>> tdelta_idx
        TimedeltaIndex(['0 days 00:00:01', '0 days 00:00:02', '0 days 00:00:03'],
                       dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.seconds
        Index([1, 2, 3], dtype='int32')
        """
        return self._wrap_field_result(self._data.seconds)

    @property
    def microseconds(self) -> Index:
        """
        Number of microseconds (>= 0 and less than 1 second) for each element.

        This attribute returns the microseconds component of each timedelta value,
        which is the number of microseconds remaining after subtracting whole
        seconds. Values range from 0 to 999999.

        See Also
        --------
        Timedelta.microseconds : Number of microseconds (>= 0 and less than
            1 second).

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="us")
        >>> tdelta_idx
        TimedeltaIndex(['0 days 00:00:00.000001', '0 days 00:00:00.000002',
                        '0 days 00:00:00.000003'],
                       dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.microseconds
        Index([1, 2, 3], dtype='int32')
        """
        return self._wrap_field_result(self._data.microseconds)

    @property
    def nanoseconds(self) -> Index:
        """
        Number of nanoseconds (>= 0 and less than 1 microsecond) for each element.

        This attribute returns the nanoseconds component of each timedelta value,
        which is the number of nanoseconds remaining after subtracting whole
        microseconds. Values range from 0 to 999.

        See Also
        --------
        Series.dt.seconds : Return number of seconds for each element.
        Series.dt.microseconds : Return number of microseconds for each element.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="ns")
        >>> tdelta_idx
        TimedeltaIndex(['0 days 00:00:00.000000001', '0 days 00:00:00.000000002',
                        '0 days 00:00:00.000000003'],
                       dtype='timedelta64[ns]', freq=None)
        >>> tdelta_idx.nanoseconds
        Index([1, 2, 3], dtype='int32')
        """
        return self._wrap_field_result(self._data.nanoseconds)

    @property
    def components(self) -> DataFrame:
        """
        Return a DataFrame of the individual resolution components of the Timedeltas.

        The components (days, hours, minutes seconds, milliseconds, microseconds,
        nanoseconds) are returned as columns in a DataFrame.

        Returns
        -------
        DataFrame

        See Also
        --------
        TimedeltaIndex.total_seconds : Return total duration expressed in seconds.
        Timedelta.components : Return a components namedtuple-like of a single
            timedelta.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta(["1 day 3 min 2 us 42 ns"])
        >>> tdelta_idx
        TimedeltaIndex(['1 days 00:03:00.000002042'],
                       dtype='timedelta64[ns]', freq=None)
        >>> tdelta_idx.components
           days  hours  minutes  seconds  milliseconds  microseconds  nanoseconds
        0     1      0        3        0             0             2           42
        """
        return self._data.components

    def to_pytimedelta(self) -> npt.NDArray[np.object_]:
        """
        Return an ndarray of datetime.timedelta objects.

        Each element of the :class:`TimedeltaIndex` is converted to the
        corresponding native Python :class:`datetime.timedelta` object.

        Returns
        -------
        numpy.ndarray
            Object-dtype array of :class:`datetime.timedelta`.

        See Also
        --------
        to_timedelta : Convert argument to timedelta format.
        Timedelta : Represents a duration between two dates or times.
        DatetimeIndex: Index of datetime64 data.
        Timedelta.components : Return a components namedtuple-like
                               of a single timedelta.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="D")
        >>> tdelta_idx
        TimedeltaIndex(['1 days', '2 days', '3 days'],
                        dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.to_pytimedelta()
        array([datetime.timedelta(days=1), datetime.timedelta(days=2),
               datetime.timedelta(days=3)], dtype=object)
        """
        return self._data.to_pytimedelta()

    def sum(
        self,
        *,
        axis: AxisInt | None = None,
        dtype: NpDtype | None = None,
        out=None,
        keepdims: bool = False,
        initial=None,
        skipna: bool = True,
        min_count: int = 0,
    ) -> Timedelta | NaTType:
        """
        Return the sum of the values over the requested axis.

        Parameters
        ----------
        axis : int, optional
            Axis for the function to be applied on.
        dtype, out, keepdims, initial
            Not implemented; kept for compatibility with :func:`numpy.sum`,
            which calls this method when ``numpy.sum(tdi)`` is used. Must be
            left at their default values.
        skipna : bool, default True
            Whether to ignore any NaT elements.
        min_count : int, default 0
            The required number of valid values to perform the operation. If fewer
            than ``min_count`` non-NaT values are present the result is NaT.

        Returns
        -------
        Timedelta
            The sum, or NaT if it cannot be computed.

        See Also
        --------
        numpy.ndarray.sum : Returns the sum of array elements along a given axis.
        Series.sum : Return the sum of the values in a Series.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="D")
        >>> tdelta_idx
        TimedeltaIndex(['1 days', '2 days', '3 days'],
                        dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.sum()
        Timedelta('6 days 00:00:00')
        """
        return self._data.sum(
            axis=axis,
            dtype=dtype,
            out=out,
            keepdims=keepdims,
            initial=initial,
            skipna=skipna,
            min_count=min_count,
        )

    def std(
        self,
        *,
        axis: AxisInt | None = None,
        dtype: NpDtype | None = None,
        out=None,
        ddof: int = 1,
        keepdims: bool = False,
        skipna: bool = True,
    ) -> Timedelta | NaTType:
        """
        Return the standard deviation of the values over the requested axis.

        Parameters
        ----------
        axis : int, optional
            Axis for the function to be applied on.
        dtype, out, keepdims
            Not implemented; kept for compatibility with :func:`numpy.std`,
            which calls this method when ``numpy.std(tdi)`` is used. Must be
            left at their default values.
        ddof : int, default 1
            Delta degrees of freedom. The divisor used in calculations is
            ``N - ddof``, where ``N`` represents the number of elements.
        skipna : bool, default True
            Whether to ignore any NaT elements.

        Returns
        -------
        Timedelta
            The standard deviation, or NaT if it cannot be computed.

        See Also
        --------
        numpy.ndarray.std : Returns the standard deviation of array elements
            along a given axis.
        Series.std : Return sample standard deviation over requested axis.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="D")
        >>> tdelta_idx
        TimedeltaIndex(['1 days', '2 days', '3 days'],
                        dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.std()
        Timedelta('1 days 00:00:00')
        """
        return self._data.std(
            axis=axis, dtype=dtype, out=out, ddof=ddof, keepdims=keepdims, skipna=skipna
        )

    def median(
        self,
        *,
        axis: AxisInt | None = None,
        skipna: bool = True,
        **kwargs,
    ) -> Timedelta | NaTType:
        """
        Return the median of the values over the requested axis.

        Parameters
        ----------
        axis : int, optional
            Axis for the function to be applied on.
        skipna : bool, default True
            Whether to ignore any NaT elements.
        **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Timedelta
            The median, or NaT if it cannot be computed.

        See Also
        --------
        numpy.median : Compute the median along the specified axis.
        Series.median : Return the median of the values in a Series.

        Examples
        --------
        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="D")
        >>> tdelta_idx
        TimedeltaIndex(['1 days', '2 days', '3 days'],
                        dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.median()
        Timedelta('2 days 00:00:00')
        """
        return self._data.median(axis=axis, skipna=skipna, **kwargs)

    # -------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data=None,
        freq=lib.no_default,
        dtype=None,
        copy: bool | None = None,
        name=None,
    ):
        name = maybe_extract_name(name, data, cls)

        # GH#63388
        data, copy = cls._maybe_copy_array_input(data, copy, dtype)

        if is_scalar(data):
            cls._raise_scalar_data_error(data)

        if dtype is not None:
            dtype = pandas_dtype(dtype)

        if (
            isinstance(data, TimedeltaArray)
            and freq is lib.no_default
            and (dtype is None or dtype == data.dtype)
        ):
            if copy:
                data = data.copy()
            return cls._simple_new(data, name=name)

        if (
            isinstance(data, TimedeltaIndex)
            and freq is lib.no_default
            and name is None
            and (dtype is None or dtype == data.dtype)
        ):
            if copy:
                return data.copy()
            else:
                return data._view()

        # - Cases checked above all return/raise before reaching here - #

        inferred_freq = data.freq if isinstance(data, TimedeltaIndex) else None

        tdarr = TimedeltaArray._from_sequence(data, dtype=dtype, copy=copy)
        refs = None
        if not copy and isinstance(data, (ABCSeries, Index)):
            refs = data._references

        result = cls._simple_new(tdarr, name=name, refs=refs)
        result._pin_freq(freq, inferred_freq)
        return result

    def __reduce__(self):
        d = {"data": self._data, "name": self.name, "freq": self._freq}
        return _new_TimedeltaIndex, (type(self), d), None

    # -------------------------------------------------------------------

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        """
        Can we compare values of the given dtype to our own?
        """
        if isinstance(dtype, ArrowDtype):
            return dtype.kind == "m"
        return lib.is_np_dtype(dtype, "m")  # aka self._data._is_recognized_dtype

    # -------------------------------------------------------------------
    # Arithmetic Methods

    def __neg__(self) -> TimedeltaIndex:
        result = self._data.__neg__()
        idx = type(self)._simple_new(result, name=self.name)
        if self.freq is not None:
            idx._freq = -self.freq
        return idx

    def __pos__(self) -> TimedeltaIndex:
        result = self._data.__pos__()
        idx = type(self)._simple_new(result, name=self.name)
        if self.freq is not None:
            idx._freq = self.freq
        return idx

    def _arith_method(self, other: object, op: Callable[..., Any]) -> Index:
        result = super()._arith_method(other, op)
        if self.freq is None or not is_scalar(other):
            return result

        new_freq = None
        if isinstance(result, type(self)):
            new_freq = self._get_arith_result_freq(other, op)
        elif (
            getattr(getattr(result, "dtype", None), "kind", None) == "M" and op is rsub
        ):
            # Timestamp/datetime - TDI produces a DatetimeIndex.
            # The array __rsub__ does (-self) + other, losing freq in
            # negation. Compute the freq at the Index level.
            new_freq = self._get_rsub_datetime_result_freq(other)

        if new_freq is not None:
            result._freq = new_freq
        return result

    def _get_arith_result_freq(
        self, other: object, op: Callable[..., Any]
    ) -> Day | Tick | None:
        """
        Compute the result freq for arithmetic operations whose result
        is also a TimedeltaIndex.

        Caller is responsible for checking self.freq is not None.
        """
        freq = self.freq
        assert freq is not None  # caller ensures this
        if op in (operator.mul, rmul):
            if bool(isna(other)):
                return None
            # error: No overload variant of "__mul__" of "BaseOffset"
            # matches argument type "object"
            new_freq = freq * other  # type: ignore[operator]
            if new_freq.n == 0:
                # GH#51575 Better to have no freq than an incorrect one
                return None
            return new_freq

        if op in (operator.truediv, operator.floordiv):
            # Note: freq gets division, not floor-division, even if op
            #  is floordiv.
            if isinstance(freq, Day):
                if freq.n % other == 0:  # type: ignore[operator]
                    new_freq = Day(freq.n // other)  # type: ignore[operator]
                else:
                    new_freq = to_offset(Timedelta(days=freq.n)) / other  # type: ignore[operator]
            else:
                new_freq = freq / other  # type: ignore[operator]
            if new_freq.nanos == 0 and freq.nanos != 0:
                # e.g. if self.freq is Nano(1) then dividing by 2
                #  rounds down to zero
                return None
            return new_freq

        if op is rsub:
            # scalar_timedelta - TDI: the array uses (-self) + other,
            # losing freq in negation. Result freq is -self.freq.
            return -freq  # type: ignore[return-value]

        return None

    def _get_rsub_datetime_result_freq(self, other: object) -> Day | Tick | None:
        """
        Compute the result freq for Timestamp/datetime - TimedeltaIndex.

        Mirrors the logic of _get_arithmetic_result_freq for the negated
        array case.
        """
        freq = self.freq
        assert freq is not None  # caller ensures this
        if isinstance(freq, Tick):
            return -freq

        # freq is a Day; only preserve with tz-naive or UTC
        if isinstance(other, Timestamp):
            tz = other.tz
        else:
            tz = Timestamp(other).tz  # type: ignore[arg-type]
        if tz is None or timezones.is_utc(tz):
            return -freq  # type: ignore[return-value]
        return None

    # -------------------------------------------------------------------
    # Indexing Methods

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int, slice, or ndarray[int]
        """
        self._check_indexing_error(key)

        try:
            key = self._data._validate_scalar(key, unbox=False)
        except TypeError as err:
            raise KeyError(key) from err

        return Index.get_loc(self, key)

    # error: Return type "tuple[Timedelta | NaTType, Resolution]" of
    # "_parse_with_reso" incompatible with return type
    # "tuple[datetime, Resolution]" in supertype
    # "pandas.core.indexes.datetimelike.DatetimeIndexOpsMixin"
    def _parse_with_reso(self, label: str) -> tuple[Timedelta | NaTType, Resolution]:  # type: ignore[override]
        # Resolution comes from the string text (GH#33603), not the value's
        # components: "720s" resolves to "s", not "min" (even though
        # 720s == 12 minutes). The parser reports both in a single pass.
        parsed, string_reso_code = parse_timedelta_string_reso(label)
        if isinstance(parsed, Timedelta):
            string_reso = Resolution(string_reso_code)
            # Fold in sub-unit precision the written unit doesn't capture,
            # e.g. "1.5min" is minute-written but second-resolution.
            value_reso = Resolution.get_reso_from_freqstr(parsed.resolution_string)
            reso = min(string_reso, value_reso)
        else:
            # i.e. pd.NaT
            reso = Resolution.RESO_SEC
        return parsed, reso

    def _parsed_string_to_bounds(self, reso: Resolution, parsed: Timedelta):
        reso_str = reso.attr_abbrev
        lbound = parsed.floor(reso_str)
        rbound = (
            lbound
            + to_offset(reso_str)
            - Timedelta(1, unit=self.unit).as_unit(self.unit)
        )
        # If reso is finer than the index unit, the window [lbound, rbound]
        # collapses to lbound alone; without the clamp rbound < lbound.
        return lbound, max(lbound, rbound)

    # -------------------------------------------------------------------

    @property
    def inferred_type(self) -> str:
        return "timedelta64"


@set_module("pandas")
def timedelta_range(
    start=None,
    end=None,
    periods: int | None = None,
    freq=None,
    name=None,
    closed=None,
    *,
    unit: TimeUnit | None = None,
) -> TimedeltaIndex:
    """
    Return a fixed frequency TimedeltaIndex with day as the default.

    This function generates a sequence of evenly spaced timedelta values
    between the specified bounds, using day as the default frequency.

    Parameters
    ----------
    start : str or timedelta-like, default None
        Left bound for generating timedeltas.
    end : str or timedelta-like, default None
        Right bound for generating timedeltas.
    periods : int, default None
        Number of periods to generate.
    freq : str, Timedelta, datetime.timedelta, or DateOffset, default 'D'
        Frequency strings can have multiples, e.g. '5h'.
    name : Hashable, default None
        Name of the resulting TimedeltaIndex.
    closed : str, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None).
    unit : {'s', 'ms', 'us', 'ns', None}, default None
        Specify the desired resolution of the result.
        If not specified, this is inferred from the 'start', 'end', and 'freq'
        using the same inference as :class:`Timedelta` taking the highest
        resolution of the three that are provided.

        .. versionadded:: 2.0.0

    Returns
    -------
    TimedeltaIndex
        Fixed frequency, with day as the default.

    See Also
    --------
    date_range : Return a fixed frequency DatetimeIndex.
    period_range : Return a fixed frequency PeriodIndex.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    a maximum of three can be specified at once. Of the three parameters
    ``start``, ``end``, and ``periods``, at least two must be specified.
    If ``freq`` is omitted, the resulting ``DatetimeIndex`` will have
    ``periods`` linearly spaced elements between ``start`` and ``end``
    (closed on both sides).

    To learn more about the frequency strings, please see
    :ref:`this link<timeseries.offset_aliases>`.

    Examples
    --------
    >>> pd.timedelta_range(start="1 day", periods=4)
    TimedeltaIndex(['1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq='D')

    The ``closed`` parameter specifies which endpoint is included.  The default
    behavior is to include both endpoints.

    >>> pd.timedelta_range(start="1 day", periods=4, closed="right")
    TimedeltaIndex(['2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq='D')

    The ``freq`` parameter specifies the frequency of the TimedeltaIndex.
    Only fixed frequencies can be passed, non-fixed frequencies such as
    'M' (month end) will raise.

    >>> pd.timedelta_range(start="1 day", end="2 days", freq="6h")
    TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 12:00:00',
                    '1 days 18:00:00', '2 days 00:00:00'],
                   dtype='timedelta64[us]', freq='6h')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.timedelta_range(start="1 day", end="5 days", periods=4)
    TimedeltaIndex(['1 days 00:00:00', '2 days 08:00:00', '3 days 16:00:00',
                    '5 days 00:00:00'],
                   dtype='timedelta64[us]', freq=None)

    **Specify a unit**

    >>> pd.timedelta_range("1 Day", periods=3, freq="100000D", unit="s")
    TimedeltaIndex(['1 days', '100001 days', '200001 days'],
                   dtype='timedelta64[s]', freq='100000D')
    """
    if freq is None and com.any_none(periods, start, end):
        freq = "D"
    freq = to_offset(freq)

    if com.count_not_none(start, end, periods, freq) != 3:
        # This check needs to come before the `unit = start.unit` line below
        raise ValueError(
            "Of the four parameters: start, end, periods, "
            "and freq, exactly three must be specified"
        )

    if unit is None:
        # Infer the unit based on the inputs

        if start is not None and end is not None:
            start = Timedelta(start)
            end = Timedelta(end)
            start = cast("Timedelta", start)
            end = cast("Timedelta", end)
            if abbrev_to_npy_unit(start.unit) > abbrev_to_npy_unit(end.unit):
                unit = cast("TimeUnit", start.unit)
            else:
                unit = cast("TimeUnit", end.unit)
        elif start is not None:
            start = Timedelta(start)
            start = cast("Timedelta", start)
            unit = cast("TimeUnit", start.unit)
        else:
            end = Timedelta(end)
            end = cast("Timedelta", end)
            unit = cast("TimeUnit", end.unit)

        # Last we need to watch out for cases where the 'freq' implies a higher
        #  unit than either start or end
        if freq is not None:
            freq = cast("Tick | Day", freq)
            creso = abbrev_to_npy_unit(unit)
            if freq._creso > creso:  # pyright: ignore[reportAttributeAccessIssue]
                unit = cast("TimeUnit", freq.base.freqstr)

    tdarr = TimedeltaArray._generate_range(
        start, end, periods, freq, closed=closed, unit=unit
    )
    result = TimedeltaIndex._simple_new(tdarr, name=name)
    result._freq = freq
    return result
