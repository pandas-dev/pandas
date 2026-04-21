"""
Base and utility classes for tseries type pandas objects.
"""

from __future__ import annotations

from abc import (
    ABC,
)
from itertools import pairwise
import operator
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
    cast,
    final,
)
import warnings

import numpy as np

from pandas._libs import (
    NaT,
    lib,
)
from pandas._libs.tslibs import (
    BaseOffset,
    Day,
    Resolution,
    Tick,
    Timedelta,
    Timestamp,
    get_resolution,
    parsing,
    timezones,
    to_offset,
)
from pandas._libs.tslibs.dtypes import abbrev_to_npy_unit
from pandas._libs.tslibs.offsets import FY5253Mixin
from pandas.compat.numpy import function as nv
from pandas.errors import (
    InvalidIndexError,
    NullFrequencyError,
    OutOfBoundsDatetime,
    OutOfBoundsTimedelta,
    Pandas4Warning,
)
from pandas.util._decorators import (
    cache_readonly,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_integer,
    is_list_like,
)
from pandas.core.dtypes.concat import concat_compat
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    PeriodDtype,
)

from pandas.core import roperator
from pandas.core.arrays import (
    DatetimeArray,
    ExtensionArray,
    PeriodArray,
    TimedeltaArray,
)
import pandas.core.common as com
from pandas.core.construction import ensure_wrapped_if_datetimelike
from pandas.core.indexers import check_array_indexer
from pandas.core.indexes.base import (
    Index,
)
from pandas.core.indexes.extension import NDArrayBackedExtensionIndex
from pandas.core.indexes.range import RangeIndex
from pandas.core.tools.timedeltas import to_timedelta

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Sequence,
    )
    from datetime import datetime

    from pandas._typing import (
        Axis,
        JoinHow,
        TimeUnit,
        npt,
    )

    from pandas import CategoricalIndex


class DatetimeIndexOpsMixin(NDArrayBackedExtensionIndex, ABC):
    """
    Common ops mixin to support a unified interface datetimelike Index.
    """

    _can_hold_strings = False
    _data: DatetimeArray | TimedeltaArray | PeriodArray

    def mean(self, *, skipna: bool = True, axis: int | None = 0):
        """
        Return the mean value of the Array.

        This method computes the arithmetic mean of the datetime or timedelta
        values in the index, optionally skipping NaT values.

        Parameters
        ----------
        skipna : bool, default True
            Whether to ignore any NaT elements.
        axis : int, optional, default 0
            Axis for the function to be applied on.

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

        Examples
        --------
        For :class:`pandas.DatetimeIndex`:

        >>> idx = pd.date_range("2001-01-01 00:00", periods=3)
        >>> idx
        DatetimeIndex(['2001-01-01', '2001-01-02', '2001-01-03'],
                      dtype='datetime64[us]', freq='D')
        >>> idx.mean()
        Timestamp('2001-01-02 00:00:00')

        For :class:`pandas.TimedeltaIndex`:

        >>> tdelta_idx = pd.to_timedelta([1, 2, 3], unit="D")
        >>> tdelta_idx
        TimedeltaIndex(['1 days', '2 days', '3 days'],
                        dtype='timedelta64[s]', freq=None)
        >>> tdelta_idx.mean()
        Timedelta('2 days 00:00:00')
        """
        return self._data.mean(skipna=skipna, axis=axis)

    @property
    def freq(self) -> BaseOffset | None:
        """
        Return the frequency object if it is set, otherwise None.

        To learn more about the frequency strings, please see
        :ref:`this link<timeseries.offset_aliases>`.

        See Also
        --------
        DatetimeIndex.freq : Return the frequency object if it is set, otherwise None.
        PeriodIndex.freq : Return the frequency object if it is set, otherwise None.

        Examples
        --------
        >>> datetimeindex = pd.date_range(
        ...     "2022-02-22 02:22:22", periods=10, tz="America/Chicago", freq="h"
        ... )
        >>> datetimeindex
        DatetimeIndex(['2022-02-22 02:22:22-06:00', '2022-02-22 03:22:22-06:00',
                       '2022-02-22 04:22:22-06:00', '2022-02-22 05:22:22-06:00',
                       '2022-02-22 06:22:22-06:00', '2022-02-22 07:22:22-06:00',
                       '2022-02-22 08:22:22-06:00', '2022-02-22 09:22:22-06:00',
                       '2022-02-22 10:22:22-06:00', '2022-02-22 11:22:22-06:00'],
                      dtype='datetime64[us, America/Chicago]', freq='h')
        >>> datetimeindex.freq
        <Hour>
        """
        return self._data.freq

    @freq.setter
    def freq(self, value) -> None:
        arr = self._data
        if not isinstance(arr, (DatetimeArray, TimedeltaArray)):
            # e.g. PeriodArray freq is derived from dtype and is read-only
            raise AttributeError(
                f"property 'freq' of {type(arr).__name__!r} object has no setter"
            )
        if value is not None:
            value = to_offset(value)
            arr._validate_frequency(arr, value)
            if arr.dtype.kind == "m" and not isinstance(value, (Tick, Day)):
                raise TypeError("TimedeltaArray/Index freq must be a Tick")

            if arr.ndim > 1:
                raise ValueError("Cannot set freq with ndim > 1")

        arr._freq = value

    @property
    def asi8(self) -> npt.NDArray[np.int64]:
        """
        Return Integer representation of the values.

        For :class:`DatetimeIndex` and :class:`TimedeltaIndex`, the
        values are the number of time units (determined by the index
        resolution) since the epoch. For :class:`PeriodIndex`, the
        values are ordinals.

        Returns
        -------
        numpy.ndarray
            An ndarray with int64 dtype.

        See Also
        --------
        Index.values : Return an array representing the data in the Index,
            using native types (datetime64, timedelta64) rather than int64.
        Index.to_numpy : Return a NumPy ndarray of the index values.

        Examples
        --------
        For :class:`DatetimeIndex` with default microsecond resolution:

        >>> idx = pd.DatetimeIndex(["2023-01-01", "2023-01-02"], dtype="datetime64[us]")
        >>> idx.asi8
        array([1672531200000000, 1672617600000000])

        For :class:`TimedeltaIndex` with millisecond resolution:

        >>> idx = pd.TimedeltaIndex(["1 day", "2 days"], dtype="timedelta64[ms]")
        >>> idx.asi8
        array([ 86400000, 172800000])

        For :class:`PeriodIndex`:

        >>> idx = pd.PeriodIndex(["2023-01", "2023-02", "2023-03"], freq="M")
        >>> idx.asi8
        array([636, 637, 638])
        """
        return self._data.asi8

    @property
    def freqstr(self) -> str:
        """
        Return the frequency object as a string if it's set, otherwise None.

        This property returns a string representation of the frequency
        (e.g., ``'D'`` for daily, ``'h'`` for hourly) when one has been set
        on the index, either explicitly or via inference.

        See Also
        --------
        DatetimeIndex.inferred_freq : Returns a string representing a frequency
            generated by infer_freq.

        Examples
        --------
        For DatetimeIndex:

        >>> idx = pd.DatetimeIndex(["1/1/2020 10:00:00+00:00"], freq="D")
        >>> idx.freqstr
        'D'

        The frequency can be inferred if there are more than 2 points:

        >>> idx = pd.DatetimeIndex(
        ...     ["2018-01-01", "2018-01-03", "2018-01-05"], freq="infer"
        ... )
        >>> idx.freqstr
        '2D'

        For PeriodIndex:

        >>> idx = pd.PeriodIndex(["2023-1", "2023-2", "2023-3"], freq="M")
        >>> idx.freqstr
        'M'
        """
        from pandas import PeriodIndex

        if self._data.freqstr is not None and isinstance(
            self._data, (PeriodArray, PeriodIndex)
        ):
            freq = PeriodDtype(self._data.freq)._freqstr
            return freq
        else:
            return self._data.freqstr  # type: ignore[return-value]

    @cache_readonly
    def _resolution_obj(self) -> Resolution:
        if isinstance(self.dtype, PeriodDtype):
            return self.dtype._resolution_obj
        elif self.dtype.kind == "M":
            return get_resolution(self.asi8, self.tz, reso=self._data._creso)  # type: ignore[attr-defined,union-attr]
        else:
            return get_resolution(self.asi8, tz=None, reso=self._data._creso)  # type: ignore[union-attr]

    @cache_readonly
    def resolution(self) -> str:
        """
        Returns day, hour, minute, second, millisecond or microsecond
        """
        return self._resolution_obj.attrname

    # ------------------------------------------------------------------------

    @cache_readonly
    def hasnans(self) -> bool:
        return self._data._hasna

    def equals(self, other: Any) -> bool:
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, Index):
            return False
        elif other.dtype.kind in "iufc":
            return False
        elif not isinstance(other, type(self)):
            if other.dtype == object:
                converted = lib.maybe_convert_objects(
                    np.asarray(other), convert_non_numeric=True
                )
                converted = ensure_wrapped_if_datetimelike(converted)
                if isinstance(converted, type(self._data)):
                    other = type(self)._simple_new(converted)
                elif self.dtype.kind == "M" and lib.infer_dtype(other) == "date":
                    # GH#65056
                    warnings.warn(
                        "Inferring datetime64 from data containing "
                        "datetime.date objects is deprecated. In a future "
                        "version, these will be left as object dtype. Use "
                        "pd.to_datetime to explicitly convert to datetime64.",
                        Pandas4Warning,
                        stacklevel=find_stack_level(),
                    )
                    try:
                        other = type(self)(other)
                    except (ValueError, TypeError, OverflowError):
                        return False
            elif isinstance(other.dtype, CategoricalDtype):
                other = cast("CategoricalIndex", other)
                cat_vals = np.asarray(other.categories)
                if cat_vals.dtype == object:
                    converted = lib.maybe_convert_objects(
                        cat_vals, convert_non_numeric=True
                    )
                    converted = ensure_wrapped_if_datetimelike(converted)
                else:
                    converted = ensure_wrapped_if_datetimelike(cat_vals)
                if isinstance(converted, type(self._data)):
                    try:
                        other = type(self)(other)
                    except (ValueError, TypeError, OverflowError):
                        return False
                elif (
                    self.dtype.kind == "M"
                    and cat_vals.dtype == object
                    and lib.infer_dtype(other.categories) == "date"
                ):
                    # GH#65056
                    warnings.warn(
                        "Inferring datetime64 from data containing "
                        "datetime.date objects is deprecated. In a "
                        "future version, these will be left as object "
                        "dtype. Use pd.to_datetime to explicitly "
                        "convert to datetime64.",
                        Pandas4Warning,
                        stacklevel=find_stack_level(),
                    )
                    try:
                        other = type(self)(other)
                    except (ValueError, TypeError, OverflowError):
                        return False

        if type(self) != type(other):
            return False
        elif self.dtype == other.dtype:
            return np.array_equal(self.asi8, other.asi8)
        elif (self.dtype.kind == "M" and self.tz == other.tz) or self.dtype.kind == "m":  # type: ignore[attr-defined]
            # different units, otherwise matching
            try:
                # TODO: do this at the EA level?
                left, right = self._data._ensure_matching_resos(other._data)  # type: ignore[union-attr]
            except (OutOfBoundsDatetime, OutOfBoundsTimedelta):
                return False
            else:
                return np.array_equal(left.view("i8"), right.view("i8"))
        return False

    def __contains__(self, key: Any) -> bool:
        """
        Return a boolean indicating whether the provided key is in the index.

        Parameters
        ----------
        key : label
            The key to check if it is present in the index.

        Returns
        -------
        bool
            Whether the key search is in the index.

        Raises
        ------
        TypeError
            If the key is not hashable.

        See Also
        --------
        Index.isin : Returns an ndarray of boolean dtype indicating whether the
            list-like key is in the index.

        Examples
        --------
        >>> idx = pd.Index([1, 2, 3, 4])
        >>> idx
        Index([1, 2, 3, 4], dtype='int64')
        >>> 2 in idx
        True
        >>> 6 in idx
        False
        """
        hash(key)
        try:
            self.get_loc(key)
        except (KeyError, TypeError, ValueError, InvalidIndexError):
            return False
        return True

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(to_timedelta(tolerance).to_numpy())
        return super()._convert_tolerance(tolerance, target)

    # --------------------------------------------------------------------
    # Rendering Methods
    _default_na_rep = "NaT"

    def _format_with_header(
        self, *, header: list[str], na_rep: str, date_format: str | None = None
    ) -> list[str]:
        # TODO: not reached in tests 2023-10-11
        # matches base class except for whitespace padding and date_format
        return header + list(
            self._get_values_for_csv(na_rep=na_rep, date_format=date_format)
        )

    def _formatter_func(self, val) -> str:
        return self._data._formatter()(val)

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value).
        """
        attrs = super()._format_attrs()
        for attrib in self._attributes:
            # iterating over _attributes prevents us from doing this for PeriodIndex
            if attrib == "freq":
                freq = self.freqstr
                if freq is not None:
                    freq = repr(freq)  # e.g. D -> 'D'
                attrs.append(("freq", freq))
        return attrs

    def _summary(self, name=None) -> str:
        """
        Return a summarized representation.

        Parameters
        ----------
        name : str
            name to use in the summary representation

        Returns
        -------
        String with a summarized representation of the index
        """
        result = super()._summary(name=name)
        if self.freq:
            result += f"\nFreq: {self.freqstr}"

        return result

    # --------------------------------------------------------------------
    # Indexing Methods

    @final
    def _can_partial_date_slice(self, reso: Resolution) -> bool:
        # e.g. test_getitem_setitem_periodindex
        # History of conversation GH#3452, GH#3931, GH#2369, GH#14826
        return reso > self._resolution_obj
        # NB: for DTI/PI, not TDI

    def _parsed_string_to_bounds(self, reso: Resolution, parsed):
        raise NotImplementedError

    def _parse_with_reso(self, label: str) -> tuple[datetime, Resolution]:
        # overridden by TimedeltaIndex
        try:
            if self.freq is None or hasattr(self.freq, "rule_code"):
                freq = self.freq
        except NotImplementedError:
            freq = getattr(self, "freqstr", getattr(self, "inferred_freq", None))

        freqstr: str | None
        if freq is not None and not isinstance(freq, str):
            freqstr = freq.rule_code
        else:
            freqstr = freq

        if isinstance(label, np.str_):
            # GH#45580
            label = str(label)

        parsed, reso_str = parsing.parse_datetime_string_with_reso(label, freqstr)
        reso = Resolution.from_attrname(reso_str)
        return parsed, reso

    def _get_string_slice(self, key: str) -> slice | npt.NDArray[np.intp]:  # type: ignore[override]
        # overridden by TimedeltaIndex
        parsed, reso = self._parse_with_reso(key)
        try:
            return self._partial_date_slice(reso, parsed)
        except KeyError as err:
            raise KeyError(key) from err

    @final
    def _partial_date_slice(
        self,
        reso: Resolution,
        parsed: datetime,
    ) -> slice | npt.NDArray[np.intp]:
        """
        Parameters
        ----------
        reso : Resolution
        parsed : datetime

        Returns
        -------
        slice or ndarray[intp]
        """
        if not self._can_partial_date_slice(reso):
            raise ValueError

        t1, t2 = self._parsed_string_to_bounds(reso, parsed)
        vals = self._data._ndarray
        unbox = self._data._unbox

        if self.is_monotonic_increasing:
            if len(self) and (
                (t1 < self[0] and t2 < self[0]) or (t1 > self[-1] and t2 > self[-1])
            ):
                # we are out of range
                raise KeyError

            # a monotonic increasing series can be sliced
            #  (searchsorted requires ascending order)
            left = vals.searchsorted(unbox(t1), side="left")
            right = vals.searchsorted(unbox(t2), side="right")
            return slice(left, right)

        elif self.is_monotonic_decreasing:
            if len(self) and (
                (t1 > self[0] and t2 > self[0]) or (t1 < self[-1] and t2 < self[-1])
            ):
                # we are out of range
                raise KeyError

            # searchsorted requires ascending order, so search the reversed
            #  array and convert the indices back
            reversed_vals = vals[::-1]
            nvals = len(vals)
            rev_left = reversed_vals.searchsorted(unbox(t1), side="left")
            rev_right = reversed_vals.searchsorted(unbox(t2), side="right")
            return slice(nvals - rev_right, nvals - rev_left)

        else:
            lhs_mask = vals >= unbox(t1)
            rhs_mask = vals <= unbox(t2)

            # try to find the dates
            return (lhs_mask & rhs_mask).nonzero()[0]

    def _maybe_cast_slice_bound(self, label, side: str):
        """
        If label is a string, cast it to scalar type according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}

        Returns
        -------
        label : object

        Notes
        -----
        Value of `side` parameter should be validated in caller.
        """
        if isinstance(label, str):
            try:
                parsed, reso = self._parse_with_reso(label)
            except ValueError as err:
                # DTI -> parsing.DateParseError
                # TDI -> 'unit abbreviation w/o a number'
                # PI -> string cannot be parsed as datetime-like
                self._raise_invalid_indexer("slice", label, err)

            lower, upper = self._parsed_string_to_bounds(reso, parsed)
            return lower if side == "left" else upper
        elif not isinstance(label, self._data._recognized_scalars):
            self._raise_invalid_indexer("slice", label)

        return label

    # --------------------------------------------------------------------
    # Arithmetic Methods

    def shift(self, periods: int = 1, freq=None) -> Self:
        """
        Shift index by desired number of time frequency increments.

        This method is for shifting the values of datetime-like indexes
        by a specified time increment a given number of times.

        Parameters
        ----------
        periods : int, default 1
            Number of periods (or increments) to shift by,
            can be positive or negative.
        freq : pandas.DateOffset, pandas.Timedelta or string, optional
            Frequency increment to shift by.
            If None, the index is shifted by its own `freq` attribute.
            Offset aliases are valid strings, e.g., 'D', 'W', 'M' etc.

        Returns
        -------
        pandas.DatetimeIndex
            Shifted index.

        See Also
        --------
        Index.shift : Shift values of Index.
        PeriodIndex.shift : Shift values of PeriodIndex.
        """
        raise NotImplementedError

    # --------------------------------------------------------------------

    def _maybe_cast_listlike_indexer(self, keyarr):
        """
        Analogue to maybe_cast_indexer for get_indexer instead of get_loc.
        """
        try:
            res = self._data._validate_listlike(keyarr, allow_object=True)
        except (ValueError, TypeError):
            if not isinstance(keyarr, ExtensionArray):
                # e.g. we don't want to cast DTA to ndarray[object]
                res = com.asarray_tuplesafe(keyarr)
                # TODO: com.asarray_tuplesafe shouldn't cast e.g. DatetimeArray
            else:
                res = keyarr
        return Index(res, dtype=res.dtype)


class DatetimeTimedeltaMixin(DatetimeIndexOpsMixin, ABC):
    """
    Mixin class for methods shared by DatetimeIndex and TimedeltaIndex,
    but not PeriodIndex
    """

    _data: DatetimeArray | TimedeltaArray
    _comparables = ["name", "freq"]
    _attributes = ["name", "freq"]

    # Compat for frequency inference, see GH#23789
    _is_monotonic_increasing = Index.is_monotonic_increasing
    _is_monotonic_decreasing = Index.is_monotonic_decreasing
    _is_unique = Index.is_unique

    def astype(self, dtype, copy: bool = True):
        result = super().astype(dtype, copy=copy)
        if isinstance(result, type(self)):
            # Preserve freq for unit conversions (e.g. datetime64[ns] -> [us])
            result._data._freq = self.freq
        return result

    def putmask(self, mask, value) -> Index:
        # GH#24555 putmask may modify values out-of-sequence; drop freq
        result = super().putmask(mask, value)
        if isinstance(result, type(self)):
            result._data._freq = None
        return result

    def _pin_freq(self, freq, validate_kwds: dict) -> None:
        """
        Constructor helper to pin the appropriate ``freq`` attribute on
        ``self._data``.  Assumes ``self._data._freq`` is currently set to any
        freq inferred from input data.
        """
        arr = self._data
        if freq is None:
            # user explicitly passed None -> override any inferred_freq
            arr._freq = None
        elif freq == "infer":
            # if arr._freq is *not* None then we already inferred a freq
            #  and there is nothing left to do
            if arr._freq is None:
                # Set _freq directly to bypass duplicative _validate_frequency
                # check.
                arr._freq = to_offset(self.inferred_freq)
        elif freq is lib.no_default:
            # user did not specify anything, keep inferred freq if the original
            #  data had one, otherwise do nothing
            pass
        elif arr._freq is None:
            # We cannot inherit a freq from the data, so we need to validate
            #  the user-passed freq
            freq = to_offset(freq)
            type(arr)._validate_frequency(self, freq, **validate_kwds)
            arr._freq = freq
        else:
            # Otherwise we just need to check that the user-passed freq
            #  doesn't conflict with the one we already have.
            freq = to_offset(freq)
            if freq != arr._freq:
                # GH#61086 freq may be equivalent but not equal (e.g.
                # QS-FEB vs QS-MAY), so validate against the actual data.
                if len(self) == 0:
                    pass
                elif len(self) == 1:
                    if not freq.is_on_offset(self[0]):
                        raise ValueError(
                            f"Inferred frequency {arr._freq} from passed "
                            "values does not conform to passed frequency "
                            f"{freq.freqstr}"
                        )
                elif self[0] + freq == self[1]:
                    # For standard offsets, the step is a deterministic
                    # function of the date, so agreement on one step proves
                    # equivalence. For Custom/FY5253 offsets, external
                    # state (holidays, 52/53-week patterns) could cause
                    # later steps to diverge, so we validate fully.
                    if hasattr(freq, "_holidays") or isinstance(freq, FY5253Mixin):
                        type(arr)._validate_frequency(self, freq, **validate_kwds)
                else:
                    raise ValueError(
                        f"Inferred frequency {arr._freq} from passed "
                        "values does not conform to passed frequency "
                        f"{freq.freqstr}"
                    )
            arr._freq = freq

    def _get_arithmetic_result_freq(self, other) -> BaseOffset | None:
        """
        Check if we can preserve self.freq in addition or subtraction.
        """
        if not lib.is_scalar(other):
            return None

        # Normalize scalar types to their pandas equivalents.
        # Array arithmetic methods internally convert Tick/timedelta/
        # np.timedelta64 → Timedelta and datetime/np.datetime64 → Timestamp.
        # Since we receive the original user argument, normalize here.
        if isinstance(other, Tick):
            other = Timedelta(other)
        elif isinstance(other, np.timedelta64):
            other = Timedelta(other)
        elif isinstance(other, np.datetime64):
            other = Timestamp(other)
        elif not isinstance(other, (Timedelta, Timestamp)):
            # Handles datetime.timedelta, datetime.datetime
            try:
                other = Timedelta(other)
            except (ValueError, TypeError, OverflowError):
                try:
                    other = Timestamp(other)
                except (ValueError, TypeError, OverflowError):
                    return None

        if isinstance(self.freq, Tick):
            return self.freq
        elif self.dtype.kind == "m" and isinstance(other, Timedelta):
            return self.freq
        elif (
            self.dtype.kind == "m"
            and isinstance(other, Timestamp)
            and (other.tz is None or timezones.is_utc(other.tz))
        ):
            return self.freq
        elif (
            lib.is_np_dtype(self.dtype, "M")
            and isinstance(self.freq, Day)
            and isinstance(other, (Timedelta, Timestamp))
        ):
            return self.freq

        return None

    def _arith_method(self, other, op):
        result = super()._arith_method(other, op)
        if result is NotImplemented:
            return result
        if op in (operator.add, roperator.radd, operator.sub, roperator.rsub):
            new_freq = self._get_arithmetic_result_freq(other)
            if new_freq is not None and isinstance(result, DatetimeTimedeltaMixin):
                if op is roperator.rsub:
                    new_freq = -new_freq
                result._data._freq = new_freq
        return result

    def factorize(
        self,
        sort: bool = False,
        use_na_sentinel: bool = True,
    ):
        if self.freq is not None:
            # We must be unique, so can short-circuit (and retain freq)
            if sort and self.freq.n < 0:
                codes = np.arange(len(self) - 1, -1, -1, dtype=np.intp)
                uniques = self[::-1]
            else:
                codes = np.arange(len(self), dtype=np.intp)
                uniques = self.copy()
            uniques._name = None
            return codes, uniques
        return super().factorize(sort=sort, use_na_sentinel=use_na_sentinel)

    @property
    def unit(self) -> TimeUnit:
        """
        The precision unit of the datetime data.

        Returns the precision unit for the dtype of the index. This is the
        smallest time frame that can be stored within this dtype.

        Returns
        -------
        str
            Unit string representation (e.g. ``"ns"``).

        See Also
        --------
        DatetimeIndex.as_unit : Convert to the given unit.
        TimedeltaIndex.as_unit : Convert to the given unit.

        Examples
        --------
        For a DatetimeIndex:

        >>> idx = pd.DatetimeIndex(["2020-01-02 01:02:03.004005006"])
        >>> idx.unit
        'ns'

        >>> idx_s = pd.DatetimeIndex(["2020-01-02 01:02:03"], dtype="datetime64[s]")
        >>> idx_s.unit
        's'

        For a TimedeltaIndex:

        >>> tdidx = pd.TimedeltaIndex(["1 day 3 min 2 us 42 ns"])
        >>> tdidx.unit
        'ns'

        >>> tdidx_s = pd.TimedeltaIndex(["1 day 3 min"], dtype="timedelta64[s]")
        >>> tdidx_s.unit
        's'
        """
        return self._data.unit

    def as_unit(self, unit: TimeUnit) -> Self:
        """
        Convert to a dtype with the given unit resolution.

        This method is for converting the dtype of a ``DatetimeIndex`` or
        ``TimedeltaIndex`` to a new dtype with the given unit
        resolution/precision.

        Parameters
        ----------
        unit : {'s', 'ms', 'us', 'ns'}

        Returns
        -------
        same type as self
            Converted to the specified unit.

        See Also
        --------
        Timestamp.as_unit : Convert to the given unit.
        Timedelta.as_unit : Convert to the given unit.
        DatetimeIndex.as_unit : Convert to the given unit.
        TimedeltaIndex.as_unit : Convert to the given unit.

        Examples
        --------
        For :class:`pandas.DatetimeIndex`:

        >>> idx = pd.DatetimeIndex(["2020-01-02 01:02:03.004005006"])
        >>> idx
        DatetimeIndex(['2020-01-02 01:02:03.004005006'],
                      dtype='datetime64[ns]', freq=None)
        >>> idx.as_unit("s")
        DatetimeIndex(['2020-01-02 01:02:03'], dtype='datetime64[s]', freq=None)

        For :class:`pandas.TimedeltaIndex`:

        >>> tdelta_idx = pd.to_timedelta(["1 day 3 min 2 us 42 ns"])
        >>> tdelta_idx
        TimedeltaIndex(['1 days 00:03:00.000002042'],
                        dtype='timedelta64[ns]', freq=None)
        >>> tdelta_idx.as_unit("s")
        TimedeltaIndex(['1 days 00:03:00'], dtype='timedelta64[s]', freq=None)
        """
        arr = self._data.as_unit(unit)
        result = type(self)._simple_new(arr, name=self.name)
        result._data._freq = self.freq
        return result

    def _with_freq(self, freq):
        # GH#29843
        if freq is None:
            pass
        elif isinstance(freq, BaseOffset):
            if self.dtype.kind == "m" and not isinstance(freq, (Tick, Day)):
                raise TypeError("TimedeltaArray/Index freq must be a Tick")
        elif freq == "infer":
            freq = to_offset(self.inferred_freq)
        else:
            raise ValueError(f"Invalid frequency: {freq!r}")

        arr = self._data.view()
        arr._freq = freq
        return type(self)._simple_new(arr, name=self._name)

    @property
    def values(self) -> np.ndarray:
        # NB: For Datetime64TZ this is lossy
        data = self._data._ndarray
        data = data.view()
        data.flags.writeable = False
        return data

    def shift(self, periods: int = 1, freq=None) -> Self:
        """
        Shift index by desired number of time frequency increments.
        This method is for shifting the values of datetime-like indexes
        by a specified time increment a given number of times.

        Parameters
        ----------
        periods : int, default 1
            Number of periods (or increments) to shift by,
            can be positive or negative.
        freq : pandas.DateOffset, pandas.Timedelta or string, optional
            Frequency increment to shift by.
            If None, the index is shifted by its own `freq` attribute.
            Offset aliases are valid strings, e.g., 'D', 'W', 'M' etc.

        Returns
        -------
        pandas.DatetimeIndex
            Shifted index.

        See Also
        --------
        Index.shift : Shift values of Index.
        PeriodIndex.shift : Shift values of PeriodIndex.
        """
        if freq is not None and freq != self.freq:
            if isinstance(freq, str):
                freq = to_offset(freq)
            offset = periods * freq
            return self + offset

        if periods == 0 or len(self) == 0:
            # GH#14811 empty case
            return self.copy()

        if self.freq is None:
            raise NullFrequencyError("Cannot shift with no freq")

        start = self[0] + periods * self.freq
        end = self[-1] + periods * self.freq

        # Note: in the DatetimeTZ case, _generate_range will infer the
        #  appropriate timezone from `start` and `end`, so tz does not need
        #  to be passed explicitly.
        result = self._data._generate_range(
            start=start, end=end, periods=None, freq=self.freq, unit=self.unit
        )
        return type(self)._simple_new(result, name=self.name)

    @cache_readonly
    def inferred_freq(self) -> str | None:
        """
        Return the inferred frequency of the index.

        Attempts to determine the frequency of the index by analyzing the
        differences between consecutive values using ``infer_freq``.

        Returns
        -------
        str or None
            A string representing a frequency generated by ``infer_freq``.
            Returns ``None`` if the frequency cannot be inferred.

        See Also
        --------
        DatetimeIndex.freqstr : Return the frequency object as a string if it's set,
            otherwise ``None``.

        Examples
        --------
        For ``DatetimeIndex``:

        >>> idx = pd.DatetimeIndex(["2018-01-01", "2018-01-03", "2018-01-05"])
        >>> idx.inferred_freq
        '2D'

        For ``TimedeltaIndex``:

        >>> tdelta_idx = pd.to_timedelta(["0 days", "10 days", "20 days"])
        >>> tdelta_idx
        TimedeltaIndex(['0 days', '10 days', '20 days'],
                       dtype='timedelta64[us]', freq=None)
        >>> tdelta_idx.inferred_freq
        '10D'
        """
        return self._data.inferred_freq

    # --------------------------------------------------------------------
    # Set Operation Methods

    @cache_readonly
    def _as_range_index(self) -> RangeIndex:
        # Convert our i8 representations to RangeIndex
        # Caller is responsible for checking isinstance(self.freq, Tick)
        freq = cast("Tick", self.freq)
        tick = Timedelta(freq).as_unit(self.unit)._value
        rng = range(self[0]._value, self[-1]._value + tick, tick)
        return RangeIndex(rng)

    def _can_range_setop(self, other) -> bool:
        return isinstance(self.freq, Tick) and isinstance(other.freq, Tick)

    def _wrap_range_setop(self, other, res_i8) -> Self:
        new_freq = None
        if not len(res_i8):
            # RangeIndex defaults to step=1, which we don't want.
            new_freq = self.freq
        elif isinstance(res_i8, RangeIndex):
            new_freq = to_offset(
                Timedelta(res_i8.step, unit=self.unit).as_unit(self.unit)
            )

        # TODO(GH#41493): we cannot just do
        #  type(self._data)(res_i8.values, dtype=self.dtype, freq=new_freq)
        # because test_setops_preserve_freq fails with _validate_frequency raising.
        # This raising is incorrect, as 'on_freq' is incorrect. This will
        # be fixed by GH#41493
        res_values = res_i8.values.view(self._data._ndarray.dtype)
        result = type(self._data)._simple_new(
            # error: Argument "dtype" to "_simple_new" of "DatetimeArray" has
            # incompatible type "Union[dtype[Any], ExtensionDtype]"; expected
            # "Union[dtype[datetime64], DatetimeTZDtype]"
            res_values,
            dtype=self.dtype,  # type: ignore[arg-type]
        )
        result._freq = new_freq
        return cast("Self", self._wrap_setop_result(other, result))

    def _range_intersect(self, other, sort) -> Self:
        # Dispatch to RangeIndex intersection logic.
        left = self._as_range_index
        right = other._as_range_index
        res_i8 = left.intersection(right, sort=sort)
        return self._wrap_range_setop(other, res_i8)

    def _range_union(self, other, sort) -> Self:
        # Dispatch to RangeIndex union logic.
        left = self._as_range_index
        right = other._as_range_index
        res_i8 = left.union(right, sort=sort)
        return self._wrap_range_setop(other, res_i8)

    def _intersection(self, other: Index, sort: bool = False) -> Index:
        """
        intersection specialized to the case with matching dtypes and both non-empty.
        """
        other = cast("DatetimeTimedeltaMixin", other)

        if self._can_range_setop(other):
            return self._range_intersect(other, sort=sort)

        if not self._can_fast_intersect(other):
            result = Index._intersection(self, other, sort=sort)
            # We need to invalidate the freq because Index._intersection
            #  uses _shallow_copy on a view of self._data, which will preserve
            #  self.freq if we're not careful.
            # At this point we should have result.dtype == self.dtype
            #  and type(result) is type(self._data)
            result = self._wrap_setop_result(other, result)
            # error: "Index" has no attribute "_with_freq"; maybe "_with_infer"?
            return result._with_freq(None)._with_freq("infer")  # type: ignore[attr-defined]

        else:
            return self._fast_intersect(other, sort)

    def _fast_intersect(self, other, sort):
        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        # after sorting, the intersection always starts with the right index
        # and ends with the index of which the last elements is smallest
        end = min(left[-1], right[-1])
        start = right[0]

        if end < start:
            result = self[:0]
        else:
            lslice = slice(*left.slice_locs(start, end))
            result = left._values[lslice]
            result._freq = self.freq  # type: ignore[union-attr]

        return result

    def _can_fast_intersect(self, other: Self) -> bool:
        # Note: we only get here with len(self) > 0 and len(other) > 0
        if self.freq is None:
            return False

        elif other.freq != self.freq:
            return False

        elif not self.is_monotonic_increasing:
            # Because freq is not None, we must then be monotonic decreasing
            return False

        # this along with matching freqs ensure that we "line up",
        #  so intersection will preserve freq
        # Note we are assuming away Ticks, as those go through _range_intersect
        # GH#42104
        return self.freq.n == 1

    def _can_fast_union(self, other: Self) -> bool:
        # Assumes that type(self) == type(other), as per the annotation
        # The ability to fast_union also implies that `freq` should be
        #  retained on union.
        freq = self.freq

        if freq is None or freq != other.freq:
            return False

        if not self.is_monotonic_increasing:
            # Because freq is not None, we must then be monotonic decreasing
            # TODO: do union on the reversed indexes?
            return False

        if len(self) == 0 or len(other) == 0:
            # only reached via union_many
            return True

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        right_start = right[0]
        left_end = left[-1]

        # Only need to "adjoin", not overlap
        return (right_start == left_end + freq) or right_start in left

    def _fast_union(self, other: Self, sort=None) -> Self:
        # Caller is responsible for ensuring self and other are non-empty

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        elif sort is False:
            # TDIs are not in the "correct" order and we don't want
            #  to sort but want to remove overlaps
            left, right = self, other
            left_start = left[0]
            loc = right.searchsorted(left_start, side="left")
            right_chunk = right._values[:loc]
            dates = concat_compat((left._values, right_chunk))
            result = type(self)._simple_new(dates, name=self.name)
            result._data._freq = self.freq
            return result
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # concatenate
        if left_end < right_end:
            loc = right.searchsorted(left_end, side="right")
            right_chunk = right._values[loc:]
            dates = concat_compat([left._values, right_chunk])
            # The can_fast_union check ensures that the result.freq
            #  should match self.freq
            assert isinstance(dates, type(self._data))
            dates._freq = self.freq  # type: ignore[union-attr]
            result = type(self)._simple_new(dates)
            return result
        else:
            return left

    def _union(self, other, sort):
        # We are called by `union`, which is responsible for this validation
        assert isinstance(other, type(self))
        assert self.dtype == other.dtype

        if self._can_range_setop(other):
            return self._range_union(other, sort=sort)

        if self._can_fast_union(other):
            # in the case with sort=None, the _can_fast_union check ensures
            #  that result.freq == self.freq
            return self._fast_union(other, sort=sort)
        else:
            # super()._union can return an ArrayLike; wrap into an Index first
            result = self._wrap_setop_result(other, super()._union(other, sort))
            return result._with_freq("infer")  # type: ignore[attr-defined]

    # --------------------------------------------------------------------
    # Join Methods

    def _get_join_freq(self, other):
        """
        Get the freq to attach to the result of a join operation.
        """
        freq = None
        if self._can_fast_union(other):
            freq = self.freq
        return freq

    def _wrap_join_result(
        self,
        joined,
        other,
        lidx: npt.NDArray[np.intp] | None,
        ridx: npt.NDArray[np.intp] | None,
        how: JoinHow,
    ) -> tuple[Self, npt.NDArray[np.intp] | None, npt.NDArray[np.intp] | None]:
        assert other.dtype == self.dtype, (other.dtype, self.dtype)
        join_index, lidx, ridx = super()._wrap_join_result(
            joined, other, lidx, ridx, how
        )
        join_index._data._freq = self._get_join_freq(other)
        return join_index, lidx, ridx

    def _get_engine_target(self) -> np.ndarray:
        # engine methods and libjoin methods need dt64/td64 values cast to i8
        return self._data._ndarray.view("i8")

    def _from_join_target(self, result: np.ndarray):
        # view e.g. i8 back to M8[ns]
        result = result.view(self._data._ndarray.dtype)
        return self._data._from_backing_data(result)

    def _searchsorted_monotonic(self, label, side: Literal["left", "right"] = "left"):
        if (
            self.is_monotonic_increasing
            and isinstance(label, (Timestamp, Timedelta))
            and abbrev_to_npy_unit(label.unit) > abbrev_to_npy_unit(self.unit)
        ):
            # For non-matching units we can safely round down (with side=right)
            # This is needed for GH#63262
            if side == "right":
                label = label.as_unit(self.unit)  # this should always be a round-down
            else:
                # round up
                label = label.ceil(self.unit).as_unit(self.unit)

        return super()._searchsorted_monotonic(label, side)

    # --------------------------------------------------------------------
    # List-like Methods

    def _getitem_slice(self, slobj: slice) -> Self:
        result = super()._getitem_slice(slobj)
        result._data._freq = self._get_getitem_freq(slobj)
        return result

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if isinstance(result, type(self)):
            result._data._freq = self._get_getitem_freq(key)
        return result

    def _get_getitem_freq(self, key) -> BaseOffset | None:
        """
        Find the `freq` attribute to assign to the result of a __getitem__ lookup.
        """
        key = check_array_indexer(self._data, key)  # maybe ndarray[bool] -> slice
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
        elif com.is_bool_indexer(key):
            new_key = lib.maybe_booleans_to_slice(key.view(np.uint8))
            if isinstance(new_key, slice):
                return self._get_getitem_freq(new_key)
        return freq

    def _concat(self, to_concat: list[Index], name: Hashable) -> Index:
        result = super()._concat(to_concat, name)

        # GH#3232: If the concat result is evenly spaced, we can retain the
        # original frequency
        obj = cast("DatetimeTimedeltaMixin", to_concat[0])
        to_concat_nonempty = cast(
            "list[DatetimeTimedeltaMixin]",
            [idx for idx in to_concat if len(idx)],
        )
        if (
            isinstance(result, type(self))
            and obj.freq is not None
            and all(idx.freq == obj.freq for idx in to_concat_nonempty)
        ):
            pairs = pairwise(to_concat_nonempty)
            if all(pair[0][-1] + obj.freq == pair[1][0] for pair in pairs):
                result._data._freq = obj.freq

        return result

    def _get_delete_freq(self, loc: int | slice | Sequence[int]):
        """
        Find the `freq` for self.delete(loc).
        """
        freq = None
        if self.freq is not None:
            if is_integer(loc):
                if loc in (0, -len(self), -1, len(self) - 1):
                    freq = self.freq
            else:
                if is_list_like(loc):
                    # error: Incompatible types in assignment (expression has
                    # type "Union[slice, ndarray]", variable has type
                    # "Union[int, slice, Sequence[int]]")
                    loc = lib.maybe_indices_to_slice(  # type: ignore[assignment]
                        np.asarray(loc, dtype=np.intp), len(self)
                    )
                if isinstance(loc, slice) and loc.step in (1, None):
                    if loc.start in (0, None) or loc.stop in (len(self), None):
                        freq = self.freq
        return freq

    def _get_insert_freq(self, loc: int, item):
        """
        Find the `freq` for self.insert(loc, item).
        """
        value = self._data._validate_scalar(item)
        item = self._data._box_func(value)

        freq = None
        if self.freq is not None:
            # freq can be preserved on edge cases
            if self.size:
                if item is NaT:
                    pass
                elif loc in (0, -len(self)) and item + self.freq == self[0]:
                    freq = self.freq
                elif (loc == len(self)) and item - self.freq == self[-1]:
                    freq = self.freq
            # Adding a single item to an empty index may preserve freq
            elif isinstance(self.freq, Tick):
                # all TimedeltaIndex cases go through here; is_on_offset
                #  would raise TypeError
                freq = self.freq
            elif self.freq.is_on_offset(item):
                freq = self.freq
        return freq

    def delete(self, loc) -> Self:
        """
        Make new Index with passed location(-s) deleted.

        Parameters
        ----------
        loc : int or list of int
            Location of item(-s) which will be deleted.
            Use a list of locations to delete more than one value at the same time.

        Returns
        -------
        Index
            Will be same type as self, except for RangeIndex.

        See Also
        --------
        numpy.delete : Delete any rows and column from NumPy array (ndarray).

        Examples
        --------
        >>> idx = pd.Index(["a", "b", "c"])
        >>> idx.delete(1)
        Index(['a', 'c'], dtype='str')
        >>> idx = pd.Index(["a", "b", "c"])
        >>> idx.delete([0, 2])
        Index(['b'], dtype='str')
        """
        result = super().delete(loc)
        result._data._freq = self._get_delete_freq(loc)
        return result

    def insert(self, loc: int, item):
        """
        Make new Index inserting new item at location.
        Follows Python numpy.insert semantics for negative values.

        Parameters
        ----------
        loc : int
            The integer location where the new item will be inserted.
        item : object
            The new item to be inserted into the Index.

        Returns
        -------
        Index
            Returns a new Index object resulting from inserting the specified item at
            the specified location within the original Index.

        See Also
        --------
        Index.append : Append a collection of Indexes together.

        Examples
        --------
        >>> idx = pd.Index(["a", "b", "c"])
        >>> idx.insert(1, "x")
        Index(['a', 'x', 'b', 'c'], dtype='str')
        """
        result = super().insert(loc, item)
        if isinstance(result, type(self)):
            # i.e. parent class method did not cast
            result._data._freq = self._get_insert_freq(loc, item)
        return result

    # --------------------------------------------------------------------
    # NDArray-Like Methods

    def take(
        self,
        indices,
        axis: Axis = 0,
        allow_fill: bool = True,
        fill_value=None,
        **kwargs,
    ) -> Self:
        """
        Return a new Index of the values selected by the indices.
        For internal compatibility with numpy arrays.

        Parameters
        ----------
        indices : array-like
            Indices to be taken.
        axis : {0 or 'index'}, optional
            The axis over which to select values, always 0 or 'index'.
        allow_fill : bool, default True
            How to handle negative values in `indices`.
            * False: negative values in `indices` indicate positional indices
                from the right (the default). This is similar to
                :func:`numpy.take`.
            * True: negative values in `indices` indicate
                missing values. These values are set to `fill_value`. Any other
                other negative values raise a ``ValueError``.
        fill_value : scalar, default None
            If allow_fill=True and fill_value is not None, indices specified by
            -1 are regarded as NA. If Index doesn't hold NA, raise ValueError.
        **kwargs
            Required for compatibility with numpy.

        Returns
        -------
        Index
            An index formed of elements at the given indices. Will be the same
            type as self, except for RangeIndex.

        See Also
        --------
        numpy.ndarray.take: Return an array formed from the
            elements of a at the given indices.

        Examples
        --------
        >>> idx = pd.Index(["a", "b", "c"])
        >>> idx.take([2, 2, 1, 2])
        Index(['c', 'c', 'b', 'c'], dtype='str')
        """
        nv.validate_take((), kwargs)
        indices = np.asarray(indices, dtype=np.intp)

        result = NDArrayBackedExtensionIndex.take(
            self, indices, axis, allow_fill, fill_value, **kwargs
        )

        maybe_slice = lib.maybe_indices_to_slice(indices, len(self))
        if isinstance(maybe_slice, slice):
            freq = self._get_getitem_freq(maybe_slice)
            result._data._freq = freq
        return result
