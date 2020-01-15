"""
Base and utility classes for tseries type pandas objects.
"""
import operator
from typing import List, Optional, Set

import numpy as np

from pandas._libs import NaT, iNaT, join as libjoin, lib
from pandas._libs.algos import unique_deltas
from pandas._libs.tslibs import timezones
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.common import (
    ensure_int64,
    is_bool_dtype,
    is_categorical_dtype,
    is_dtype_equal,
    is_float,
    is_integer,
    is_list_like,
    is_period_dtype,
    is_scalar,
    needs_i8_conversion,
)
from pandas.core.dtypes.concat import concat_compat
from pandas.core.dtypes.generic import ABCIndex, ABCIndexClass, ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core import algorithms
from pandas.core.accessor import PandasDelegate
from pandas.core.arrays import DatetimeArray, ExtensionArray, TimedeltaArray
from pandas.core.arrays.datetimelike import DatetimeLikeArrayMixin
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import Index, _index_shared_docs
from pandas.core.indexes.extension import (
    ExtensionIndex,
    inherit_names,
    make_wrapped_arith_op,
)
from pandas.core.indexes.numeric import Int64Index
from pandas.core.ops import get_op_result_name
from pandas.core.tools.timedeltas import to_timedelta

from pandas.tseries.frequencies import DateOffset, to_offset

_index_doc_kwargs = dict(ibase._index_doc_kwargs)


def _join_i8_wrapper(joinf, with_indexers: bool = True):
    """
    Create the join wrapper methods.
    """

    @staticmethod  # type: ignore
    def wrapper(left, right):
        if isinstance(left, (np.ndarray, ABCIndex, ABCSeries, DatetimeLikeArrayMixin)):
            left = left.view("i8")
        if isinstance(right, (np.ndarray, ABCIndex, ABCSeries, DatetimeLikeArrayMixin)):
            right = right.view("i8")

        results = joinf(left, right)
        if with_indexers:
            # dtype should be timedelta64[ns] for TimedeltaIndex
            #  and datetime64[ns] for DatetimeIndex
            dtype = left.dtype.base

            join_index, left_indexer, right_indexer = results
            join_index = join_index.view(dtype)
            return join_index, left_indexer, right_indexer
        return results

    return wrapper


@inherit_names(
    ["inferred_freq", "_isnan", "_resolution", "resolution"],
    DatetimeLikeArrayMixin,
    cache=True,
)
@inherit_names(
    ["__iter__", "mean", "freq", "freqstr", "_ndarray_values", "asi8", "_box_values"],
    DatetimeLikeArrayMixin,
)
class DatetimeIndexOpsMixin(ExtensionIndex):
    """
    Common ops mixin to support a unified interface datetimelike Index.
    """

    _data: ExtensionArray
    freq: Optional[DateOffset]
    freqstr: Optional[str]
    _resolution: int
    _bool_ops: List[str] = []
    _field_ops: List[str] = []

    hasnans = cache_readonly(DatetimeLikeArrayMixin._hasnans.fget)  # type: ignore
    _hasnans = hasnans  # for index / array -agnostic code

    @property
    def is_all_dates(self) -> bool:
        return True

    # ------------------------------------------------------------------------
    # Abstract data attributes

    @property
    def values(self):
        # Note: PeriodArray overrides this to return an ndarray of objects.
        return self._data._data

    def __array_wrap__(self, result, context=None):
        """
        Gets called after a ufunc.
        """
        result = lib.item_from_zerodim(result)
        if is_bool_dtype(result) or lib.is_scalar(result):
            return result

        attrs = self._get_attributes_dict()
        if not is_period_dtype(self) and attrs["freq"]:
            # no need to infer if freq is None
            attrs["freq"] = "infer"
        return Index(result, **attrs)

    # ------------------------------------------------------------------------

    def equals(self, other) -> bool:
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, ABCIndexClass):
            return False
        elif not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except (ValueError, TypeError, OverflowError):
                # e.g.
                #  ValueError -> cannot parse str entry, or OutOfBoundsDatetime
                #  TypeError  -> trying to convert IntervalIndex to DatetimeIndex
                #  OverflowError -> Index([very_large_timedeltas])
                return False

        if not is_dtype_equal(self.dtype, other.dtype):
            # have different timezone
            return False

        return np.array_equal(self.asi8, other.asi8)

    @Appender(_index_shared_docs["contains"] % _index_doc_kwargs)
    def __contains__(self, key):
        try:
            res = self.get_loc(key)
        except (KeyError, TypeError, ValueError):
            return False
        return bool(
            is_scalar(res) or isinstance(res, slice) or (is_list_like(res) and len(res))
        )

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, mapper, na_action=None):
        try:
            result = mapper(self)

            # Try to use this result if we can
            if isinstance(result, np.ndarray):
                result = Index(result)

            if not isinstance(result, Index):
                raise TypeError("The map function must return an Index object")
            return result
        except Exception:
            return self.astype(object).map(mapper)

    def sort_values(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index.
        """
        if return_indexer:
            _as = self.argsort()
            if not ascending:
                _as = _as[::-1]
            sorted_index = self.take(_as)
            return sorted_index, _as
        else:
            # NB: using asi8 instead of _ndarray_values matters in numpy 1.18
            #  because the treatment of NaT has been changed to put NaT last
            #  instead of first.
            sorted_values = np.sort(self.asi8)
            attribs = self._get_attributes_dict()
            freq = attribs["freq"]

            if freq is not None and not is_period_dtype(self):
                if freq.n > 0 and not ascending:
                    freq = freq * -1
                elif freq.n < 0 and ascending:
                    freq = freq * -1
            attribs["freq"] = freq

            if not ascending:
                sorted_values = sorted_values[::-1]

            return self._simple_new(sorted_values, **attribs)

    @Appender(_index_shared_docs["take"] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = ensure_int64(indices)

        maybe_slice = lib.maybe_indices_to_slice(indices, len(self))
        if isinstance(maybe_slice, slice):
            return self[maybe_slice]

        return ExtensionIndex.take(
            self, indices, axis, allow_fill, fill_value, **kwargs
        )

    _can_hold_na = True

    _na_value = NaT
    """The expected NA value to use with this index."""

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(to_timedelta(tolerance).to_numpy())

        if target.size != tolerance.size and tolerance.size > 1:
            raise ValueError("list-like tolerance size must match target index size")
        return tolerance

    def tolist(self) -> List:
        """
        Return a list of the underlying data.
        """
        return list(self.astype(object))

    def min(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the minimum value of the Index or minimum along
        an axis.

        See Also
        --------
        numpy.ndarray.min
        Series.min : Return the minimum value in a Series.
        """
        nv.validate_min(args, kwargs)
        nv.validate_minmax_axis(axis)

        if not len(self):
            return self._na_value

        i8 = self.asi8
        try:
            # quick check
            if len(i8) and self.is_monotonic:
                if i8[0] != iNaT:
                    return self._box_func(i8[0])

            if self.hasnans:
                if skipna:
                    min_stamp = self[~self._isnan].asi8.min()
                else:
                    return self._na_value
            else:
                min_stamp = i8.min()
            return self._box_func(min_stamp)
        except ValueError:
            return self._na_value

    def argmin(self, axis=None, skipna=True, *args, **kwargs):
        """
        Returns the indices of the minimum values along an axis.

        See `numpy.ndarray.argmin` for more information on the
        `axis` parameter.

        See Also
        --------
        numpy.ndarray.argmin
        """
        nv.validate_argmin(args, kwargs)
        nv.validate_minmax_axis(axis)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all() or not skipna:
                return -1
            i8 = i8.copy()
            i8[mask] = np.iinfo("int64").max
        return i8.argmin()

    def max(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the maximum value of the Index or maximum along
        an axis.

        See Also
        --------
        numpy.ndarray.max
        Series.max : Return the maximum value in a Series.
        """
        nv.validate_max(args, kwargs)
        nv.validate_minmax_axis(axis)

        if not len(self):
            return self._na_value

        i8 = self.asi8
        try:
            # quick check
            if len(i8) and self.is_monotonic:
                if i8[-1] != iNaT:
                    return self._box_func(i8[-1])

            if self.hasnans:
                if skipna:
                    max_stamp = self[~self._isnan].asi8.max()
                else:
                    return self._na_value
            else:
                max_stamp = i8.max()
            return self._box_func(max_stamp)
        except ValueError:
            return self._na_value

    def argmax(self, axis=None, skipna=True, *args, **kwargs):
        """
        Returns the indices of the maximum values along an axis.

        See `numpy.ndarray.argmax` for more information on the
        `axis` parameter.

        See Also
        --------
        numpy.ndarray.argmax
        """
        nv.validate_argmax(args, kwargs)
        nv.validate_minmax_axis(axis)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all() or not skipna:
                return -1
            i8 = i8.copy()
            i8[mask] = 0
        return i8.argmax()

    # --------------------------------------------------------------------
    # Rendering Methods

    def _format_with_header(self, header, na_rep="NaT", **kwargs):
        return header + list(self._format_native_types(na_rep, **kwargs))

    @property
    def _formatter_func(self):
        raise AbstractMethodError(self)

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value).
        """
        attrs = super()._format_attrs()
        for attrib in self._attributes:
            if attrib == "freq":
                freq = self.freqstr
                if freq is not None:
                    freq = repr(freq)
                attrs.append(("freq", freq))
        return attrs

    # --------------------------------------------------------------------

    def _convert_scalar_indexer(self, key, kind=None):
        """
        We don't allow integer or float indexing on datetime-like when using
        loc.

        Parameters
        ----------
        key : label of the slice bound
        kind : {'loc', 'getitem', 'iloc'} or None
        """

        assert kind in ["loc", "getitem", "iloc", None]

        # we don't allow integer/float indexing for loc
        # we don't allow float indexing for ix/getitem
        if is_scalar(key):
            is_int = is_integer(key)
            is_flt = is_float(key)
            if kind in ["loc"] and (is_int or is_flt):
                self._invalid_indexer("index", key)
            elif kind in ["getitem"] and is_flt:
                self._invalid_indexer("index", key)

        return super()._convert_scalar_indexer(key, kind=kind)

    __add__ = make_wrapped_arith_op("__add__")
    __radd__ = make_wrapped_arith_op("__radd__")
    __sub__ = make_wrapped_arith_op("__sub__")
    __rsub__ = make_wrapped_arith_op("__rsub__")
    __pow__ = make_wrapped_arith_op("__pow__")
    __rpow__ = make_wrapped_arith_op("__rpow__")
    __mul__ = make_wrapped_arith_op("__mul__")
    __rmul__ = make_wrapped_arith_op("__rmul__")
    __floordiv__ = make_wrapped_arith_op("__floordiv__")
    __rfloordiv__ = make_wrapped_arith_op("__rfloordiv__")
    __mod__ = make_wrapped_arith_op("__mod__")
    __rmod__ = make_wrapped_arith_op("__rmod__")
    __divmod__ = make_wrapped_arith_op("__divmod__")
    __rdivmod__ = make_wrapped_arith_op("__rdivmod__")
    __truediv__ = make_wrapped_arith_op("__truediv__")
    __rtruediv__ = make_wrapped_arith_op("__rtruediv__")

    def isin(self, values, level=None):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values.

        Parameters
        ----------
        values : set or sequence of values

        Returns
        -------
        is_contained : ndarray (boolean dtype)
        """
        if level is not None:
            self._validate_index_level(level)

        if not isinstance(values, type(self)):
            try:
                values = type(self)(values)
            except ValueError:
                return self.astype(object).isin(values)

        return algorithms.isin(self.asi8, values.asi8)

    @Appender(_index_shared_docs["where"] % _index_doc_kwargs)
    def where(self, cond, other=None):
        values = self.view("i8")

        if is_scalar(other) and isna(other):
            other = NaT.value

        else:
            # Do type inference if necessary up front
            # e.g. we passed PeriodIndex.values and got an ndarray of Periods
            other = Index(other)

            if is_categorical_dtype(other):
                # e.g. we have a Categorical holding self.dtype
                if needs_i8_conversion(other.categories):
                    other = other._internal_get_values()

            if not is_dtype_equal(self.dtype, other.dtype):
                raise TypeError(f"Where requires matching dtype, not {other.dtype}")

            other = other.view("i8")

        result = np.where(cond, values, other).astype("i8")
        return self._shallow_copy(result)

    def _summary(self, name=None):
        """
        Return a summarized representation.

        Parameters
        ----------
        name : str
            Name to use in the summary representation.

        Returns
        -------
        str
            Summarized representation of the index.
        """
        formatter = self._formatter_func
        if len(self) > 0:
            index_summary = f", {formatter(self[0])} to {formatter(self[-1])}"
        else:
            index_summary = ""

        if name is None:
            name = type(self).__name__
        result = f"{name}: {len(self)} entries{index_summary}"
        if self.freq:
            result += f"\nFreq: {self.freqstr}"

        # display as values, not quoted
        result = result.replace("'", "")
        return result

    def _concat_same_dtype(self, to_concat, name):
        """
        Concatenate to_concat which has the same class.
        """
        attribs = self._get_attributes_dict()
        attribs["name"] = name
        # do not pass tz to set because tzlocal cannot be hashed
        if len({str(x.dtype) for x in to_concat}) != 1:
            raise ValueError("to_concat must have the same tz")

        new_data = type(self._values)._concat_same_type(to_concat).asi8

        # GH 3232: If the concat result is evenly spaced, we can retain the
        # original frequency
        is_diff_evenly_spaced = len(unique_deltas(new_data)) == 1
        if not is_period_dtype(self) and not is_diff_evenly_spaced:
            # reset freq
            attribs["freq"] = None

        return self._simple_new(new_data, **attribs)

    def shift(self, periods=1, freq=None):
        """
        Shift index by desired number of time frequency increments.

        This method is for shifting the values of datetime-like indexes
        by a specified time increment a given number of times.

        Parameters
        ----------
        periods : int, default 1
            Number of periods (or increments) to shift by,
            can be positive or negative.

            .. versionchanged:: 0.24.0

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
        result = self._data._time_shift(periods, freq=freq)
        return type(self)(result, name=self.name)

    # --------------------------------------------------------------------
    # List-like Methods

    def delete(self, loc):
        new_i8s = np.delete(self.asi8, loc)

        freq = None
        if is_period_dtype(self):
            freq = self.freq
        elif is_integer(loc):
            if loc in (0, -len(self), -1, len(self) - 1):
                freq = self.freq
        else:
            if is_list_like(loc):
                loc = lib.maybe_indices_to_slice(ensure_int64(np.array(loc)), len(self))
            if isinstance(loc, slice) and loc.step in (1, None):
                if loc.start in (0, None) or loc.stop in (len(self), None):
                    freq = self.freq

        return self._shallow_copy(new_i8s, freq=freq)


class DatetimeTimedeltaMixin(DatetimeIndexOpsMixin, Int64Index):
    """
    Mixin class for methods shared by DatetimeIndex and TimedeltaIndex,
    but not PeriodIndex
    """

    # Compat for frequency inference, see GH#23789
    _is_monotonic_increasing = Index.is_monotonic_increasing
    _is_monotonic_decreasing = Index.is_monotonic_decreasing
    _is_unique = Index.is_unique

    def _set_freq(self, freq):
        """
        Set the _freq attribute on our underlying DatetimeArray.

        Parameters
        ----------
        freq : DateOffset, None, or "infer"
        """
        # GH#29843
        if freq is None:
            # Always valid
            pass
        elif len(self) == 0 and isinstance(freq, DateOffset):
            # Always valid.  In the TimedeltaIndex case, we assume this
            #  is a Tick offset.
            pass
        else:
            # As an internal method, we can ensure this assertion always holds
            assert freq == "infer"
            freq = to_offset(self.inferred_freq)

        self._data._freq = freq

    def _shallow_copy(self, values=None, **kwargs):
        if values is None:
            values = self._data
        if isinstance(values, type(self)):
            values = values._data

        attributes = self._get_attributes_dict()

        if "freq" not in kwargs and self.freq is not None:
            if isinstance(values, (DatetimeArray, TimedeltaArray)):
                if values.freq is None:
                    del attributes["freq"]

        attributes.update(kwargs)
        return self._simple_new(values, **attributes)

    # --------------------------------------------------------------------
    # Set Operation Methods

    @Appender(Index.difference.__doc__)
    def difference(self, other, sort=None):
        new_idx = super().difference(other, sort=sort)
        new_idx._set_freq(None)
        return new_idx

    def intersection(self, other, sort=False):
        """
        Specialized intersection for DatetimeIndex/TimedeltaIndex.

        May be much faster than Index.intersection

        Parameters
        ----------
        other : Same type as self or array-like
        sort : False or None, default False
            Sort the resulting index if possible.

            .. versionadded:: 0.24.0

            .. versionchanged:: 0.24.1

               Changed the default to ``False`` to match the behaviour
               from before 0.24.0.

            .. versionchanged:: 0.25.0

               The `sort` keyword is added

        Returns
        -------
        y : Index or same type as self
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)

        if self.equals(other):
            return self._get_reconciled_name_object(other)

        if len(self) == 0:
            return self.copy()
        if len(other) == 0:
            return other.copy()

        if not isinstance(other, type(self)):
            result = Index.intersection(self, other, sort=sort)
            if isinstance(result, type(self)):
                if result.freq is None:
                    result._set_freq("infer")
            return result

        elif (
            other.freq is None
            or self.freq is None
            or other.freq != self.freq
            or not other.freq.is_anchored()
            or (not self.is_monotonic or not other.is_monotonic)
        ):
            result = Index.intersection(self, other, sort=sort)

            # Invalidate the freq of `result`, which may not be correct at
            # this point, depending on the values.

            result._set_freq(None)
            result = self._shallow_copy(
                result._data, name=result.name, dtype=result.dtype, freq=None
            )
            if result.freq is None:
                result._set_freq("infer")
            return result

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
            return type(self)(data=[])
        else:
            lslice = slice(*left.slice_locs(start, end))
            left_chunk = left.values[lslice]
            return self._shallow_copy(left_chunk)

    def _can_fast_union(self, other) -> bool:
        if not isinstance(other, type(self)):
            return False

        freq = self.freq

        if freq is None or freq != other.freq:
            return False

        if not self.is_monotonic or not other.is_monotonic:
            return False

        if len(self) == 0 or len(other) == 0:
            return True

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        right_start = right[0]
        left_end = left[-1]

        # Only need to "adjoin", not overlap
        try:
            return (right_start == left_end + freq) or right_start in left
        except ValueError:
            # if we are comparing a freq that does not propagate timezones
            # this will raise
            return False

    def _fast_union(self, other, sort=None):
        if len(other) == 0:
            return self.view(type(self))

        if len(self) == 0:
            return other.view(type(self))

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        elif sort is False:
            # TDIs are not in the "correct" order and we don't want
            #  to sort but want to remove overlaps
            left, right = self, other
            left_start = left[0]
            loc = right.searchsorted(left_start, side="left")
            right_chunk = right.values[:loc]
            dates = concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # concatenate
        if left_end < right_end:
            loc = right.searchsorted(left_end, side="right")
            right_chunk = right.values[loc:]
            dates = concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        else:
            return left

    def _union(self, other, sort):
        if not len(other) or self.equals(other) or not len(self):
            return super()._union(other, sort=sort)

        # We are called by `union`, which is responsible for this validation
        assert isinstance(other, type(self))

        this, other = self._maybe_utc_convert(other)

        if this._can_fast_union(other):
            return this._fast_union(other, sort=sort)
        else:
            result = Index._union(this, other, sort=sort)
            if isinstance(result, type(self)):
                assert result._data.dtype == this.dtype
                if result.freq is None:
                    result._set_freq("infer")
            return result

    # --------------------------------------------------------------------
    # Join Methods
    _join_precedence = 10

    _inner_indexer = _join_i8_wrapper(libjoin.inner_join_indexer)
    _outer_indexer = _join_i8_wrapper(libjoin.outer_join_indexer)
    _left_indexer = _join_i8_wrapper(libjoin.left_join_indexer)
    _left_indexer_unique = _join_i8_wrapper(
        libjoin.left_join_indexer_unique, with_indexers=False
    )

    def join(
        self, other, how: str = "left", level=None, return_indexers=False, sort=False
    ):
        """
        See Index.join
        """
        if self._is_convertible_to_index_for_join(other):
            try:
                other = type(self)(other)
            except (TypeError, ValueError):
                pass

        this, other = self._maybe_utc_convert(other)
        return Index.join(
            this,
            other,
            how=how,
            level=level,
            return_indexers=return_indexers,
            sort=sort,
        )

    def _maybe_utc_convert(self, other):
        this = self
        if not hasattr(self, "tz"):
            return this, other

        if isinstance(other, type(self)):
            if self.tz is not None:
                if other.tz is None:
                    raise TypeError("Cannot join tz-naive with tz-aware DatetimeIndex")
            elif other.tz is not None:
                raise TypeError("Cannot join tz-naive with tz-aware DatetimeIndex")

            if not timezones.tz_compare(self.tz, other.tz):
                this = self.tz_convert("UTC")
                other = other.tz_convert("UTC")
        return this, other

    @classmethod
    def _is_convertible_to_index_for_join(cls, other: Index) -> bool:
        """
        return a boolean whether I can attempt conversion to a
        DatetimeIndex/TimedeltaIndex
        """
        if isinstance(other, cls):
            return False
        elif len(other) > 0 and other.inferred_type not in (
            "floating",
            "mixed-integer",
            "integer",
            "integer-na",
            "mixed-integer-float",
            "mixed",
        ):
            return True
        return False

    def _wrap_joined_index(self, joined, other):
        name = get_op_result_name(self, other)
        if self._can_fast_union(other):
            joined = self._shallow_copy(joined)
            joined.name = name
            return joined
        else:
            kwargs = {}
            if hasattr(self, "tz"):
                kwargs["tz"] = getattr(other, "tz", None)
            return self._simple_new(joined, name, **kwargs)


class DatetimelikeDelegateMixin(PandasDelegate):
    """
    Delegation mechanism, specific for Datetime, Timedelta, and Period types.

    Functionality is delegated from the Index class to an Array class. A
    few things can be customized

    * _delegated_methods, delegated_properties : List
        The list of property / method names being delagated.
    * raw_methods : Set
        The set of methods whose results should should *not* be
        boxed in an index, after being returned from the array
    * raw_properties : Set
        The set of properties whose results should should *not* be
        boxed in an index, after being returned from the array
    """

    # raw_methods : dispatch methods that shouldn't be boxed in an Index
    _raw_methods: Set[str] = set()
    # raw_properties : dispatch properties that shouldn't be boxed in an Index
    _raw_properties: Set[str] = set()
    _data: ExtensionArray

    def _delegate_property_get(self, name, *args, **kwargs):
        result = getattr(self._data, name)
        if name not in self._raw_properties:
            result = Index(result, name=self.name)
        return result

    def _delegate_property_set(self, name, value, *args, **kwargs):
        setattr(self._data, name, value)

    def _delegate_method(self, name, *args, **kwargs):
        result = operator.methodcaller(name, *args, **kwargs)(self._data)
        if name not in self._raw_methods:
            result = Index(result, name=self.name)
        return result
