# pylint: disable=E1101
from datetime import datetime, time, timedelta
import operator
import warnings

import numpy as np

from pandas._libs import (
    Timestamp, index as libindex, join as libjoin, lib, tslib as libts)
from pandas._libs.tslibs import ccalendar, fields, parsing, timezones
from pandas.util._decorators import Appender, Substitution, cache_readonly

from pandas.core.dtypes.common import (
    _NS_DTYPE, ensure_int64, is_float, is_integer, is_list_like, is_scalar,
    is_string_like)
import pandas.core.dtypes.concat as _concat
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import isna

from pandas.core.accessor import delegate_names
from pandas.core.arrays.datetimes import (
    DatetimeArray, _to_M8, tz_to_dtype, validate_tz_from_dtype)
from pandas.core.base import _shared_docs
import pandas.core.common as com
from pandas.core.indexes.base import Index
from pandas.core.indexes.datetimelike import (
    DatetimeIndexOpsMixin, DatetimelikeDelegateMixin, ea_passthrough)
from pandas.core.indexes.numeric import Int64Index
from pandas.core.ops import get_op_result_name
import pandas.core.tools.datetimes as tools

from pandas.tseries.frequencies import Resolution, to_offset
from pandas.tseries.offsets import CDay, Nano, prefix_mapping


def _new_DatetimeIndex(cls, d):
    """ This is called upon unpickling, rather than the default which doesn't
    have arguments and breaks __new__ """

    if "data" in d and not isinstance(d["data"], DatetimeIndex):
        # Avoid need to verify integrity by calling simple_new directly
        data = d.pop("data")
        result = cls._simple_new(data, **d)
    else:
        with warnings.catch_warnings():
            # we ignore warnings from passing verify_integrity=False
            # TODO: If we knew what was going in to **d, we might be able to
            #  go through _simple_new instead
            warnings.simplefilter("ignore")
            result = cls.__new__(cls, verify_integrity=False, **d)

    return result


class DatetimeDelegateMixin(DatetimelikeDelegateMixin):
    # Most attrs are dispatched via datetimelike_{ops,methods}
    # Some are "raw" methods, the result is not not re-boxed in an Index
    # We also have a few "extra" attrs, which may or may not be raw,
    # which we we dont' want to expose in the .dt accessor.
    _extra_methods = [
        'to_period',
        'to_perioddelta',
        'to_julian_date',
    ]
    _extra_raw_methods = [
        'to_pydatetime',
        '_local_timestamps',
        '_has_same_tz',
    ]
    _extra_raw_properties = [
        '_box_func',
        'tz', 'tzinfo',
    ]
    _delegated_properties = (
        DatetimeArray._datetimelike_ops + _extra_raw_properties
    )
    _delegated_methods = (
        DatetimeArray._datetimelike_methods + _extra_methods +
        _extra_raw_methods
    )
    _raw_properties = {
        'date',
        'time',
        'timetz',
    } | set(DatetimeArray._bool_ops) | set(_extra_raw_properties)
    _raw_methods = set(_extra_raw_methods)
    _delegate_class = DatetimeArray


@delegate_names(DatetimeArray,
                DatetimeDelegateMixin._delegated_properties,
                typ="property")
@delegate_names(DatetimeArray,
                DatetimeDelegateMixin._delegated_methods,
                typ="method", overwrite=False)
class DatetimeIndex(DatetimeIndexOpsMixin, Int64Index, DatetimeDelegateMixin):
    """
    Immutable ndarray of datetime64 data, represented internally as int64, and
    which can be boxed to Timestamp objects that are subclasses of datetime and
    carry metadata such as frequency information.

    Parameters
    ----------
    data  : array-like (1-dimensional), optional
        Optional datetime-like data to construct index with
    copy  : bool
        Make a copy of input ndarray
    freq : string or pandas offset object, optional
        One of pandas date offset strings or corresponding objects. The string
        'infer' can be passed in order to set the frequency of the index as the
        inferred frequency upon creation

    start : starting value, datetime-like, optional
        If data is None, start is used as the start point in generating regular
        timestamp data.

        .. deprecated:: 0.24.0

    periods  : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument

        .. deprecated:: 0.24.0

    end : end time, datetime-like, optional
        If periods is none, generated index will extend to first conforming
        time on or just past end argument

        .. deprecated:: 0.24.0

    closed : string or None, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)

        .. deprecated:: 0.24. 0

    tz : pytz.timezone or dateutil.tz.tzfile
    ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
        When clocks moved backward due to DST, ambiguous times may arise.
        For example in Central European Time (UTC+01), when going from 03:00
        DST to 02:00 non-DST, 02:30:00 local time occurs both at 00:30:00 UTC
        and at 01:30:00 UTC. In such a situation, the `ambiguous` parameter
        dictates how ambiguous times should be handled.

        - 'infer' will attempt to infer fall dst-transition hours based on
          order
        - bool-ndarray where True signifies a DST time, False signifies a
          non-DST time (note that this flag is only applicable for ambiguous
          times)
        - 'NaT' will return NaT where there are ambiguous times
        - 'raise' will raise an AmbiguousTimeError if there are ambiguous times
    name : object
        Name to be stored in the index
    dayfirst : bool, default False
        If True, parse dates in `data` with the day first order
    yearfirst : bool, default False
        If True parse dates in `data` with the year first order

    Attributes
    ----------
    year
    month
    day
    hour
    minute
    second
    microsecond
    nanosecond
    date
    time
    timetz
    dayofyear
    weekofyear
    week
    dayofweek
    weekday
    quarter
    tz
    freq
    freqstr
    is_month_start
    is_month_end
    is_quarter_start
    is_quarter_end
    is_year_start
    is_year_end
    is_leap_year
    inferred_freq

    Methods
    -------
    normalize
    strftime
    snap
    tz_convert
    tz_localize
    round
    floor
    ceil
    to_period
    to_perioddelta
    to_pydatetime
    to_series
    to_frame
    month_name
    day_name

    Notes
    -----
    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Creating a DatetimeIndex based on `start`, `periods`, and `end` has
    been deprecated in favor of :func:`date_range`.

    See Also
    ---------
    Index : The base pandas Index type.
    TimedeltaIndex : Index of timedelta64 data.
    PeriodIndex : Index of Period data.
    to_datetime : Convert argument to datetime.
    date_range : Create a fixed-frequency DatetimeIndex.
    """
    _typ = 'datetimeindex'
    _join_precedence = 10

    def _join_i8_wrapper(joinf, **kwargs):
        return DatetimeIndexOpsMixin._join_i8_wrapper(joinf, dtype='M8[ns]',
                                                      **kwargs)

    _inner_indexer = _join_i8_wrapper(libjoin.inner_join_indexer_int64)
    _outer_indexer = _join_i8_wrapper(libjoin.outer_join_indexer_int64)
    _left_indexer = _join_i8_wrapper(libjoin.left_join_indexer_int64)
    _left_indexer_unique = _join_i8_wrapper(
        libjoin.left_join_indexer_unique_int64, with_indexers=False)

    _engine_type = libindex.DatetimeEngine

    _tz = None
    _freq = None
    _comparables = ['name', 'freqstr', 'tz']
    _attributes = ['name', 'tz', 'freq']

    # dummy attribute so that datetime.__eq__(DatetimeArray) defers
    # by returning NotImplemented
    timetuple = None

    _is_numeric_dtype = False
    _infer_as_myclass = True

    # some things like freq inference make use of these attributes.
    _bool_ops = DatetimeArray._bool_ops
    _object_ops = DatetimeArray._object_ops
    _field_ops = DatetimeArray._field_ops
    _datetimelike_ops = DatetimeArray._datetimelike_ops
    _datetimelike_methods = DatetimeArray._datetimelike_methods

    # --------------------------------------------------------------------
    # Constructors

    def __new__(cls, data=None,
                freq=None, start=None, end=None, periods=None, tz=None,
                normalize=False, closed=None, ambiguous='raise',
                dayfirst=False, yearfirst=False, dtype=None,
                copy=False, name=None, verify_integrity=None):

        if verify_integrity is not None:
            warnings.warn("The 'verify_integrity' argument is deprecated, "
                          "will be removed in a future version.",
                          FutureWarning, stacklevel=2)
        else:
            verify_integrity = True

        if data is None:
            dtarr = DatetimeArray._generate_range(
                start, end, periods,
                freq=freq, tz=tz, normalize=normalize,
                closed=closed, ambiguous=ambiguous)
            warnings.warn("Creating a DatetimeIndex by passing range "
                          "endpoints is deprecated.  Use "
                          "`pandas.date_range` instead.",
                          FutureWarning, stacklevel=2)
            return cls._simple_new(
                dtarr._data, freq=dtarr.freq, tz=dtarr.tz, name=name)

        if is_scalar(data):
            raise TypeError("{cls}() must be called with a "
                            "collection of some kind, {data} was passed"
                            .format(cls=cls.__name__, data=repr(data)))

        # - Cases checked above all return/raise before reaching here - #

        if name is None and hasattr(data, 'name'):
            name = data.name

        dtarr = DatetimeArray._from_sequence(
            data, dtype=dtype, copy=copy, tz=tz, freq=freq,
            dayfirst=dayfirst, yearfirst=yearfirst, ambiguous=ambiguous,
            int_as_wall_time=True)

        subarr = cls._simple_new(dtarr, name=name,
                                 freq=dtarr.freq, tz=dtarr.tz)
        return subarr

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, tz=None, dtype=None):
        """
        we require the we have a dtype compat for the values
        if we are passed a non-dtype compat, then coerce using the constructor
        """
        if isinstance(values, DatetimeArray):
            if tz:
                tz = validate_tz_from_dtype(dtype, tz)
                dtype = DatetimeTZDtype(tz=tz)
            elif dtype is None:
                dtype = _NS_DTYPE

            values = DatetimeArray(values, freq=freq, dtype=dtype)
            tz = values.tz
            freq = values.freq
            values = values._data

        # DatetimeArray._simple_new will accept either i8 or M8[ns] dtypes
        if isinstance(values, DatetimeIndex):
            values = values._data

        dtype = tz_to_dtype(tz)
        dtarr = DatetimeArray._simple_new(values, freq=freq, dtype=dtype)
        assert isinstance(dtarr, DatetimeArray)

        result = object.__new__(cls)
        result._data = dtarr
        result.name = name
        # For groupby perf. See note in indexes/base about _index_data
        result._index_data = dtarr._data
        result._reset_identity()
        return result

    # --------------------------------------------------------------------

    def __array__(self, dtype=None):
        if (dtype is None and isinstance(self._data, DatetimeArray)
                and getattr(self.dtype, 'tz', None)):
            msg = (
                "Converting timezone-aware DatetimeArray to timezone-naive "
                "ndarray with 'datetime64[ns]' dtype. In the future, this "
                "will return an ndarray with 'object' dtype where each "
                "element is a 'pandas.Timestamp' with the correct 'tz'.\n\t"
                "To accept the future behavior, pass 'dtype=object'.\n\t"
                "To keep the old behavior, pass 'dtype=\"datetime64[ns]\"'."
            )
            warnings.warn(msg, FutureWarning, stacklevel=3)
            dtype = 'M8[ns]'
        return np.asarray(self._data, dtype=dtype)

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def tz(self):
        # GH 18595
        return self._data.tz

    @tz.setter
    def tz(self, value):
        # GH 3746: Prevent localizing or converting the index by setting tz
        raise AttributeError("Cannot directly set timezone. Use tz_localize() "
                             "or tz_convert() as appropriate")

    tzinfo = tz

    @cache_readonly
    def _is_dates_only(self):
        """Return a boolean if we are only dates (and don't have a timezone)"""
        from pandas.io.formats.format import _is_dates_only
        return _is_dates_only(self.values) and self.tz is None

    def __reduce__(self):

        # we use a special reudce here because we need
        # to simply set the .tz (and not reinterpret it)

        d = dict(data=self._data)
        d.update(self._get_attributes_dict())
        return _new_DatetimeIndex, (self.__class__, d), None

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if isinstance(state, dict):
            super(DatetimeIndex, self).__setstate__(state)

        elif isinstance(state, tuple):

            # < 0.15 compat
            if len(state) == 2:
                nd_state, own_state = state
                data = np.empty(nd_state[1], dtype=nd_state[2])
                np.ndarray.__setstate__(data, nd_state)

                freq = own_state[1]
                tz = timezones.tz_standardize(own_state[2])
                dtype = tz_to_dtype(tz)
                dtarr = DatetimeArray._simple_new(data, freq=freq, dtype=dtype)

                self.name = own_state[0]

            else:  # pragma: no cover
                data = np.empty(state)
                np.ndarray.__setstate__(data, state)
                dtarr = DatetimeArray(data)

            self._data = dtarr
            self._reset_identity()

        else:
            raise Exception("invalid pickle state")
    _unpickle_compat = __setstate__

    def _convert_for_op(self, value):
        """ Convert value to be insertable to ndarray """
        if self._has_same_tz(value):
            return _to_M8(value)
        raise ValueError('Passed item and index have different timezone')

    def _maybe_update_attributes(self, attrs):
        """ Update Index attributes (e.g. freq) depending on op """
        freq = attrs.get('freq', None)
        if freq is not None:
            # no need to infer if freq is None
            attrs['freq'] = 'infer'
        return attrs

    # --------------------------------------------------------------------
    # Rendering Methods

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return libts.ints_to_pydatetime(self.asi8, self.tz)

    def _format_native_types(self, na_rep='NaT', date_format=None, **kwargs):
        from pandas.io.formats.format import _get_format_datetime64_from_values
        fmt = _get_format_datetime64_from_values(self, date_format)

        return libts.format_array_from_datetime(self.asi8,
                                                tz=self.tz,
                                                format=fmt,
                                                na_rep=na_rep)

    @property
    def _formatter_func(self):
        from pandas.io.formats.format import _get_format_datetime64
        formatter = _get_format_datetime64(is_dates_only=self._is_dates_only)
        return lambda x: "'%s'" % formatter(x, tz=self.tz)

    # --------------------------------------------------------------------
    # Set Operation Methods

    def union(self, other, sort=None):
        """
        Specialized union for DatetimeIndex objects. If combine
        overlapping ranges with the same DateOffset, will be much
        faster than Index.union

        Parameters
        ----------
        other : DatetimeIndex or array-like
        sort : bool or None, default None
            Whether to sort the resulting Index.

            * None : Sort the result, except when

              1. `self` and `other` are equal.
              2. `self` or `other` has length 0.
              3. Some values in `self` or `other` cannot be compared.
                 A RuntimeWarning is issued in this case.

            * False : do not sort the result

            .. versionadded:: 0.25.0

        Returns
        -------
        y : Index or DatetimeIndex
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)

        if len(other) == 0 or self.equals(other) or len(self) == 0:
            return super(DatetimeIndex, self).union(other, sort=sort)

        if not isinstance(other, DatetimeIndex):
            try:
                other = DatetimeIndex(other)
            except TypeError:
                pass

        this, other = self._maybe_utc_convert(other)

        if this._can_fast_union(other):
            return this._fast_union(other, sort=sort)
        else:
            result = Index.union(this, other, sort=sort)
            if isinstance(result, DatetimeIndex):
                # TODO: we shouldn't be setting attributes like this;
                #  in all the tests this equality already holds
                result._data._dtype = this.dtype
                if (result.freq is None and
                        (this.freq is not None or other.freq is not None)):
                    result.freq = to_offset(result.inferred_freq)
            return result

    def union_many(self, others):
        """
        A bit of a hack to accelerate unioning a collection of indexes
        """
        this = self

        for other in others:
            if not isinstance(this, DatetimeIndex):
                this = Index.union(this, other)
                continue

            if not isinstance(other, DatetimeIndex):
                try:
                    other = DatetimeIndex(other)
                except TypeError:
                    pass

            this, other = this._maybe_utc_convert(other)

            if this._can_fast_union(other):
                this = this._fast_union(other)
            else:
                dtype = this.dtype
                this = Index.union(this, other)
                if isinstance(this, DatetimeIndex):
                    # TODO: we shouldn't be setting attributes like this;
                    #  in all the tests this equality already holds
                    this._data._dtype = dtype
        return this

    def _can_fast_union(self, other):
        if not isinstance(other, DatetimeIndex):
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
        except (ValueError):

            # if we are comparing a freq that does not propagate timezones
            # this will raise
            return False

    def _fast_union(self, other, sort=None):
        if len(other) == 0:
            return self.view(type(self))

        if len(self) == 0:
            return other.view(type(self))

        # Both DTIs are monotonic. Check if they are already
        # in the "correct" order
        if self[0] <= other[0]:
            left, right = self, other
        # DTIs are not in the "correct" order and we don't want
        # to sort but want to remove overlaps
        elif sort is False:
            left, right = self, other
            left_start = left[0]
            loc = right.searchsorted(left_start, side='left')
            right_chunk = right.values[:loc]
            dates = _concat._concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        # DTIs are not in the "correct" order and we want
        # to sort
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # TODO: consider re-implementing freq._should_cache for fastpath

        # concatenate dates
        if left_end < right_end:
            loc = right.searchsorted(left_end, side='right')
            right_chunk = right.values[loc:]
            dates = _concat._concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        else:
            return left

    def _wrap_setop_result(self, other, result):
        name = get_op_result_name(self, other)
        return self._shallow_copy(result, name=name, freq=None, tz=self.tz)

    def intersection(self, other, sort=False):
        """
        Specialized intersection for DatetimeIndex objects. May be much faster
        than Index.intersection

        Parameters
        ----------
        other : DatetimeIndex or array-like
        sort : False or None, default False
            Sort the resulting index if possible.

            .. versionadded:: 0.24.0

            .. versionchanged:: 0.24.1

               Changed the default to ``False`` to match the behaviour
               from before 0.24.0.

        Returns
        -------
        y : Index or DatetimeIndex
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)

        if self.equals(other):
            return self._get_reconciled_name_object(other)

        if not isinstance(other, DatetimeIndex):
            try:
                other = DatetimeIndex(other)
            except (TypeError, ValueError):
                pass
            result = Index.intersection(self, other, sort=sort)
            if isinstance(result, DatetimeIndex):
                if result.freq is None:
                    result.freq = to_offset(result.inferred_freq)
            return result

        elif (other.freq is None or self.freq is None or
              other.freq != self.freq or
              not other.freq.isAnchored() or
              (not self.is_monotonic or not other.is_monotonic)):
            result = Index.intersection(self, other, sort=sort)
            # Invalidate the freq of `result`, which may not be correct at
            # this point, depending on the values.
            result.freq = None
            result = self._shallow_copy(result._values, name=result.name,
                                        tz=result.tz, freq=None)
            if result.freq is None:
                result.freq = to_offset(result.inferred_freq)
            return result

        if len(self) == 0:
            return self
        if len(other) == 0:
            return other
        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        end = min(left[-1], right[-1])
        start = right[0]

        if end < start:
            return type(self)(data=[])
        else:
            lslice = slice(*left.slice_locs(start, end))
            left_chunk = left.values[lslice]
            return self._shallow_copy(left_chunk)

    # --------------------------------------------------------------------

    def _get_time_micros(self):
        values = self.asi8
        if self.tz is not None and not timezones.is_utc(self.tz):
            values = self._data._local_timestamps()
        return fields.get_time_micros(values)

    def to_series(self, keep_tz=None, index=None, name=None):
        """
        Create a Series with both index and values equal to the index keys
        useful with map for returning an indexer based on an index

        Parameters
        ----------
        keep_tz : optional, defaults False
            Return the data keeping the timezone.

            If keep_tz is True:

              If the timezone is not set, the resulting
              Series will have a datetime64[ns] dtype.

              Otherwise the Series will have an datetime64[ns, tz] dtype; the
              tz will be preserved.

            If keep_tz is False:

              Series will have a datetime64[ns] dtype. TZ aware
              objects will have the tz removed.

            .. versionchanged:: 0.24
                The default value will change to True in a future release.
                You can set ``keep_tz=True`` to already obtain the future
                behaviour and silence the warning.

        index : Index, optional
            index of resulting Series. If None, defaults to original index
        name : string, optional
            name of resulting Series. If None, defaults to name of original
            index

        Returns
        -------
        Series
        """
        from pandas import Series

        if index is None:
            index = self._shallow_copy()
        if name is None:
            name = self.name

        if keep_tz is None and self.tz is not None:
            warnings.warn("The default of the 'keep_tz' keyword will change "
                          "to True in a future release. You can set "
                          "'keep_tz=True' to obtain the future behaviour and "
                          "silence this warning.", FutureWarning, stacklevel=2)
            keep_tz = False
        elif keep_tz is False:
            warnings.warn("Specifying 'keep_tz=False' is deprecated and this "
                          "option will be removed in a future release. If "
                          "you want to remove the timezone information, you "
                          "can do 'idx.tz_convert(None)' before calling "
                          "'to_series'.", FutureWarning, stacklevel=2)

        if keep_tz and self.tz is not None:
            # preserve the tz & copy
            values = self.copy(deep=True)
        else:
            values = self.values.copy()

        return Series(values, index=index, name=name)

    def snap(self, freq='S'):
        """
        Snap time stamps to nearest occurring frequency
        """
        # Superdumb, punting on any optimizing
        freq = to_offset(freq)

        snapped = np.empty(len(self), dtype=_NS_DTYPE)

        for i, v in enumerate(self):
            s = v
            if not freq.onOffset(s):
                t0 = freq.rollback(s)
                t1 = freq.rollforward(s)
                if abs(s - t0) < abs(t1 - s):
                    s = t0
                else:
                    s = t1
            snapped[i] = s

        # we know it conforms; skip check
        return DatetimeIndex._simple_new(snapped, name=self.name, tz=self.tz,
                                         freq=freq)

    def join(self, other, how='left', level=None, return_indexers=False,
             sort=False):
        """
        See Index.join
        """
        if (not isinstance(other, DatetimeIndex) and len(other) > 0 and
            other.inferred_type not in ('floating', 'integer', 'mixed-integer',
                                        'mixed-integer-float', 'mixed')):
            try:
                other = DatetimeIndex(other)
            except (TypeError, ValueError):
                pass

        this, other = self._maybe_utc_convert(other)
        return Index.join(this, other, how=how, level=level,
                          return_indexers=return_indexers, sort=sort)

    def _maybe_utc_convert(self, other):
        this = self
        if isinstance(other, DatetimeIndex):
            if self.tz is not None:
                if other.tz is None:
                    raise TypeError('Cannot join tz-naive with tz-aware '
                                    'DatetimeIndex')
            elif other.tz is not None:
                raise TypeError('Cannot join tz-naive with tz-aware '
                                'DatetimeIndex')

            if not timezones.tz_compare(self.tz, other.tz):
                this = self.tz_convert('UTC')
                other = other.tz_convert('UTC')
        return this, other

    def _wrap_joined_index(self, joined, other):
        name = get_op_result_name(self, other)
        if (isinstance(other, DatetimeIndex) and
                self.freq == other.freq and
                self._can_fast_union(other)):
            joined = self._shallow_copy(joined)
            joined.name = name
            return joined
        else:
            tz = getattr(other, 'tz', None)
            return self._simple_new(joined, name, tz=tz)

    def _parsed_string_to_bounds(self, reso, parsed):
        """
        Calculate datetime bounds for parsed time string and its resolution.

        Parameters
        ----------
        reso : Resolution
            Resolution provided by parsed string.
        parsed : datetime
            Datetime from parsed string.

        Returns
        -------
        lower, upper: pd.Timestamp

        """
        valid_resos = {'year', 'month', 'quarter', 'day', 'hour', 'minute',
                       'second', 'minute', 'second', 'microsecond'}
        if reso not in valid_resos:
            raise KeyError
        if reso == 'year':
            start = Timestamp(parsed.year, 1, 1)
            end = Timestamp(parsed.year, 12, 31, 23, 59, 59, 999999)
        elif reso == 'month':
            d = ccalendar.get_days_in_month(parsed.year, parsed.month)
            start = Timestamp(parsed.year, parsed.month, 1)
            end = Timestamp(parsed.year, parsed.month, d, 23, 59, 59, 999999)
        elif reso == 'quarter':
            qe = (((parsed.month - 1) + 2) % 12) + 1  # two months ahead
            d = ccalendar.get_days_in_month(parsed.year, qe)  # at end of month
            start = Timestamp(parsed.year, parsed.month, 1)
            end = Timestamp(parsed.year, qe, d, 23, 59, 59, 999999)
        elif reso == 'day':
            start = Timestamp(parsed.year, parsed.month, parsed.day)
            end = start + timedelta(days=1) - Nano(1)
        elif reso == 'hour':
            start = Timestamp(parsed.year, parsed.month, parsed.day,
                              parsed.hour)
            end = start + timedelta(hours=1) - Nano(1)
        elif reso == 'minute':
            start = Timestamp(parsed.year, parsed.month, parsed.day,
                              parsed.hour, parsed.minute)
            end = start + timedelta(minutes=1) - Nano(1)
        elif reso == 'second':
            start = Timestamp(parsed.year, parsed.month, parsed.day,
                              parsed.hour, parsed.minute, parsed.second)
            end = start + timedelta(seconds=1) - Nano(1)
        elif reso == 'microsecond':
            start = Timestamp(parsed.year, parsed.month, parsed.day,
                              parsed.hour, parsed.minute, parsed.second,
                              parsed.microsecond)
            end = start + timedelta(microseconds=1) - Nano(1)
        # GH 24076
        # If an incoming date string contained a UTC offset, need to localize
        # the parsed date to this offset first before aligning with the index's
        # timezone
        if parsed.tzinfo is not None:
            if self.tz is None:
                raise ValueError("The index must be timezone aware "
                                 "when indexing with a date string with a "
                                 "UTC offset")
            start = start.tz_localize(parsed.tzinfo).tz_convert(self.tz)
            end = end.tz_localize(parsed.tzinfo).tz_convert(self.tz)
        elif self.tz is not None:
            start = start.tz_localize(self.tz)
            end = end.tz_localize(self.tz)
        return start, end

    def _partial_date_slice(self, reso, parsed, use_lhs=True, use_rhs=True):
        is_monotonic = self.is_monotonic
        if (is_monotonic and reso in ['day', 'hour', 'minute', 'second'] and
                self._resolution >= Resolution.get_reso(reso)):
            # These resolution/monotonicity validations came from GH3931,
            # GH3452 and GH2369.

            # See also GH14826
            raise KeyError

        if reso == 'microsecond':
            # _partial_date_slice doesn't allow microsecond resolution, but
            # _parsed_string_to_bounds allows it.
            raise KeyError

        t1, t2 = self._parsed_string_to_bounds(reso, parsed)
        stamps = self.asi8

        if is_monotonic:

            # we are out of range
            if (len(stamps) and ((use_lhs and t1.value < stamps[0] and
                                  t2.value < stamps[0]) or
                                 ((use_rhs and t1.value > stamps[-1] and
                                   t2.value > stamps[-1])))):
                raise KeyError

            # a monotonic (sorted) series can be sliced
            left = stamps.searchsorted(
                t1.value, side='left') if use_lhs else None
            right = stamps.searchsorted(
                t2.value, side='right') if use_rhs else None

            return slice(left, right)

        lhs_mask = (stamps >= t1.value) if use_lhs else True
        rhs_mask = (stamps <= t2.value) if use_rhs else True

        # try to find a the dates
        return (lhs_mask & rhs_mask).nonzero()[0]

    def _maybe_promote(self, other):
        if other.inferred_type == 'date':
            other = DatetimeIndex(other)
        return self, other

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """

        if isinstance(key, datetime):

            # needed to localize naive datetimes
            if self.tz is not None:
                if key.tzinfo is not None:
                    key = Timestamp(key).tz_convert(self.tz)
                else:
                    key = Timestamp(key).tz_localize(self.tz)

            return self.get_value_maybe_box(series, key)

        if isinstance(key, time):
            locs = self.indexer_at_time(key)
            return series.take(locs)

        try:
            return com.maybe_box(self, Index.get_value(self, series, key),
                                 series, key)
        except KeyError:
            try:
                loc = self._get_string_slice(key)
                return series[loc]
            except (TypeError, ValueError, KeyError):
                pass

            try:
                return self.get_value_maybe_box(series, key)
            except (TypeError, ValueError, KeyError):
                raise KeyError(key)

    def get_value_maybe_box(self, series, key):
        # needed to localize naive datetimes
        if self.tz is not None:
            key = Timestamp(key)
            if key.tzinfo is not None:
                key = key.tz_convert(self.tz)
            else:
                key = key.tz_localize(self.tz)
        elif not isinstance(key, Timestamp):
            key = Timestamp(key)
        values = self._engine.get_value(com.values_from_object(series),
                                        key, tz=self.tz)
        return com.maybe_box(self, values, series, key)

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """

        if tolerance is not None:
            # try converting tolerance now, so errors don't get swallowed by
            # the try/except clauses below
            tolerance = self._convert_tolerance(tolerance, np.asarray(key))

        if isinstance(key, datetime):
            # needed to localize naive datetimes
            if key.tzinfo is None:
                key = Timestamp(key, tz=self.tz)
            else:
                key = Timestamp(key).tz_convert(self.tz)
            return Index.get_loc(self, key, method, tolerance)

        elif isinstance(key, timedelta):
            # GH#20464
            raise TypeError("Cannot index {cls} with {other}"
                            .format(cls=type(self).__name__,
                                    other=type(key).__name__))

        if isinstance(key, time):
            if method is not None:
                raise NotImplementedError('cannot yet lookup inexact labels '
                                          'when key is a time object')
            return self.indexer_at_time(key)

        try:
            return Index.get_loc(self, key, method, tolerance)
        except (KeyError, ValueError, TypeError):
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError, ValueError, OverflowError):
                pass

            try:
                stamp = Timestamp(key)
                if stamp.tzinfo is not None and self.tz is not None:
                    stamp = stamp.tz_convert(self.tz)
                else:
                    stamp = stamp.tz_localize(self.tz)
                return Index.get_loc(self, stamp, method, tolerance)
            except KeyError:
                raise KeyError(key)
            except ValueError as e:
                # list-like tolerance size must match target index size
                if 'list-like' in str(e):
                    raise e
                raise KeyError(key)

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        If label is a string, cast it to datetime according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'ix', 'loc', 'getitem'}

        Returns
        -------
        label :  object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """
        assert kind in ['ix', 'loc', 'getitem', None]

        if is_float(label) or isinstance(label, time) or is_integer(label):
            self._invalid_indexer('slice', label)

        if isinstance(label, str):
            freq = getattr(self, 'freqstr',
                           getattr(self, 'inferred_freq', None))
            _, parsed, reso = parsing.parse_time_string(label, freq)
            lower, upper = self._parsed_string_to_bounds(reso, parsed)
            # lower, upper form the half-open interval:
            #   [parsed, parsed + 1 freq)
            # because label may be passed to searchsorted
            # the bounds need swapped if index is reverse sorted and has a
            # length > 1 (is_monotonic_decreasing gives True for empty
            # and length 1 index)
            if self._is_strictly_monotonic_decreasing and len(self) > 1:
                return upper if side == 'left' else lower
            return lower if side == 'left' else upper
        else:
            return label

    def _get_string_slice(self, key, use_lhs=True, use_rhs=True):
        freq = getattr(self, 'freqstr',
                       getattr(self, 'inferred_freq', None))
        _, parsed, reso = parsing.parse_time_string(key, freq)
        loc = self._partial_date_slice(reso, parsed, use_lhs=use_lhs,
                                       use_rhs=use_rhs)
        return loc

    def slice_indexer(self, start=None, end=None, step=None, kind=None):
        """
        Return indexer for specified label slice.
        Index.slice_indexer, customized to handle time slicing.

        In addition to functionality provided by Index.slice_indexer, does the
        following:

        - if both `start` and `end` are instances of `datetime.time`, it
          invokes `indexer_between_time`
        - if `start` and `end` are both either string or None perform
          value-based selection in non-monotonic cases.

        """
        # For historical reasons DatetimeIndex supports slices between two
        # instances of datetime.time as if it were applying a slice mask to
        # an array of (self.hour, self.minute, self.seconds, self.microsecond).
        if isinstance(start, time) and isinstance(end, time):
            if step is not None and step != 1:
                raise ValueError('Must have step size of 1 with time slices')
            return self.indexer_between_time(start, end)

        if isinstance(start, time) or isinstance(end, time):
            raise KeyError('Cannot mix time and non-time slice keys')

        try:
            return Index.slice_indexer(self, start, end, step, kind=kind)
        except KeyError:
            # For historical reasons DatetimeIndex by default supports
            # value-based partial (aka string) slices on non-monotonic arrays,
            # let's try that.
            if ((start is None or isinstance(start, str)) and
                    (end is None or isinstance(end, str))):
                mask = True
                if start is not None:
                    start_casted = self._maybe_cast_slice_bound(
                        start, 'left', kind)
                    mask = start_casted <= self

                if end is not None:
                    end_casted = self._maybe_cast_slice_bound(
                        end, 'right', kind)
                    mask = (self <= end_casted) & mask

                indexer = mask.nonzero()[0][::step]
                if len(indexer) == len(self):
                    return slice(None)
                else:
                    return indexer
            else:
                raise

    # --------------------------------------------------------------------
    # Wrapping DatetimeArray

    # Compat for frequency inference, see GH#23789
    _is_monotonic_increasing = Index.is_monotonic_increasing
    _is_monotonic_decreasing = Index.is_monotonic_decreasing
    _is_unique = Index.is_unique

    _timezone = cache_readonly(DatetimeArray._timezone.fget)
    is_normalized = cache_readonly(DatetimeArray.is_normalized.fget)
    _resolution = cache_readonly(DatetimeArray._resolution.fget)

    strftime = ea_passthrough(DatetimeArray.strftime)
    _has_same_tz = ea_passthrough(DatetimeArray._has_same_tz)

    @property
    def offset(self):
        """
        get/set the frequency of the instance
        """
        msg = ('{cls}.offset has been deprecated and will be removed '
               'in a future version; use {cls}.freq instead.'
               .format(cls=type(self).__name__))
        warnings.warn(msg, FutureWarning, stacklevel=2)
        return self.freq

    @offset.setter
    def offset(self, value):
        """
        get/set the frequency of the instance
        """
        msg = ('{cls}.offset has been deprecated and will be removed '
               'in a future version; use {cls}.freq instead.'
               .format(cls=type(self).__name__))
        warnings.warn(msg, FutureWarning, stacklevel=2)
        self.freq = value

    def __getitem__(self, key):
        result = self._data.__getitem__(key)
        if is_scalar(result):
            return result
        elif result.ndim > 1:
            # To support MPL which performs slicing with 2 dim
            # even though it only has 1 dim by definition
            assert isinstance(result, np.ndarray), result
            return result
        return type(self)(result, name=self.name)

    @property
    def _box_func(self):
        return lambda x: Timestamp(x, tz=self.tz)

    # --------------------------------------------------------------------

    @Substitution(klass='DatetimeIndex')
    @Appender(_shared_docs['searchsorted'])
    def searchsorted(self, value, side='left', sorter=None):
        if isinstance(value, (np.ndarray, Index)):
            value = np.array(value, dtype=_NS_DTYPE, copy=False)
        else:
            value = _to_M8(value, tz=self.tz)

        return self.values.searchsorted(value, side=side)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type or typ == 'datetime'

    @property
    def inferred_type(self):
        # b/c datetime is represented as microseconds since the epoch, make
        # sure we can't have ambiguous indexing
        return 'datetime64'

    @property
    def is_all_dates(self):
        return True

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location

        Parameters
        ----------
        loc : int
        item : object
            if not either a Python datetime or a numpy integer-like, returned
            Index dtype will be object rather than datetime.

        Returns
        -------
        new_index : Index
        """
        if is_scalar(item) and isna(item):
            # GH 18295
            item = self._na_value

        freq = None

        if isinstance(item, (datetime, np.datetime64)):
            self._assert_can_do_op(item)
            if not self._has_same_tz(item) and not isna(item):
                raise ValueError(
                    'Passed item and index have different timezone')
            # check freq can be preserved on edge cases
            if self.size and self.freq is not None:
                if ((loc == 0 or loc == -len(self)) and
                        item + self.freq == self[0]):
                    freq = self.freq
                elif (loc == len(self)) and item - self.freq == self[-1]:
                    freq = self.freq
            item = _to_M8(item, tz=self.tz)

        try:
            new_dates = np.concatenate((self[:loc].asi8, [item.view(np.int64)],
                                        self[loc:].asi8))
            return self._shallow_copy(new_dates, freq=freq)
        except (AttributeError, TypeError):

            # fall back to object index
            if isinstance(item, str):
                return self.astype(object).insert(loc, item)
            raise TypeError(
                "cannot insert DatetimeIndex with incompatible label")

    def delete(self, loc):
        """
        Make a new DatetimeIndex with passed location(s) deleted.

        Parameters
        ----------
        loc: int, slice or array of ints
            Indicate which sub-arrays to remove.

        Returns
        -------
        new_index : DatetimeIndex
        """
        new_dates = np.delete(self.asi8, loc)

        freq = None
        if is_integer(loc):
            if loc in (0, -len(self), -1, len(self) - 1):
                freq = self.freq
        else:
            if is_list_like(loc):
                loc = lib.maybe_indices_to_slice(
                    ensure_int64(np.array(loc)), len(self))
            if isinstance(loc, slice) and loc.step in (1, None):
                if (loc.start in (0, None) or loc.stop in (len(self), None)):
                    freq = self.freq

        return self._shallow_copy(new_dates, freq=freq)

    def indexer_at_time(self, time, asof=False):
        """
        Return index locations of index values at particular time of day
        (e.g. 9:30AM).

        Parameters
        ----------
        time : datetime.time or string
            datetime.time or string in appropriate format ("%H:%M", "%H%M",
            "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
            "%I%M%S%p").

        Returns
        -------
        values_at_time : array of integers

        See Also
        --------
        indexer_between_time, DataFrame.at_time
        """
        if asof:
            raise NotImplementedError("'asof' argument is not supported")

        if isinstance(time, str):
            from dateutil.parser import parse
            time = parse(time).time()

        if time.tzinfo:
            if self.tz is None:
                raise ValueError("Index must be timezone aware.")
            time_micros = self.tz_convert(time.tzinfo)._get_time_micros()
        else:
            time_micros = self._get_time_micros()
        micros = _time_to_micros(time)
        return (micros == time_micros).nonzero()[0]

    def indexer_between_time(self, start_time, end_time, include_start=True,
                             include_end=True):
        """
        Return index locations of values between particular times of day
        (e.g., 9:00-9:30AM).

        Parameters
        ----------
        start_time, end_time : datetime.time, str
            datetime.time or string in appropriate format ("%H:%M", "%H%M",
            "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
            "%I%M%S%p").
        include_start : boolean, default True
        include_end : boolean, default True

        Returns
        -------
        values_between_time : array of integers

        See Also
        --------
        indexer_at_time, DataFrame.between_time
        """
        start_time = tools.to_time(start_time)
        end_time = tools.to_time(end_time)
        time_micros = self._get_time_micros()
        start_micros = _time_to_micros(start_time)
        end_micros = _time_to_micros(end_time)

        if include_start and include_end:
            lop = rop = operator.le
        elif include_start:
            lop = operator.le
            rop = operator.lt
        elif include_end:
            lop = operator.lt
            rop = operator.le
        else:
            lop = rop = operator.lt

        if start_time <= end_time:
            join_op = operator.and_
        else:
            join_op = operator.or_

        mask = join_op(lop(start_micros, time_micros),
                       rop(time_micros, end_micros))

        return mask.nonzero()[0]


DatetimeIndex._add_comparison_ops()
DatetimeIndex._add_numeric_methods_disabled()
DatetimeIndex._add_logical_methods_disabled()
DatetimeIndex._add_datetimelike_methods()


def date_range(start=None, end=None, periods=None, freq=None, tz=None,
               normalize=False, name=None, closed=None, **kwargs):
    """
    Return a fixed frequency DatetimeIndex.

    Parameters
    ----------
    start : str or datetime-like, optional
        Left bound for generating dates.
    end : str or datetime-like, optional
        Right bound for generating dates.
    periods : integer, optional
        Number of periods to generate.
    freq : str or DateOffset, default 'D'
        Frequency strings can have multiples, e.g. '5H'. See
        :ref:`here <timeseries.offset_aliases>` for a list of
        frequency aliases.
    tz : str or tzinfo, optional
        Time zone name for returning localized DatetimeIndex, for example
        'Asia/Hong_Kong'. By default, the resulting DatetimeIndex is
        timezone-naive.
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    name : str, default None
        Name of the resulting DatetimeIndex.
    closed : {None, 'left', 'right'}, optional
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None, the default).
    **kwargs
        For compatibility. Has no effect on the result.

    Returns
    -------
    rng : DatetimeIndex

    See Also
    --------
    DatetimeIndex : An immutable container for datetimes.
    timedelta_range : Return a fixed frequency TimedeltaIndex.
    period_range : Return a fixed frequency PeriodIndex.
    interval_range : Return a fixed frequency IntervalIndex.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``DatetimeIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end`` (closed on both sides).

    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Examples
    --------
    **Specifying the values**

    The next four examples generate the same `DatetimeIndex`, but vary
    the combination of `start`, `end` and `periods`.

    Specify `start` and `end`, with the default daily frequency.

    >>> pd.date_range(start='1/1/2018', end='1/08/2018')
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
                   '2018-01-05', '2018-01-06', '2018-01-07', '2018-01-08'],
                  dtype='datetime64[ns]', freq='D')

    Specify `start` and `periods`, the number of periods (days).

    >>> pd.date_range(start='1/1/2018', periods=8)
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
                   '2018-01-05', '2018-01-06', '2018-01-07', '2018-01-08'],
                  dtype='datetime64[ns]', freq='D')

    Specify `end` and `periods`, the number of periods (days).

    >>> pd.date_range(end='1/1/2018', periods=8)
    DatetimeIndex(['2017-12-25', '2017-12-26', '2017-12-27', '2017-12-28',
                   '2017-12-29', '2017-12-30', '2017-12-31', '2018-01-01'],
                  dtype='datetime64[ns]', freq='D')

    Specify `start`, `end`, and `periods`; the frequency is generated
    automatically (linearly spaced).

    >>> pd.date_range(start='2018-04-24', end='2018-04-27', periods=3)
    DatetimeIndex(['2018-04-24 00:00:00', '2018-04-25 12:00:00',
                   '2018-04-27 00:00:00'],
                  dtype='datetime64[ns]', freq=None)

    **Other Parameters**

    Changed the `freq` (frequency) to ``'M'`` (month end frequency).

    >>> pd.date_range(start='1/1/2018', periods=5, freq='M')
    DatetimeIndex(['2018-01-31', '2018-02-28', '2018-03-31', '2018-04-30',
                   '2018-05-31'],
                  dtype='datetime64[ns]', freq='M')

    Multiples are allowed

    >>> pd.date_range(start='1/1/2018', periods=5, freq='3M')
    DatetimeIndex(['2018-01-31', '2018-04-30', '2018-07-31', '2018-10-31',
                   '2019-01-31'],
                  dtype='datetime64[ns]', freq='3M')

    `freq` can also be specified as an Offset object.

    >>> pd.date_range(start='1/1/2018', periods=5, freq=pd.offsets.MonthEnd(3))
    DatetimeIndex(['2018-01-31', '2018-04-30', '2018-07-31', '2018-10-31',
                   '2019-01-31'],
                  dtype='datetime64[ns]', freq='3M')

    Specify `tz` to set the timezone.

    >>> pd.date_range(start='1/1/2018', periods=5, tz='Asia/Tokyo')
    DatetimeIndex(['2018-01-01 00:00:00+09:00', '2018-01-02 00:00:00+09:00',
                   '2018-01-03 00:00:00+09:00', '2018-01-04 00:00:00+09:00',
                   '2018-01-05 00:00:00+09:00'],
                  dtype='datetime64[ns, Asia/Tokyo]', freq='D')

    `closed` controls whether to include `start` and `end` that are on the
    boundary. The default includes boundary points on either end.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed=None)
    DatetimeIndex(['2017-01-01', '2017-01-02', '2017-01-03', '2017-01-04'],
                  dtype='datetime64[ns]', freq='D')

    Use ``closed='left'`` to exclude `end` if it falls on the boundary.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed='left')
    DatetimeIndex(['2017-01-01', '2017-01-02', '2017-01-03'],
                  dtype='datetime64[ns]', freq='D')

    Use ``closed='right'`` to exclude `start` if it falls on the boundary.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed='right')
    DatetimeIndex(['2017-01-02', '2017-01-03', '2017-01-04'],
                  dtype='datetime64[ns]', freq='D')
    """

    if freq is None and com._any_none(periods, start, end):
        freq = 'D'

    dtarr = DatetimeArray._generate_range(
        start=start, end=end, periods=periods,
        freq=freq, tz=tz, normalize=normalize,
        closed=closed, **kwargs)
    return DatetimeIndex._simple_new(
        dtarr, tz=dtarr.tz, freq=dtarr.freq, name=name)


def bdate_range(start=None, end=None, periods=None, freq='B', tz=None,
                normalize=True, name=None, weekmask=None, holidays=None,
                closed=None, **kwargs):
    """
    Return a fixed frequency DatetimeIndex, with business day as the default
    frequency

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates.
    end : string or datetime-like, default None
        Right bound for generating dates.
    periods : integer, default None
        Number of periods to generate.
    freq : string or DateOffset, default 'B' (business daily)
        Frequency strings can have multiples, e.g. '5H'.
    tz : string or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Beijing.
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    name : string, default None
        Name of the resulting DatetimeIndex.
    weekmask : string or None, default None
        Weekmask of valid business days, passed to ``numpy.busdaycalendar``,
        only used when custom frequency strings are passed.  The default
        value None is equivalent to 'Mon Tue Wed Thu Fri'.

        .. versionadded:: 0.21.0

    holidays : list-like or None, default None
        Dates to exclude from the set of valid business days, passed to
        ``numpy.busdaycalendar``, only used when custom frequency strings
        are passed.

        .. versionadded:: 0.21.0

    closed : string, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None).
    **kwargs
        For compatibility. Has no effect on the result.

    Returns
    -------
    DatetimeIndex

    Notes
    -----
    Of the four parameters: ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified.  Specifying ``freq`` is a requirement
    for ``bdate_range``.  Use ``date_range`` if specifying ``freq`` is not
    desired.

    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Examples
    --------
    Note how the two weekend days are skipped in the result.

    >>> pd.bdate_range(start='1/1/2018', end='1/08/2018')
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
               '2018-01-05', '2018-01-08'],
              dtype='datetime64[ns]', freq='B')
    """
    if freq is None:
        msg = 'freq must be specified for bdate_range; use date_range instead'
        raise TypeError(msg)

    if is_string_like(freq) and freq.startswith('C'):
        try:
            weekmask = weekmask or 'Mon Tue Wed Thu Fri'
            freq = prefix_mapping[freq](holidays=holidays, weekmask=weekmask)
        except (KeyError, TypeError):
            msg = 'invalid custom frequency string: {freq}'.format(freq=freq)
            raise ValueError(msg)
    elif holidays or weekmask:
        msg = ('a custom frequency string is required when holidays or '
               'weekmask are passed, got frequency {freq}').format(freq=freq)
        raise ValueError(msg)

    return date_range(start=start, end=end, periods=periods,
                      freq=freq, tz=tz, normalize=normalize, name=name,
                      closed=closed, **kwargs)


def cdate_range(start=None, end=None, periods=None, freq='C', tz=None,
                normalize=True, name=None, closed=None, **kwargs):
    """
    Return a fixed frequency DatetimeIndex, with CustomBusinessDay as the
    default frequency

    .. deprecated:: 0.21.0

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer, default None
        Number of periods to generate
    freq : string or DateOffset, default 'C' (CustomBusinessDay)
        Frequency strings can have multiples, e.g. '5H'
    tz : string, default None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Beijing
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    name : string, default None
        Name of the resulting DatetimeIndex
    weekmask : string, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    closed : string, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)

    Notes
    -----
    Of the three parameters: ``start``, ``end``, and ``periods``, exactly two
    must be specified.

    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Returns
    -------
    rng : DatetimeIndex
    """
    warnings.warn("cdate_range is deprecated and will be removed in a future "
                  "version, instead use pd.bdate_range(..., freq='{freq}')"
                  .format(freq=freq), FutureWarning, stacklevel=2)

    if freq == 'C':
        holidays = kwargs.pop('holidays', [])
        weekmask = kwargs.pop('weekmask', 'Mon Tue Wed Thu Fri')
        freq = CDay(holidays=holidays, weekmask=weekmask)

    return date_range(start=start, end=end, periods=periods, freq=freq,
                      tz=tz, normalize=normalize, name=name,
                      closed=closed, **kwargs)


def _time_to_micros(time):
    seconds = time.hour * 60 * 60 + 60 * time.minute + time.second
    return 1000000 * seconds + time.microsecond
