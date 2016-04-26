""" implement the TimedeltaIndex """

from datetime import timedelta
import numpy as np
from pandas.core.common import (ABCSeries, _TD_DTYPE, _INT64_DTYPE,
                                _maybe_box,
                                _values_from_object, isnull,
                                is_integer, is_float)
from pandas.core.index import Index, Int64Index
import pandas.compat as compat
from pandas.compat import u
from pandas.tseries.frequencies import to_offset
from pandas.core.base import _shared_docs
import pandas.core.common as com
import pandas.types.concat as _concat
from pandas.util.decorators import Appender, Substitution
from pandas.tseries.base import TimelikeOps, DatetimeIndexOpsMixin
from pandas.tseries.timedeltas import (to_timedelta,
                                       _coerce_scalar_to_timedelta_type)
from pandas.tseries.offsets import Tick, DateOffset

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos
import pandas.index as _index

Timedelta = tslib.Timedelta


def _td_index_cmp(opname, nat_result=False):
    """
    Wrap comparison operations to convert timedelta-like to timedelta64
    """

    def wrapper(self, other):
        func = getattr(super(TimedeltaIndex, self), opname)
        if _is_convertible_to_td(other):
            other = _to_m8(other)
            result = func(other)
            if com.isnull(other):
                result.fill(nat_result)
        else:
            if not com.is_list_like(other):
                raise TypeError("cannot compare a TimedeltaIndex with type "
                                "{0}".format(type(other)))

            other = TimedeltaIndex(other).values
            result = func(other)
            result = _values_from_object(result)

            if isinstance(other, Index):
                o_mask = other.values.view('i8') == tslib.iNaT
            else:
                o_mask = other.view('i8') == tslib.iNaT

            if o_mask.any():
                result[o_mask] = nat_result

        if self.hasnans:
            result[self._isnan] = nat_result

        # support of bool dtype indexers
        if com.is_bool_dtype(result):
            return result
        return Index(result)

    return wrapper


class TimedeltaIndex(DatetimeIndexOpsMixin, TimelikeOps, Int64Index):
    """
    Immutable ndarray of timedelta64 data, represented internally as int64, and
    which can be boxed to timedelta objects

    Parameters
    ----------
    data  : array-like (1-dimensional), optional
        Optional timedelta-like data to construct index with
    unit: unit of the arg (D,h,m,s,ms,us,ns) denote the unit, optional
        which is an integer/float number
    freq: a frequency for the index, optional
    copy  : bool
        Make a copy of input ndarray
    start : starting value, timedelta-like, optional
        If data is None, start is used as the start point in generating regular
        timedelta data.
    periods  : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument
    end   : end time, timedelta-like, optional
        If periods is none, generated index will extend to first conforming
        time on or just past end argument
    closed : string or None, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)
    name : object
        Name to be stored in the index
    """

    _typ = 'timedeltaindex'
    _join_precedence = 10

    def _join_i8_wrapper(joinf, **kwargs):
        return DatetimeIndexOpsMixin._join_i8_wrapper(
            joinf, dtype='m8[ns]', **kwargs)

    _inner_indexer = _join_i8_wrapper(_algos.inner_join_indexer_int64)
    _outer_indexer = _join_i8_wrapper(_algos.outer_join_indexer_int64)
    _left_indexer = _join_i8_wrapper(_algos.left_join_indexer_int64)
    _left_indexer_unique = _join_i8_wrapper(
        _algos.left_join_indexer_unique_int64, with_indexers=False)
    _arrmap = None
    _datetimelike_ops = ['days', 'seconds', 'microseconds', 'nanoseconds',
                         'freq', 'components']

    __eq__ = _td_index_cmp('__eq__')
    __ne__ = _td_index_cmp('__ne__', nat_result=True)
    __lt__ = _td_index_cmp('__lt__')
    __gt__ = _td_index_cmp('__gt__')
    __le__ = _td_index_cmp('__le__')
    __ge__ = _td_index_cmp('__ge__')

    _engine_type = _index.TimedeltaEngine

    _comparables = ['name', 'freq']
    _attributes = ['name', 'freq']
    _is_numeric_dtype = True
    _infer_as_myclass = True

    freq = None

    def __new__(cls, data=None, unit=None,
                freq=None, start=None, end=None, periods=None,
                copy=False, name=None,
                closed=None, verify_integrity=True, **kwargs):

        if isinstance(data, TimedeltaIndex) and freq is None and name is None:
            if copy:
                data = data.copy()
            return data

        freq_infer = False
        if not isinstance(freq, DateOffset):

            # if a passed freq is None, don't infer automatically
            if freq != 'infer':
                freq = to_offset(freq)
            else:
                freq_infer = True
                freq = None

        if periods is not None:
            if is_float(periods):
                periods = int(periods)
            elif not is_integer(periods):
                raise ValueError('Periods must be a number, got %s' %
                                 str(periods))

        if data is None and freq is None:
            raise ValueError("Must provide freq argument if no data is "
                             "supplied")

        if data is None:
            return cls._generate(start, end, periods, name, freq,
                                 closed=closed)

        if unit is not None:
            data = to_timedelta(data, unit=unit, box=False)

        if not isinstance(data, (np.ndarray, Index, ABCSeries)):
            if lib.isscalar(data):
                raise ValueError('TimedeltaIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

        # convert if not already
        if getattr(data, 'dtype', None) != _TD_DTYPE:
            data = to_timedelta(data, unit=unit, box=False)
        elif copy:
            data = np.array(data, copy=True)

        # check that we are matching freqs
        if verify_integrity and len(data) > 0:
            if freq is not None and not freq_infer:
                index = cls._simple_new(data, name=name)
                inferred = index.inferred_freq
                if inferred != freq.freqstr:
                    on_freq = cls._generate(
                        index[0], None, len(index), name, freq)
                    if not np.array_equal(index.asi8, on_freq.asi8):
                        raise ValueError('Inferred frequency {0} from passed '
                                         'timedeltas does not conform to '
                                         'passed frequency {1}'
                                         .format(inferred, freq.freqstr))
                index.freq = freq
                return index

        if freq_infer:
            index = cls._simple_new(data, name=name)
            inferred = index.inferred_freq
            if inferred:
                index.freq = to_offset(inferred)
            return index

        return cls._simple_new(data, name=name, freq=freq)

    @classmethod
    def _generate(cls, start, end, periods, name, offset, closed=None):
        if com._count_not_none(start, end, periods) != 2:
            raise ValueError('Must specify two of start, end, or periods')

        if start is not None:
            start = Timedelta(start)

        if end is not None:
            end = Timedelta(end)

        left_closed = False
        right_closed = False

        if start is None and end is None:
            if closed is not None:
                raise ValueError("Closed has to be None if not both of start"
                                 "and end are defined")

        if closed is None:
            left_closed = True
            right_closed = True
        elif closed == "left":
            left_closed = True
        elif closed == "right":
            right_closed = True
        else:
            raise ValueError("Closed has to be either 'left', 'right' or None")

        index = _generate_regular_range(start, end, periods, offset)
        index = cls._simple_new(index, name=name, freq=offset)

        if not left_closed:
            index = index[1:]
        if not right_closed:
            index = index[:-1]

        return index

    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit='ns')

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, **kwargs):
        if not getattr(values, 'dtype', None):
            values = np.array(values, copy=False)
        if values.dtype == np.object_:
            values = tslib.array_to_timedelta64(values)
        if values.dtype != _TD_DTYPE:
            values = com._ensure_int64(values).view(_TD_DTYPE)

        result = object.__new__(cls)
        result._data = values
        result.name = name
        result.freq = freq
        result._reset_identity()
        return result

    _na_value = tslib.NaT
    """The expected NA value to use with this index."""

    @property
    def _formatter_func(self):
        from pandas.formats.format import _get_format_timedelta64
        return _get_format_timedelta64(self, box=True)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if isinstance(state, dict):
            super(TimedeltaIndex, self).__setstate__(state)
        else:
            raise Exception("invalid pickle state")
    _unpickle_compat = __setstate__

    def _maybe_update_attributes(self, attrs):
        """ Update Index attributes (e.g. freq) depending on op """
        freq = attrs.get('freq', None)
        if freq is not None:
            # no need to infer if freq is None
            attrs['freq'] = 'infer'
        return attrs

    def _add_delta(self, delta):
        if isinstance(delta, (Tick, timedelta, np.timedelta64)):
            new_values = self._add_delta_td(delta)
            name = self.name
        elif isinstance(delta, TimedeltaIndex):
            new_values = self._add_delta_tdi(delta)
            # update name when delta is index
            name = com._maybe_match_name(self, delta)
        else:
            raise ValueError("cannot add the type {0} to a TimedeltaIndex"
                             .format(type(delta)))

        result = TimedeltaIndex(new_values, freq='infer', name=name)
        return result

    def _evaluate_with_timedelta_like(self, other, op, opstr):

        # allow division by a timedelta
        if opstr in ['__div__', '__truediv__']:
            if _is_convertible_to_td(other):
                other = Timedelta(other)
                if isnull(other):
                    raise NotImplementedError(
                        "division by pd.NaT not implemented")

                i8 = self.asi8
                result = i8 / float(other.value)
                result = self._maybe_mask_results(result, convert='float64')
                return Index(result, name=self.name, copy=False)

        return NotImplemented

    def _add_datelike(self, other):
        # adding a timedeltaindex to a datetimelike
        from pandas import Timestamp, DatetimeIndex
        if other is tslib.NaT:
            result = self._nat_new(box=False)
        else:
            other = Timestamp(other)
            i8 = self.asi8
            result = i8 + other.value
            result = self._maybe_mask_results(result, fill_value=tslib.iNaT)
        return DatetimeIndex(result, name=self.name, copy=False)

    def _sub_datelike(self, other):
        from pandas import DatetimeIndex
        if other is tslib.NaT:
            result = self._nat_new(box=False)
        else:
            raise TypeError("cannot subtract a datelike from a TimedeltaIndex")
        return DatetimeIndex(result, name=self.name, copy=False)

    def _format_native_types(self, na_rep=u('NaT'),
                             date_format=None, **kwargs):
        from pandas.formats.format import Timedelta64Formatter
        return Timedelta64Formatter(values=self,
                                    nat_rep=na_rep,
                                    justify='all').get_result()

    def _get_field(self, m):

        values = self.asi8
        hasnans = self.hasnans
        if hasnans:
            result = np.empty(len(self), dtype='float64')
            mask = self._isnan
            imask = ~mask
            result.flat[imask] = np.array(
                [getattr(Timedelta(val), m) for val in values[imask]])
            result[mask] = np.nan
        else:
            result = np.array([getattr(Timedelta(val), m)
                               for val in values], dtype='int64')
        return result

    @property
    def days(self):
        """ Number of days for each element. """
        return self._get_field('days')

    @property
    def seconds(self):
        """ Number of seconds (>= 0 and less than 1 day) for each element. """
        return self._get_field('seconds')

    @property
    def microseconds(self):
        """
        Number of microseconds (>= 0 and less than 1 second) for each
        element. """
        return self._get_field('microseconds')

    @property
    def nanoseconds(self):
        """
        Number of nanoseconds (>= 0 and less than 1 microsecond) for each
        element.
        """
        return self._get_field('nanoseconds')

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

        columns = ['days', 'hours', 'minutes', 'seconds',
                   'milliseconds', 'microseconds', 'nanoseconds']
        hasnans = self.hasnans
        if hasnans:
            def f(x):
                if isnull(x):
                    return [np.nan] * len(columns)
                return x.components
        else:
            def f(x):
                return x.components

        result = DataFrame([f(x) for x in self])
        result.columns = columns
        if not hasnans:
            result = result.astype('int64')
        return result

    def total_seconds(self):
        """
        Total duration of each element expressed in seconds.

        .. versionadded:: 0.17.0
        """
        return self._maybe_mask_results(1e-9 * self.asi8)

    def to_pytimedelta(self):
        """
        Return TimedeltaIndex as object ndarray of datetime.timedelta objects

        Returns
        -------
        datetimes : ndarray
        """
        return tslib.ints_to_pytimedelta(self.asi8)

    def astype(self, dtype):
        dtype = np.dtype(dtype)

        if dtype == np.object_:
            return self.asobject
        elif dtype == _INT64_DTYPE:
            return self.asi8.copy()
        elif dtype == _TD_DTYPE:
            return self
        elif dtype.kind == 'm':

            # return an index (essentially this is division)
            result = self.values.astype(dtype)
            if self.hasnans:
                return Index(self._maybe_mask_results(result,
                                                      convert='float64'),
                             name=self.name)

            return Index(result.astype('i8'), name=self.name)

        else:  # pragma: no cover
            raise ValueError('Cannot cast TimedeltaIndex to dtype %s' % dtype)

    def union(self, other):
        """
        Specialized union for TimedeltaIndex objects. If combine
        overlapping ranges with the same DateOffset, will be much
        faster than Index.union

        Parameters
        ----------
        other : TimedeltaIndex or array-like

        Returns
        -------
        y : Index or TimedeltaIndex
        """
        self._assert_can_do_setop(other)
        if not isinstance(other, TimedeltaIndex):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass
        this, other = self, other

        if this._can_fast_union(other):
            return this._fast_union(other)
        else:
            result = Index.union(this, other)
            if isinstance(result, TimedeltaIndex):
                if result.freq is None:
                    result.freq = to_offset(result.inferred_freq)
            return result

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        name = self.name
        to_concat = [self]

        if isinstance(other, (list, tuple)):
            to_concat = to_concat + list(other)
        else:
            to_concat.append(other)

        for obj in to_concat:
            if isinstance(obj, Index) and obj.name != name:
                name = None
                break

        to_concat = self._ensure_compat_concat(to_concat)
        return Index(_concat._concat_compat(to_concat), name=name)

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        if _is_convertible_to_index(other):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass

        return Index.join(self, other, how=how, level=level,
                          return_indexers=return_indexers)

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        if (isinstance(other, TimedeltaIndex) and self.freq == other.freq and
                self._can_fast_union(other)):
            joined = self._shallow_copy(joined, name=name)
            return joined
        else:
            return self._simple_new(joined, name)

    def _can_fast_union(self, other):
        if not isinstance(other, TimedeltaIndex):
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
        return (right_start == left_end + freq) or right_start in left

    def _fast_union(self, other):
        if len(other) == 0:
            return self.view(type(self))

        if len(self) == 0:
            return other.view(type(self))

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # concatenate
        if left_end < right_end:
            loc = right.searchsorted(left_end, side='right')
            right_chunk = right.values[loc:]
            dates = _concat._concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        else:
            return left

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        return self._simple_new(result, name=name, freq=None)

    def intersection(self, other):
        """
        Specialized intersection for TimedeltaIndex objects. May be much faster
        than Index.intersection

        Parameters
        ----------
        other : TimedeltaIndex or array-like

        Returns
        -------
        y : Index or TimedeltaIndex
        """
        self._assert_can_do_setop(other)
        if not isinstance(other, TimedeltaIndex):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass
            result = Index.intersection(self, other)
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

    def _possibly_promote(self, other):
        if other.inferred_type == 'timedelta':
            other = TimedeltaIndex(other)
        return self, other

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """

        if _is_convertible_to_td(key):
            key = Timedelta(key)
            return self.get_value_maybe_box(series, key)

        try:
            return _maybe_box(self, Index.get_value(self, series, key),
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
        if not isinstance(key, Timedelta):
            key = Timedelta(key)
        values = self._engine.get_value(_values_from_object(series), key)
        return _maybe_box(self, values, series, key)

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
            tolerance = self._convert_tolerance(tolerance)

        if _is_convertible_to_td(key):
            key = Timedelta(key)
            return Index.get_loc(self, key, method, tolerance)

        try:
            return Index.get_loc(self, key, method, tolerance)
        except (KeyError, ValueError, TypeError):
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError, ValueError):
                pass

            try:
                stamp = Timedelta(key)
                return Index.get_loc(self, stamp, method, tolerance)
            except (KeyError, ValueError):
                raise KeyError(key)

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        If label is a string, cast it to timedelta according to resolution.


        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'ix', 'loc', 'getitem'}

        Returns
        -------
        label :  object

        """
        assert kind in ['ix', 'loc', 'getitem', None]

        if isinstance(label, compat.string_types):
            parsed = _coerce_scalar_to_timedelta_type(label, box=True)
            lbound = parsed.round(parsed.resolution)
            if side == 'left':
                return lbound
            else:
                return (lbound + to_offset(parsed.resolution) -
                        Timedelta(1, 'ns'))
        elif is_integer(label) or is_float(label):
            self._invalid_indexer('slice', label)

        return label

    def _get_string_slice(self, key, use_lhs=True, use_rhs=True):
        freq = getattr(self, 'freqstr',
                       getattr(self, 'inferred_freq', None))
        if is_integer(key) or is_float(key):
            self._invalid_indexer('slice', key)
        loc = self._partial_td_slice(key, freq, use_lhs=use_lhs,
                                     use_rhs=use_rhs)
        return loc

    def _partial_td_slice(self, key, freq, use_lhs=True, use_rhs=True):

        # given a key, try to figure out a location for a partial slice
        if not isinstance(key, compat.string_types):
            return key

        raise NotImplementedError

        # TODO(wesm): dead code
        # parsed = _coerce_scalar_to_timedelta_type(key, box=True)

        # is_monotonic = self.is_monotonic

        # # figure out the resolution of the passed td
        # # and round to it

        # # t1 = parsed.round(reso)

        # t2 = t1 + to_offset(parsed.resolution) - Timedelta(1, 'ns')

        # stamps = self.asi8

        # if is_monotonic:

        #     # we are out of range
        #     if (len(stamps) and ((use_lhs and t1.value < stamps[0] and
        #                           t2.value < stamps[0]) or
        #                          ((use_rhs and t1.value > stamps[-1] and
        #                            t2.value > stamps[-1])))):
        #         raise KeyError

        #     # a monotonic (sorted) series can be sliced
        #     left = (stamps.searchsorted(t1.value, side='left')
        #             if use_lhs else None)
        #     right = (stamps.searchsorted(t2.value, side='right')
        #              if use_rhs else None)

        #     return slice(left, right)

        # lhs_mask = (stamps >= t1.value) if use_lhs else True
        # rhs_mask = (stamps <= t2.value) if use_rhs else True

        # # try to find a the dates
        # return (lhs_mask & rhs_mask).nonzero()[0]

    @Substitution(klass='TimedeltaIndex', value='key')
    @Appender(_shared_docs['searchsorted'])
    def searchsorted(self, key, side='left', sorter=None):
        if isinstance(key, (np.ndarray, Index)):
            key = np.array(key, dtype=_TD_DTYPE, copy=False)
        else:
            key = _to_m8(key)

        return self.values.searchsorted(key, side=side, sorter=sorter)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type or typ == 'timedelta'

    @property
    def inferred_type(self):
        return 'timedelta64'

    @property
    def dtype(self):
        return _TD_DTYPE

    @property
    def is_all_dates(self):
        return True

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if (not hasattr(other, 'inferred_type') or
                other.inferred_type != 'timedelta64'):
            try:
                other = TimedeltaIndex(other)
            except:
                return False

        return np.array_equal(self.asi8, other.asi8)

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

        # try to convert if possible
        if _is_convertible_to_td(item):
            try:
                item = Timedelta(item)
            except:
                pass

        freq = None
        if isinstance(item, (Timedelta, tslib.NaTType)):

            # check freq can be preserved on edge cases
            if self.freq is not None:
                if ((loc == 0 or loc == -len(self)) and
                        item + self.freq == self[0]):
                    freq = self.freq
                elif (loc == len(self)) and item - self.freq == self[-1]:
                    freq = self.freq
            item = _to_m8(item)

        try:
            new_tds = np.concatenate((self[:loc].asi8, [item.view(np.int64)],
                                      self[loc:].asi8))
            return TimedeltaIndex(new_tds, name=self.name, freq=freq)

        except (AttributeError, TypeError):

            # fall back to object index
            if isinstance(item, compat.string_types):
                return self.asobject.insert(loc, item)
            raise TypeError(
                "cannot insert TimedeltaIndex with incompatible label")

    def delete(self, loc):
        """
        Make a new DatetimeIndex with passed location(s) deleted.

        Parameters
        ----------
        loc: int, slice or array of ints
            Indicate which sub-arrays to remove.

        Returns
        -------
        new_index : TimedeltaIndex
        """
        new_tds = np.delete(self.asi8, loc)

        freq = 'infer'
        if is_integer(loc):
            if loc in (0, -len(self), -1, len(self) - 1):
                freq = self.freq
        else:
            if com.is_list_like(loc):
                loc = lib.maybe_indices_to_slice(
                    com._ensure_int64(np.array(loc)), len(self))
            if isinstance(loc, slice) and loc.step in (1, None):
                if (loc.start in (0, None) or loc.stop in (len(self), None)):
                    freq = self.freq

        return TimedeltaIndex(new_tds, name=self.name, freq=freq)


TimedeltaIndex._add_numeric_methods()
TimedeltaIndex._add_logical_methods_disabled()
TimedeltaIndex._add_datetimelike_methods()


def _is_convertible_to_index(other):
    """
    return a boolean whether I can attempt conversion to a TimedeltaIndex
    """
    if isinstance(other, TimedeltaIndex):
        return True
    elif (len(other) > 0 and
          other.inferred_type not in ('floating', 'mixed-integer', 'integer',
                                      'mixed-integer-float', 'mixed')):
        return True
    return False


def _is_convertible_to_td(key):
    return isinstance(key, (DateOffset, timedelta, Timedelta,
                            np.timedelta64, compat.string_types))


def _to_m8(key):
    """
    Timedelta-like => dt64
    """
    if not isinstance(key, Timedelta):
        # this also converts strings
        key = Timedelta(key)

    # return an type that can be compared
    return np.int64(key.value).view(_TD_DTYPE)


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
        raise ValueError("at least 'start' or 'end' should be specified "
                         "if a 'period' is given.")

    data = np.arange(b, e, stride, dtype=np.int64)
    data = TimedeltaIndex._simple_new(data, None)

    return data


def timedelta_range(start=None, end=None, periods=None, freq='D',
                    name=None, closed=None):
    """
    Return a fixed frequency timedelta index, with day as the default
    frequency

    Parameters
    ----------
    start : string or timedelta-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer or None, default None
        If None, must specify start and end
    freq : string or DateOffset, default 'D' (calendar daily)
        Frequency strings can have multiples, e.g. '5H'
    name : str, default None
        Name of the resulting index
    closed : string or None, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)

    Notes
    -----
    2 of start, end, or periods must be specified

    Returns
    -------
    rng : TimedeltaIndex
    """
    return TimedeltaIndex(start=start, end=end, periods=periods,
                          freq=freq, name=name,
                          closed=closed)
