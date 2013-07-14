# pylint: disable=E1101
import operator

from datetime import time, datetime
from datetime import timedelta

import numpy as np

from pandas.core.common import isnull, _NS_DTYPE, _INT64_DTYPE
from pandas.core.index import Index, Int64Index
from pandas.tseries.frequencies import (
    infer_freq, to_offset, get_period_alias,
    Resolution, get_reso_string)
from pandas.tseries.offsets import DateOffset, generate_range, Tick, CDay
from pandas.tseries.tools import parse_time_string, normalize_date
from pandas.util.decorators import cache_readonly
import pandas.core.common as com
import pandas.tseries.offsets as offsets
import pandas.tseries.tools as tools

from pandas.lib import Timestamp
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos
import pandas.index as _index


def _utc():
    import pytz
    return pytz.utc

# -------- some conversion wrapper functions


def _field_accessor(name, field):
    def f(self):
        values = self.asi8
        if self.tz is not None:
            utc = _utc()
            if self.tz is not utc:
                values = self._local_timestamps()
        return tslib.get_date_field(values, field)
    f.__name__ = name
    return property(f)


def _join_i8_wrapper(joinf, with_indexers=True):
    @staticmethod
    def wrapper(left, right):
        if isinstance(left, np.ndarray):
            left = left.view('i8', type=np.ndarray)
        if isinstance(right, np.ndarray):
            right = right.view('i8', type=np.ndarray)
        results = joinf(left, right)
        if with_indexers:
            join_index, left_indexer, right_indexer = results
            join_index = join_index.view('M8[ns]')
            return join_index, left_indexer, right_indexer
        return results
    return wrapper


def _dt_index_cmp(opname):
    """
    Wrap comparison operations to convert datetime-like to datetime64
    """
    def wrapper(self, other):
        func = getattr(super(DatetimeIndex, self), opname)
        if isinstance(other, datetime):
            other = _to_m8(other, tz=self.tz)
        elif isinstance(other, list):
            other = DatetimeIndex(other)
        elif isinstance(other, basestring):
            other = _to_m8(other, tz=self.tz)
        elif not isinstance(other, np.ndarray):
            other = _ensure_datetime64(other)
        result = func(other)

        return result.view(np.ndarray)

    return wrapper


def _ensure_datetime64(other):
    if isinstance(other, np.datetime64):
        return other
    raise TypeError('%s type object %s' % (type(other), str(other)))


_midnight = time(0, 0)

class DatetimeIndex(Int64Index):
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
        One of pandas date offset strings or corresponding objects
    start : starting value, datetime-like, optional
        If data is None, start is used as the start point in generating regular
        timestamp data.
    periods  : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument
    end   : end time, datetime-like, optional
        If periods is none, generated index will extend to first conforming
        time on or just past end argument
    """
    _join_precedence = 10

    _inner_indexer = _join_i8_wrapper(_algos.inner_join_indexer_int64)
    _outer_indexer = _join_i8_wrapper(_algos.outer_join_indexer_int64)
    _left_indexer = _join_i8_wrapper(_algos.left_join_indexer_int64)
    _left_indexer_unique = _join_i8_wrapper(
        _algos.left_join_indexer_unique_int64, with_indexers=False)
    _arrmap = None

    __eq__ = _dt_index_cmp('__eq__')
    __ne__ = _dt_index_cmp('__ne__')
    __lt__ = _dt_index_cmp('__lt__')
    __gt__ = _dt_index_cmp('__gt__')
    __le__ = _dt_index_cmp('__le__')
    __ge__ = _dt_index_cmp('__ge__')

    # structured array cache for datetime fields
    _sarr_cache = None

    _engine_type = _index.DatetimeEngine

    offset = None

    def __new__(cls, data=None,
                freq=None, start=None, end=None, periods=None,
                copy=False, name=None, tz=None,
                verify_integrity=True, normalize=False, **kwds):

        dayfirst = kwds.pop('dayfirst', None)
        yearfirst = kwds.pop('yearfirst', None)
        warn = False
        if 'offset' in kwds and kwds['offset']:
            freq = kwds['offset']
            warn = True

        freq_infer = False
        if not isinstance(freq, DateOffset):
            if freq != 'infer':
                freq = to_offset(freq)
            else:
                freq_infer = True
                freq = None

        if warn:
            import warnings
            warnings.warn("parameter 'offset' is deprecated, "
                          "please use 'freq' instead",
                          FutureWarning)

        offset = freq

        if periods is not None:
            if com.is_float(periods):
                periods = int(periods)
            elif not com.is_integer(periods):
                raise ValueError('Periods must be a number, got %s' %
                                 str(periods))

        if data is None and offset is None:
            raise ValueError("Must provide freq argument if no data is "
                             "supplied")

        if data is None:
            return cls._generate(start, end, periods, name, offset,
                                 tz=tz, normalize=normalize)

        if not isinstance(data, np.ndarray):
            if np.isscalar(data):
                raise ValueError('DatetimeIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            data = np.asarray(data, dtype='O')

            # try a few ways to make it datetime64
            if lib.is_string_array(data):
                data = _str_to_dt_array(data, offset, dayfirst=dayfirst,
                                        yearfirst=yearfirst)
            else:
                data = tools.to_datetime(data)
                data.offset = offset
                if isinstance(data, DatetimeIndex):
                    if name is not None:
                        data.name = name

                    if tz is not None:
                        return data.tz_localize(tz)

                    return data

        if issubclass(data.dtype.type, basestring):
            data = _str_to_dt_array(data, offset, dayfirst=dayfirst,
                                      yearfirst=yearfirst)

        if issubclass(data.dtype.type, np.datetime64):
            if isinstance(data, DatetimeIndex):
                if tz is None:
                    tz = data.tz

                subarr = data.values

                if offset is None:
                    offset = data.offset
                    verify_integrity = False
            else:
                if data.dtype != _NS_DTYPE:
                    subarr = tslib.cast_to_nanoseconds(data)
                else:
                    subarr = data
        elif data.dtype == _INT64_DTYPE:
            if isinstance(data, Int64Index):
                raise TypeError('cannot convert Int64Index->DatetimeIndex')
            if copy:
                subarr = np.asarray(data, dtype=_NS_DTYPE)
            else:
                subarr = data.view(_NS_DTYPE)
        else:
            try:
                subarr = tools.to_datetime(data)
            except ValueError:
                # tz aware
                subarr = tools.to_datetime(data, utc=True)

            if not np.issubdtype(subarr.dtype, np.datetime64):
                raise TypeError('Unable to convert %s to datetime dtype'
                                % str(data))

        if isinstance(subarr, DatetimeIndex):
            if tz is None:
                tz = subarr.tz
        else:
            if tz is not None:
                tz = tools._maybe_get_tz(tz)

                if (not isinstance(data, DatetimeIndex) or
                        getattr(data, 'tz', None) is None):
                    # Convert tz-naive to UTC
                    ints = subarr.view('i8')
                    subarr = tslib.tz_localize_to_utc(ints, tz)

                subarr = subarr.view(_NS_DTYPE)

        subarr = subarr.view(cls)
        subarr.name = name
        subarr.offset = offset
        subarr.tz = tz

        if verify_integrity and len(subarr) > 0:
            if offset is not None and not freq_infer:
                inferred = subarr.inferred_freq
                if inferred != offset.freqstr:
                    raise ValueError('Dates do not conform to passed '
                                     'frequency')

        if freq_infer:
            inferred = subarr.inferred_freq
            if inferred:
                subarr.offset = to_offset(inferred)

        return subarr

    @classmethod
    def _generate(cls, start, end, periods, name, offset,
                  tz=None, normalize=False):
        if com._count_not_none(start, end, periods) != 2:
            raise ValueError('Must specify two of start, end, or periods')

        _normalized = True

        if start is not None:
            start = Timestamp(start)

        if end is not None:
            end = Timestamp(end)

        try:
            inferred_tz = tools._infer_tzinfo(start, end)
        except:
            raise ValueError('Start and end cannot both be tz-aware with '
                             'different timezones')

        inferred_tz = tools._maybe_get_tz(inferred_tz)
        tz = tools._maybe_get_tz(tz)

        if tz is not None and inferred_tz is not None:
            if not inferred_tz == tz:
                raise AssertionError()

        elif inferred_tz is not None:
            tz = inferred_tz


        if start is not None:
            if normalize:
                start = normalize_date(start)
                _normalized = True
            else:
                _normalized = _normalized and start.time() == _midnight

        if end is not None:
            if normalize:
                end = normalize_date(end)
                _normalized = True
            else:
                _normalized = _normalized and end.time() == _midnight

        if hasattr(offset, 'delta') and offset != offsets.Day():
            if inferred_tz is None and tz is not None:
                # naive dates
                if start is not None and start.tz is None:
                    start = start.tz_localize(tz)

                if end is not None and end.tz is None:
                    end = end.tz_localize(tz)

            if start and end:
                if start.tz is None and end.tz is not None:
                    start = start.tz_localize(end.tz)

                if end.tz is None and start.tz is not None:
                    end = end.tz_localize(start.tz)

            if (offset._should_cache() and
                not (offset._normalize_cache and not _normalized) and
                    _naive_in_cache_range(start, end)):
                index = cls._cached_range(start, end, periods=periods,
                                          offset=offset, name=name)
            else:
                index = _generate_regular_range(start, end, periods, offset)

        else:

            if inferred_tz is None and tz is not None:
                # naive dates
                if start is not None and start.tz is not None:
                    start = start.replace(tzinfo=None)

                if end is not None and end.tz is not None:
                    end = end.replace(tzinfo=None)

            if start and end:
                if start.tz is None and end.tz is not None:
                    end = end.replace(tzinfo=None)

                if end.tz is None and start.tz is not None:
                    start = start.replace(tzinfo=None)

            if (offset._should_cache() and
                not (offset._normalize_cache and not _normalized) and
                    _naive_in_cache_range(start, end)):
                index = cls._cached_range(start, end, periods=periods,
                                          offset=offset, name=name)
            else:
                index = _generate_regular_range(start, end, periods, offset)

            if tz is not None and getattr(index, 'tz', None) is None:
                index = tslib.tz_localize_to_utc(com._ensure_int64(index), tz)
                index = index.view(_NS_DTYPE)

        index = index.view(cls)
        index.name = name
        index.offset = offset
        index.tz = tz

        return index

    def _box_values(self, values):
        return lib.map_infer(values, lib.Timestamp)

    def _local_timestamps(self):
        utc = _utc()

        if self.is_monotonic:
            return tslib.tz_convert(self.asi8, utc, self.tz)
        else:
            values = self.asi8
            indexer = values.argsort()
            result = tslib.tz_convert(values.take(indexer), utc, self.tz)

            n = len(indexer)
            reverse = np.empty(n, dtype=np.int_)
            reverse.put(indexer, np.arange(n))
            return result.take(reverse)

    @classmethod
    def _simple_new(cls, values, name, freq=None, tz=None):
        if values.dtype != _NS_DTYPE:
            values = com._ensure_int64(values).view(_NS_DTYPE)

        result = values.view(cls)
        result.name = name
        result.offset = freq
        result.tz = tools._maybe_get_tz(tz)

        return result

    @property
    def tzinfo(self):
        """
        Alias for tz attribute
        """
        return self.tz

    @classmethod
    def _cached_range(cls, start=None, end=None, periods=None, offset=None,
                      name=None):
        if start is None and end is None:
            # I somewhat believe this should never be raised externally and therefore
            # should be a `PandasError` but whatever...
            raise TypeError('Must specify either start or end.')
        if start is not None:
            start = Timestamp(start)
        if end is not None:
            end = Timestamp(end)
        if (start is None or end is None) and periods is None:
            raise TypeError('Must either specify period or provide both start and end.')

        if offset is None:
            # This can't happen with external-facing code, therefore PandasError
            raise TypeError('Must provide offset.')

        drc = _daterange_cache
        if offset not in _daterange_cache:
            xdr = generate_range(offset=offset, start=_CACHE_START,
                                 end=_CACHE_END)

            arr = tools.to_datetime(list(xdr), box=False)

            cachedRange = arr.view(DatetimeIndex)
            cachedRange.offset = offset
            cachedRange.tz = None
            cachedRange.name = None
            drc[offset] = cachedRange
        else:
            cachedRange = drc[offset]

        if start is None:
            if not (isinstance(end, Timestamp)):
                raise AssertionError()

            end = offset.rollback(end)

            endLoc = cachedRange.get_loc(end) + 1
            startLoc = endLoc - periods
        elif end is None:
            if not (isinstance(start, Timestamp)):
                raise AssertionError()

            start = offset.rollforward(start)

            startLoc = cachedRange.get_loc(start)
            endLoc = startLoc + periods
        else:
            if not offset.onOffset(start):
                start = offset.rollforward(start)

            if not offset.onOffset(end):
                end = offset.rollback(end)

            startLoc = cachedRange.get_loc(start)
            endLoc = cachedRange.get_loc(end) + 1

        indexSlice = cachedRange[startLoc:endLoc]
        indexSlice.name = name
        indexSlice.offset = offset

        return indexSlice

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return tslib.ints_to_pydatetime(self.asi8, self.tz)

    def __unicode__(self):
        from pandas.core.format import _format_datetime64
        values = self.values

        freq = None
        if self.offset is not None:
            freq = self.offset.freqstr

        summary = str(self.__class__)
        if len(self) == 1:
            first = _format_datetime64(values[0], tz=self.tz)
            summary += '\n[%s]' % first
        elif len(self) == 2:
            first = _format_datetime64(values[0], tz=self.tz)
            last = _format_datetime64(values[-1], tz=self.tz)
            summary += '\n[%s, %s]' % (first, last)
        elif len(self) > 2:
            first = _format_datetime64(values[0], tz=self.tz)
            last = _format_datetime64(values[-1], tz=self.tz)
            summary += '\n[%s, ..., %s]' % (first, last)

        tagline = '\nLength: %d, Freq: %s, Timezone: %s'
        summary += tagline % (len(self), freq, self.tz)

        return summary

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = self.name, self.offset, self.tz
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if len(state) == 2:
            nd_state, own_state = state
            self.name = own_state[0]
            self.offset = own_state[1]
            self.tz = own_state[2]
            np.ndarray.__setstate__(self, nd_state)
        else:  # pragma: no cover
            np.ndarray.__setstate__(self, state)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        elif isinstance(other, (DateOffset, timedelta)):
            return self._add_delta(other)
        elif isinstance(other, np.timedelta64):
            raise NotImplementedError
        elif com.is_integer(other):
            return self.shift(other)
        else:  # pragma: no cover
            raise TypeError(other)

    def __sub__(self, other):
        if isinstance(other, Index):
            return self.diff(other)
        elif isinstance(other, (DateOffset, timedelta)):
            return self._add_delta(-other)
        elif isinstance(other, np.timedelta64):
            raise NotImplementedError
        elif com.is_integer(other):
            return self.shift(-other)
        else:  # pragma: no cover
            raise TypeError(other)

    def _add_delta(self, delta):
        if isinstance(delta, (Tick, timedelta)):
            inc = offsets._delta_to_nanoseconds(delta)
            new_values = (self.asi8 + inc).view(_NS_DTYPE)
            tz = 'UTC' if self.tz is not None else None
            result = DatetimeIndex(new_values, tz=tz, freq='infer')
            utc = _utc()
            if self.tz is not None and self.tz is not utc:
                result = result.tz_convert(self.tz)
        else:
            new_values = self.astype('O') + delta
            result = DatetimeIndex(new_values, tz=self.tz, freq='infer')
        return result

    def __contains__(self, key):
        try:
            res = self.get_loc(key)
            return np.isscalar(res) or type(res) == slice
        except (KeyError, TypeError):
            return False

    def _format_with_header(self, header, **kwargs):
        return header + self._format_native_types(**kwargs)

    def _format_native_types(self, na_rep=u'NaT', **kwargs):
        data = list(self)

        # tz formatter or time formatter
        zero_time = time(0, 0)
        for d in data:
            if d.time() != zero_time or d.tzinfo is not None:
                return [u'%s' % x for x in data ]

        values = np.array(data,dtype=object)
        mask = isnull(self.values)
        values[mask] = na_rep

        imask = -mask
        values[imask] = np.array([ u'%d-%.2d-%.2d' % (dt.year, dt.month, dt.day) for dt in values[imask] ])
        return values.tolist()

    def isin(self, values):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values

        Parameters
        ----------
        values : set or sequence of values

        Returns
        -------
        is_contained : ndarray (boolean dtype)
        """
        if not isinstance(values, DatetimeIndex):
            try:
                values = DatetimeIndex(values)
            except ValueError:
                return self.asobject.isin(values)

        value_set = set(values.asi8)
        return lib.ismember(self.asi8, value_set)

    def to_datetime(self, dayfirst=False):
        return self.copy()

    def groupby(self, f):
        objs = self.asobject
        return _algos.groupby_object(objs, f)

    def summary(self, name=None):
        if len(self) > 0:
            index_summary = ', %s to %s' % (com.pprint_thing(self[0]),
                                            com.pprint_thing(self[-1]))
        else:
            index_summary = ''

        if name is None:
            name = type(self).__name__
        result = '%s: %s entries%s' % (com.pprint_thing(name),
                                       len(self), index_summary)
        if self.freq:
            result += '\nFreq: %s' % self.freqstr

        return result

    def get_duplicates(self):
        values = Index.get_duplicates(self)
        return DatetimeIndex(values)

    def astype(self, dtype):
        dtype = np.dtype(dtype)

        if dtype == np.object_:
            return self.asobject
        elif dtype == _INT64_DTYPE:
            return self.asi8.copy()
        else:  # pragma: no cover
            raise ValueError('Cannot cast DatetimeIndex to dtype %s' % dtype)

    def _get_time_micros(self):
        utc = _utc()
        values = self.asi8
        if self.tz is not None and self.tz is not utc:
            values = self._local_timestamps()
        return tslib.get_time_micros(values)

    @property
    def asobject(self):
        """
        Convert to Index of datetime objects
        """
        if isnull(self).any():
            msg = 'DatetimeIndex with NaT cannot be converted to object'
            raise ValueError(msg)
        return self._get_object_index()

    def tolist(self):
        """
        See ndarray.tolist
        """
        return list(self.asobject)

    def _get_object_index(self):
        boxfunc = lambda x: Timestamp(x, offset=self.offset, tz=self.tz)
        boxed_values = lib.map_infer(self.asi8, boxfunc)
        return Index(boxed_values, dtype=object)

    def to_pydatetime(self):
        """
        Return DatetimeIndex as object ndarray of datetime.datetime objects

        Returns
        -------
        datetimes : ndarray
        """
        return tslib.ints_to_pydatetime(self.asi8, tz=self.tz)

    def to_period(self, freq=None):
        """
        Cast to PeriodIndex at a particular frequency
        """
        from pandas.tseries.period import PeriodIndex

        if self.freq is None and freq is None:
            msg = "You must pass a freq argument as current index has none."
            raise ValueError(msg)

        if freq is None:
            freq = get_period_alias(self.freqstr)

        return PeriodIndex(self.values, freq=freq, tz=self.tz)

    def order(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index
        """
        if return_indexer:
            _as = self.argsort()
            if not ascending:
                _as = _as[::-1]
            sorted_index = self.take(_as)
            return sorted_index, _as
        else:
            sorted_values = np.sort(self.values)
            if not ascending:
                sorted_values = sorted_values[::-1]
            return self._simple_new(sorted_values, self.name, None,
                                    self.tz)

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
        return DatetimeIndex(snapped, freq=freq, verify_integrity=False)

    def shift(self, n, freq=None):
        """
        Specialized shift which produces a DatetimeIndex

        Parameters
        ----------
        n : int
            Periods to shift by
        freq : DateOffset or timedelta-like, optional

        Returns
        -------
        shifted : DatetimeIndex
        """
        if freq is not None and freq != self.offset:
            if isinstance(freq, basestring):
                freq = to_offset(freq)
            result = Index.shift(self, n, freq)
            result.tz = self.tz

            return result

        if n == 0:
            # immutable so OK
            return self

        if self.offset is None:
            raise ValueError("Cannot shift with no offset")

        start = self[0] + n * self.offset
        end = self[-1] + n * self.offset
        return DatetimeIndex(start=start, end=end, freq=self.offset,
                             name=self.name, tz=self.tz)

    def repeat(self, repeats, axis=None):
        """
        Analogous to ndarray.repeat
        """
        return DatetimeIndex(self.values.repeat(repeats),
                             name=self.name)

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take
        """
        maybe_slice = lib.maybe_indices_to_slice(com._ensure_int64(indices))
        if isinstance(maybe_slice, slice):
            return self[maybe_slice]
        indices = com._ensure_platform_int(indices)
        taken = self.values.take(indices, axis=axis)
        return self._simple_new(taken, self.name, None, self.tz)

    def unique(self):
        """
        Index.unique with handling for DatetimeIndex metadata

        Returns
        -------
        result : DatetimeIndex
        """
        result = Int64Index.unique(self)
        return DatetimeIndex._simple_new(result, tz=self.tz,
                                         name=self.name)

    def union(self, other):
        """
        Specialized union for DatetimeIndex objects. If combine
        overlapping ranges with the same DateOffset, will be much
        faster than Index.union

        Parameters
        ----------
        other : DatetimeIndex or array-like

        Returns
        -------
        y : Index or DatetimeIndex
        """
        if not isinstance(other, DatetimeIndex):
            try:
                other = DatetimeIndex(other)
            except TypeError:
                pass

        this, other = self._maybe_utc_convert(other)

        if this._can_fast_union(other):
            return this._fast_union(other)
        else:
            result = Index.union(this, other)
            if isinstance(result, DatetimeIndex):
                result.tz = this.tz
                if result.freq is None:
                    result.offset = to_offset(result.inferred_freq)
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
                tz = this.tz
                this = Index.union(this, other)
                if isinstance(this, DatetimeIndex):
                    this.tz = tz

        if this.freq is None:
            this.offset = to_offset(this.inferred_freq)
        return this

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
        to_concat, factory = _process_concat_data(to_concat, name)

        return factory(to_concat)

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        if (not isinstance(other, DatetimeIndex) and len(other) > 0 and
            other.inferred_type != 'mixed-integer'):
            try:
                other = DatetimeIndex(other)
            except TypeError:
                pass

        this, other = self._maybe_utc_convert(other)
        return Index.join(this, other, how=how, level=level,
                          return_indexers=return_indexers)

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

            if self.tz != other.tz:
                this = self.tz_convert('UTC')
                other = other.tz_convert('UTC')
        return this, other

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        if (isinstance(other, DatetimeIndex)
            and self.offset == other.offset
                and self._can_fast_union(other)):
            joined = self._view_like(joined)
            joined.name = name
            return joined
        else:
            tz = getattr(other, 'tz', None)
            return self._simple_new(joined, name, tz=tz)

    def _can_fast_union(self, other):
        if not isinstance(other, DatetimeIndex):
            return False

        offset = self.offset

        if offset is None or offset != other.offset:
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

        left_end = left[-1]
        right_start = right[0]

        # Only need to "adjoin", not overlap
        return (left_end + offset) >= right_start

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

        left_start, left_end = left[0], left[-1]
        right_end = right[-1]

        if not self.offset._should_cache():
            # concatenate dates
            if left_end < right_end:
                loc = right.searchsorted(left_end, side='right')
                right_chunk = right.values[loc:]
                dates = com._concat_compat((left.values, right_chunk))
                return self._view_like(dates)
            else:
                return left
        else:
            return type(self)(start=left_start,
                              end=max(left_end, right_end),
                              freq=left.offset)

    def __array_finalize__(self, obj):
        if self.ndim == 0:  # pragma: no cover
            return self.item()

        self.offset = getattr(obj, 'offset', None)
        self.tz = getattr(obj, 'tz', None)
        self.name = getattr(obj, 'name', None)

    def intersection(self, other):
        """
        Specialized intersection for DatetimeIndex objects. May be much faster
        than Index.intersection

        Parameters
        ----------
        other : DatetimeIndex or array-like

        Returns
        -------
        y : Index or DatetimeIndex
        """
        if not isinstance(other, DatetimeIndex):
            try:
                other = DatetimeIndex(other)
            except TypeError:
                pass
            result = Index.intersection(self, other)
            if isinstance(result, DatetimeIndex):
                if result.freq is None:
                    result.offset = to_offset(result.inferred_freq)
            return result

        elif (other.offset is None or self.offset is None or
              other.offset != self.offset or
              not other.offset.isAnchored() or
              (not self.is_monotonic or not other.is_monotonic)):
            result = Index.intersection(self, other)
            if isinstance(result, DatetimeIndex):
                if result.freq is None:
                    result.offset = to_offset(result.inferred_freq)
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
            return self._view_like(left_chunk)

    def _partial_date_slice(self, reso, parsed, use_lhs=True, use_rhs=True):

        is_monotonic = self.is_monotonic

        if reso == 'year':
            t1 = Timestamp(datetime(parsed.year, 1, 1), tz=self.tz)
            t2 = Timestamp(datetime(parsed.year, 12, 31, 23, 59, 59, 999999), tz=self.tz)
        elif reso == 'month':
            d = tslib.monthrange(parsed.year, parsed.month)[1]
            t1 = Timestamp(datetime(parsed.year, parsed.month, 1), tz=self.tz)
            t2 = Timestamp(datetime(parsed.year, parsed.month, d, 23, 59, 59, 999999), tz=self.tz)
        elif reso == 'quarter':
            qe = (((parsed.month - 1) + 2) % 12) + 1  # two months ahead
            d = tslib.monthrange(parsed.year, qe)[1]   # at end of month
            t1 = Timestamp(datetime(parsed.year, parsed.month, 1), tz=self.tz)
            t2 = Timestamp(datetime(parsed.year, qe, d, 23, 59, 59, 999999), tz=self.tz)
        elif (reso == 'day' and (self._resolution < Resolution.RESO_DAY or not is_monotonic)):
            st = datetime(parsed.year, parsed.month, parsed.day)
            t1 = Timestamp(st, tz=self.tz)
            t2 = st + offsets.Day()
            t2 = Timestamp(Timestamp(t2, tz=self.tz).value - 1)
        elif (reso == 'hour' and (
                self._resolution < Resolution.RESO_HR or not is_monotonic)):
            st = datetime(parsed.year, parsed.month, parsed.day,
                          hour=parsed.hour)
            t1 = Timestamp(st, tz=self.tz)
            t2 = Timestamp(Timestamp(st + offsets.Hour(),
                                     tz=self.tz).value - 1)
        elif (reso == 'minute' and (
                self._resolution < Resolution.RESO_MIN or not is_monotonic)):
            st = datetime(parsed.year, parsed.month, parsed.day,
                          hour=parsed.hour, minute=parsed.minute)
            t1 = Timestamp(st, tz=self.tz)
            t2 = Timestamp(Timestamp(st + offsets.Minute(),
                                     tz=self.tz).value - 1)
        elif (reso == 'second' and (
                self._resolution == Resolution.RESO_SEC or not is_monotonic)):
            st = datetime(parsed.year, parsed.month, parsed.day,
                          hour=parsed.hour, minute=parsed.minute, second=parsed.second)
            t1 = Timestamp(st, tz=self.tz)
            t2 = Timestamp(Timestamp(st + offsets.Second(),
                                     tz=self.tz).value - 1)
        else:
            raise KeyError


        stamps = self.asi8

        if is_monotonic:

            # we are out of range
            if len(stamps) and (
                (use_lhs and t1.value < stamps[0] and t2.value < stamps[0]) or (
                (use_rhs and t1.value > stamps[-1] and t2.value > stamps[-1]))):
                raise KeyError

            # a monotonic (sorted) series can be sliced
            left = stamps.searchsorted(t1.value, side='left') if use_lhs else None
            right = stamps.searchsorted(t2.value, side='right') if use_rhs else None

            return slice(left, right)

        lhs_mask = (stamps>=t1.value) if use_lhs else True
        rhs_mask = (stamps<=t2.value) if use_rhs else True

        # try to find a the dates
        return (lhs_mask & rhs_mask).nonzero()[0]

    def _possibly_promote(self, other):
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
            stamp = Timestamp(key, tz=self.tz)
            return self._engine.get_value(series, stamp)

        try:
            return Index.get_value(self, series, key)
        except KeyError:
            try:
                loc = self._get_string_slice(key)
                return series[loc]
            except (TypeError, ValueError, KeyError):
                pass

            if isinstance(key, time):
                locs = self.indexer_at_time(key)
                return series.take(locs)

            try:
                stamp = Timestamp(key, tz=self.tz)
                return self._engine.get_value(series, stamp)
            except (KeyError, ValueError):
                raise KeyError(key)

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        if isinstance(key, datetime):
            # needed to localize naive datetimes
            stamp = Timestamp(key, tz=self.tz)
            return self._engine.get_loc(stamp)

        try:
            return Index.get_loc(self, key)
        except (KeyError, ValueError):
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError, ValueError):
                pass

            if isinstance(key, time):
                return self.indexer_at_time(key)

            try:
                stamp = Timestamp(key, tz=self.tz)
                return self._engine.get_loc(stamp)
            except (KeyError, ValueError):
                raise KeyError(key)

    def _get_string_slice(self, key, use_lhs=True, use_rhs=True):
        freq = getattr(self, 'freqstr',
                       getattr(self, 'inferred_freq', None))
        _, parsed, reso = parse_time_string(key, freq)
        loc = self._partial_date_slice(reso, parsed, use_lhs=use_lhs, use_rhs=use_rhs)
        return loc

    def slice_indexer(self, start=None, end=None, step=None):
        """
        Index.slice_indexer, customized to handle time slicing
        """
        if isinstance(start, time) and isinstance(end, time):
            if step is not None and step != 1:
                raise ValueError('Must have step size of 1 with time slices')
            return self.indexer_between_time(start, end)

        if isinstance(start, time) or isinstance(end, time):
            raise KeyError('Cannot mix time and non-time slice keys')

        if isinstance(start, float) or isinstance(end, float):
            raise TypeError('Cannot index datetime64 with float keys')

        return Index.slice_indexer(self, start, end, step)

    def slice_locs(self, start=None, end=None):
        """
        Index.slice_locs, customized to handle partial ISO-8601 string slicing
        """
        if isinstance(start, basestring) or isinstance(end, basestring):

            if self.is_monotonic:
                try:
                    if start:
                        start_loc = self._get_string_slice(start).start
                    else:
                        start_loc = 0

                    if end:
                        end_loc = self._get_string_slice(end).stop
                    else:
                        end_loc = len(self)

                    return start_loc, end_loc
                except KeyError:
                    pass

            else:
                # can't use a slice indexer because we are not sorted!
                # so create an indexer directly
                try:
                    if start:
                        start_loc = self._get_string_slice(start,use_rhs=False)
                    else:
                        start_loc = np.arange(len(self))

                    if end:
                        end_loc = self._get_string_slice(end,use_lhs=False)
                    else:
                        end_loc = np.arange(len(self))

                    return start_loc, end_loc
                except KeyError:
                    pass

        if isinstance(start, time) or isinstance(end, time):
            raise KeyError('Cannot use slice_locs with time slice keys')

        return Index.slice_locs(self, start, end)

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            val = arr_idx[key]
            return Timestamp(val, offset=self.offset, tz=self.tz)
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)
                key = lib.maybe_booleans_to_slice(key.view(np.uint8))

            new_offset = None
            if isinstance(key, slice):
                if self.offset is not None and key.step is not None:
                    new_offset = key.step * self.offset
                else:
                    new_offset = self.offset

            result = arr_idx[key]
            if result.ndim > 1:
                return result

            return self._simple_new(result, self.name, new_offset, self.tz)

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, f):
        try:
            result = f(self)
            if not isinstance(result, np.ndarray):
                raise TypeError
            return result
        except Exception:
            return _algos.arrmap_object(self.asobject, f)

    # alias to offset
    @property
    def freq(self):
        return self.offset

    @cache_readonly
    def inferred_freq(self):
        try:
            return infer_freq(self)
        except ValueError:
            return None

    @property
    def freqstr(self):
        return self.offset.freqstr

    year = _field_accessor('year', 'Y')
    month = _field_accessor('month', 'M')
    day = _field_accessor('day', 'D')
    hour = _field_accessor('hour', 'h')
    minute = _field_accessor('minute', 'm')
    second = _field_accessor('second', 's')
    microsecond = _field_accessor('microsecond', 'us')
    nanosecond = _field_accessor('nanosecond', 'ns')
    weekofyear = _field_accessor('weekofyear', 'woy')
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 'dow')
    weekday = dayofweek
    dayofyear = _field_accessor('dayofyear', 'doy')
    quarter = _field_accessor('quarter', 'q')

    @property
    def time(self):
        """
        Returns numpy array of datetime.time. The time part of the Timestamps.
        """
        # can't call self.map() which tries to treat func as ufunc
        # and causes recursion warnings on python 2.6
        return _algos.arrmap_object(self.asobject, lambda x: x.time())

    @property
    def date(self):
        """
        Returns numpy array of datetime.date. The date part of the Timestamps.
        """
        return _algos.arrmap_object(self.asobject, lambda x: x.date())


    def normalize(self):
        """
        Return DatetimeIndex with times to midnight. Length is unaltered

        Returns
        -------
        normalized : DatetimeIndex
        """
        new_values = tslib.date_normalize(self.asi8, self.tz)
        return DatetimeIndex(new_values, freq='infer', name=self.name,
                             tz=self.tz)

    def __iter__(self):
        return iter(self._get_object_index())

    def searchsorted(self, key, side='left'):
        if isinstance(key, np.ndarray):
            key = np.array(key, dtype=_NS_DTYPE, copy=False)
        else:
            key = _to_m8(key, tz=self.tz)

        return self.values.searchsorted(key, side=side)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type or typ == 'datetime'

    def argmin(self):
        # hack to workaround argmin failure
        try:
            return self.values.argmin()
        except Exception:  # pragma: no cover
            return self.asi8.argmin()

    @property
    def inferred_type(self):
        # b/c datetime is represented as microseconds since the epoch, make
        # sure we can't have ambiguous indexing
        return 'datetime64'

    @property
    def dtype(self):
        return _NS_DTYPE

    @property
    def is_all_dates(self):
        return True

    @cache_readonly
    def is_normalized(self):
        """
        Returns True if all of the dates are at midnight ("no time")
        """
        return tslib.dates_normalized(self.asi8, self.tz)

    @cache_readonly
    def resolution(self):
        """
        Returns day, hour, minute, second, or microsecond
        """
        reso = self._resolution
        return get_reso_string(reso)

    @cache_readonly
    def _resolution(self):
        return tslib.resolution(self.asi8, self.tz)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self is other:
            return True

        if (not hasattr(other, 'inferred_type') or
                other.inferred_type != 'datetime64'):
            if self.offset is not None:
                return False
            try:
                other = DatetimeIndex(other)
            except:
                return False

        if self.tz is not None:
            if other.tz is None:
                return False
            same_zone = tslib.get_timezone(
                self.tz) == tslib.get_timezone(other.tz)
        else:
            if other.tz is not None:
                return False
            same_zone = True

        return same_zone and np.array_equal(self.asi8, other.asi8)

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : Index
        """
        if isinstance(item, datetime):
            item = _to_m8(item, tz=self.tz)

        new_index = np.concatenate((self[:loc].asi8,
                                    [item.view(np.int64)],
                                    self[loc:].asi8))
        return DatetimeIndex(new_index, freq='infer')

    def delete(self, loc):
        """
        Make new DatetimeIndex with passed location deleted

        Returns
        -------
        new_index : DatetimeIndex
        """
        arr = np.delete(self.values, loc)
        return DatetimeIndex(arr, tz=self.tz)

    def _view_like(self, ndarray):
        result = ndarray.view(type(self))
        result.offset = self.offset
        result.tz = self.tz
        result.name = self.name
        return result

    def tz_convert(self, tz):
        """
        Convert DatetimeIndex from one time zone to another (using pytz)

        Returns
        -------
        normalized : DatetimeIndex
        """
        tz = tools._maybe_get_tz(tz)

        if self.tz is None:
            # tz naive, use tz_localize
            raise TypeError('Cannot convert tz-naive timestamps, use '
                            'tz_localize to localize')

        # No conversion since timestamps are all UTC to begin with
        return self._simple_new(self.values, self.name, self.offset, tz)

    def tz_localize(self, tz):
        """
        Localize tz-naive DatetimeIndex to given time zone (using pytz)

        Returns
        -------
        localized : DatetimeIndex
        """
        if self.tz is not None:
            raise TypeError("Already tz-aware, use tz_convert to convert.")
        tz = tools._maybe_get_tz(tz)

        # Convert to UTC
        new_dates = tslib.tz_localize_to_utc(self.asi8, tz)
        new_dates = new_dates.view(_NS_DTYPE)

        return self._simple_new(new_dates, self.name, self.offset, tz)

    def indexer_at_time(self, time, asof=False):
        """
        Select values at particular time of day (e.g. 9:30AM)

        Parameters
        ----------
        time : datetime.time or string
        tz : string or pytz.timezone
            Time zone for time. Corresponding timestamps would be converted to
            time zone of the TimeSeries

        Returns
        -------
        values_at_time : TimeSeries
        """
        from dateutil.parser import parse

        if asof:
            raise NotImplementedError

        if isinstance(time, basestring):
            time = parse(time).time()

        if time.tzinfo:
            # TODO
            raise NotImplementedError

        time_micros = self._get_time_micros()
        micros = _time_to_micros(time)
        return (micros == time_micros).nonzero()[0]

    def indexer_between_time(self, start_time, end_time, include_start=True,
                             include_end=True):
        """
        Select values between particular times of day (e.g., 9:00-9:30AM)

        Parameters
        ----------
        start_time : datetime.time or string
        end_time : datetime.time or string
        include_start : boolean, default True
        include_end : boolean, default True
        tz : string or pytz.timezone, default None

        Returns
        -------
        values_between_time : TimeSeries
        """
        from dateutil.parser import parse

        if isinstance(start_time, basestring):
            start_time = parse(start_time).time()

        if isinstance(end_time, basestring):
            end_time = parse(end_time).time()

        if start_time.tzinfo or end_time.tzinfo:
            raise NotImplementedError

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

    def min(self, axis=None):
        """
        Overridden ndarray.min to return a Timestamp
        """
        if self.is_monotonic:
            return self[0]
        else:
            min_stamp = self.asi8.min()
            return Timestamp(min_stamp, tz=self.tz)

    def max(self, axis=None):
        """
        Overridden ndarray.max to return a Timestamp
        """
        if self.is_monotonic:
            return self[-1]
        else:
            max_stamp = self.asi8.max()
            return Timestamp(max_stamp, tz=self.tz)


def _generate_regular_range(start, end, periods, offset):
    if isinstance(offset, Tick):
        stride = offset.nanos
        if periods is None:
            b = Timestamp(start).value
            e = Timestamp(end).value
            e += stride - e % stride
            # end.tz == start.tz by this point due to _generate implementation
            tz = start.tz
        elif start is not None:
            b = Timestamp(start).value
            e = b + periods * stride
            tz = start.tz
        elif end is not None:
            e = Timestamp(end).value + stride
            b = e - periods * stride
            tz = end.tz
        else:
            raise NotImplementedError

        data = np.arange(b, e, stride, dtype=np.int64)
        data = DatetimeIndex._simple_new(data, None, tz=tz)
    else:
        if isinstance(start, Timestamp):
            start = start.to_pydatetime()

        if isinstance(end, Timestamp):
            end = end.to_pydatetime()

        xdr = generate_range(start=start, end=end,
                             periods=periods, offset=offset)

        dates = list(xdr)
        # utc = len(dates) > 0 and dates[0].tzinfo is not None
        data = tools.to_datetime(dates)

    return data


def date_range(start=None, end=None, periods=None, freq='D', tz=None,
               normalize=False, name=None):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer or None, default None
        If None, must specify start and end
    freq : string or DateOffset, default 'D' (calendar daily)
        Frequency strings can have multiples, e.g. '5H'
    tz : string or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Hong_Kong
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    name : str, default None
        Name of the resulting index

    Notes
    -----
    2 of start, end, or periods must be specified

    Returns
    -------
    rng : DatetimeIndex
    """
    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz, normalize=normalize, name=name)


def bdate_range(start=None, end=None, periods=None, freq='B', tz=None,
                normalize=True, name=None):
    """
    Return a fixed frequency datetime index, with business day as the default
    frequency

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer or None, default None
        If None, must specify start and end
    freq : string or DateOffset, default 'B' (business daily)
        Frequency strings can have multiples, e.g. '5H'
    tz : string or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Beijing
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    name : str, default None
        Name for the resulting index

    Notes
    -----
    2 of start, end, or periods must be specified

    Returns
    -------
    rng : DatetimeIndex
    """

    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz, normalize=normalize, name=name)


def cdate_range(start=None, end=None, periods=None, freq='C', tz=None,
                normalize=True, name=None, **kwargs):
    """
    **EXPERIMENTAL** Return a fixed frequency datetime index, with
    CustomBusinessDay as the default frequency

    .. warning:: EXPERIMENTAL

        The CustomBusinessDay class is not officially supported and the API is
        likely to change in future versions. Use this at your own risk.

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer or None, default None
        If None, must specify start and end
    freq : string or DateOffset, default 'C' (CustomBusinessDay)
        Frequency strings can have multiples, e.g. '5H'
    tz : string or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Beijing
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    name : str, default None
        Name for the resulting index
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``

    Notes
    -----
    2 of start, end, or periods must be specified

    Returns
    -------
    rng : DatetimeIndex
    """

    if freq=='C':
        holidays = kwargs.pop('holidays', [])
        weekmask = kwargs.pop('weekmask', 'Mon Tue Wed Thu Fri')
        freq = CDay(holidays=holidays, weekmask=weekmask)
    return DatetimeIndex(start=start, end=end, periods=periods, freq=freq,
                         tz=tz, normalize=normalize, name=name, **kwargs)


def _to_m8(key, tz=None):
    '''
    Timestamp-like => dt64
    '''
    if not isinstance(key, Timestamp):
        # this also converts strings
        key = Timestamp(key, tz=tz)

    return np.int64(tslib.pydt_to_i8(key)).view(_NS_DTYPE)


def _str_to_dt_array(arr, offset=None, dayfirst=None, yearfirst=None):
    def parser(x):
        result = parse_time_string(x, offset, dayfirst=dayfirst,
                                   yearfirst=yearfirst)
        return result[0]

    arr = np.asarray(arr, dtype=object)
    data = _algos.arrmap_object(arr, parser)
    return tools.to_datetime(data)


_CACHE_START = Timestamp(datetime(1950, 1, 1))
_CACHE_END = Timestamp(datetime(2030, 1, 1))

_daterange_cache = {}


def _naive_in_cache_range(start, end):
    if start is None or end is None:
        return False
    else:
        if start.tzinfo is not None or end.tzinfo is not None:
            return False
        return _in_range(start, end, _CACHE_START, _CACHE_END)


def _in_range(start, end, rng_start, rng_end):
    return start > rng_start and end < rng_end


def _time_to_micros(time):
    seconds = time.hour * 60 * 60 + 60 * time.minute + time.second
    return 1000000 * seconds + time.microsecond


def _process_concat_data(to_concat, name):
    klass = Index
    kwargs = {}
    concat = np.concatenate

    all_dti = True
    need_utc_convert = False
    has_naive = False
    tz = None

    for x in to_concat:
        if not isinstance(x, DatetimeIndex):
            all_dti = False
        else:
            if tz is None:
                tz = x.tz

            if x.tz is None:
                has_naive = True

            if x.tz != tz:
                need_utc_convert = True
                tz = 'UTC'

    if all_dti:
        need_obj_convert = False
        if has_naive and tz is not None:
            need_obj_convert = True

        if need_obj_convert:
            to_concat = [x.asobject.values for x in to_concat]

        else:
            if need_utc_convert:
                to_concat = [x.tz_convert('UTC').values for x in to_concat]
            else:
                to_concat = [x.values for x in to_concat]

            # well, technically not a "class" anymore...oh well
            klass = DatetimeIndex._simple_new
            kwargs = {'tz': tz}
            concat = com._concat_compat
    else:
        for i, x in enumerate(to_concat):
            if isinstance(x, DatetimeIndex):
                to_concat[i] = x.asobject.values
            elif isinstance(x, Index):
                to_concat[i] = x.values

    factory_func = lambda x: klass(concat(x), name=name, **kwargs)
    return to_concat, factory_func
