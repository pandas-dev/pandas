# pylint: disable=E1101

from datetime import time, datetime
from datetime import timedelta

import numpy as np

from pandas.core.index import Index, Int64Index
from pandas.tseries.frequencies import infer_freq, to_offset
from pandas.tseries.offsets import DateOffset, generate_range, Tick
from pandas.tseries.tools import parse_time_string, normalize_date
from pandas.util.decorators import cache_readonly
import pandas.core.common as com
import pandas.tseries.offsets as offsets
import pandas.tseries.tools as tools

from pandas._tseries import Timestamp
import pandas._tseries as lib
import pandas._algos as _algos

def _utc():
    import pytz
    return pytz.utc

# -------- some conversion wrapper functions

def _as_i8(arg):
    if isinstance(arg, np.ndarray) and arg.dtype == np.datetime64:
        return arg.view('i8', type=np.ndarray)
    else:
        return arg


def _field_accessor(name, field):
    def f(self):
        values = self.asi8
        if self.tz is not None:
            utc = _utc()
            if self.tz is not utc:
                values = lib.tz_convert(values, utc, self.tz)
        return lib.fast_field_accessor(values, field)
    f.__name__ = name
    return property(f)

def _wrap_i8_function(f):
    @staticmethod
    def wrapper(*args, **kwargs):
        view_args = [_as_i8(arg) for arg in args]
        return f(*view_args, **kwargs)
    return wrapper

def _wrap_dt_function(f):
    @staticmethod
    def wrapper(*args, **kwargs):
        view_args = [_dt_box_array(_as_i8(arg)) for arg in args]
        return f(*view_args, **kwargs)
    return wrapper

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
        if isinstance(other, datetime):
            func = getattr(self, opname)
            result = func(_to_m8(other))
        elif isinstance(other, np.ndarray):
            func = getattr(super(DatetimeIndex, self), opname)
            result = func(other)
        else:
            other = _ensure_datetime64(other)
            func = getattr(super(DatetimeIndex, self), opname)
            result = func(other)
        try:
            return result.view(np.ndarray)
        except:
            return result
    return wrapper

def _ensure_datetime64(other):
    if isinstance(other, np.datetime64):
        return other
    elif com.is_integer(other):
        return np.int64(other).view('M8[us]')
    else:
        raise TypeError(other)

def _dt_index_op(opname):
    """
    Wrap arithmetic operations to convert timedelta to a timedelta64.
    """
    def wrapper(self, other):
        if isinstance(other, timedelta):
            func = getattr(self, opname)
            return func(np.timedelta64(other))
        else:
            func = getattr(super(DatetimeIndex, self), opname)
            return func(other)
    return wrapper


class TimeSeriesError(Exception):
    pass


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
    _left_indexer  = _join_i8_wrapper(_algos.left_join_indexer_int64,
                                      with_indexers=False)
    _groupby = lib.groupby_arrays # _wrap_i8_function(lib.groupby_int64)

    _arrmap = _wrap_dt_function(_algos.arrmap_object)

    __eq__ = _dt_index_cmp('__eq__')
    __ne__ = _dt_index_cmp('__ne__')
    __lt__ = _dt_index_cmp('__lt__')
    __gt__ = _dt_index_cmp('__gt__')
    __le__ = _dt_index_cmp('__le__')
    __ge__ = _dt_index_cmp('__ge__')

    # structured array cache for datetime fields
    _sarr_cache = None

    _engine_type = lib.DatetimeEngine

    offset = None

    def __new__(cls, data=None,
                freq=None, start=None, end=None, periods=None,
                copy=False, name=None, tz=None,
                verify_integrity=True, normalize=False, **kwds):

        warn = False
        if 'offset' in kwds and kwds['offset']:
            freq = kwds['offset']
            warn = True

        infer_freq = False
        if not isinstance(freq, DateOffset):
            if freq != 'infer':
                freq = to_offset(freq)
            else:
                infer_freq = True
                freq = None

        if warn:
            import warnings
            warnings.warn("parameter 'offset' is deprecated, "
                          "please use 'freq' instead",
                          FutureWarning)
            if isinstance(freq, basestring):
                freq = to_offset(freq)
        else:
            if isinstance(freq, basestring):
                freq = to_offset(freq)

        offset = freq

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

            if isinstance(data, datetime):
                data = [data]

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            data = np.asarray(data, dtype='O')

            # try a few ways to make it datetime64
            if lib.is_string_array(data):
                data = _str_to_dt_array(data, offset)
            else:
                data = tools.to_datetime(data)
                data = np.asarray(data, dtype='M8[ns]')

        if issubclass(data.dtype.type, basestring):
            subarr = _str_to_dt_array(data, offset)
        elif issubclass(data.dtype.type, np.datetime64):
            if isinstance(data, DatetimeIndex):
                subarr = data.values
                offset = data.offset
                verify_integrity = False
            else:
                subarr = np.array(data, dtype='M8[ns]', copy=copy)
        elif issubclass(data.dtype.type, np.integer):
            subarr = np.array(data, dtype='M8[ns]', copy=copy)
        else:
            subarr = tools.to_datetime(data)
            if not np.issubdtype(subarr.dtype, np.datetime64):
                raise TypeError('Unable to convert %s to datetime dtype'
                                % str(data))

        if tz is not None:
            tz = tools._maybe_get_tz(tz)
            # Convert local to UTC
            ints = subarr.view('i8')
            lib.tz_localize_check(ints, tz)
            subarr = lib.tz_convert(ints, tz, _utc())
            subarr = subarr.view('M8[ns]')

        subarr = subarr.view(cls)
        subarr.name = name
        subarr.offset = offset
        subarr.tz = tz

        if verify_integrity and len(subarr) > 0:
            if offset is not None and not infer_freq:
                inferred = subarr.inferred_freq
                if inferred != offset.freqstr:
                    raise ValueError('Dates do not conform to passed '
                                     'frequency')

        if infer_freq:
            inferred = subarr.inferred_freq
            if inferred:
                subarr.offset = to_offset(inferred)

        return subarr

    @classmethod
    def _generate(cls, start, end, periods, name, offset,
                  tz=None, normalize=False):
        _normalized = True

        if start is not None:
            start = Timestamp(start)
            if not isinstance(start, Timestamp):
                raise ValueError('Failed to convert %s to timestamp'
                                 % start)

            if normalize:
                start = normalize_date(start)
                _normalized = True
            else:
                _normalized = _normalized and start.time() == _midnight

        if end is not None:
            end = Timestamp(end)
            if not isinstance(end, Timestamp):
                raise ValueError('Failed to convert %s to timestamp'
                                 % end)

            if normalize:
                end = normalize_date(end)
                _normalized = True
            else:
                _normalized = _normalized and end.time() == _midnight

        start, end, tz = tools._figure_out_timezone(start, end, tz)

        if (offset._should_cache() and
            not (offset._normalize_cache and not _normalized) and
            _naive_in_cache_range(start, end)):
            index = cls._cached_range(start, end, periods=periods,
                                      offset=offset, name=name)
        else:
            index = _generate_regular_range(start, end, periods, offset)

        if tz is not None:
            # Convert local to UTC
            ints = index.view('i8')
            lib.tz_localize_check(ints, tz)
            index = lib.tz_convert(ints, tz, _utc())
            index = index.view('M8[ns]')

        index = index.view(cls)
        index.name = name
        index.offset = offset
        index.tz = tz

        return index

    @classmethod
    def _simple_new(cls, values, name, freq=None, tz=None):
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
        if start is not None:
            start = Timestamp(start)
        if end is not None:
            end = Timestamp(end)

        if offset is None:
            raise Exception('Must provide a DateOffset!')

        drc = _daterange_cache
        if offset not in _daterange_cache:
            xdr = generate_range(offset=offset, start=_CACHE_START,
                                 end=_CACHE_END)

            arr = np.array(_to_m8_array(list(xdr)),
                           dtype='M8[ns]', copy=False)

            cachedRange = arr.view(DatetimeIndex)
            cachedRange.offset = offset
            cachedRange.tz = None
            cachedRange.name = None
            drc[offset] = cachedRange
        else:
            cachedRange = drc[offset]

        if start is None:
            if end is None:
                raise Exception('Must provide start or end date!')
            if periods is None:
                raise Exception('Must provide number of periods!')

            assert(isinstance(end, Timestamp))

            end = offset.rollback(end)

            endLoc = cachedRange.get_loc(end) + 1
            startLoc = endLoc - periods
        elif end is None:
            assert(isinstance(start, Timestamp))
            start = offset.rollforward(start)

            startLoc = cachedRange.get_loc(start)
            if periods is None:
                raise Exception('Must provide number of periods!')

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
        return self.values.astype('O')

    def __repr__(self):
        if self.offset is not None:
            output = str(self.__class__)
            if len(self) > 0:
                output += '\n[%s, ..., %s]' % (self[0], self[-1])
            tagline = '\nLength: %d, Freq: %s, Timezone: %s'
            output += tagline % (len(self), self.offset.freqstr, self.tz)
            return output
        else:
            return Index.__repr__(self)

    __str__ = __repr__

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
        elif len(state) == 3:
            # legacy format: daterange
            offset = state[1]

            if len(state) > 2:
                tzinfo = state[2]
            else: # pragma: no cover
                tzinfo = None

            self.offset = offset
            self.tzinfo = tzinfo

            # extract the raw datetime data, turn into datetime64
            index_state = state[0]
            raw_data = index_state[0][4]
            raw_data = np.array(raw_data, dtype='M8[ns]')
            new_state = raw_data.__reduce__()
            np.ndarray.__setstate__(self, new_state[2])
        else:  # pragma: no cover
            np.ndarray.__setstate__(self, state)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        elif isinstance(other, (DateOffset, timedelta)):
            return self._add_delta(other)
        elif com.is_integer(other):
            return self.shift(other)
        else:
            return Index(self.view(np.ndarray) + other)

    def __sub__(self, other):
        if isinstance(other, Index):
            return self.diff(other)
        elif isinstance(other, (DateOffset, timedelta)):
            return self._add_delta(-other)
        elif com.is_integer(other):
            return self.shift(-other)
        else:
            return Index(self.view(np.ndarray) - other)

    def _add_delta(self, delta):
        if isinstance(delta, (Tick, timedelta)):
            inc = offsets._delta_to_nanoseconds(delta)
            new_values = (self.asi8 + inc).view('M8[ns]')
        else:
            new_values = self.astype('O') + delta
        return DatetimeIndex(new_values, tz=self.tz, freq='infer')

    def summary(self, name=None):
        if len(self) > 0:
            index_summary = ', %s to %s' % (str(self[0]), str(self[-1]))
        else:
            index_summary = ''

        if name is None:
            name = type(self).__name__
        result = '%s: %s entries%s' % (name, len(self), index_summary)
        if self.freq:
            result += '\nFreq: %s' % self.freqstr

        return result

    def astype(self, dtype):
        dtype = np.dtype(dtype)

        if dtype == np.object_:
            return self.asobject
        return Index.astype(self, dtype)

    @property
    def asi8(self):
        # do not cache or you'll create a memory leak
        return self.values.view('i8')

    @property
    def asstruct(self):
        if self._sarr_cache is None:
            self._sarr_cache = lib.build_field_sarray(self.asi8)
        return self._sarr_cache

    @property
    def asobject(self):
        """
        Convert to Index of datetime objects
        """
        boxed_values = _dt_box_array(self.asi8, self.offset, self.tz)
        return Index(boxed_values, dtype=object)

    def to_period(self, freq=None):
        """
        Cast to PeriodIndex at a particular frequency
        """
        from pandas.tseries.period import PeriodIndex

        if self.freq is None and freq is None:
            msg = "You must pass a freq argument as current index has none."
            raise ValueError(msg)

        if freq is None:
            freq = self.freqstr

        return PeriodIndex(self.values, freq=freq)

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
            return self._simple_new(sorted_values, self.name, None,
                                    self.tz)

    def snap(self, freq='S'):
        """
        Snap time stamps to nearest occuring frequency

        """
        # Superdumb, punting on any optimizing
        freq = to_offset(freq)

        snapped = np.empty(len(self), dtype='M8[ns]')

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
            return Index.shift(self, n, freq)

        if n == 0:
            # immutable so OK
            return self

        if self.offset is None:
            raise ValueError("Cannot shift with no offset")

        start = self[0] + n * self.offset
        end = self[-1] + n * self.offset
        return DatetimeIndex(start=start, end=end, freq=self.offset,
                             name=self.name)

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
        return DatetimeIndex(taken, tz=self.tz, name=self.name)

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
                result.tz = self.tz
                if result.freq is None:
                    result.offset = to_offset(result.inferred_freq)
            return result

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        if not isinstance(other, DatetimeIndex) and len(other) > 0:
            try:
                other = DatetimeIndex(other)
            except ValueError:
                pass

        this, other = self._maybe_utc_convert(other)
        return Index.join(this, other, how=how, level=level,
                          return_indexers=return_indexers)

    def _maybe_utc_convert(self, other):
        this = self
        if isinstance(other, DatetimeIndex):
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
            return DatetimeIndex(joined, name=name)

    def _can_fast_union(self, other):
        if not isinstance(other, DatetimeIndex):
            return False

        offset = self.offset

        if offset is None:
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
                dates = np.concatenate((left.values, right_chunk))
                return self._view_like(dates)
            else:
                return left
        else:
            return type(self)(start=left_start,
                              end=max(left_end, right_end),
                              freq=left.offset)

    def __array_finalize__(self, obj):
        if self.ndim == 0: # pragma: no cover
            return self.item()

        self.offset = getattr(obj, 'offset', None)
        self.tz = getattr(obj, 'tz', None)

    def intersection(self, other):
        """
        Specialized intersection for DatetimeIndex objects. May be much faster
        than Index.union

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

        elif other.offset != self.offset or (not self.is_monotonic or
                                             not other.is_monotonic):
            result = Index.intersection(self, other)
            if isinstance(result, DatetimeIndex):
                if result.freq is None:
                    result.offset = to_offset(result.inferred_freq)
            return result

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

    def _partial_date_slice(self, reso, parsed):
        if not self.is_monotonic:
            raise TimeSeriesError('Partial indexing only valid for ordered time'
                                  ' series')

        if reso == 'year':
            t1 = Timestamp(datetime(parsed.year, 1, 1))
            t2 = Timestamp(datetime(parsed.year, 12, 31))
        elif reso == 'month':
            d = lib.monthrange(parsed.year, parsed.month)[1]
            t1 = Timestamp(datetime(parsed.year, parsed.month, 1))
            t2 = Timestamp(datetime(parsed.year, parsed.month, d))
        elif reso == 'quarter':
            qe = (((parsed.month - 1) + 2) % 12) + 1 # two months ahead
            d = lib.monthrange(parsed.year, qe)[1]   # at end of month
            t1 = Timestamp(datetime(parsed.year, parsed.month, 1))
            t2 = Timestamp(datetime(parsed.year, qe, d))
        else:
            raise KeyError

        stamps = self.asi8
        left = stamps.searchsorted(t1.value, side='left')
        right = stamps.searchsorted(t2.value, side='right')
        return slice(left, right)

    def _possibly_promote(self, other):
        if other.inferred_type == 'date':
            other = DatetimeIndex(other)
        return self, other

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        try:
            return Index.get_value(self, series, key)
        except KeyError:

            try:
                freq = getattr(self, 'freq', getattr(self, 'inferred_freq',
                                                     None))
                asdt, parsed, reso = parse_time_string(key, freq)
                key = asdt
                loc = self._partial_date_slice(reso, parsed)
                return series[loc]
            except (TypeError, ValueError, KeyError):
                pass

            if isinstance(key, time):
                locs = self._indices_at_time(key)
                return series.take(locs)

            stamp = Timestamp(key)
            try:
                return self._engine.get_value(series, stamp)
            except KeyError:
                raise KeyError(stamp)

    def get_loc(self, key):
        """y
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        try:
            return self._engine.get_loc(key)
        except KeyError:
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError):
                pass

            if isinstance(key, time):
                return self._indices_at_time(key)

            stamp = Timestamp(key)
            try:
                return self._engine.get_loc(stamp)
            except KeyError:
                raise KeyError(stamp)

    def _indices_at_time(self, key):
        from dateutil.parser import parse

        # TODO: time object with tzinfo?

        nanos = _time_to_nanosecond(key)
        indexer = lib.values_at_time(self.asi8, nanos)
        return com._ensure_platform_int(indexer)

    def _get_string_slice(self, key):
        freq = getattr(self, 'freq', getattr(self, 'inferred_freq', None))
        asdt, parsed, reso = parse_time_string(key, freq)
        key = asdt
        loc = self._partial_date_slice(reso, parsed)
        return loc

    def slice_locs(self, start=None, end=None):
        """
        Index.slice_locs, customized to handle partial ISO-8601 string slicing
        """
        if isinstance(start, basestring) or isinstance(end, basestring):
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
            return f(self)
        except:
            return Index.map(self, f)

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
    weekofyear = _field_accessor('weekofyear', 'woy')
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 'dow')
    weekday = dayofweek
    dayofyear = _field_accessor('dayofyear', 'doy')
    quarter = _field_accessor('quarter', 'q')

    def normalize(self):
        """
        Return DatetimeIndex with times to midnight. Length is unaltered

        Returns
        -------
        normalized : DatetimeIndex
        """
        new_values = lib.date_normalize(self.asi8)
        return DatetimeIndex(new_values, freq='infer', name=self.name)

    def __iter__(self):
        return iter(self.asobject)

    def searchsorted(self, key, side='left'):
        if isinstance(key, np.ndarray):
            key = np.array(key, dtype='M8[ns]', copy=False)
        else:
            key = _to_m8(key)

        return self.values.searchsorted(key, side=side)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type or typ == 'datetime'

    # hack to workaround argmin failure
    def argmin(self):
        return (-self).argmax()

    @property
    def inferred_type(self):
        # b/c datetime is represented as microseconds since the epoch, make
        # sure we can't have ambiguous indexing
        return 'datetime64'

    @property
    def _constructor(self):
        return DatetimeIndex

    @property
    def dtype(self):
        return np.dtype('M8[ns]')

    @property
    def is_all_dates(self):
        return True

    @cache_readonly
    def is_normalized(self):
        """
        Returns True if all of the dates are at midnight ("no time")
        """
        return lib.dates_normalized(self.asi8)

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

        return self.tz == other.tz and np.array_equal(self.asi8, other.asi8)

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
        if type(item) == datetime:
            item = _to_m8(item)

        if self.offset is not None and not self.offset.onOffset(item):
            raise ValueError("Cannot insert value at non-conforming time")

        return super(DatetimeIndex, self).insert(loc, item)

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
            return self.tz_localize(tz)

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
            raise ValueError("Already have timezone info, "
                             "use tz_convert to convert.")
        tz = tools._maybe_get_tz(tz)

        lib.tz_localize_check(self.asi8, tz)

        # Convert to UTC
        new_dates = lib.tz_convert(self.asi8, tz, _utc())
        new_dates = new_dates.view('M8[ns]')
        return self._simple_new(new_dates, self.name, self.offset, tz)

    def tz_validate(self):
        """
        For a localized time zone, verify that there are no DST ambiguities
        (using pytz)

        Returns
        -------
        result : boolean
            True if there are no DST ambiguities
        """
        import pytz

        if self.tz is None or self.tz is pytz.utc:
            return True

        # See if there are any DST resolution problems
        try:
            lib.tz_localize_check(self.asi8, self.tz)
        except:
            return False

        return True

def _generate_regular_range(start, end, periods, offset):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify two of start, end, or periods')

    if isinstance(offset, Tick):
        stride = offset.nanos
        if periods is None:
            b = Timestamp(start).value
            e = Timestamp(end).value
            e += stride - e % stride
        elif start is not None:
            b = Timestamp(start).value
            e = b + periods * stride
        elif end is not None:
            e = Timestamp(end).value + stride
            b = e - periods * stride
        else:
            raise NotImplementedError

        data = np.arange(b, e, stride, dtype=np.int64)
        data = data.view('M8[ns]')
    else:
        xdr = generate_range(start=start, end=end,
            periods=periods, offset=offset)

        data = np.array(list(xdr), dtype='M8[ns]')

    return data


def date_range(start=None, end=None, periods=None, freq='D', tz=None,
               normalize=False):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency


    Parameters
    ----------
    start :
    end :
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range

    Returns
    -------

    """
    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz, normalize=normalize)


def bdate_range(start=None, end=None, periods=None, freq='B', tz=None,
                normalize=True):
    """
    Return a fixed frequency datetime index, with business day as the default
    frequency

    Parameters
    ----------

    normalize : bool, default False
        Normalize start/end dates to midnight before generating date
        range. Defaults to True for legacy reasons

    Returns
    -------
    date_range : DatetimeIndex

    """

    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz, normalize=normalize)


def _dt_box_array(arr, offset=None, tz=None):
    if arr is None:
        return arr

    if not isinstance(arr, np.ndarray):
        return arr

    boxfunc = lambda x: Timestamp(x, offset=offset, tz=tz)
    return lib.map_infer(arr, boxfunc)


def _to_m8(key):
    '''
    Timestamp-like => dt64
    '''
    if not isinstance(key, datetime):
        # this also converts strings
        key = Timestamp(key)

    return np.int64(lib.pydt_to_i8(key)).view('M8[ns]')


def _to_m8_array(arr):
    if arr is None:
        return arr
    return np.frompyfunc(_to_m8, 1, 1)(arr)


def _str_to_dt_array(arr, offset=None):
    def parser(x):
        result = parse_time_string(x, offset)
        return result[0]

    p_ufunc = np.frompyfunc(parser, 1, 1)
    data = p_ufunc(arr)
    return np.array(data, dtype='M8[ns]')


_CACHE_START = Timestamp(datetime(1950, 1, 1))
_CACHE_END   = Timestamp(datetime(2030, 1, 1))

_daterange_cache = {}


def _naive_in_cache_range(start, end):
    if start is None or end is None:
        return False
    else:
        return _in_range(start, end, _CACHE_START, _CACHE_END)

def _in_range(start, end, rng_start, rng_end):
    return start > rng_start and end < rng_end

def _time_to_nanosecond(time):
    seconds = time.hour * 60 * 60 + 60 * time.minute + time.second
    return (1000000 * seconds + time.microsecond) * 1000
