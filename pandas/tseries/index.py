from datetime import time, datetime
from datetime import timedelta

import numpy as np

from pandas.core.index import Index, Int64Index
from pandas.tseries.tools import parse_time_string
import pandas.core.common as com
import pandas.core.datetools as datetools
import pandas.tseries.tools as tools

from pandas._engines import DatetimeEngine
from pandas._tseries import Timestamp
import pandas._tseries as lib

# -------- some conversion wrapper functions

def _as_i8(arg):
    if isinstance(arg, np.ndarray) and arg.dtype == np.datetime64:
        return arg.view('i8', type=np.ndarray)
    else:
        return arg

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
            join_index = join_index.view('M8')
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
        return np.datetime64(other)
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
    dtype : NumPy dtype (default: M8[us])
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

    _inner_indexer = _join_i8_wrapper(lib.inner_join_indexer_int64)
    _outer_indexer = _join_i8_wrapper(lib.outer_join_indexer_int64)
    _left_indexer  = _join_i8_wrapper(lib.left_join_indexer_int64,
                                      with_indexers=False)
    _groupby = lib.groupby_arrays # _wrap_i8_function(lib.groupby_int64)

    _arrmap = _wrap_dt_function(lib.arrmap_object)

    __eq__ = _dt_index_cmp('__eq__')
    __ne__ = _dt_index_cmp('__ne__')
    __lt__ = _dt_index_cmp('__lt__')
    __gt__ = _dt_index_cmp('__gt__')
    __le__ = _dt_index_cmp('__le__')
    __ge__ = _dt_index_cmp('__ge__')

    __add__ = _dt_index_op('__add__')
    __sub__ = _dt_index_op('__sub__')

    # structured array cache for datetime fields
    _sarr_cache = None

    _engine_type = DatetimeEngine

    offset = None

    def __new__(cls, data=None,
                freq=None, start=None, end=None, periods=None,
                dtype=None, copy=False, name=None, tz=None,
                verify_integrity=True, normalize=False, **kwds):

        warn = False
        if 'offset' in kwds and kwds['offset']:
            freq = kwds['offset']
            warn = True

        if not isinstance(freq, datetools.DateOffset):
            freq = datetools.to_offset(freq)

        if warn:
            import warnings
            warnings.warn("parameter 'offset' is deprecated, "
                          "please use 'freq' instead",
                          FutureWarning)
            if isinstance(freq, basestring):
                freq = datetools.get_offset(freq)
        else:
            if isinstance(freq, basestring):
                freq = datetools.to_offset(freq)

        offset = freq

        if data is None and offset is None:
            raise ValueError("Must provide freq argument if no data is "
                             "supplied")

        if data is None:
            _normalized = True

            if start is not None:
                start = Timestamp(start)
                if not isinstance(start, Timestamp):
                    raise ValueError('Failed to convert %s to timestamp'
                                     % start)

                if normalize:
                    start = datetools.normalize_date(start)
                    _normalized = True
                else:
                    _normalized = _normalized and start.time() == _midnight

            if end is not None:
                end = Timestamp(end)
                if not isinstance(end, Timestamp):
                    raise ValueError('Failed to convert %s to timestamp'
                                     % end)

                if normalize:
                    end = datetools.normalize_date(end)
                    _normalized = True
                else:
                    _normalized = _normalized and end.time() == _midnight

            start, end, tz = tools._figure_out_timezone(start, end, tz)

            if (offset._should_cache() and
                not (offset._normalize_cache and not _normalized) and
                datetools._naive_in_cache_range(start, end)):
                index = cls._cached_range(start, end, periods=periods,
                                          offset=offset, name=name)
            else:
                index = _generate_regular_range(start, end, periods, offset)

            index = index.view(cls)
            index.name = name
            index.offset = offset
            index.tz = tz

            return index

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
                data = _str_to_dt_array(data)
            else:
                data = np.asarray(data, dtype='M8[us]')

        if issubclass(data.dtype.type, basestring):
            subarr = _str_to_dt_array(data)
        elif issubclass(data.dtype.type, np.integer):
            subarr = np.array(data, dtype='M8[us]', copy=copy)
        elif issubclass(data.dtype.type, np.datetime64):
            subarr = np.array(data, dtype='M8[us]', copy=copy)
        else:
            subarr = np.array(data, dtype='M8[us]', copy=copy)

        # TODO: this is horribly inefficient. If user passes data + offset, we
        # need to make sure data points conform. Punting on this

        if verify_integrity:
            if offset is not None:
                for i, ts in enumerate(subarr):
                    if not offset.onOffset(Timestamp(ts)):
                        val = Timestamp(offset.rollforward(ts)).value
                        subarr[i] = val

        subarr = subarr.view(cls)
        subarr.name = name
        subarr.offset = offset
        subarr.tz = tz

        return subarr

    @classmethod
    def _simple_new(cls, values, name, offset, tz):
        result = values.view(cls)
        result.name = name
        result.offset = offset
        result.tz = tz

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
            xdr = datetools.generate_range(offset=offset,
                    start=_CACHE_START, end=_CACHE_END)

            arr = np.array(_to_m8_array(list(xdr)),
                           dtype='M8[us]', copy=False)

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
            output = str(self.__class__) + '\n'
            output += 'freq: %s, timezone: %s\n' % (self.offset, self.tz)
            if len(self) > 0:
                output += '[%s, ..., %s]\n' % (self[0], self[-1])
            output += 'length: %d' % len(self)
            return output
        else:
            return super(DatetimeIndex, self).__repr__()

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
            raw_data = np.array(raw_data, dtype='M8[us]')
            new_state = raw_data.__reduce__()
            np.ndarray.__setstate__(self, new_state[2])
        else:  # pragma: no cover
            np.ndarray.__setstate__(self, state)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        elif isinstance(other, (datetools.DateOffset, timedelta)):
            new_values = self.astype('O') + other
            return DatetimeIndex(new_values, tz=self.tz)
        else:
            return Index(self.view(np.ndarray) + other)

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
            freq = self.freq

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
        freq = datetools.to_offset(freq)

        snapped = np.empty(len(self), dtype='M8[us]')

        for i, v in enumerate(self):
            s = v
            if not freq.onOffset(s):
                t0 = freq.rollback(s)
                t1 = freq.rollforward(s)
                if abs(s - t0) < abs(t1 - s):
                    s = t0
                else:
                    s = t1
            snapped[i] = np.datetime64(s)

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
                freq = datetools.to_offset(freq)
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
        if self._can_fast_union(other):
            return self._fast_union(other)
        else:
            return Index.union(self, other)

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        if not isinstance(other, DatetimeIndex) and len(other) > 0:
            try:
                other = DatetimeIndex(other)
            except ValueError:
                pass

        return Index.join(self, other, how=how, level=level,
                          return_indexers=return_indexers)


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

        if len(other) == 0:
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
            return Index.intersection(self, other)
        elif other.offset != self.offset:
            return Index.intersection(self, other)

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

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        try:
            return Index.get_value(self, series, key)
        except KeyError:

            try:
                asdt, parsed, reso = datetools.parse_time_string(key)
                key = asdt
                loc = self._partial_date_slice(reso, parsed)
                return series[loc]
            except (TypeError, ValueError, KeyError):
                pass

            stamp = Timestamp(key)
            try:
                return self._engine.get_value(series, stamp)
            except KeyError:
                raise KeyError(stamp)

    def get_loc(self, key):
        """
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

            stamp = Timestamp(key)
            try:
                return self._engine.get_loc(stamp)
            except KeyError:
                raise KeyError(stamp)

    def _get_string_slice(self, key):
        asdt, parsed, reso = datetools.parse_time_string(key)
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
                key = lib.maybe_booleans_to_slice(key)

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

    # Fast field accessors for periods of datetime index
    # --------------------------------------------------------------

    @property
    def year(self):
        return lib.fast_field_accessor(self.asi8, 'Y')

    @property
    def month(self):
        return lib.fast_field_accessor(self.asi8, 'M')

    @property
    def day(self):
        return lib.fast_field_accessor(self.asi8, 'D')

    @property
    def hour(self):
        return lib.fast_field_accessor(self.asi8, 'h')

    @property
    def minute(self):
        return lib.fast_field_accessor(self.asi8, 'm')

    @property
    def second(self):
        return lib.fast_field_accessor(self.asi8, 's')

    @property
    def microsecond(self):
        return lib.fast_field_accessor(self.asi8, 'us')

    @property
    def weekofyear(self):
        return lib.fast_field_accessor(self.asi8, 'woy')

    @property
    def dayofweek(self):
        return lib.fast_field_accessor(self.asi8, 'dow')

    @property
    def dayofyear(self):
        return lib.fast_field_accessor(self.asi8, 'doy')

    @property
    def quarter(self):
        return lib.fast_field_accessor(self.asi8, 'q')

    def __iter__(self):
        return iter(self.asobject)

    def searchsorted(self, key, side='left'):
        if isinstance(key, np.ndarray):
            key = np.array(key, dtype='M8[us]', copy=False)
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
        return np.dtype('M8')

    @property
    def is_all_dates(self):
        return True

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

        return np.array_equal(self.asi8, other.asi8)

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

    def tz_normalize(self, tz):
        """
        Convert DatetimeIndex from one time zone to another (using pytz)

        Returns
        -------
        normalized : DatetimeIndex
        """
        new_dates = lib.tz_normalize_array(self.asi8, self.tz, tz)
        new_dates = new_dates.view('M8[us]')
        new_dates = new_dates.view(self.__class__)
        new_dates.offset = self.offset
        new_dates.tz = tz
        new_dates.name = self.name
        return new_dates

    def tz_localize(self, tz):
        """
        Localize tz-naive DatetimeIndex to given time zone (using pytz)

        Returns
        -------
        localized : DatetimeIndex
        """
        if self.tz is not None:
            raise ValueError("Already have timezone info, "
                             "use tz_normalize to convert.")

        new_dates = lib.tz_localize_array(self.asi8, tz)
        new_dates = new_dates.view('M8[us]')
        new_dates = new_dates.view(self.__class__)
        new_dates.offset = self.offset
        new_dates.tz = tz
        new_dates.name = self.name
        return new_dates

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
            lib.tz_localize_array(self.asi8, self.tz)
        except:
            return False

        return True


def _generate_regular_range(start, end, periods, offset):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify two of start, end, or periods')

    if isinstance(offset, datetools.Tick):
        stride = offset.us_stride()
        if periods is None:
            b = Timestamp(start).value
            e = Timestamp(end).value + stride
        elif start is not None:
            b = Timestamp(start).value
            e = b + periods * stride
        else:
            e = Timestamp(end).value + stride
            b = e - periods * stride

        data = np.arange(b, e, stride, dtype=np.int64)
        data = data.view('M8[us]')
    else:
        xdr = datetools.generate_range(start=start, end=end,
            periods=periods, offset=offset)

        data = np.array(list(xdr), dtype='M8[us]')

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
    boxer = np.frompyfunc(boxfunc, 1, 1)
    return boxer(arr)


def _to_m8(key):
    '''
    Timestamp-like => dt64
    '''
    if not isinstance(key, datetime):
        # this also converts strings
        key = Timestamp(key)

    return np.datetime64(lib.pydt_to_i8(key))


def _to_m8_array(arr):
    if arr is None:
        return arr
    return np.frompyfunc(_to_m8, 1, 1)(arr)


def _str_to_dt_array(arr):
    def parser(x):
        result = parse_time_string(x)
        return result[0]

    p_ufunc = np.frompyfunc(parser, 1, 1)
    data = p_ufunc(arr)
    return np.array(data, dtype='M8[us]')


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
