# pylint: disable=E1101,E1103

from datetime import datetime
import operator

import numpy as np

from pandas.core.index import Index
import pandas.core.datetools as datetools

__all__ = ['DateRange']

#-------------------------------------------------------------------------------
# DateRange class

def _bin_op(op):
    def f(self, other):
        return op(self.view(np.ndarray), other)

    return f

_CACHE_START = datetime(1950, 1, 1)
_CACHE_END   = datetime(2030, 1, 1)

class DateRange(Index):
    """
    Fixed frequency date range according to input parameters.

    Input dates satisfy:
        begin <= d <= end, where d lies on the given offset

    Parameters
    ----------
    start : {datetime, None}
        left boundary for range
    end : {datetime, None}
        right boundary for range
    periods : int
        Number of periods to generate.
    offset : DateOffset, default is 1 BusinessDay
        Used to determine the dates returned
    time_rule : time_rule to use
    tzinfo : pytz.timezone
        To endow DateRange with time zone information
    """
    _cache = {}
    def __new__(cls, start=None, end=None, periods=None,
                offset=datetools.bday, time_rule=None,
                tzinfo=None, **kwds):

        time_rule = kwds.get('timeRule', time_rule)
        if time_rule is not None:
            offset = datetools.getOffset(time_rule)

        if time_rule is None:
            if offset in datetools._offsetNames:
                time_rule = datetools._offsetNames[offset]

        # Cachable
        if not start:
            start = kwds.get('begin')
        if not periods:
            periods = kwds.get('nPeriods')

        start = datetools.to_datetime(start)
        end = datetools.to_datetime(end)

        # inside cache range. Handle UTC case

        useCache = (offset.isAnchored() and
                    isinstance(offset, datetools.CacheableOffset))

        start, end, tzinfo = _figure_out_timezone(start, end, tzinfo)
        useCache = useCache and _naive_in_cache_range(start, end)

        if useCache:
            index = cls._cached_range(start, end, periods=periods,
                                      offset=offset, time_rule=time_rule)
            if tzinfo is None:
                return index
        else:
            xdr = generate_range(start=start, end=end, periods=periods,
                                 offset=offset, time_rule=time_rule)
            index = list(xdr)

        if tzinfo is not None:
            index = [d.replace(tzinfo=tzinfo) for d in index]

        index = np.array(index, dtype=object, copy=False)
        index = index.view(cls)
        index.offset = offset
        index.tzinfo = tzinfo
        return index

    def __reduce__(self):
        """Necessary for making this object picklable"""
        a, b, state = Index.__reduce__(self)
        aug_state = state, self.offset, self.tzinfo

        return a, b, aug_state

    def __setstate__(self, aug_state):
        """Necessary for making this object picklable"""
        index_state = aug_state[:1]
        offset = aug_state[1]

        # for backwards compatibility
        if len(aug_state) > 2:
            tzinfo = aug_state[2]
        else: # pragma: no cover
            tzinfo = None

        self.offset = offset
        self.tzinfo = tzinfo
        Index.__setstate__(self, *index_state)

    def equals(self, other):
        if self is other:
            return True

        if not isinstance(other, Index):
            return False

        return Index.equals(self.view(Index), other)

    def is_all_dates(self):
        return True

    @classmethod
    def _cached_range(cls, start=None, end=None, periods=None, offset=None,
                      time_rule=None):

        # HACK: fix this dependency later
        if time_rule is not None:
            offset = datetools.getOffset(time_rule)

        if offset is None:
            raise Exception('Must provide a DateOffset!')

        if offset not in cls._cache:
            xdr = generate_range(_CACHE_START, _CACHE_END, offset=offset)
            arr = np.array(list(xdr), dtype=object, copy=False)

            cachedRange = arr.view(DateRange)
            cachedRange.offset = offset
            cachedRange.tzinfo = None
            cls._cache[offset] = cachedRange
        else:
            cachedRange = cls._cache[offset]

        if start is None:
            if end is None:
                raise Exception('Must provide start or end date!')
            if periods is None:
                raise Exception('Must provide number of periods!')

            assert(isinstance(end, datetime))

            end = offset.rollback(end)

            endLoc = cachedRange.indexMap[end] + 1
            startLoc = endLoc - periods
        elif end is None:
            assert(isinstance(start, datetime))
            start = offset.rollforward(start)

            startLoc = cachedRange.indexMap[start]
            if periods is None:
                raise Exception('Must provide number of periods!')

            endLoc = startLoc + periods
        else:
            start = offset.rollforward(start)
            end = offset.rollback(end)

            startLoc = cachedRange.indexMap[start]
            endLoc = cachedRange.indexMap[end] + 1

        indexSlice = cachedRange[startLoc:endLoc]

        return indexSlice

    def __array_finalize__(self, obj):
        if self.ndim == 0: # pragma: no cover
            return self.item()

        self.offset = getattr(obj, 'offset', None)

    __lt__ = _bin_op(operator.lt)
    __le__ = _bin_op(operator.le)
    __gt__ = _bin_op(operator.gt)
    __ge__ = _bin_op(operator.ge)
    __eq__ = _bin_op(operator.eq)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        result = self.view(np.ndarray)[key]

        if isinstance(key, (int, np.integer)):
            return result
        elif isinstance(key, slice):
            new_index = result.view(DateRange)
            if key.step is not None:
                new_index.offset = key.step * self.offset
            else:
                new_index.offset = self.offset

            new_index.tzinfo = self.tzinfo
            return new_index
        else:
            return Index(result)

    def summary(self):
        if len(self) > 0:
            index_summary = ', %s to %s' % (self[0], self[-1])
        else:
            index_summary = ''
        sum_line = 'DateRange: %s entries%s' % (len(self), index_summary)
        sum_line += '\noffset: %s' % self.offset
        if self.tzinfo is not None:
            sum_line += ', tzinfo: %s' % self.tzinfo

        return sum_line

    def __repr__(self):
        output = str(self.__class__) + '\n'
        output += 'offset: %s, tzinfo: %s\n' % (self.offset, self.tzinfo)
        if len(self) > 0:
            output += '[%s, ..., %s]\n' % (self[0], self[-1])
        output += 'length: %d' % len(self)
        return output

    __str__ = __repr__

    def shift(self, n, offset=None):
        """
        Specialized shift which produces a DateRange

        Parameters
        ----------
        n : int
            Periods to shift by
        offset : DateOffset or timedelta-like, optional

        Returns
        -------
        shifted : DateRange
        """
        if offset is not None and offset != self.offset:
            return Index.shift(self, n, offset)

        if n == 0:
            # immutable so OK
            return self

        start = self[0] + n * self.offset
        end = self[-1] + n * self.offset
        return DateRange(start, end, offset=self.offset)

    def union(self, other):
        """
        Specialized union for DateRange objects. If combine
        overlapping ranges with the same DateOffset, will be much
        faster than Index.union

        Parameters
        ----------
        other : DateRange or array-like

        Returns
        -------
        y : Index or DateRange
        """
        if not isinstance(other, DateRange) or other.offset != self.offset:
            return Index.union(self.view(Index), other)

        offset = self.offset

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        left_start, left_end = left[0], left[-1]
        right_start, right_end = right[0], right[-1]

        # Only need to "adjoin", not overlap
        if (left_end + offset) >= right_start:
            return DateRange(left_start, max(left_end, right_end),
                             offset=offset)
        else:
            return Index.union(self, other)

    def tz_normalize(self, tz):
        """
        Convert DateRange from one time zone to another (using pytz)

        Returns
        -------
        normalized : DateRange
        """
        new_dates = np.array([tz.normalize(x) for x in self])
        new_dates = new_dates.view(DateRange)
        new_dates.offset = self.offset
        new_dates.tzinfo = tz
        return new_dates

    def tz_localize(self, tz):
        """
        Localize tzinfo-naive DateRange to given time zone (using pytz)

        Returns
        -------
        localized : DateRange
        """
        new_dates = np.array([tz.localize(x) for x in self])
        new_dates = new_dates.view(DateRange)
        new_dates.offset = self.offset
        new_dates.tzinfo = tz
        return new_dates

    def tz_validate(self):
        """
        For a localized time zone, verify that there are no DST ambiguities

        Returns
        -------
        result : boolean
            True if there are no DST ambiguities
        """
        import pytz

        tz = self.tzinfo
        if tz is None or tz is pytz.utc:
            return True

        # See if there are any DST resolution problems
        for date in self:
            try:
                tz.utcoffset(date.replace(tzinfo=None))
            except pytz.InvalidTimeError:
                return False

        return True

def generate_range(start=None, end=None, periods=None,
                   offset=datetools.BDay(), time_rule=None):
    """
    Generates a sequence of dates corresponding to the specified time
    offset. Similar to dateutil.rrule except uses pandas DateOffset
    objects to represent time increments

    Parameters
    ----------
    start : datetime (default None)
    end : datetime (default None)
    periods : int, optional

    Note
    ----
    * This method is faster for generating weekdays than dateutil.rrule
    * At least two of (start, end, periods) must be specified.
    * If both start and end are specified, the returned dates will
    satisfy start <= date <= end.

    Returns
    -------
    dates : generator object

    See also
    --------
    DateRange, dateutil.rrule
    """

    if time_rule is not None:
        offset = datetools.getOffset(time_rule)

    if time_rule is None:
        if offset in datetools._offsetNames:
            time_rule = datetools._offsetNames[offset]

    start = datetools.to_datetime(start)
    end = datetools.to_datetime(end)

    if start and not offset.onOffset(start):
        start = offset.rollforward(start)

    if end and not offset.onOffset(end):
        end = offset.rollback(end)

        if periods is None and end < start:
            end = None
            periods = 0

    if end is None:
        end = start + (periods - 1) * offset

    if start is None:
        start = end - (periods - 1) * offset

    cur = start
    if offset._normalizeFirst:
        cur = datetools.normalize_date(cur)

    while cur <= end:
        yield cur

        # faster than cur + offset
        cur = offset.apply(cur)

# Do I want to cache UTC dates? Can't decide...

# def _utc_in_cache_range(start, end):
#     import pytz
#     if start is None or end is None:
#         return False

#     _CACHE_START = datetime(1950, 1, 1, tzinfo=pytz.utc)
#     _CACHE_END   = datetime(2030, 1, 1, tzinfo=pytz.utc)

#     try:
#         assert(_isutc(start))
#         assert(_isutc(end))
#     except AssertionError:
#         raise Exception('To use localized time zone, create '
#                         'DateRange with pytz.UTC then call '
#                         'tz_normalize')
#     return _in_range(start, end, _CACHE_START, _CACHE_END)

# def _isutc(dt):
#     import pytz
#     return dt.tzinfo is pytz.utc

# def _hastz(dt):
#     return dt is not None and dt.tzinfo is not None

# def _have_pytz():
#     try:
#         import pytz
#         return True
#     except ImportError:
#         return False

def _in_range(start, end, rng_start, rng_end):
    return start > rng_start and end < rng_end

def _naive_in_cache_range(start, end):
    if start is None or end is None:
        return False
    else:
        return _in_range(start, end, _CACHE_START, _CACHE_END)

def _figure_out_timezone(start, end, tzinfo):
    inferred_tz = _infer_tzinfo(start, end)
    tz = inferred_tz
    if inferred_tz is None and tzinfo is not None:
        tz = tzinfo
    elif tzinfo is not None:
        assert(inferred_tz == tzinfo)
        # make tz naive for now

    start = start if start is None else start.replace(tzinfo=None)
    end = end if end is None else end.replace(tzinfo=None)

    return start, end, tz

def _infer_tzinfo(start, end):
    def _infer(a, b):
        tz = a.tzinfo
        if b and b.tzinfo:
            assert(tz == b.tzinfo)
        return tz
    tz = None
    if start is not None:
        tz = _infer(start, end)
    elif end is not None:
        tz = _infer(end, start)
    return tz

if __name__ == '__main__':
    import pytz
    # just want it to work
    tz = pytz.timezone('US/Eastern')
    dr = DateRange(datetime(2011, 3, 12, tzinfo=pytz.utc),
                   periods=50, offset=datetools.Hour())
    dr2 = dr.tz_normalize(tz)
