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
    timeRule : timeRule to use
    """
    _cache = {}
    def __new__(cls, start=None, end=None, periods=None,
                offset=datetools.bday, timeRule=None, **kwds):

        # Allow us to circumvent hitting the cache
        index = kwds.get('index')
        if index is None:
            if timeRule is not None:
                offset = datetools.getOffset(timeRule)

            if timeRule is None:
                if offset in datetools._offsetNames:
                    timeRule = datetools._offsetNames[offset]

            # Cachable
            if not start:
                start = kwds.get('begin')
            if not end:
                end = kwds.get('end')
            if not periods:
                periods = kwds.get('nPeriods')

            start = datetools.to_datetime(start)
            end = datetools.to_datetime(end)

            # inside cache range
            fromInside = start is not None and start > _CACHE_START
            toInside = end is not None and end < _CACHE_END

            useCache = fromInside and toInside

            if (useCache and offset.isAnchored() and
                isinstance(offset, datetools.CacheableOffset)):

                index = cls._cached_range(start, end, periods=periods,
                                          offset=offset, timeRule=timeRule)

            else:
                xdr = generate_range(start=start, end=end, periods=periods,
                                     offset=offset, timeRule=timeRule)

                index = np.array(list(xdr), dtype=object, copy=False)

                index = index.view(cls)
                index.offset = offset
        else:
            index = index.view(cls)

        return index

    def __reduce__(self):
        """Necessary for making this object picklable"""
        a, b, state = Index.__reduce__(self)
        aug_state = state, self.offset

        return a, b, aug_state

    def __setstate__(self, aug_state):
        """Necessary for making this object picklable"""
        state, offset = aug_state[:-1], aug_state[-1]

        self.offset = offset
        Index.__setstate__(self, *state)

    @property
    def _allDates(self):
        return True

    @classmethod
    def _cached_range(cls, start=None, end=None, periods=None, offset=None,
                       timeRule=None):

        # HACK: fix this dependency later
        if timeRule is not None:
            offset = datetools.getOffset(timeRule)

        if offset is None:
            raise Exception('Must provide a DateOffset!')

        if offset not in cls._cache:
            xdr = generate_range(_CACHE_START, _CACHE_END, offset=offset)
            arr = np.array(list(xdr), dtype=object, copy=False)

            cachedRange = DateRange.fromIndex(arr)
            cachedRange.offset = offset

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

    @classmethod
    def fromIndex(cls, index):
        index = cls(index=index)
        return index

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
            newIndex = result.view(DateRange)

            if key.step is not None:
                newIndex.offset = key.step * self.offset
            else:
                newIndex.offset = self.offset

            return newIndex
        else:
            return Index(result)

    def __repr__(self):
        output = str(self.__class__) + '\n'
        output += 'offset: %s\n' % self.offset
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

def generate_range(start=None, end=None, periods=None,
                   offset=datetools.BDay(), timeRule=None):
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

    if timeRule is not None:
        offset = datetools.getOffset(timeRule)

    if timeRule is None:
        if offset in datetools._offsetNames:
            timeRule = datetools._offsetNames[offset]

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
