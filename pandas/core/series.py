"""
Data structure for 1-dimensional cross-sectional and time series data
"""

# pylint: disable-msg=E1101
# pylint: disable-msg=E1103

from datetime import datetime
from itertools import izip

from numpy import array, NaN, ndarray
import numpy as np

from pandas.core.daterange import DateRange
from pandas.core.index import Index, NULL_INDEX
from pandas.core.mixins import Picklable, Groupable
import pandas.core.datetools as datetools

from pandas.lib.tseries import isnull, notnull
import pandas.lib.tseries as tseries

#-------------------------------------------------------------------------------
# Wrapper function for Series arithmetic methods

def _seriesOpWrap(opname, comp=False):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    MIRROR_OPS = {
        '__add__' : '__radd__',
        '__sub__' : '__rsub__',
        '__div__' : '__rdiv__',
        '__mul__' : '__rmul__'
    }
    def wrapper(self, other):
        from pandas.core.frame import DataFrame

        func = getattr(self.view(ndarray), opname)
        cls = self.__class__
        if isinstance(other, Series):
            if self.index.equals(other.index):
                return cls(func(other.view(ndarray)), index=self.index)
            if len(self.index) + len(other.index) > 0:
                newIndex = self.index + other.index
            else:
                newIndex = NULL_INDEX
            try:
                arr = tseries.combineFunc(opname, newIndex, self, other,
                                          self.index.indexMap,
                                          other.index.indexMap)
            except Exception, e:
                arr = Series.combineFunc(self, other,
                                         getattr(type(self[0]), opname))
            result = cls(arr, index=newIndex)
            if comp:
                result[isnull(result)] = 0
                return result.astype(np.bool)
            else:
                return result
        elif isinstance(other, DataFrame):
            reverse_op = MIRROR_OPS.get(opname)
            if reverse_op is None:
                raise Exception('Cannot do %s op, sorry!')
            return getattr(other, reverse_op)(self)
        else:
            return cls(func(other), index=self.index)
    return wrapper

#-------------------------------------------------------------------------------
# Series class

class Series(np.ndarray, Picklable, Groupable):
    """Generic indexed series (time series or otherwise) object.

    Contains values in a numpy-ndarray with an optional bound index
    (also an array of dates, strings, or whatever you want the 'row
    names' of your series to be)

    Rows can be retrieved by index value (date, string, etc.) or
    relative position in the underlying array.

    Operations between Series (+, -, /, *, **) objects are
    *index-safe*, meaning that values will be combined by their
    respective index positions rather than relative positions in the
    underlying ndarray. In other words, there is no 'matching' or
    'aligning' to do, it's all taken care of for you.

    Note
    ----
    If you combine two series, all values for an index position must
    be present or the value for that index position will be nan. The
    new index is the sorted union of the two Series indices.

    ALSO NOTE: There is currently no restriction on what can be in the
    index.

    Parameters
    ----------
    data:  array-like
        Underlying values of Series, preferably as numpy ndarray
    index: array-like, optional
        Index object (or other iterable of same length as data)


    Example usage:
        >>> s = Series(arr, index=Index(dates))
        >>> t = Series(otherArr, index=Index(otherDates))
        >>> s / t # --> new Series resulting from by-index division of elements
        >>> d = s.index[5]
        >>> s[5]
        >>> s[d]    # Valid
    """
    def __new__(cls, data, index=None, dtype=None, copy=False):
        if index is None and isinstance(data, Series):
            index = data.index

        if index is None:
            raise Exception('Index cannot be None!')

        # Make a copy of the data, infer type
        subarr = array(data, dtype=dtype, copy=copy)

        if subarr.ndim == 0:
            return subarr.item()

        # This is to prevent mixed-type Series getting all casted to
        # NumPy string type, e.g. NaN --> '-1#IND'.

        if issubclass(subarr.dtype.type, basestring):
            subarr = array(data, dtype=object, copy=copy)

        # Change the class of the array to be the subclass type.
        subarr = subarr.view(cls)
        subarr.index = index

        if subarr.index._allDates:
            subarr = subarr.view(TimeSeries)

        return subarr

    def _get_index(self):
        return self._index

    def _set_index(self, index):
        indexTypes = ndarray, Index, list, tuple
        if not isinstance(index, indexTypes):
            raise TypeError("Expected index to be in %s; was %s."
                            % (indexTypes, type(index)))

        if len(self) != len(index):
            raise AssertionError('Lengths of index and values did not match!')

        if not isinstance(index, Index):
            index = Index(index)

        self._index = index

    _index = None
    index = property(fget=_get_index, fset=_set_index)

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)

    @classmethod
    def fromDict(cls, input={}, castFloat=True, **kwds):
        """
        Analogous to asDataFrame, but turns dict into Series

        Parameters
        ----------
        input: dict object
            Keys become indices of returned Series
        kwds: optionally provide arguments as keywords

        Returns
        -------
        Series
        """
        input = input.copy()
        input.update(kwds)

        index = Index(sorted(input.keys()))
        if castFloat:
            try:
                useData = [float(input[idx]) for idx in index]
            except Exception, e:
                useData = [input[idx] for idx in index]
        else:
            useData = [input[idx] for idx in index]
        return Series(useData, index=index)

    @classmethod
    def fromValue(cls, value=np.NaN, index=None, dtype=None):
        """
        Analogous to asDataFrame, but turns dict into Series

        Parameters
        ----------
        input: dict object
            Keys become indices of returned Series
        kwds: optionally provide arguments as keywords

        Returns
        -------
        Series
        """
        if not isinstance(index, Index):
            index = Index(index)

        # If we create an empty array using a string to infer
        # the dtype, NumPy will only allocate one character per entry
        # so this is kind of bad. Alternately we could use np.repeat
        # instead of np.empty (but then you still don't want things
        # coming out as np.str_!
        if isinstance(value, basestring):
            dtype = np.object_

        if dtype is None:
            arr = np.empty(len(index), dtype=type(value))
        else:
            arr = np.empty(len(index), dtype=dtype)
        arr.fill(value)

        return Series(arr, index=index)

    @classmethod
    def load(cls, baseFile):
        """
        Load Series from file.

        Parameters
        ----------
        baseFile: string
            Filename base where index/values are stored.
            e.g. baseFile='myfile' --> 'myfile_index.npy', 'myfile_values.npy'

        Returns
        -------
        Series or TimeSeries
        """
        indexFile = baseFile + '_index.npy'
        valuesFile = baseFile + '_values.npy'
        index = np.load(indexFile)
        values = np.load(valuesFile)

        return cls(values, index=index)

    def save(self, baseFile):
        """
        Save Series to file.

        Parameters
        ----------
        baseFile: string
            Filename base where index/values are stored.
            e.g. baseFile='myfile' --> 'myfile_index.npy', 'myfile_values.npy'
        """
        indexFile = baseFile + '_index'
        valuesFile = baseFile + '_values'

        np.save(indexFile, self.index)
        np.save(valuesFile, self)

    def __contains__(self, key):
        return key in self.index

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))
        subclass_state = (self.index, )
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)
        index, = own_state
        self.index = index

    def __getitem__(self, key):
        """
        Returns item(s) for requested index/sequence, overrides default behavior
        for series[key].

        Logic is as follows:
            - If key is in the index, return the value corresponding
              to that index
            - Otherwise, use key (presumably one integer or a sequence
              of integers) to obtain values from the series. In the case
              of a sequence, a 'slice' of the series (with corresponding dates)
              will be returned, otherwise a single value.
        """
        if key is None and key not in self.index:
            raise Exception('None/Null object requested of Series!')
        if not hasattr(key, '__iter__'):
            try:
                # Check that we can even look for this in the index
                return ndarray.__getitem__(self, self.index.indexMap[key])
            except KeyError:
                if isinstance(key, int):
                    return ndarray.__getitem__(self, key)
                raise Exception('Requested index not in this series!')
            except TypeError:
                # Could not hash item
                pass
        dataSlice = self.view(ndarray)[key]
        if self.index is not None:
            indices = self.index.view(ndarray)[key]
            if isinstance(indices, ndarray):
                indexSlice = Index(indices)
                return self.__class__(dataSlice, index=indexSlice)
            else:
                return dataSlice
        else:
            if isinstance(dataSlice, ndarray):
                return self.__class__(dataSlice)
            else:
                return dataSlice    # Just one element

    def __getslice__(self, i, j):
        """
        Returns a slice of the Series.

        Note that the underlying values are COPIES.

        The reason that the getslice returns copies is that otherwise you
        will have a reference to the original series which could be
        inadvertently changed if the slice were altered (made mutable).
        """
        newArr = self.view(ndarray)[i:j].copy()

        if self.index is not None:
            newIndex = self.index[i:j]
            return self.__class__(newArr, index = newIndex)
        else:
            return self.__class__(newArr)

    def __setitem__(self, key, value):
        """
        If this series is mutable, set specified indices equal to given values.
        """
        try:
            loc = self.index.indexMap[key]
            ndarray.__setitem__(self, loc, value)
        except Exception, e:
            ndarray.__setitem__(self, key, value)

    def __setslice__(self, i, j, value):
        """Set slice equal to given value(s)"""
        ndarray.__setslice__(self, i, j, value)

    def __repr__(self):
        """Clean string representation of a Series"""
        vals = self.view(ndarray)
        index = self.index
        if index is not None and len(index) > 0:
            if len(index) > 500:
                head = _seriesRepr(index[:50], vals[:50])
                tail = _seriesRepr(index[-50:], vals[-50:])
                return head + '\n...\n' + tail + '\nlength: %d' % len(vals)
            else:
                return _seriesRepr(index, vals)
        else:
            return 'No index!\n' + ndarray.__repr__(self)

    def __str__(self):
        return self.__repr__()

    def __iter__(self):
        return iter(self.view(ndarray))

#-------------------------------------------------------------------------------
#   Arithmetic operators

    __add__ = _seriesOpWrap('__add__')
    __sub__ = _seriesOpWrap('__sub__')
    __mul__ = _seriesOpWrap('__mul__')
    __div__ = _seriesOpWrap('__div__')
    __pow__ = _seriesOpWrap('__pow__')

    # Inplace operators
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __idiv__ = __div__
    __ipow__ = __pow__

#-------------------------------------------------------------------------------
# Overridden ndarray methods

    def sum(self, axis=None, dtype=None, out=None):
        """
        Overridden version of ndarray.sum for Series which excludes
        NaN automatically
        """
        arr = self.view(ndarray)
        retVal = arr.sum(axis, dtype, out)

        if isnull(retVal):
            arr = remove_na(arr)
            retVal = arr.sum(axis, dtype, out)

        return retVal

    def mean(self, axis=None, dtype=None, out=None):
        """
        Overridden version of ndarray.mean for Series which excludes
        NaN automatically
        """
        arr = self.view(ndarray)
        retVal = arr.mean(axis, dtype, out)

        if isnull(retVal):
            arr = remove_na(arr)
            retVal = arr.mean(axis, dtype, out)

        return retVal

    def min(self, axis=None, out=None):
        """
        Overridden version of ndarray.min for Series which excludes
        NaN automatically
        """
        arr = self.view(ndarray)
        retVal = arr.min(axis, out)

        if isnull(retVal):
            arr = remove_na(arr)
            retVal = arr.min(axis, out)

        return retVal

    def max(self, axis=None, out=None):
        """
        Overridden version of ndarray.max for Series which excludes
        NaN automatically
        """
        arr = self.view(ndarray)
        retVal = arr.max(axis, out)

        if isnull(retVal):
            arr = remove_na(arr)
            retVal = arr.max(axis, out)

        return retVal

    def std(self, axis=None, dtype=None, out=None, ddof=1):
        """
        Overridden version of ndarray.std for Series which excludes
        NaN automatically
        """
        nona = remove_na(self.view(ndarray))
        if len(nona) < 2:
            return NaN
        return ndarray.std(nona, axis, dtype, out, ddof)

    def var(self, axis=None, dtype=None, out=None, ddof=1):
        """
        Overridden version of ndarray.var for Series which excludes
        NaN automatically
        """
        nona = remove_na(self.view(ndarray))
        if len(nona) < 2:
            return NaN
        return ndarray.var(nona, axis, dtype, out, ddof)

    def skew(self, bias=False):
        """Computes the skewness of the Series

        For normally distributed data, the skewness should be about 0.
        A skewness value > 0 means that there is more weight in the left
        tail of the distribution. The function skewtest() can be used to
        determine if the skewness value is close enough to 0, statistically
        speaking.

        Parameters
        ----------
        bias : bool
        If False, then the calculations are corrected for statistical bias.
        """

        from scipy.stats import skew
        nona = remove_na(self.view(ndarray))
        if len(nona) < 3:
            return NaN
        theSkew = skew(nona, bias=bias)

        if isinstance(theSkew, ndarray):
            theSkew = theSkew.item()

        return theSkew

    def keys(self):
        """
        Return Series index

        Analogous to dict.keys()
        """
        return self.index

    def values(self):
        """
        Return Series values

        Analogous to dict.values()
        """
        return self.view(ndarray)

    def iteritems(self):
        """
        Iterate over (index, value) tuples
        """
        if self.index is not None:
            return izip(iter(self.index), iter(self))
        else:
            raise Exception('This series has no index!')

    def get(self, key, missingVal=None):
        """
        Returns value occupying requested index, and
        return missingVal if not in Index

        Parameters
        ----------
        key: object
            Index value looking for
        missingVal: object, optional
            Value to return if key not in index
        """
        if key in self.index:
            return ndarray.__getitem__(self, self.index.indexMap[key])
        else:
            return missingVal

    def append(self, other):
        """
        Concatenate two Series
        """
        newIndex = np.concatenate((self.index, other.index))
        newValues = np.concatenate((self, other))
        return self.__class__(newValues, index = newIndex)

    def merge(self, other):
        """
        If self is {A}->{B} and other is another mapping of {B}->{C}
        then returns a new Series that is {A}->{C}

        Parameters
        ----------
        other: dict or Series

        Returns
        -------
        Series having same index as calling instance, with values from
        input Series
        """
        if isinstance(other, dict):
            other = Series.fromDict(other)
        if not isinstance(other, Series):
            raise Exception('Argument must be a Series!')
        fillVec, mask = tseries.getMergeVec(self, other.index.indexMap)

        newValues = other.view(np.ndarray).take(fillVec)
        newValues[-mask] = np.nan

        newSer = Series(newValues, index=self.index)
        return newSer

    def combineFunc(self, other, func):
        """
        Combines this Series with another Series index by index using
        the given function.
        """
        if self.index.equals(other.index):
            newIndex = self.index
        else:
            newIndex = self.index + other.index
        newArr = np.empty(len(newIndex), dtype = self.dtype)
        for i, idx in enumerate(newIndex):
            newArr[i] = func(self.get(idx, NaN), other.get(idx, NaN))
        return self.__class__(newArr, index=newIndex)

    def combineFirst(self, other):
        """
        Combine Series values, choosing calling Series's values first.

        Parameters
        ----------
        other: Series

        Returns
        -------
        Series formed as union of
        """
        if self.index.equals(other.index):
            newIndex = self.index
        else:
            newIndex = self.index + other.index

        this = self.reindex(newIndex)
        other = other.reindex(newIndex)
        result = Series(np.where(isnull(this), other, this), index=newIndex)

        return result

    def argsort(self, axis = 0, kind='quicksort', order=None):
        """
        Overriding numpy's built-in cumsum functionality
        """
        arr = self.view(ndarray).copy()
        okLocs = notnull(arr)
        arr[okLocs] = np.argsort(arr[okLocs])
        return self.__class__(arr, index=self.index)

    def cumsum(self, axis = 0, dtype = None, out = None):
        """
        Overriding numpy's built-in cumsum functionality
        """
        arr = self.copy()
        okLocs = notnull(arr)
        result = np.cumsum(arr.view(ndarray)[okLocs])
        arr = arr.astype(result.dtype)
        arr[okLocs] = result
        return arr

    def cumprod(self, axis = 0, dtype = None, out = None):
        """
        Overriding numpy's built-in cumprod functionality
        """
        arr = self.copy()
        okLocs = notnull(arr)
        arr[okLocs] = np.cumprod(arr.view(ndarray)[okLocs])
        return arr

    def copy(self):
        return self.__class__(self.view(ndarray).copy(), index=self.index)

    def corr(self, other):
        """
        Correlation of this Series with another Series, NaN excluded

        Parameters
        ----------
        other: Series object

        Returns
        -------
        float (the correlation coefficient)
        """
        commonIdx = list(set(remove_na(self).index) &
                         set(remove_na(other).index))

        if len(commonIdx) == 0:
            return NaN

        this = self.reindex(commonIdx)
        that = other.reindex(commonIdx)

        return np.corrcoef(this, that)[0, 1]

    def count(self):
        """
        Return number of observations of Series.

        Returns
        -------
        int (# obs)
        """
        return np.isfinite(self.view(ndarray)).sum()

    def median(self):
        """
        Return median value of Series
        """
        selfExNaN = remove_na(self.view(ndarray))
        med = np.median(selfExNaN)
        return med

    def sort(self, axis=0, kind='quicksort', order=None):
        sortedSeries = self.order(missingAtEnd=True)
        self[:] = sortedSeries
        self.index = sortedSeries.index

    def order(self, missingAtEnd = True):
        """
        Sorts Series object, by value, maintaining index-value object

        Parameters
        ----------
        missingAtEnd: boolean (optional, default=True)
            Put NaN's at beginning or end

        In general, AVOID sorting Series unless you absolutely need to.

        Returns
        -------
        SORTED series by values (indices correspond to the appropriate values)
        """
        arr = self.view(ndarray)
        sortedIdx = np.empty(len(self), dtype=np.int32)

        bad = isnull(arr)

        good = -bad
        idx = np.arange(len(self))
        if missingAtEnd:
            n = sum(good)
            sortedIdx[:n] = idx[good][arr[good].argsort()]
            sortedIdx[n:] = idx[bad]
        else:
            n = sum(bad)
            sortedIdx[n:] = idx[good][arr[good].argsort()]
            sortedIdx[:n] = idx[bad]

        return Series(arr[sortedIdx], index=self.index[sortedIdx])

    def map(self, func):
        """
        Apply input Python function element-wise to each element of
        Series.

        Parameters
        ----------
        func: function
            Element-wise function to apply

        Returns
        -------
        Series with same index
        """
        return Series([func(x) for x in self], index = self.index)

    def plot(self, label=None, kind='line', **kwds):
        """
        Plot the input series with the index on the x-axis using
        matplotlib / pylab.

        Params
        ------
        label: label argument to provide to plot

        kind: {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series

        kwds: other plotting keyword arguments

        Default plot-types:
            TimeSeries: line chart
            Series: histogram
                Also support for bar charts

        Type show() (make sure to do 'from pylab import *') to see graph if you
        do not.

        Note
        ----
        See matplotlib documentation online for more on this subject
        """
        import pylab

        if label is not None:
            kwds = kwds.copy()
            kwds['label'] = label

        # I can't get this to work

        #fig = pylab.gcf()
        #fig.autofmt_xdate(bottom=0.1)

        #ax = fig.gca()
        #if not ax.has_data():
            #ax = fig.add_subplot(111)

        #ax.plot(self.index, self, **kwds)

        pylab.plot(self.index, self, **kwds)

    def unstack(self):
        """
        Inverse operator for *stack*
        """
        from pandas.core.frame import DataFrame
        data = {}
        for idx, value in self.iteritems():
            row, col = idx.split(';')
            try:
                row = datetime.fromordinal(int(row))
            except Exception, e:
                pass
            data.setdefault(row, {})[col] = value
        return DataFrame.fromDict(data)

    def toCSV(self, path=None):
        """
        Write the Series to a CSV file

        Parameters
        ----------
        path: string or None
            Output filepath. If None, write to stdout
        """
        if not path:
            import sys
            f = sys.stdout
        else:
            f = open(path, 'wb')
        for idx, value in self.iteritems():
            f.write(str(idx) + ',' + str(value) + ',\n')
        if path is not None:
            f.close()

    def toDict(self):
        return dict(self.iteritems())

    def cap(self, value):
        """Return copy of series with values above given value truncated"""
        myCopy = self.copy()
        myCopy[notnull(myCopy) & (myCopy > value)] = value
        return myCopy

    def floor(self, value):
        """Return copy of series with values BELOW given value truncated"""
        myCopy = self.copy()
        myCopy[notnull(myCopy) & (myCopy < value)] = value
        return myCopy

    def valid(self):
        """
        Return Series without NaN values

        Returns
        -------
        Series
        """
        return remove_na(self)

#-------------------------------------------------------------------------------
# TimeSeries methods

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataMatrix and Series objects within
        by given number (positive or negative) of business/weekdays.

        Parameters
        ----------
        periods: int (+ or -)
            Number of periods to move
        offset: DateOffset, optional
            Increment to use from datetools module
        timeRule: string
            time rule name to use by name (e.g. 'WEEKDAY')

        Returns
        -------
        TimeSeries
        """
        if periods == 0:
            return self

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        if offset is None:
            if periods > 0:
                newIndex = self.index[periods:]
                newValues = np.array(self)[:-periods]
            else:
                newIndex = self.index[:periods]
                newValues = np.array(self)[-periods:]
            return self.__class__(newValues, index=newIndex)
        else:
            offset = periods * offset
            newIndex = Index([idx + offset for idx in self.index])
            return self.__class__(self, index=newIndex)

    def slice(self, before, after):
        """

        """

        import bisect

        if before is not None:
            binsearch = bisect.bisect_left(self.index, before)
            cur = self.index[binsearch]
            next = self.index[min(binsearch + 1, len(self.index) - 1)]
            leftDate = cur if cur >= before else next
        else:
            leftDate = self.index[0]

        if after is not None:
            if after < self.index[-1]:
                binsearch = bisect.bisect_right(self.index, after)
                cur = self.index[binsearch]
                prior = self.index[max(binsearch - 1, 0)]
                rightDate = cur if cur <= after else prior
            else:
                rightDate = self.index[-1]
        else:
            rightDate = self.index[-1]

        beg_slice = max(0, self.index.indexMap[leftDate])
        end_slice = min(len(self.index), self.index.indexMap[rightDate] + 1)

        return self[beg_slice:end_slice]

    def asOf(self, date):
        """
        Return last good (non-NaN) value in TimeSeries if value is NaN for
        requested date.

        If there is no good value, NaN is returned.

        Returns
        -------
        value or NaN
        """
        if isinstance(date, basestring):
            date = datetools.to_datetime(date)

        v = self.get(date)

        if isnull(v):
            candidates = self.index[notnull(self)]
            candidates = candidates[candidates <= date]

            if any(candidates):
                asOfDate = max(candidates)
            else:
                return NaN

            return self.get(asOfDate)
        else:
            return v

    def fill(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Parameters
        ----------
        value: any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        method: {'backfill', 'pad', None}
            Method to use for filling holes in new inde

        Returns
        -------
        TimeSeries with NaN's filled

        See also
        --------
        reindex, asfreq
        """
        if value is not None:
            newSeries = self.copy()
            newSeries[isnull(newSeries)] = value
            return newSeries
        else: # Using reindex to pad / backfill
            withoutna = remove_na(self)
            return withoutna.reindex(self.index, fillMethod=method)

    def asfreq(self, freq, fillMethod=None):
        """
        Convert this TimeSeries to the provided frequency using DateOffset
        objects. Optionally provide fill method to pad/backfill/interpolate
        missing values.

        Parameters
        ----------
        offset: DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)

        fillMethod: {'backfill', 'pad', 'interpolate', None}
                    Method to use for filling holes in new inde

        Returns
        -------
        TimeSeries
        """
        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], timeRule=freq)

        return self.reindex(dateRange, fillMethod=fillMethod)

    def reindex(self, newIndex, fillMethod=None):
        """Overloaded version of reindex for TimeSeries. Supports filling
        with values based on new index.

        See analogous method for DataFrame, will be faster for multiple
        TimeSeries

        Parameters
        ----------
        newIndex:   array-like, preferably an Index object (to avoid
                    duplicating data)
        fillMethod: {'backfill', 'pad', 'interpolate', None}
                    Method to use for filling holes in reindexed Series

        Returns
        -------
        TimeSeries
        """
        if fillMethod is None:
            if self.index is newIndex:
                return self.copy()

            if not isinstance(newIndex, Index):
                newIndex = Index(newIndex)

            idxMap = self.index.indexMap

            if self.dtype == float:
                return self.__class__(tseries.reindex(newIndex, self, idxMap),
                                      index=newIndex)
            elif self.dtype == int:
                # This could be unsafe, but NaN will not work in int arrays.
                reindexed = tseries.reindex(newIndex, self.astype(float),
                                            idxMap)
                return self.__class__(reindexed, index=newIndex)

            else:
                if self.dtype.type == np.object_:
                    result = tseries.reindexObj(newIndex, self, idxMap)
                    return self.__class__(result, index=newIndex)
                else:
                    thisVals = self.view(np.ndarray).astype(object)
                    vals = tseries.reindexObj(newIndex, thisVals, idxMap)

                    if not isnull(vals).any():
                        vals = vals.astype(self.dtype)

                    return self.__class__(vals, index=newIndex)

        if not isinstance(newIndex, Index):
            newIndex = Index(newIndex)

        oldMap = self.index.indexMap
        newMap = newIndex.indexMap

        if not fillMethod:
            fillMethod = ''

        fillMethod = fillMethod.upper()

        if fillMethod not in ['BACKFILL', 'PAD', '']:
            raise Exception("Don't recognize fillMethod: %s" % fillMethod)

        # Cython for blazing speed
        fillVec, mask = tseries.getFillVec(self.index, newIndex, oldMap,
                                           newMap, kind=fillMethod)

        newValues = self.view(ndarray).take(fillVec)
        newValues[-mask] = NaN

        return self.__class__(newValues, index = newIndex)

    @property
    def weekday(self):
        return self.__class__([d.weekday() for d in self.index],
                              index = self.index)

    def truncate(self, before=None, after=None):
        """Function truncate a TimeSeries before and/or after some
        particular dates.

        Parameters
        ----------
        before: date
            Truncate before date
        after: date
            Truncate after date

        Note
        ----
        If TimeSeries is contained in a DataFrame, consider using the version
        of the function there.

        Returns
        -------
        TimeSeries
        """
        before = datetools.to_datetime(before)
        after = datetools.to_datetime(after)

        if before is None:
            before = min(self.index)
        if after is None:
            after = max(self.index)
        return self.slice(before, after)

    def diff(self):
        """
        1st discrete difference of object

        Returns
        -------
        TimeSeries
        """
        return (self - self.shift(1))

    def autocorr(self):
        """
        1st period autocorrelation coefficient

        Returns
        -------
        TimeSeries
        """
        return self.corr(self.shift(1))

    def _firstTimeWithValue(self):
        noNA = remove_na(self)

        if len(noNA) > 0:
            return noNA.index[0]
        else:
            return None

    def _lastTimeWithValue(self):
        noNA = remove_na(self)

        if len(noNA) > 0:
            return noNA.index[-1]
        else:
            return None


class TimeSeries(Series):
    pass

#-------------------------------------------------------------------------------
# Supplementary functions

def remove_na(arr):
    """
    Return array containing only true/non-NaN values, possibly empty.
    """
    return arr[notnull(arr)]

def _seriesRepr(index, vals):
    maxlen = max([len(str(idx)) for idx in index])
    padSpace = min(maxlen, 60)
    return '\n'.join([str(x).ljust(padSpace) + '\t' + str(v)
                      for x, v in izip(index, vals)])
