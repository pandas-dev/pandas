"""
Data structure for 1-dimensional cross-sectional and time series data
"""

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622

import itertools
import sys
import warnings

from numpy import NaN, ndarray
import numpy as np

from pandas.core.common import isnull, notnull
from pandas.core.daterange import DateRange
from pandas.core.index import Index, NULL_INDEX
from pandas.core.mixins import Picklable, Groupable
import pandas.core.datetools as datetools
import pandas.lib.tseries as tseries

__all__ = ['Series', 'TimeSeries']

#-------------------------------------------------------------------------------
# Wrapper function for Series arithmetic methods

def _seriesOpWrap(opname):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    MIRROR_OPS = {
        '__add__' : '__radd__',
        '__sub__' : '__rsub__',
        '__div__' : '__rdiv__',
        '__mul__' : '__rmul__',
    }
    def wrapper(self, other):
        from pandas.core.frame import DataFrame

        func = getattr(self.values, opname)
        if isinstance(other, Series):
            if self.index.equals(other.index):
                return Series(func(other.values), index=self.index)

            newIndex = self.index + other.index

            try:
                if self.dtype != np.float_:
                    this = self.astype(float)
                else:
                    this = self

                if other.dtype != np.float_:
                    other = other.astype(float)

                # buffered Cython function expects double type

                arr = tseries.combineFunc(opname, newIndex,
                                          this, other,
                                          self.index.indexMap,
                                          other.index.indexMap)
            except Exception:
                arr = Series.combineFunc(self, other,
                                         getattr(type(self[0]), opname))
            result = Series(arr, index=newIndex)
            return result

        elif isinstance(other, DataFrame):
            reverse_op = MIRROR_OPS.get(opname)

            if reverse_op is None:
                raise Exception('Cannot do %s op, sorry!')

            return getattr(other, reverse_op)(self)
        else:
            return Series(func(other), index=self.index)
    return wrapper

#-------------------------------------------------------------------------------
# Series class

class Series(np.ndarray, Picklable, Groupable):
    """
    Generic indexed (labeled) vector (time series or cross-section)

    Contains values in a numpy-ndarray with an optional bound index
    (also an array of dates, strings, or whatever you want the 'row
    names' of your series to be)

    Rows can be retrieved by index value (date, string, etc.) or
    relative position in the underlying array.

    Operations between Series (+, -, /, *, **) align values based on
    their associated index values-- they need not be the same length.

    Parameters
    ----------
    data : array-like, dict, or scalar value
        Contains data stored in Series
    index : array-like
        Index object (or other iterable of same length as data)
        Must be input if first argument is not a dict. If both a dict
        and index sequence are used, the index will override the keys
        found in the dict.
    dtype : numpy.dtype or None
        If None, dtype will be inferred
    copy : boolean, default False
        Copy input data

    Notes
    -----
    If you combine two series, all values for an index position must
    be present or the value for that index position will be nan. The
    new index is the sorted union of the two Series indices.

    Data is *not* copied from input arrays by default
    """
    def __new__(cls, data, index=None, dtype=None, copy=False):
        if isinstance(data, Series):
            if index is None:
                index = data.index
        elif isinstance(data, dict):
            if index is None:
                index = Index(sorted(data.keys()))
            data = [data[idx] for idx in index]

        # Create array, do *not* copy data by default, infer type
        try:
            subarr = np.array(data, dtype=dtype, copy=copy)
        except ValueError:
            if dtype:
                raise

            subarr = np.array(data, dtype=object)

        if subarr.ndim == 0:
            if isinstance(data, list): # pragma: no cover
                subarr = np.array(data, dtype=object)
            elif index is not None:
                value = data

                # If we create an empty array using a string to infer
                # the dtype, NumPy will only allocate one character per entry
                # so this is kind of bad. Alternately we could use np.repeat
                # instead of np.empty (but then you still don't want things
                # coming out as np.str_!
                if isinstance(value, basestring) and dtype is None:
                    dtype = np.object_

                if dtype is None:
                    subarr = np.empty(len(index), dtype=type(value))
                else:
                    subarr = np.empty(len(index), dtype=dtype)
                subarr.fill(value)
            else:
                return subarr.item()

        elif subarr.ndim > 1:
            raise Exception('Data must be 1-dimensional')

        if index is None:
            raise Exception('Index cannot be None!')

        # This is to prevent mixed-type Series getting all casted to
        # NumPy string type, e.g. NaN --> '-1#IND'.
        if issubclass(subarr.dtype.type, basestring):
            subarr = np.array(data, dtype=object, copy=copy)

        # Change the class of the array to be the subclass type.
        subarr = subarr.view(cls)
        subarr.index = index

        if subarr.index._allDates:
            subarr = subarr.view(TimeSeries)

        return subarr

    def __hash__(self):
        raise TypeError('unhashable type')

    _index = None
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

    index = property(fget=_get_index, fset=_set_index)

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)

    def toDict(self):
        return dict(self.iteritems())

    @classmethod
    def fromValue(cls, value=np.NaN, index=None, dtype=None): # pragma: no cover
        warnings.warn("'fromValue', can call Series(value, index=index) now",
                      FutureWarning)

        return Series(value, index=index, dtype=dtype)

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
        values = self.values

        try:
            # Check that we can even look for this in the index
            return values[self.index.indexMap[key]]
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return values[key]
            raise Exception('Requested index not in this series!')
        except TypeError:
            # Could not hash item
            pass

        # is there a case where this would NOT be an ndarray?
        # need to find an example, I took out the case for now

        dataSlice = values[key]
        indices = Index(self.index.view(ndarray)[key])
        return Series(dataSlice, index=indices)

    def get(self, key, default=None):
        """
        Returns value occupying requested index, default to specified
        missing value if not present

        Parameters
        ----------
        key : object
            Index value looking for
        default : object, optional
            Value to return if key not in index

        Returns
        -------
        y : scalar
        """
        if key in self.index:
            return ndarray.__getitem__(self, self.index.indexMap[key])
        else:
            return default

    def __getslice__(self, i, j):
        """
        Returns a slice of the Series.

        Note that the underlying values are COPIES.

        The reason that the getslice returns copies is that otherwise you
        will have a reference to the original series which could be
        inadvertently changed if the slice were altered (made mutable).
        """
        newArr = self.values[i:j].copy()
        newIndex = self.index[i:j]

        return Series(newArr, index=newIndex)

    def __setitem__(self, key, value):
        """
        If this series is mutable, set specified indices equal to given values.
        """
        try:
            loc = self.index.indexMap[key]
            ndarray.__setitem__(self, loc, value)
        except Exception:
            values = self.values
            values[key] = value

    def __setslice__(self, i, j, value):
        """Set slice equal to given value(s)"""
        ndarray.__setslice__(self, i, j, value)

    def __repr__(self):
        """Clean string representation of a Series"""
        vals = self.values
        index = self.index

        if len(index) > 500:
            head = _seriesRepr(index[:50], vals[:50])
            tail = _seriesRepr(index[-50:], vals[-50:])
            return head + '\n...\n' + tail + '\nlength: %d' % len(vals)
        elif len(index) > 0:
            return _seriesRepr(index, vals)
        else:
            return '%s' % ndarray.__repr__(self)

    def toString(self, buffer=sys.stdout, nanRep='NaN'):
        print >> buffer, _seriesRepr(self.index, self.values,
                                     nanRep=nanRep)

    def __str__(self):
        return repr(self)

    def __iter__(self):
        return iter(self.values)

    def copy(self):
        return Series(self.values.copy(), index=self.index)

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
# Statistics, overridden ndarray methods

    def count(self):
        """
        Return number of observations of Series.

        Returns
        -------
        nobs : int
        """
        return notnull(self.values).sum()

    def sum(self, axis=None, dtype=None, out=None):
        """
        Sum of non-null values
        """
        return self._ndarray_statistic('sum')

    def mean(self, axis=None, dtype=None, out=None):
        """
        Mean of non-null values
        """
        return self._ndarray_statistic('mean')

    def _ndarray_statistic(self, funcname):
        arr = self.values
        retVal = getattr(arr, funcname)()

        if isnull(retVal):
            arr = remove_na(arr)
            retVal = getattr(arr, funcname)()

        return retVal

    def min(self, axis=None, out=None):
        """
        Minimum of non-null values
        """
        arr = self.values.copy()
        if not issubclass(arr.dtype.type, np.int_):
            arr[isnull(arr)] = np.inf
        return arr.min()

    def max(self, axis=None, out=None):
        """
        Maximum of non-null values
        """
        arr = self.values.copy()
        if not issubclass(arr.dtype.type, np.int_):
            arr[isnull(arr)] = -np.inf
        return arr.max()

    def std(self, axis=None, dtype=None, out=None, ddof=1):
        """
        Unbiased standard deviation of non-null values
        """
        nona = remove_na(self.values)
        if len(nona) < 2:
            return NaN
        return ndarray.std(nona, axis, dtype, out, ddof)

    def var(self, axis=None, dtype=None, out=None, ddof=1):
        """
        Unbiased variance of non-null values
        """
        nona = remove_na(self.values)
        if len(nona) < 2:
            return NaN
        return ndarray.var(nona, axis, dtype, out, ddof)

    def skew(self):
        """
        Unbiased skewness of the non-null values

        Returns
        -------
        skew : float
        """
        y = np.array(self.values)
        mask = notnull(y)
        count = mask.sum()
        np.putmask(y, -mask, 0)

        A = y.sum() / count
        B = (y**2).sum() / count  - A**2
        C = (y**3).sum() / count - A**3 - 3*A*B

        return (np.sqrt((count**2-count))*C) / ((count-2)*np.sqrt(B)**3)

    def cumsum(self, axis=0, dtype=None, out=None):
        """
        Cumulative sum of values. Preserves NaN values

        Extra parameters are to preserve ndarray interface.

        Returns
        -------

        """
        arr = self.copy()
        okLocs = notnull(arr)
        result = np.cumsum(arr.view(ndarray)[okLocs])
        arr = arr.astype(result.dtype)
        arr[okLocs] = result
        return arr

    def cumprod(self, axis=0, dtype=None, out=None):
        """
        Overriding numpy's built-in cumprod functionality
        """
        arr = self.copy()
        okLocs = notnull(arr)
        arr[okLocs] = np.cumprod(arr.view(ndarray)[okLocs])
        return arr

    def median(self):
        """
        Compute median value of non-null values
        """
        arr = self.values

        if arr.dtype != np.float_:
            arr = arr.astype(float)

        arr = arr[notnull(arr)]
        return tseries.median(arr)

    def corr(self, other):
        """
        Compute correlation two Series, excluding missing values

        Parameters
        ----------
        other : Series object

        Returns
        -------
        correlation : float
        """
        commonIdx = self.valid().index.intersection(other.valid().index)

        if len(commonIdx) == 0:
            return NaN

        this = self.reindex(commonIdx)
        that = other.reindex(commonIdx)

        return np.corrcoef(this, that)[0, 1]

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
        Lag-1 autocorrelation

        Returns
        -------
        TimeSeries
        """
        return self.corr(self.shift(1))

    def clip(self, upper=None, lower=None):
        """
        Trim values at input threshold(s)

        Parameters
        ----------
        lower : float, default None
        upper : float, default None

        Returns
        -------
        y : Series
        """
        result = self
        if lower is not None:
            result = result.clip_lower(lower)
        if upper is not None:
            result = result.clip_upper(upper)

        return result

    def clip_upper(self, threshold):
        """Return copy of series with values above given value truncated"""
        return np.where(self > threshold, threshold, self)

    def clip_lower(self, threshold):
        """Return copy of series with values below given value truncated"""
        return np.where(self < threshold, threshold, self)

#-------------------------------------------------------------------------------
# Iteration

    def keys(self):
        "Alias for Series index"
        return self.index

    @property
    def values(self):
        """
        Return Series as ndarray

        Returns
        -------
        arr : numpy.ndarray
        """
        return self.view(ndarray)

    def iteritems(self):
        """
        Lazily iterate over (index, value) tuples
        """
        return itertools.izip(iter(self.index), iter(self))

#-------------------------------------------------------------------------------
# Combination

    def append(self, other):
        """
        Concatenate two Series. The indices should not overlap

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
        """
        newIndex = np.concatenate((self.index, other.index))

        # Force overlap check
        try:
            newIndex = Index(newIndex)
        except Exception:
            raise

        newValues = np.concatenate((self, other))
        return Series(newValues, index=newIndex)

    def combineFunc(self, other, func):
        """
        Combines this Series using the given function with either
          * another Series index by index
          * a scalar value
          * DataFrame

        Parameters
        ----------
        other : {Series, DataFrame, scalar value}

        Returns
        -------
        y : {Series or DataFrame}
            Output depends on input. If a DataFrame is inputted, that
            will be the return type.
        """
        if isinstance(other, Series):
            newIndex = self.index + other.index

            newArr = np.empty(len(newIndex), dtype=self.dtype)
            for i, idx in enumerate(newIndex):
                newArr[i] = func(self.get(idx, NaN), other.get(idx, NaN))
        else:
            newIndex = self.index
            newArr = func(self.values, other)

        return Series(newArr, index=newIndex)

    def combineFirst(self, other):
        """
        Combine Series values, choosing calling Series's values first.

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
            formed as union of two Series
        """
        if self.index.equals(other.index):
            newIndex = self.index
            # save ourselves the copying in this case
            this = self
        else:
            newIndex = self.index + other.index

            this = self.reindex(newIndex)
            other = other.reindex(newIndex)

        result = Series(np.where(isnull(this), other, this), index=newIndex)
        return result

#-------------------------------------------------------------------------------
# Reindexing, sorting

    def sort(self, axis=0, kind='quicksort', order=None):
        """
        Overridden NumPy sort, taking care with missing values
        """
        sortedSeries = self.order(missingAtEnd=True)
        self[:] = sortedSeries
        self.index = sortedSeries.index

    def argsort(self, axis=0, kind='quicksort', order=None):
        """
        Overriding numpy's built-in cumsum functionality
        """
        values = self.values
        mask = isnull(values)

        if mask.any():
            result = values.copy()
            notmask = -mask
            result[notmask] = np.argsort(values[notmask])
            return Series(result, index=self.index)
        else:
            return Series(np.argsort(values), index=self.index)

    def order(self, missingAtEnd=True):
        """
        Sorts Series object, by value, maintaining index-value object

        Parameters
        ----------
        missingAtEnd : boolean (optional, default=True)
            Put NaN's at beginning or end

        In general, AVOID sorting Series unless you absolutely need to.

        Returns
        -------
        y : Series
            sorted by values
        """
        def _try_mergesort(arr):
            # easier to ask forgiveness than permission
            try:
                return arr.argsort(kind='mergesort')
            except TypeError:
                # stable sort not available for object dtype
                return arr.argsort()

        arr = self.values
        sortedIdx = np.empty(len(self), dtype=np.int32)

        bad = isnull(arr)

        good = -bad
        idx = np.arange(len(self))
        if missingAtEnd:
            n = sum(good)
            sortedIdx[:n] = idx[good][_try_mergesort(arr[good])]
            sortedIdx[n:] = idx[bad]
        else:
            n = sum(bad)
            sortedIdx[n:] = idx[good][_try_mergesort(arr[good])]
            sortedIdx[:n] = idx[bad]

        return Series(arr[sortedIdx], index=self.index[sortedIdx])

    def map(self, arg):
        """
        Map values of Series using input correspondence (which can be
        a dict, Series, or function).

        Parameters
        ----------
        arg : function, dict, or Series

        Returns
        -------
        y : Series
            same index as caller
        """
        if isinstance(arg, (dict, Series)):
            if isinstance(arg, dict):
                arg = Series(arg)

            indexer, mask = tseries.getMergeVec(self, arg.index.indexMap)
            notmask = -mask

            newValues = arg.view(np.ndarray).take(indexer)

            if notmask.any():
                if issubclass(newValues.dtype.type, np.integer):
                    newValues = newValues.astype(float)

                np.putmask(newValues, notmask, np.nan)

            newSer = Series(newValues, index=self.index)
            return newSer
        else:
            return Series([arg(x) for x in self], index=self.index)

    merge = map

    def apply(self, func):
        """
        Call function on elements on array. Can be ufunc or Python function
        expecting only single values

        Parameters
        ----------
        func : function

        Returns
        -------
        y : Series
        """
        try:
            return func(self)
        except Exception:
            return Series([func(x) for x in self], index=self.index)

    applymap = apply

    def reindex(self, new_index, method=None, fillMethod=None):
        """Conform Series to new Index

        Parameters
        ----------
        new_index :   array-like
            Preferably an Index object (to avoid duplicating data)
        method : {'backfill', 'pad', None}
            Method to use for filling holes in reindexed Series

            pad : propagate last valid observation forward to next valid
            backfill : use NEXT valid observation to fill gap

        Returns
        -------
        reindexed : Series
        """
        if fillMethod is not None: # pragma: no cover
            warnings.warn("'fillMethod' is deprecated. Use 'method' instead",
                          FutureWarning)

            method = fillMethod

        if self.index.equals(new_index):
            return self.copy()

        if not isinstance(new_index, Index):
            new_index = Index(new_index)

        if len(self.index) == 0:
            return Series(NaN, index=new_index)

        if method is not None:
            method = method.upper()

        # Cython for blazing speed
        fillVec, mask = tseries.getFillVec(self.index, new_index,
                                           self.index.indexMap,
                                           new_index.indexMap,
                                           kind=method)

        newValues = self.values.take(fillVec)

        notmask = -mask
        if notmask.any():
            if issubclass(newValues.dtype.type, np.int_):
                newValues = newValues.astype(float)
            elif issubclass(newValues.dtype.type, np.bool_):
                newValues = newValues.astype(object)

            np.putmask(newValues, notmask, NaN)

        return Series(newValues, index=new_index)

    def reindex_like(self, other, method=None):
        """
        Reindex Series to match index of another Series

        Parameters
        ----------
        other : Series
        method : string or None
            See Series.reindex docstring

        Notes
        -----
        Like calling s.reindex(other.index)

        Returns
        -------
        reindexed : Series
        """
        return self.reindex(other.index, method=method)

    def fill(self, value=None, method='pad'): # pragma: no cover
        warnings.warn("fill is being replaced by fillna, and the fill function "
                      "behavior will disappear in the next release: please "
                      "modify your code accordingly",
                      FutureWarning)
        return self.fillna(value=value, method=method)

    def fillna(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Parameters
        ----------
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        method : {'backfill', 'pad', None}
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
            return withoutna.reindex(self.index, method=method)

#-------------------------------------------------------------------------------
# Miscellaneous

    def plot(self, label=None, kind='line', rot=30, axes=None, style='-',
             **kwds): # pragma: no cover
        """
        Plot the input series with the index on the x-axis using
        matplotlib / pylab.

        Parameters
        ----------
        label : label argument to provide to plot
        kind : {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series
        kwds : other plotting keyword arguments

        Notes
        -----
        See matplotlib documentation online for more on this subject

        Default plot-types: TimeSeries (line), Series (bar)

        Intended to be used in ipython -pylab mode
        """
        import matplotlib.pyplot as plt

        if label is not None:
            kwds = kwds.copy()
            kwds['label'] = label

        N = len(self)

        if axes is None:
            axes = plt.gca()

        if kind == 'line':
            axes.plot(self.index, self.values, style, **kwds)
        elif kind == 'bar':
            xinds = np.arange(N) + 0.25
            axes.bar(xinds, self.values, 0.5, bottom=np.zeros(N), linewidth=1)

            if N < 10:
                fontsize = 12
            else:
                fontsize = 10

            axes.set_xticks(xinds + 0.25)
            axes.set_xticklabels(self.index, rotation=rot, fontsize=fontsize)

        plt.draw_if_interactive()

    def toCSV(self, path):
        """
        Write the Series to a CSV file

        Parameters
        ----------
        path : string or None
            Output filepath. If None, write to stdout
        """
        f = open(path, 'wb')

        for idx, value in self.iteritems():
            f.write(str(idx) + ',' + str(value) + '\n')

        f.close()

    def valid(self):
        """
        Return Series without NaN values

        Returns
        -------
        Series
        """
        return remove_na(self)

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

#-------------------------------------------------------------------------------
# Time series-oriented methods

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataMatrix and Series objects within
        by given number (positive or negative) of business/weekdays.

        Parameters
        ----------
        periods : int (+ or -)
            Number of periods to move
        offset : DateOffset, optional
            Increment to use from datetools module
        timeRule : string
            time rule name to use by name (e.g. 'WEEKDAY')

        Returns
        -------
        TimeSeries
        """
        if periods == 0:
            return self.copy()

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        if offset is None:
            newValues = np.empty(len(self), dtype=self.dtype)

            if periods > 0:
                newValues[periods:] = self.values[:-periods]
                newValues[:periods] = np.NaN
            elif periods < 0:
                newValues[:periods] = self.values[-periods:]
                newValues[periods:] = np.NaN

            return Series(newValues, index=self.index)
        else:
            return Series(self, index=self.index.shift(periods, offset))

    def truncate(self, before=None, after=None):
        """Function truncate a sorted TimeSeries before and/or after
        some particular dates.

        Parameters
        ----------
        before : date
            Truncate before date
        after : date
            Truncate after date

        Notes
        -----
        If TimeSeries is contained in a DataFrame, consider using the version
        of the function there.

        Returns
        -------
        TimeSeries
        """
        before = datetools.to_datetime(before)
        after = datetools.to_datetime(after)

        if before is None:
            beg_slice = 0
        elif before in self.index:
            beg_slice = self.index.indexMap[before]
        elif before < self.index[-1]:
            beg_slice = self.index.searchsorted(before, side='left')
        else:
            return Series([], index=NULL_INDEX)

        if after is None:
            end_slice = len(self)
        elif after in self.index:
            end_slice = self.index.indexMap[after] + 1
        elif after > self.index[0]:
            end_slice = self.index.searchsorted(after, side='right')
        else:
            return Series([], index=NULL_INDEX)

        return self[beg_slice:end_slice]

    def asOf(self, date):
        """
        Return last good (non-NaN) value in TimeSeries if value is NaN for
        requested date.

        If there is no good value, NaN is returned.

        Parameters
        ----------
        date : datetime or similar value

        Notes
        -----
        Dates are assumed to be sorted

        Returns
        -------
        value or NaN
        """
        if isinstance(date, basestring):
            date = datetools.to_datetime(date)

        v = self.get(date)

        if isnull(v):
            candidates = self.index[notnull(self)]
            index = candidates.searchsorted(date)

            if index > 0:
                asOfDate = candidates[index - 1]
            else:
                return NaN

            return self.get(asOfDate)
        else:
            return v

    def asfreq(self, freq, method=None, fillMethod=None):
        """
        Convert this TimeSeries to the provided frequency using DateOffset
        objects. Optionally provide fill method to pad/backfill/interpolate
        missing values.

        Parameters
        ----------
        offset : DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)
        method : {'backfill', 'pad', None}
            Method to use for filling holes in new index

        Returns
        -------
        TimeSeries
        """
        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], timeRule=freq)

        return self.reindex(dateRange, method=method, fillMethod=fillMethod)

    def interpolate(self, method='linear'):
        """
        Interpolate missing values (after the first valid value)

        Parameters
        ----------
        method : {'linear', 'time'}
            Interpolation method.

            Time interpolation works on daily and higher resolution
            data to interpolate given length of interval

        Returns
        -------
        Series with values interpolated
        """
        if method == 'time':
            if not isinstance(self, TimeSeries):
                raise Exception('time-weighted interpolation only works'
                                'on TimeSeries')
            inds = np.array([d.toordinal() for d in self.index])
        else:
            inds = np.arange(len(self))

        values = self.values

        invalid = isnull(values)
        valid = -invalid

        firstIndex = valid.argmax()
        valid = valid[firstIndex:]
        invalid = invalid[firstIndex:]
        inds = inds[firstIndex:]

        result = values.copy()
        result[firstIndex:][invalid] = np.interp(inds[invalid], inds[valid],
                                                 values[firstIndex:][valid])

        return Series(result, index=self.index)

    def rename(self, mapper):
        """
        Alter Series index using dict or function

        Parameters
        ----------
        mapper : dict-like or function
            Transformation to apply to each index

        Notes
        -----
        Function / dict values must be unique (1-to-1)

        Returns
        -------
        y : Series (new object)
        """
        if isinstance(mapper, (dict, Series)):
            mapper = mapper.__getitem__

        result = self.copy()
        result.index = [mapper(x) for x in self.index]

        return result

    @property
    def weekday(self):
        return Series([d.weekday() for d in self.index],
                      index=self.index)


class TimeSeries(Series):
    pass

#-------------------------------------------------------------------------------
# Supplementary functions

def remove_na(arr):
    """
    Return array containing only true/non-NaN values, possibly empty.
    """
    return arr[notnull(arr)]

def _seriesRepr(index, vals, nanRep='NaN'):
    string_index = [str(x) for x in index]
    maxlen = max(len(x) for x in string_index)
    padSpace = min(maxlen, 60)

    if vals.dtype == np.object_:
        def _format(k, v):
            return '%s    %s' % (str(k).ljust(padSpace), v)
    elif vals.dtype == np.float_:
        def _format(k, v):
            if np.isnan(v):
                v = nanRep
            else:
                v = str(v)

            return '%s    %s' % (str(k).ljust(padSpace), v)
    else:
        def _format(k, v):
            return '%s    %s' % (str(k).ljust(padSpace), v)

    it = itertools.starmap(_format,
                           itertools.izip(string_index, vals))

    return '\n'.join(it)
