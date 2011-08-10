"""
Data structure for 1-dimensional cross-sectional and time series data
"""

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import itertools
import operator
import sys
import warnings

from numpy import nan, ndarray
import numpy as np

from pandas.core.common import (isnull, notnull, _ensure_index,
                                _is_bool_indexer, _default_index)
from pandas.core.daterange import DateRange
from pandas.core.generic import PandasObject
from pandas.core.index import Index, MultiIndex
import pandas.core.datetools as datetools
import pandas._tseries as _tseries

__all__ = ['Series', 'TimeSeries']

def _numpy_lt_151():
    return np.__version__ < '1.5.1'

#-------------------------------------------------------------------------------
# Wrapper function for Series arithmetic methods

def _arith_method(op, name):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        from pandas.core.frame import DataFrame

        if isinstance(other, Series):
            if self.index.equals(other.index):
                return Series(op(self.values, other.values), index=self.index)

            new_index = self.index + other.index
            this_reindexed = self.reindex(new_index)
            other_reindexed = other.reindex(new_index)
            arr = op(this_reindexed.values, other_reindexed.values)
            return Series(arr, index=new_index)
        elif isinstance(other, DataFrame):
            return NotImplemented
        else:
            # scalars
            return Series(op(self.values, other), index=self.index)
    return wrapper

def _flex_method(op, name):
    def f(self, other, fill_value=None):
        return self._binop(other, op, fill_value=fill_value)

    f.__doc__ = """
    Binary operator %s with support to substitute a fill_value for missing data
    in one of the inputs

    Parameters
    ----------
    other: Series or scalar value
    fill_value : None or float value, default None
        Fill missing (NaN) values with this value. If both Series are
        missing, the result will be missing

    Returns
    -------
    result : Series
    """ % name
    f.__name__ = name
    return f

#-------------------------------------------------------------------------------
# Series class

class Series(np.ndarray, PandasObject):
    """
    Generic indexed (labeled) vector, including time series

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
    _AXIS_NUMBERS = {
        'index' : 0
    }

    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    def __new__(cls, data, index=None, dtype=None, name=None, copy=False):
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
            index = _default_index(len(subarr))

        # This is to prevent mixed-type Series getting all casted to
        # NumPy string type, e.g. NaN --> '-1#IND'.
        if issubclass(subarr.dtype.type, basestring):
            subarr = np.array(data, dtype=object, copy=copy)

        # Change the class of the array to be the subclass type.
        subarr = subarr.view(cls)
        subarr.index = index
        subarr.name = name

        if subarr.index.is_all_dates():
            subarr = subarr.view(TimeSeries)

        return subarr

    def __init__(self, *args, **kwargs):
        pass

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

        self._index = _ensure_index(index)

    index = property(fget=_get_index, fset=_set_index)

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)

    def toDict(self):
        return dict(self.iteritems())

    def to_sparse(self, kind='block', fill_value=None):
        """
        Convert Series to SparseSeries

        Parameters
        ----------
        kind : {'block', 'integer'}
        fill_value : float, defaults to NaN (missing)

        Returns
        -------
        sp : SparseSeries
        """
        from pandas.core.sparse import SparseSeries
        return SparseSeries(self, kind=kind, fill_value=fill_value)

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
        try:
            if isinstance(key, int):
                return self._regular_index(key)
            elif isinstance(self.index, MultiIndex):
                return self._multilevel_index(key)
            else:
                return self._regular_index(key)
        except Exception:
            pass

        self._check_bool_indexer(key)

        def _index_with(indexer):
            return Series(self.values[indexer],
                          index=self.index[indexer])

        # special handling of boolean data with NAs stored in object
        # arrays. Sort of an elaborate hack since we can't represent boolean
        # NA. Hmm
        if _is_bool_indexer(key):
            key = np.asarray(key, dtype=bool)
            return _index_with(key)

        # TODO: [slice(0, 5, None)] will break if you convert to ndarray,
        # e.g. as requested by np.median

        try:
            return _index_with(key)
        except Exception:
            key = np.asarray(key)
            return _index_with(key)

    def _regular_index(self, key):
        values = self.values

        try:
            return values[self.index.get_loc(key)]
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return values[key]
            raise Exception('Requested index not in this series!')

    def _multilevel_index(self, key):
        values = self.values
        try:
            loc = self.index.get_loc(key)
            if isinstance(loc, slice):
                return Series(values[loc], index=self.index[loc])
            else:
                return values[loc]
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return values[key]
            raise Exception('Requested index not in this series!')

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
            return self._get_val_at(self.index.get_loc(key))
        else:
            return default

    # help out SparseSeries
    _get_val_at = ndarray.__getitem__

    def __getslice__(self, i, j):
        """
        Returns a slice of the Series.

        Note that the underlying values are COPIES.

        The reason that the getslice returns copies is that otherwise you
        will have a reference to the original series which could be
        inadvertently changed
        """
        return Series(self.values[i:j].copy(), index=self.index[i:j])

    def __setitem__(self, key, value):
        values = self.values
        try:
            loc = self.index.get_loc(key)
            values[loc] = value
            return
        except KeyError:
            if isinstance(key, (int, np.integer)):
                values[key] = value
                return
            raise Exception('Requested index not in this series!')
        except TypeError:
            # Could not hash item
            pass

        self._check_bool_indexer(key)

        # special handling of boolean data with NAs stored in object
        # arrays. Sort of an elaborate hack since we can't represent boolean
        # NA. Hmm
        if isinstance(key, np.ndarray) and key.dtype == np.object_:
            mask = isnull(key)
            if mask.any():
                raise ValueError('cannot index with vector containing '
                                 'NA / NaN values')

            if set([True, False]).issubset(set(key)):
                key = np.asarray(key, dtype=bool)
                values[key] = value
                return

        values[key] = value

    def _check_bool_indexer(self, key):
        # boolean indexing, need to check that the data are aligned, otherwise
        # disallowed
        if isinstance(key, Series) and key.dtype == np.bool_:
            if not key.index.equals(self.index):
                raise Exception('can only boolean index with like-indexed '
                                'Series or raw ndarrays')

    def __setslice__(self, i, j, value):
        """Set slice equal to given value(s)"""
        ndarray.__setslice__(self, i, j, value)

    def __repr__(self):
        """Clean string representation of a Series"""
        if len(self.index) > 500:
            return self._make_repr(50)
        elif len(self.index) > 0:
            return _seriesRepr(self.index, self.values)
        else:
            return '%s' % ndarray.__repr__(self)

    def _make_repr(self, max_vals=50):
        vals = self.values
        index = self.index

        num = max_vals // 2
        head = _seriesRepr(index[:num], vals[:num])
        tail = _seriesRepr(index[-(max_vals - num):], vals[-(max_vals - num):])
        return head + '\n...\n' + tail + '\nlength: %d' % len(vals)

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

    __add__ = _arith_method(operator.add, '__add__')
    __sub__ = _arith_method(operator.sub, '__sub__')
    __mul__ = _arith_method(operator.mul, '__mul__')
    __div__ = _arith_method(operator.div, '__div__')
    __truediv__ = _arith_method(operator.truediv, '__truediv__')
    __pow__ = _arith_method(operator.pow, '__pow__')
    __truediv__ = _arith_method(operator.truediv, '__truediv__')

    __radd__ = _arith_method(operator.add, '__add__')
    __rmul__ = _arith_method(operator.mul, '__mul__')
    __rsub__ = _arith_method(lambda x, y: y - x, '__sub__')
    __rdiv__ = _arith_method(lambda x, y: y / x, '__div__')
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__truediv__')
    __rpow__ = _arith_method(lambda x, y: y ** x, '__pow__')

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

    def quantile(self, q=0.5):
        """
        Return value at the given quantile

        Parameters
        ----------
        q : quantile
            0 <= q <= 1

        Returns
        -------
        q : float
        """
        from scipy.stats import scoreatpercentile
        return scoreatpercentile(self.valid().values, q * 100)

    def describe(self):
        """
        Generate various summary statistics of columns, excluding NaN values

        Returns
        -------
        DataFrame
        """
        names = ['count', 'mean', 'std', 'min',
                 '10%', '50%', '90%', 'max']

        data = [self.count(), self.mean(), self.std(), self.min(),
                self.quantile(.1), self.median(), self.quantile(.9),
                self.max()]

        return Series(data, index=names)

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
            return nan
        return ndarray.std(nona, axis, dtype, out, ddof)

    def var(self, axis=None, dtype=None, out=None, ddof=1):
        """
        Unbiased variance of non-null values
        """
        nona = remove_na(self.values)
        if len(nona) < 2:
            return nan
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
        arr = self.values.copy()

        do_mask = not issubclass(self.dtype.type, np.int_)
        if do_mask:
            mask = isnull(arr)
            np.putmask(arr, mask, 0.)

        result = arr.cumsum()

        if do_mask:
            np.putmask(result, mask, np.nan)

        return Series(result, index=self.index)

    def cumprod(self, axis=0, dtype=None, out=None):
        """
        Overriding numpy's built-in cumprod functionality
        """
        arr = self.values.copy()

        do_mask = not issubclass(self.dtype.type, np.int_)
        if do_mask:
            mask = isnull(arr)
            np.putmask(arr, mask, 1.)

        result = arr.cumprod()

        if do_mask:
            np.putmask(result, mask, np.nan)

        return Series(result, index=self.index)

    def median(self):
        """
        Compute median value of non-null values
        """
        arr = self.values

        if arr.dtype != np.float_:
            arr = arr.astype(float)

        arr = arr[notnull(arr)]
        return _tseries.median(arr)

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
            return nan

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
        new_index = np.concatenate((self.index, other.index))
        new_index = Index(new_index)
        new_index._verify_integrity()

        new_values = np.concatenate((self, other))
        return Series(new_values, index=new_index)

    def _binop(self, other, func, fill_value=None):
        """
        Parameters
        ----------
        other : Series

        Returns
        -------
        combined : Series
        """
        # TODO: docstring

        assert(isinstance(other, Series))

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            new_index = self.index + other.index
            this = self.reindex(new_index)
            other = other.reindex(new_index)

        this_vals = this.values
        other_vals = other.values

        if fill_value is not None:
            this_mask = isnull(this_vals)
            other_mask = isnull(other_vals)
            this_vals = this_vals.copy()
            other_vals = other_vals.copy()

            # one but not both
            mask = this_mask ^ other_mask
            this_vals[this_mask & mask] = fill_value
            other_vals[other_mask & mask] = fill_value

        result = func(this_vals, other_vals)
        return Series(result, index=new_index)

    add = _flex_method(operator.add, 'add')
    sub = _flex_method(operator.sub, 'subtract')
    mul = _flex_method(operator.mul, 'multiply')
    div = _flex_method(operator.div, 'divide')

    def combine(self, other, func, fill_value=nan):
        """
        Perform elementwise binary operation on two Series using given function
        with optional fill value when an index is missing from one Series or the
        other

        Parameters
        ----------
        other : Series or scalar value
        func : function
        fill_value : scalar value

        Returns
        -------
        result : Series
        """
        if isinstance(other, Series):
            new_index = self.index + other.index

            new_values = np.empty(len(new_index), dtype=self.dtype)
            for i, idx in enumerate(new_index):
                new_values[i] = func(self.get(idx, fill_value),
                                 other.get(idx, fill_value))
        else:
            new_index = self.index
            new_values = func(self.values, other)

        return Series(new_values, index=new_index)

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
            new_index = self.index
            # save ourselves the copying in this case
            this = self
        else:
            new_index = self.index + other.index

            this = self.reindex(new_index)
            other = other.reindex(new_index)

        result = Series(np.where(isnull(this), other, this), index=new_index)
        return result

    #----------------------------------------------------------------------
    # Reindexing, sorting

    def sort(self, axis=0, kind='quicksort', order=None):
        """
        Overridden NumPy sort, taking care with missing values
        """
        sortedSeries = self.order(na_last=True)
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

    def order(self, na_last=True, ascending=True, **kwds):
        """
        Sorts Series object, by value, maintaining index-value object

        Parameters
        ----------
        na_last : boolean (optional, default=True)
            Put NaN's at beginning or end

        Returns
        -------
        y : Series
        """
        def _try_mergesort(arr):
            # easier to ask forgiveness than permission
            try:
                return arr.argsort(kind='mergesort')
            except TypeError:
                # stable sort not available for object dtype
                return arr.argsort()

        if 'missingAtEnd' in kwds: # pragma: no cover
            warnings.warn("missingAtEnd is deprecated, use na_last",
                          FutureWarning)
            na_last = kwds['missingAtEnd']

        arr = self.values
        sortedIdx = np.empty(len(self), dtype=np.int32)

        bad = isnull(arr)

        good = -bad
        idx = np.arange(len(self))

        argsorted = _try_mergesort(arr[good])

        if not ascending:
            argsorted = argsorted[::-1]

        if na_last:
            n = sum(good)
            sortedIdx[:n] = idx[good][argsorted]
            sortedIdx[n:] = idx[bad]
        else:
            n = sum(bad)
            sortedIdx[n:] = idx[good][argsorted]
            sortedIdx[:n] = idx[bad]

        return Series(arr[sortedIdx], index=self.index[sortedIdx])

    def sortlevel(self, level=0, ascending=True):
        """
        Sort multilevel index by chosen level. Data will be lexicographically
        sorted by the chosen level followed by the other levels (in order)

        Parameters
        ----------
        level : int
        ascending : bool, default True

        Returns
        -------
        sorted : Series
        """
        if not isinstance(self.index, MultiIndex):
            raise Exception('can only sort by level with a hierarchical index')

        new_index, indexer = self.index.sortlevel(level, ascending=ascending)
        new_values = self.values.take(indexer)
        return Series(new_values, index=new_index)

    def unstack(self, level=-1):
        """
        "Unstack" Series with multi-level index to produce DataFrame

        Parameters
        ----------
        level : int, default last level
            Level to "unstack"

        Examples
        --------
        >>> s
        one  a   1.
        one  b   2.
        two  a   3.
        two  b   4.

        >>> s.unstack(level=-1)
             a   b
        one  1.  2.
        two  3.  4.

        >>> s.unstack(level=0)
           one  two
        a  1.   2.
        b  3.   4.

        Returns
        -------
        unstacked : DataFrame
        """
        unstacker = _Unstacker(self.values, self.index, level=level)
        return unstacker.get_result()

    #----------------------------------------------------------------------
    # function application

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

            indexer, mask = _tseries.getMergeVec(self, arg.index.indexMap)
            notmask = -mask

            new_values = arg.view(np.ndarray).take(indexer)

            if notmask.any():
                if issubclass(new_values.dtype.type, np.integer):
                    new_values = new_values.astype(float)

                np.putmask(new_values, notmask, np.nan)

            newSer = Series(new_values, index=self.index)
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

    def reindex(self, index=None, method=None):
        """Conform Series to new Index

        Parameters
        ----------
        index : array-like
            Preferably an Index object (to avoid duplicating data)
        method : {'backfill', 'bfill', 'pad', 'ffill', None}
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Returns
        -------
        reindexed : Series
        """
        if self.index.equals(index):
            return self.copy()

        index = _ensure_index(index)
        if len(self.index) == 0:
            return Series(nan, index=index)

        fill_vec, mask = self.index.get_indexer(index, method=method)
        new_values = self.values.take(fill_vec)

        notmask = -mask
        if notmask.any():
            if issubclass(new_values.dtype.type, np.int_):
                new_values = new_values.astype(float)
            elif issubclass(new_values.dtype.type, np.bool_):
                new_values = new_values.astype(object)

            np.putmask(new_values, notmask, nan)

        return Series(new_values, index=index)

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

    def fillna(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Parameters
        ----------
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

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
            if method is None: # pragma: no cover
                raise ValueError('must specify a fill method')

            method = method.lower()

            if method == 'ffill':
                method = 'pad'
            if method == 'bfill':
                method = 'backfill'

            mask = isnull(self.values)

            if _numpy_lt_151(): # pragma: no cover
                mask = mask.astype(np.uint8)

            if method == 'pad':
                indexer = _tseries.get_pad_indexer(mask)
            elif method == 'backfill':
                indexer = _tseries.get_backfill_indexer(mask)

            new_values = self.values.take(indexer)
            return Series(new_values, index=self.index)

#-------------------------------------------------------------------------------
# Miscellaneous

    def plot(self, label=None, kind='line', use_index=True, rot=30, ax=None,
             style='-', **kwds): # pragma: no cover
        """
        Plot the input series with the index on the x-axis using
        matplotlib / pylab.

        Parameters
        ----------
        label : label argument to provide to plot
        kind : {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series
        auto_x : if True, it will use range(len(self)) as x-axis
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

        if ax is None:
            ax = plt.gca()

        if kind == 'line':
            if use_index:
                x = np.asarray(self.index)
            else:
                x = range(len(self))

            ax.plot(x, self.values.astype(float), style, **kwds)
        elif kind == 'bar':
            xinds = np.arange(N) + 0.25
            ax.bar(xinds, self.values.astype(float), 0.5,
                   bottom=np.zeros(N), linewidth=1)

            if N < 10:
                fontsize = 12
            else:
                fontsize = 10

            ax.set_xticks(xinds + 0.25)
            ax.set_xticklabels(self.index, rotation=rot, fontsize=fontsize)

        # kludge
        try:
            fig = plt.gcf()
            fig.autofmt_xdate()
        except Exception:
            pass

        plt.draw_if_interactive()

    def hist(self, ax=None): # pragma: no cover
        """
        Draw histogram of the input series using matplotlib / pylab.

        Parameters
        ----------

        Notes
        -----
        See matplotlib documentation online for more on this subject

        Default plot-types: TimeSeries (line), Series (bar)

        Intended to be used in ipython -pylab mode
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        ax.hist(self.values)

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

    def first_valid_index(self):
        if len(self) == 0:
            return None

        mask = isnull(self.values)
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[i]

    def last_valid_index(self):
        if len(self) == 0:
            return None

        mask = isnull(self.values[::-1])
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[len(self) - i - 1]

#-------------------------------------------------------------------------------
# Time series-oriented methods

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataFrame and Series objects within
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
            new_values = np.empty(len(self), dtype=self.dtype)

            if periods > 0:
                new_values[periods:] = self.values[:-periods]
                new_values[:periods] = nan
            elif periods < 0:
                new_values[:periods] = self.values[-periods:]
                new_values[periods:] = nan

            return Series(new_values, index=self.index)
        else:
            return Series(self, index=self.index.shift(periods, offset))

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
                return nan

            return self.get(asOfDate)
        else:
            return v

    def asfreq(self, freq, method=None):
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

        return self.reindex(dateRange, method=method)

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

    #----------------------------------------------------------------------
    # Deprecated stuff

    @classmethod
    def fromValue(cls, value=nan, index=None, dtype=None): # pragma: no cover
        warnings.warn("'fromValue', can call Series(value, index=index) now",
                      FutureWarning)

        return Series(value, index=index, dtype=dtype)

    def _firstTimeWithValue(self): # pragma: no cover
        warnings.warn("_firstTimeWithValue is deprecated. Use "
                      "first_valid_index instead", FutureWarning)
        return self.first_valid_index()

    def _lastTimeWithValue(self): # pragma: no cover
        warnings.warn("_firstTimeWithValue is deprecated. Use "
                      "last_valid_index instead", FutureWarning)
        return self.last_valid_index()

    _ix = None
    @property
    def ix(self):
        from pandas.core.indexing import _SeriesIndexer

        if self._ix is None:
            self._ix = _SeriesIndexer(self)

        return self._ix

class TimeSeries(Series):
    pass

#-------------------------------------------------------------------------------
# Supplementary functions

class _Unstacker(object):
    """
    Helper class to unstack data / pivot with multi-level index

    Parameters
    ----------
    level : int, default last level
        Level to "unstack"

    Examples
    --------
    >>> s
    one  a   1.
    one  b   2.
    two  a   3.
    two  b   4.

    >>> s.unstack(level=-1)
         a   b
    one  1.  2.
    two  3.  4.

    >>> s.unstack(level=0)
       one  two
    a  1.   2.
    b  3.   4.

    Returns
    -------
    unstacked : DataFrame
    """
    def __init__(self, values, index, level=-1, value_columns=None):
        if values.ndim == 1:
            values = values[:, np.newaxis]
        self.values = values
        self.value_columns = value_columns

        if value_columns is None and values.shape[1] != 1:
            raise ValueError('must pass column labels for multi-column data')

        self.index = index

        if level < 0:
            level += index.nlevels
        self.level = level

        self.new_index_levels = list(index.levels)
        self.removed_level = self.new_index_levels.pop(level)

        v = self.level
        lshape = self.index.levshape
        self.full_shape = np.prod(lshape[:v] + lshape[v+1:]), lshape[v]

        self._make_sorted_values_labels()
        self._make_selectors()

    def _make_sorted_values_labels(self):
        v = self.level

        labs = self.index.labels
        to_sort = labs[:v] + labs[v+1:] + [labs[v]]
        indexer = np.lexsort(to_sort[::-1])

        self.sorted_values = self.values.take(indexer, axis=0)
        self.sorted_labels = [l.take(indexer) for l in to_sort]

    def _make_selectors(self):
        new_levels = self.new_index_levels

        # make the mask
        group_index = self.sorted_labels[0]
        prev_stride = np.prod([len(x) for x in new_levels[1:]])

        for lev, lab in zip(new_levels[1:], self.sorted_labels[1:-1]):
            group_index = group_index * prev_stride + lab
            prev_stride /= len(lev)

        group_mask = np.zeros(self.full_shape[0], dtype=bool)
        group_mask.put(group_index, True)

        stride = self.index.levshape[self.level]
        selector = self.sorted_labels[-1] + stride * group_index
        mask = np.zeros(np.prod(self.full_shape), dtype=bool)
        mask.put(selector, True)

        # compress labels
        unique_groups = np.arange(self.full_shape[0])[group_mask]
        compressor = group_index.searchsorted(unique_groups)

        self.group_mask = group_mask
        self.group_index = group_index
        self.mask = mask
        self.unique_groups = unique_groups
        self.compressor = compressor

    def get_result(self):
        from pandas.core.frame import DataFrame
        return DataFrame(self.get_new_values(),
                         index=self.get_new_index(),
                         columns=self.get_new_columns())

    def get_new_values(self):
        # place the values
        length, width = self.full_shape
        stride = self.values.shape[1]
        result_width = width * stride

        new_values = np.empty((length, result_width), dtype=self.values.dtype)
        new_values.fill(np.nan)

        # is there a simpler / faster way of doing this?
        for i in xrange(self.values.shape[1]):
            chunk = new_values[:, i * width : (i + 1) * width]
            chunk.flat[self.mask] = self.sorted_values[:, i]

        new_values = new_values.take(self.unique_groups, axis=0)
        return new_values

    def get_new_columns(self):
        if self.value_columns is None:
            return self.removed_level

        stride = len(self.removed_level)
        width = len(self.value_columns)
        propagator = np.repeat(np.arange(width), stride)
        if isinstance(self.value_columns, MultiIndex):
            new_levels = self.value_columns.levels + [self.removed_level]
            new_labels = [lab.take(propagator)
                          for lab in self.value_columns.labels]
            new_labels.append(np.tile(np.arange(stride), width))
        else:
            new_levels = [self.value_columns, self.removed_level]
            new_labels = []

            new_labels.append(propagator)
            new_labels.append(np.tile(np.arange(stride), width))

        return MultiIndex(levels=new_levels, labels=new_labels)

    def get_new_index(self):
        result_labels = []
        for cur in self.sorted_labels[:-1]:
            result_labels.append(cur.take(self.compressor))

        # construct the new index
        if len(self.new_index_levels) == 1:
            new_index = Index(self.new_index_levels[0])
        else:
            new_index = MultiIndex(levels=self.new_index_levels,
                                   labels=result_labels)

        return new_index


def remove_na(arr):
    """
    Return array containing only true/non-NaN values, possibly empty.
    """
    return arr[notnull(arr)]

def _seriesRepr(index, vals, nanRep='NaN'):
    string_index = index.format()
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
