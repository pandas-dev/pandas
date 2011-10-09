"""
Data structure for 1-dimensional cross-sectional and time series data
"""

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import csv
import itertools
import operator
import sys
import warnings

from numpy import nan, ndarray
import numpy as np

from pandas.core.common import (isnull, notnull, _is_bool_indexer,
                                _default_index, _maybe_upcast)
from pandas.core.daterange import DateRange
from pandas.core.generic import PandasObject
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.indexing import _SeriesIndexer, _maybe_droplevels
from pandas.util.decorators import deprecate
from pandas.util import py3compat
import pandas.core.common as common
import pandas.core.datetools as datetools
import pandas._tseries as lib

__all__ = ['Series', 'TimeSeries']

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
                name = _maybe_match_name(self, other)
                return Series(op(self.values, other.values), index=self.index,
                              name=name)

            this_reindexed, other_reindexed = self.align(other, join='outer',
                                                         copy=False)
            arr = op(this_reindexed.values, other_reindexed.values)

            name = _maybe_match_name(self, other)
            return Series(arr, index=this_reindexed.index, name=name)
        elif isinstance(other, DataFrame):
            return NotImplemented
        else:
            # scalars
            return Series(op(self.values, other), index=self.index,
                          name=self.name)
    return wrapper

def _maybe_match_name(a, b):
    name = None
    if a.name == b.name:
        name = a.name
    return name

def _flex_method(op, name):
    def f(self, other, fill_value=None):
        return self._binop(other, op, fill_value=fill_value)

    f.__doc__ = """
    Binary operator %s with support to substitute a fill_value for missing data
    in one of the inputs

    Parameters
    ----------
    other: Series or scalar value
    fill_value : None or float value, default None (NaN)
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
            data = [data.get(idx, np.nan) for idx in index]

        # Create array, do *not* copy data by default, infer type
        try:
            subarr = np.array(data, dtype=dtype, copy=copy)
        except ValueError:
            if dtype:
                raise
            else:  # pragma: no cover
                subarr = np.array(data, dtype=object)

        if subarr.ndim == 0:
            if isinstance(data, list):  # pragma: no cover
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

    def __init__(self, data, index=None, dtype=None, name=None, copy=False):
        """One-dimensional ndarray with axis labels (including time
series). Labels must be unique and can any hashable type. The object supports
both integer- and label-based indexing and provides a host of methods for
performing operations involving the index. Statistical methods from ndarray have
been overridden to automatically exclude missing data (currently represented as
NaN)

Operations between Series (+, -, /, *, **) align values based on their
associated index values-- they need not be the same length. The result
index will be the sorted union of the two indexes.

Parameters
----------
data : array-like, dict, or scalar value
    Contains data stored in Series
index : array-like or Index (1d)

    Values must be unique and hashable, same length as data. Index object
    (or other iterable of same length as data) Will default to
    np.arange(len(data)) if not provided. If both a dict and index sequence
    are used, the index will override the keys found in the dict.

dtype : numpy.dtype or None
    If None, dtype will be inferred copy : boolean, default False Copy
    input data
copy : boolean, default False
        """
        pass

    @property
    def _constructor(self):
        return Series

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
        self.name = getattr(obj, 'name', None)

    def __contains__(self, key):
        return key in self.index

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))
        subclass_state = (self.index, self.name)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)

        # backwards compat
        index, name = own_state[0], None
        if len(own_state) > 1:
            name = own_state[1]

        self.index = index
        self.name = name

    _ix = None

    @property
    def ix(self):
        if self._ix is None:
            self._ix = _SeriesIndexer(self)

        return self._ix

    def __getitem__(self, key):
        # Label-based
        try:
            if isinstance(self.index, MultiIndex):
                return self._multilevel_index(key)
            else:
                values = self.values
                try:
                    return values[self.index.get_loc(key)]
                except KeyError, e1:
                    try:
                        return values[key]
                    except Exception, _:
                        pass
                    raise e1
        except TypeError:
            pass

        def _index_with(indexer):
            return Series(self.values[indexer], index=self.index[indexer],
                          name=self.name)

        # boolean

        # special handling of boolean data with NAs stored in object
        # arrays. Since we can't represent NA with dtype=bool
        if _is_bool_indexer(key):
            self._check_bool_indexer(key)
            key = np.asarray(key, dtype=bool)
            return _index_with(key)

        # other: fancy integer or otherwise

        # [slice(0, 5, None)] will break if you convert to ndarray,
        # e.g. as requested by np.median

        return _index_with(key)

    def _multilevel_index(self, key):
        values = self.values
        try:
            loc = self.index.get_loc(key)
            if isinstance(loc, (slice, np.ndarray)):
                # TODO: what if a level contains tuples??
                new_index = self.index[loc]
                new_index = _maybe_droplevels(new_index, key)
                return Series(values[loc], index=new_index, name=self.name)
            else:
                return values[loc]
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return values[key]
            raise KeyError('%s not in this series!' % str(key))

    # help out SparseSeries
    _get_val_at = ndarray.__getitem__

    def __getslice__(self, i, j):
        return self._constructor(self.values[i:j], index=self.index[i:j],
                                 name=self.name)

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
            raise KeyError('%s not in this series!' % str(key))
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
            return self._tidy_repr(30)
        elif len(self.index) > 0:
            return self._get_repr(name=True)
        else:
            return '%s' % ndarray.__repr__(self)

    def _tidy_repr(self, max_vals=20):
        num = max_vals // 2
        head = self[:num]._get_repr(name=False)
        tail = self[-(max_vals - num):]._get_repr(name=False)
        result = head + '\n...\n' + tail
        result = '%s\nName: %s, Length: %d' % (result, self.name, len(self))
        return result

    def to_string(self, buffer=sys.stdout, nanRep='NaN'):
        print >> buffer, self._get_repr(nanRep=nanRep)

    def _get_repr(self, name=False, nanRep='NaN'):
        vals = self.values
        index = self.index

        string_index = index.format()
        maxlen = max(len(x) for x in string_index)
        padSpace = min(maxlen, 60)

        def _format_float(k, v):
            if np.isnan(v):
                v = nanRep
            else:
                v = str(v)
            return '%s    %s' % (str(k).ljust(padSpace), v)

        def _format_nonfloat(k, v):
            return '%s    %s' % (str(k).ljust(padSpace), v)

        if vals.dtype == np.float_:
            _format = _format_float
        else:
            _format = _format_nonfloat

        it = itertools.starmap(_format,
                               itertools.izip(string_index, vals))
        it = list(it)
        if name:
            it.append('Name: %s, Length: %d' % (str(self.name), len(self)))
        return '\n'.join(it)

    def __str__(self):
        return repr(self)

    def __iter__(self):
        return iter(self.values)

    def iteritems(self):
        """
        Lazily iterate over (index, value) tuples
        """
        return itertools.izip(iter(self.index), iter(self))

    iterkv = iteritems
    if py3compat.PY3:
        items = iteritems

    #----------------------------------------------------------------------
    #   Arithmetic operators

    __add__ = _arith_method(operator.add, '__add__')
    __sub__ = _arith_method(operator.sub, '__sub__')
    __mul__ = _arith_method(operator.mul, '__mul__')
    __truediv__ = _arith_method(operator.truediv, '__truediv__')
    __floordiv__ = _arith_method(operator.floordiv, '__floordiv__')
    __pow__ = _arith_method(operator.pow, '__pow__')

    __radd__ = _arith_method(operator.add, '__add__')
    __rmul__ = _arith_method(operator.mul, '__mul__')
    __rsub__ = _arith_method(lambda x, y: y - x, '__sub__')
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__truediv__')
    __rfloordiv__ = _arith_method(lambda x, y: y // x, '__floordiv__')
    __rpow__ = _arith_method(lambda x, y: y ** x, '__pow__')

    # Inplace operators
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __itruediv__ = __truediv__
    __ifloordiv__ = __floordiv__
    __ipow__ = __pow__

    # Python 2 division operators
    if not py3compat.PY3:
        __div__ = _arith_method(operator.div, '__div__')
        __rdiv__ = _arith_method(lambda x, y: y / x, '__div__')
        __idiv__ = __div__

    #----------------------------------------------------------------------
    # Misc public methods

    def keys(self):
        "Alias for index"
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

    def copy(self):
        """
        Return new Series with copy of underlying values

        Returns
        -------
        cp : Series
        """
        return Series(self.values.copy(), index=self.index, name=self.name)

    def to_dict(self):
        """
        Convert Series to {label -> value} dict

        Returns
        -------
        value_dict : dict
        """
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
        return SparseSeries(self, kind=kind, fill_value=fill_value,
                            name=self.name)

    def get(self, key, default=None):
        """
        Returns value occupying requested index, default to specified
        missing value if not present. Analogous to dict.get

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

    #----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck

    def count(self, level=None):
        """
        Return number of non-NA/null observations in the Series

        Returns
        -------
        nobs : int
        """
        if level is not None:
            return self._count_level(level)

        return notnull(self.values).sum()

    def _count_level(self, level):
        # TODO: GENERALIZE CODE OVERLAP WITH DATAFRAME
        # TODO: deal with sortedness??
        obj = self.sortlevel(level)
        mask = notnull(obj.values)

        level_index = obj.index.levels[level]

        if len(self) == 0:
            return Series(0, index=level_index)

        n = len(level_index)
        locs = obj.index.labels[level].searchsorted(np.arange(n))

        # WORKAROUND: reduceat fusses about the endpoints. should file ticket?
        start = locs.searchsorted(0, side='right') - 1
        end = locs.searchsorted(len(mask), side='left')

        result = np.zeros((n), dtype=int)
        out = result[start:end]
        np.add.reduceat(mask, locs[start:end], out=out)

        # WORKAROUND: to see why, try this
        # arr = np.ones((10, 4), dtype=bool)
        # np.add.reduceat(arr, [0, 3, 3, 7, 9], axis=0)

        # this stinks
        if len(locs) > 1:
            workaround_mask = locs[:-1] == locs[1:]
            result[:-1][workaround_mask] = 0

        return Series(result, index=level_index)

    def sum(self, axis=0, dtype=None, out=None, skipna=True):
        """
        Sum of values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        sum : float
        """
        values = self.values.copy()

        if skipna:
            mask = isnull(values)
            if mask.all():
                return np.nan
            np.putmask(values, mask, 0)

        return values.sum()

    def mean(self, axis=0, dtype=None, out=None, skipna=True):
        """
        Mean of values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        mean : float
        """
        return self._ndarray_statistic('mean', dtype=dtype, skipna=skipna)

    def median(self, skipna=True):
        """
        Compute median of values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        median : float
        """
        arr = self.values
        if arr.dtype != np.float_:
            arr = arr.astype(float)
        mask = notnull(arr)

        if skipna:
            arr = arr[mask]
        else:
            if not mask.all():
                return np.nan

        return lib.median(arr)

    def prod(self, axis=0, dtype=None, out=None, skipna=True):
        """
        Product of all values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        product : float
        """
        return self._ndarray_statistic('prod', dtype=dtype, skipna=skipna)

    def min(self, axis=None, out=None, skipna=True):
        """
        Minimum of values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        min : float
        """
        arr = self.values.copy()
        if skipna:
            if not issubclass(arr.dtype.type, np.int_):
                np.putmask(arr, isnull(arr), np.inf)
        return arr.min()

    def max(self, axis=None, out=None, skipna=True):
        """
        Maximum of values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        max : float
        """
        arr = self.values.copy()
        if skipna:
            if not issubclass(arr.dtype.type, np.int_):
                np.putmask(arr, isnull(arr), -np.inf)
        return arr.max()

    def std(self, axis=None, dtype=None, out=None, ddof=1, skipna=True):
        """
        Unbiased standard deviation of values

        Extra parameters are to preserve ndarray interface.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        stdev : float
        """
        if skipna:
            nona = remove_na(self.values)
            if len(nona) < 2:
                return nan
            return ndarray.std(nona, axis, dtype, out, ddof)
        else:
            return self.values.std(axis, dtype, out, ddof)

    def var(self, axis=None, dtype=None, out=None, ddof=1, skipna=True):
        """
        Unbiased variance of non-NA/null values

        Extra parameters are to preserve ndarray interface.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        var : float
        """
        if skipna:
            nona = remove_na(self.values)
            if len(nona) < 2:
                return nan
            return ndarray.var(nona, axis, dtype, out, ddof)
        else:
            return self.values.var(axis, dtype, out, ddof)

    def skew(self, skipna=True):
        """
        Unbiased skewness of the non-NA/null values

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        skew : float
        """
        y = np.array(self.values)
        mask = notnull(y)
        count = mask.sum()

        if count < len(self) and not skipna:
            return np.nan

        np.putmask(y, -mask, 0)
        A = y.sum() / count
        B = (y**2).sum() / count  - A**2
        C = (y**3).sum() / count - A**3 - 3*A*B

        return (np.sqrt((count**2-count))*C) / ((count-2)*np.sqrt(B)**3)

    def cumsum(self, axis=0, dtype=None, out=None, skipna=True):
        """
        Cumulative sum of values. Preserves locations of NaN values

        Extra parameters are to preserve ndarray interface.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        cumsum : Series
        """
        arr = self.values.copy()

        do_mask = skipna and not issubclass(self.dtype.type, np.int_)
        if do_mask:
            mask = isnull(arr)
            np.putmask(arr, mask, 0.)

        result = arr.cumsum()

        if do_mask:
            np.putmask(result, mask, np.nan)

        return Series(result, index=self.index)

    def cumprod(self, axis=0, dtype=None, out=None, skipna=True):
        """
        Cumulative product of values. Preserves locations of NaN values

        Extra parameters are to preserve ndarray interface.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values

        Returns
        -------
        cumprod : Series
        """
        arr = self.values.copy()

        do_mask = skipna and not issubclass(self.dtype.type, np.int_)
        if do_mask:
            mask = isnull(arr)
            np.putmask(arr, mask, 1.)

        result = arr.cumprod()

        if do_mask:
            np.putmask(result, mask, np.nan)

        return Series(result, index=self.index)

    def _ndarray_statistic(self, funcname, dtype=None, skipna=True):
        arr = self.values
        retVal = getattr(arr, funcname)(dtype=dtype)

        if skipna and isnull(retVal):
            arr = remove_na(arr)
            if len(arr) == 0:
                return np.nan
            retVal = getattr(arr, funcname)(dtype=dtype)

        return retVal

    def quantile(self, q=0.5):
        """
        Return value at the given quantile, a la scoreatpercentile in
        scipy.stats

        Parameters
        ----------
        q : quantile
            0 <= q <= 1

        Returns
        -------
        quantile : float
        """
        from scipy.stats import scoreatpercentile
        return scoreatpercentile(self.dropna().values, q * 100)

    def describe(self):
        """
        Generate various summary statistics of Series, excluding NaN
        values. These include: count, mean, std, min, max, and 10%/50%/90%
        quantiles

        Returns
        -------
        desc : Series
        """
        names = ['count', 'mean', 'std', 'min',
                 '25%', '50%', '75%', 'max']

        data = [self.count(), self.mean(), self.std(), self.min(),
                self.quantile(.25), self.median(), self.quantile(.75),
                self.max()]

        return Series(data, index=names)

    def corr(self, other):
        """
        Compute correlation two Series, excluding missing values

        Parameters
        ----------
        other : Series

        Returns
        -------
        correlation : float
        """
        commonIdx = self.dropna().index.intersection(other.dropna().index)

        if len(commonIdx) == 0:
            return nan

        this = self.reindex(commonIdx)
        that = other.reindex(commonIdx)

        return np.corrcoef(this, that)[0, 1]

    def diff(self, periods=1):
        """
        1st discrete difference of object

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for forming difference

        Returns
        -------
        diffed : Series
        """
        return (self - self.shift(periods))

    def autocorr(self):
        """
        Lag-1 autocorrelation

        Returns
        -------
        autocorr : float
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
        clipped : Series
        """
        result = self
        if lower is not None:
            result = result.clip_lower(lower)
        if upper is not None:
            result = result.clip_upper(upper)

        return result

    def clip_upper(self, threshold):
        """
        Return copy of series with values above given value truncated

        See also
        --------
        clip

        Returns
        -------
        clipped : Series
        """
        return np.where(self > threshold, threshold, self)

    def clip_lower(self, threshold):
        """
        Return copy of series with values below given value truncated

        See also
        --------
        clip

        Returns
        -------
        clipped : Series
        """
        return np.where(self < threshold, threshold, self)

#-------------------------------------------------------------------------------
# Combination

    def append(self, other):
        """
        Concatenate two Series. The indexes must not overlap

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
        """
        new_index = self.index.append(other.index)
        new_index._verify_integrity()

        new_values = np.concatenate((self.values, other.values))
        name = _maybe_match_name(self, other)
        return self._constructor(new_values, index=new_index, name=name)

    def _binop(self, other, func, fill_value=None):
        """
        Perform generic binary operation with optional fill value

        Parameters
        ----------
        other : Series
        func : binary operator
        fill_value : float or object
            Value to substitute for NA/null values. If both Series are NA in a
            location, the result will be NA regardless of the passed fill value

        Returns
        -------
        combined : Series
        """
        assert(isinstance(other, Series))

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            this, other = self.align(other, join='outer')
            new_index = this.index

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
        name = _maybe_match_name(self, other)
        return Series(result, index=new_index, name=name)

    add = _flex_method(operator.add, 'add')
    sub = _flex_method(operator.sub, 'subtract')
    mul = _flex_method(operator.mul, 'multiply')
    try:
        div = _flex_method(operator.div, 'divide')
    except AttributeError:    # Python 3
        div = _flex_method(operator.truediv, 'divide')

    def combine(self, other, func, fill_value=nan):
        """
        Perform elementwise binary operation on two Series using given function
        with optional fill value when an index is missing from one Series or
        the other

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
            new_name = _maybe_match_name(self, other)
            new_values = np.empty(len(new_index), dtype=self.dtype)
            for i, idx in enumerate(new_index):
                new_values[i] = func(self.get(idx, fill_value),
                                     other.get(idx, fill_value))
        else:
            new_index = self.index
            new_values = func(self.values, other)
            new_name = self.name
        return Series(new_values, index=new_index, name=new_name)

    def combine_first(self, other):
        """
        Combine Series values, choosing the calling Series's values
        first. Result index will be the union of the two indexes

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
        """
        new_index = self.index + other.index
        this = self.reindex(new_index, copy=False)
        other = other.reindex(new_index, copy=False)
        name = _maybe_match_name(self, other)
        return Series(np.where(isnull(this), other, this), index=new_index,
                      name=name)

    #----------------------------------------------------------------------
    # Reindexing, sorting

    def sort(self, axis=0, kind='quicksort', order=None):
        """
        Sort values and index labels in place, for compatibility with
        ndarray. No return value
        """
        sortedSeries = self.order(na_last=True)
        self[:] = sortedSeries
        self.index = sortedSeries.index

    def sort_index(self, ascending=True):
        """
        Sort object by labels (along an axis)

        Parameters
        ----------
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted_obj : Series
        """
        labels = self.index
        sort_index = labels.argsort()
        if not ascending:
            sort_index = sort_index[::-1]
        new_labels = labels.take(sort_index)
        new_values = self.values.take(sort_index)
        return Series(new_values, new_labels, name=self.name)

    def argsort(self, axis=0, kind='quicksort', order=None):
        """
        Overrides ndarray.argsort. Argsorts the value, omitting NA/null values,
        and places the result in the same locations as the non-NA values

        Returns
        -------
        argsorted : Series
        """
        values = self.values
        mask = isnull(values)

        if mask.any():
            result = values.copy()
            notmask = -mask
            result[notmask] = np.argsort(values[notmask])
            return Series(result, index=self.index, name=self.name)
        else:
            return Series(np.argsort(values), index=self.index, name=self.name)

    def order(self, na_last=True, ascending=True, **kwds):
        """
        Sorts Series object, by value, maintaining index-value link

        Parameters
        ----------
        na_last : boolean (optional, default=True)
            Put NaN's at beginning or end
        ascending : boolean, default True
            Sort ascending. Passing False sorts descending

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

        if 'missingAtEnd' in kwds:  # pragma: no cover
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

        return Series(arr[sortedIdx], index=self.index[sortedIdx],
                      name=self.name)

    def sortlevel(self, level=0, ascending=True):
        """
        Sort Series with MultiIndex by chosen level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

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
        return Series(new_values, index=new_index, name=self.name)

    def swaplevel(self, i, j, copy=True):
        """
        Swap levels i and j in a MultiIndex

        Returns
        -------
        swapped : Series
        """
        new_index = self.index.swaplevel(i, j)
        return Series(self.values, index=new_index, copy=copy, name=self.name)

    def unstack(self, level=-1):
        """
        Unstack, a.k.a. pivot, Series with MultiIndex to produce DataFrame

        Parameters
        ----------
        level : int, default last level
            Level to unstack

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
        from pandas.core.reshape import _Unstacker
        unstacker = _Unstacker(self.values, self.index, level=level)
        return unstacker.get_result()

    #----------------------------------------------------------------------
    # function application

    def map(self, arg):
        """
        Map values of Series using input correspondence (which can be
        a dict, Series, or function)

        Parameters
        ----------
        arg : function, dict, or Series

        Examples
        --------
        >>> x
        one   1
        two   2
        three 3

        >>> y
        1  foo
        2  bar
        3  baz

        >>> x.map(y)
        one   foo
        two   bar
        three baz

        Returns
        -------
        y : Series
            same index as caller
        """
        if isinstance(arg, (dict, Series)):
            if isinstance(arg, dict):
                arg = Series(arg)

            indexer = lib.merge_indexer_object(self.values.astype(object),
                                               arg.index.indexMap)

            new_values = common.take_1d(np.asarray(arg), indexer)
            return Series(new_values, index=self.index, name=self.name)
        else:
            return Series([arg(x) for x in self], index=self.index,
                          name=self.name)

    def apply(self, func):
        """
        Invoke function on values of Series. Can be ufunc or Python function
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
            return Series([func(x) for x in self], index=self.index,
                          name=self.name)

    def align(self, other, join='outer', copy=True):
        """
        Align two Series object with the specified join method

        Parameters
        ----------
        other : Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'

        Returns
        -------
        (left, right) : (Series, Series)
            Aligned Series
        """
        join_index, lidx, ridx = self.index.join(other.index, how=join,
                                                 return_indexers=True)

        if lidx is not None:
            left = Series(common.take_1d(self.values, lidx), join_index,
                          name=self.name)
        else:
            if copy:
                new_values = self.values.copy()
            else:
                new_values = self.values
            left = Series(new_values, join_index, name=self.name)

        if ridx is not None:
            right = Series(common.take_1d(other.values, ridx), join_index,
                           name=other.name)
        else:
            if copy:
                new_values = other.values.copy()
            else:
                new_values = other.values
            right = Series(new_values, join_index, name=other.name)

        return left, right

    def reindex(self, index=None, method=None, copy=True):
        """Conform Series to new index with optional filling logic, placing
        NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        copy=False

        Parameters
        ----------
        index : array-like or Index
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        method : {'backfill', 'bfill', 'pad', 'ffill', None}
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate LAST valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same

        Returns
        -------
        reindexed : Series
        """
        if self.index.equals(index):
            if copy:
                return self.copy()
            else:
                return self

        index = _ensure_index(index)
        if len(self.index) == 0:
            return Series(nan, index=index, name=self.name)

        new_index, fill_vec = self.index.reindex(index, method=method)
        new_values = common.take_1d(self.values, fill_vec)
        return Series(new_values, index=new_index, name=self.name)

    def reindex_like(self, other, method=None):
        """
        Reindex Series to match index of another Series, optionally with
        filling logic

        Parameters
        ----------
        other : Series
        method : string or None
            See Series.reindex docstring

        Notes
        -----
        Like calling s.reindex(other.index, method=...)

        Returns
        -------
        reindexed : Series
        """
        return self.reindex(other.index, method=method)

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take, return Series corresponding to requested
        indices

        Parameters
        ----------
        indices : list / array of ints

        Returns
        -------
        taken : Series
        """
        new_index = self.index.take(indices)
        new_values = self.values.take(indices)
        return Series(new_values, index=new_index, name=self.name)

    def fillna(self, value=None, method='pad'):
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : Series
        """
        if value is not None:
            newSeries = self.copy()
            newSeries[isnull(newSeries)] = value
            return newSeries
        else:
            if method is None:  # pragma: no cover
                raise ValueError('must specify a fill method')

            method = method.lower()

            if method == 'ffill':
                method = 'pad'
            if method == 'bfill':
                method = 'backfill'

            mask = isnull(self.values)

            # sadness. for Python 2.5 compatibility
            mask = mask.astype(np.uint8)

            if method == 'pad':
                indexer = lib.get_pad_indexer(mask)
            elif method == 'backfill':
                indexer = lib.get_backfill_indexer(mask)

            new_values = self.values.take(indexer)
            return Series(new_values, index=self.index, name=self.name)

#-------------------------------------------------------------------------------
# Miscellaneous

    def plot(self, label=None, kind='line', use_index=True, rot=30, ax=None,
             style='-', grid=True, **kwds):  # pragma: no cover
        """
        Plot the input series with the index on the x-axis using matplotlib

        Parameters
        ----------
        label : label argument to provide to plot
        kind : {'line', 'bar'}
        rot : int, default 30
            Rotation for tick labels
        use_index : boolean, default True
            Plot index as axis tick labels
        ax : matplotlib axis object
            If not passed, uses gca()
        style : string, default '-'
            matplotlib line style to use
        kwds : keywords
            To be passed to the actual plotting function

        Notes
        -----
        See matplotlib documentation online for more on this subject
        Intended to be used in ipython --pylab mode
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
                   bottom=np.zeros(N), linewidth=1, **kwds)

            if N < 10:
                fontsize = 12
            else:
                fontsize = 10

            ax.set_xticks(xinds + 0.25)
            ax.set_xticklabels(self.index, rotation=rot, fontsize=fontsize)

        ax.grid(grid)

        # try to make things prettier
        try:
            fig = plt.gcf()
            fig.autofmt_xdate()
        except Exception:
            pass

        plt.draw_if_interactive()

    def hist(self, ax=None, grid=True, **kwds):  # pragma: no cover
        """
        Draw histogram of the input series using matplotlib

        Parameters
        ----------
        ax : matplotlib axis object
            If not passed, uses gca()
        kwds : keywords
            To be passed to the actual plotting function

        Notes
        -----
        See matplotlib documentation online for more on this

        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        values = self.dropna().values

        ax.hist(values)
        ax.grid(grid)

    def to_csv(self, path):
        """
        Write the Series to a CSV file

        Parameters
        ----------
        path : string or None
            Output filepath. If None, write to stdout
        """
        f = open(path, 'w')
        csvout = csv.writer(f, lineterminator='\n')
        csvout.writerows(self.iteritems())
        f.close()

    def dropna(self):
        """
        Return Series without null values

        Returns
        -------
        valid : Series
        """
        return remove_na(self)

    valid = dropna

    isnull = isnull
    notnull = notnull

    def first_valid_index(self):
        """
        Return label for first non-NA/null value
        """
        if len(self) == 0:
            return None

        mask = isnull(self.values)
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[i]

    def last_valid_index(self):
        """
        Return label for last non-NA/null value
        """
        if len(self) == 0:
            return None

        mask = isnull(self.values[::-1])
        i = mask.argmin()
        if mask[i]:
            return None
        else:
            return self.index[len(self) - i - 1]

    #----------------------------------------------------------------------
    # Time series-oriented methods

    def shift(self, periods, offset=None, **kwds):
        """
        Shift the index of the Series by desired number of periods with an
        optional time offset

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        offset : DateOffset, timedelta, or time rule string, optional
            Increment to use from datetools module or time rule (e.g. 'EOM')

        Returns
        -------
        shifted : Series
        """
        if periods == 0:
            return self.copy()

        offset = kwds.get('timeRule', offset)
        if isinstance(offset, basestring):
            offset = datetools.getOffset(offset)

        if offset is None:
            new_values = np.empty(len(self), dtype=self.dtype)
            new_values = _maybe_upcast(new_values)

            if periods > 0:
                new_values[periods:] = self.values[:-periods]
                new_values[:periods] = nan
            elif periods < 0:
                new_values[:periods] = self.values[-periods:]
                new_values[periods:] = nan

            return Series(new_values, index=self.index, name=self.name)
        else:
            return Series(self, index=self.index.shift(periods, offset),
                          name=self.name)

    def asof(self, date):
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
        object or time rule. Optionally provide fill method to pad/backfill
        missing values.

        Parameters
        ----------
        offset : DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)
        method : {'backfill', 'pad', None}
            Method to use for filling holes in new index

        Returns
        -------
        converted : TimeSeries
        """
        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], time_rule=freq)

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
        interpolated : Series
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

        return Series(result, index=self.index, name=self.name)

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

        Examples
        --------
        >>> x
        foo 1
        bar 2
        baz 3

        >>> x.rename(str.upper)
        FOO 1
        BAR 2
        BAZ 3

        >>> x.rename({'foo' : 'a', 'bar' : 'b', 'baz' : 'c'})
        a 1
        b 2
        c 3

        Returns
        -------
        renamed : Series (new object)
        """
        if isinstance(mapper, (dict, Series)):
            def mapper_f(x):
                if x in mapper:
                    return mapper[x]
                else:
                    return x
        else:
            mapper_f = mapper

        result = self.copy()
        result.index = [mapper_f(x) for x in self.index]

        return result

    @property
    def weekday(self):
        return Series([d.weekday() for d in self.index], index=self.index)

    #----------------------------------------------------------------------
    # Deprecated stuff

    @classmethod
    def fromValue(cls, value=nan, index=None, dtype=None):  # pragma: no cover
        warnings.warn("'fromValue', can call Series(value, index=index) now",
                      FutureWarning)
        return Series(value, index=index, dtype=dtype)

    asOf = deprecate('asOf', asof)
    toDict = deprecate('toDict', to_dict)
    toString = deprecate('toString', to_string)
    merge = deprecate('merge', map)
    applymap = deprecate('applymap', apply)
    combineFirst = deprecate('combineFirst', combine_first)
    _firstTimeWithValue = deprecate('_firstTimeWithValue', first_valid_index)
    _lastTimeWithValue = deprecate('_lastTimeWithValue', last_valid_index)
    toCSV = deprecate('toCSV', to_csv)

class TimeSeries(Series):
    pass

#-------------------------------------------------------------------------------
# Supplementary functions

def remove_na(arr):
    """
    Return array containing only true/non-NaN values, possibly empty.
    """
    return arr[notnull(arr)]
