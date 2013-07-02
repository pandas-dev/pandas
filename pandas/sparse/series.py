"""
Data structures for sparse float data. Life is made simpler by dealing only
with float64 data
"""

# pylint: disable=E1101,E1103,W0231

from numpy import nan, ndarray
import numpy as np

import operator

from pandas.core.common import isnull
from pandas.core.index import Index, _ensure_index
from pandas.core.series import Series, TimeSeries, _maybe_match_name
from pandas.core.frame import DataFrame
import pandas.core.common as com
import pandas.core.datetools as datetools

from pandas.util import py3compat

from pandas.sparse.array import (make_sparse, _sparse_array_op, SparseArray)
from pandas._sparse import BlockIndex, IntIndex
import pandas._sparse as splib

from pandas.util.decorators import Appender

#------------------------------------------------------------------------------
# Wrapper function for Series arithmetic methods


def _sparse_op_wrap(op, name):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        if isinstance(other, Series):
            if not isinstance(other, SparseSeries):
                other = other.to_sparse(fill_value=self.fill_value)
            return _sparse_series_op(self, other, op, name)
        elif isinstance(other, DataFrame):
            return NotImplemented
        elif np.isscalar(other):
            if isnull(other) or isnull(self.fill_value):
                new_fill_value = np.nan
            else:
                new_fill_value = op(np.float64(self.fill_value),
                                    np.float64(other))

            return SparseSeries(op(self.sp_values, other),
                                index=self.index,
                                sparse_index=self.sp_index,
                                fill_value=new_fill_value,
                                name=self.name)
        else:  # pragma: no cover
            raise TypeError('operation with %s not supported' % type(other))

    wrapper.__name__ = name
    return wrapper


def _sparse_series_op(left, right, op, name):
    left, right = left.align(right, join='outer', copy=False)
    new_index = left.index
    new_name = _maybe_match_name(left, right)

    result = _sparse_array_op(left, right, op, name)
    result = result.view(SparseSeries)
    result.index = new_index
    result.name = new_name

    return result


class SparseSeries(SparseArray, Series):
    """Data structure for labeled, sparse floating point data

    Parameters
    ----------
    data : {array-like, Series, SparseSeries, dict}
    kind : {'block', 'integer'}
    fill_value : float
        Defaults to NaN (code for missing)
    sparse_index : {BlockIndex, IntIndex}, optional
        Only if you have one. Mainly used internally

    Notes
    -----
    SparseSeries objects are immutable via the typical Python means. If you
    must change values, convert to dense, make your changes, then convert back
    to sparse
    """
    __array_priority__ = 15

    sp_index = None
    fill_value = None

    def __new__(cls, data, index=None, sparse_index=None, kind='block',
                fill_value=None, name=None, copy=False):

        is_sparse_array = isinstance(data, SparseArray)
        if fill_value is None:
            if is_sparse_array:
                fill_value = data.fill_value
            else:
                fill_value = nan

        if is_sparse_array:
            if isinstance(data, SparseSeries) and index is None:
                index = data.index
            elif index is not None:
                if not (len(index) == len(data)):
                    raise AssertionError()

            sparse_index = data.sp_index
            values = np.asarray(data)
        elif isinstance(data, (Series, dict)):
            if index is None:
                index = data.index

            data = Series(data)
            values, sparse_index = make_sparse(data, kind=kind,
                                               fill_value=fill_value)
        elif isinstance(data, (tuple, list, np.ndarray)):
            # array-like
            if sparse_index is None:
                values, sparse_index = make_sparse(data, kind=kind,
                                                   fill_value=fill_value)
            else:
                values = data
                if not (len(values) == sparse_index.npoints):
                    raise AssertionError()
        else:
            if index is None:
                raise TypeError('must pass index!')

            length = len(index)

            if data == fill_value or (isnull(data)
                                      and isnull(fill_value)):
                if kind == 'block':
                    sparse_index = BlockIndex(length, [], [])
                else:
                    sparse_index = IntIndex(length, [])
                values = np.array([])
            else:
                if kind == 'block':
                    locs, lens = ([0], [length]) if length else ([], [])
                    sparse_index = BlockIndex(length, locs, lens)
                else:
                    sparse_index = IntIndex(length, index)
                values = np.empty(length)
                values.fill(data)

        if index is None:
            index = com._default_index(sparse_index.length)
        index = _ensure_index(index)

        # Create array, do *not* copy data by default
        if copy:
            subarr = np.array(values, dtype=np.float64, copy=True)
        else:
            subarr = np.asarray(values, dtype=np.float64)

        if index.is_all_dates:
            cls = SparseTimeSeries

        # Change the class of the array to be the subclass type.
        output = subarr.view(cls)
        output.sp_index = sparse_index
        output.fill_value = np.float64(fill_value)
        output.index = index
        output.name = name
        return output

    def _make_time_series(self):
        # oh boy #2139
        self.__class__ = SparseTimeSeries

    @classmethod
    def from_array(cls, arr, index=None, name=None, copy=False, fill_value=None):
        """
        Simplified alternate constructor
        """
        return SparseSeries(arr, index=index, name=name, copy=copy, fill_value=fill_value)

    def __init__(self, data, index=None, sparse_index=None, kind='block',
                 fill_value=None, name=None, copy=False):
        pass

    @property
    def _constructor(self):
        def make_sp_series(data, index=None, name=None):
            return SparseSeries(data, index=index, fill_value=self.fill_value,
                                kind=self.kind, name=name)

        return make_sp_series

    @property
    def kind(self):
        if isinstance(self.sp_index, BlockIndex):
            return 'block'
        elif isinstance(self.sp_index, IntIndex):
            return 'integer'

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)
        self.name = getattr(obj, 'name', None)
        self.sp_index = getattr(obj, 'sp_index', None)
        self.fill_value = getattr(obj, 'fill_value', None)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))

        subclass_state = (self.index, self.fill_value, self.sp_index,
                          self.name)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)

        index, fill_value, sp_index = own_state[:3]
        name = None
        if len(own_state) > 3:
            name = own_state[3]

        self.sp_index = sp_index
        self.fill_value = fill_value
        self.index = index
        self.name = name

    def __len__(self):
        return self.sp_index.length

    def __unicode__(self):
        # currently, unicode is same as repr...fixes infinite loop
        series_rep = Series.__unicode__(self)
        rep = '%s\n%s' % (series_rep, repr(self.sp_index))
        return rep

    # Arithmetic operators

    __add__ = _sparse_op_wrap(operator.add, 'add')
    __sub__ = _sparse_op_wrap(operator.sub, 'sub')
    __mul__ = _sparse_op_wrap(operator.mul, 'mul')
    __truediv__ = _sparse_op_wrap(operator.truediv, 'truediv')
    __floordiv__ = _sparse_op_wrap(operator.floordiv, 'floordiv')
    __pow__ = _sparse_op_wrap(operator.pow, 'pow')

    # reverse operators
    __radd__ = _sparse_op_wrap(operator.add, '__radd__')
    __rsub__ = _sparse_op_wrap(lambda x, y: y - x, '__rsub__')
    __rmul__ = _sparse_op_wrap(operator.mul, '__rmul__')
    __rtruediv__ = _sparse_op_wrap(lambda x, y: y / x, '__rtruediv__')
    __rfloordiv__ = _sparse_op_wrap(lambda x, y: y // x, 'floordiv')
    __rpow__ = _sparse_op_wrap(lambda x, y: y ** x, '__rpow__')

    # Python 2 division operators
    if not py3compat.PY3:
        __div__ = _sparse_op_wrap(operator.div, 'div')
        __rdiv__ = _sparse_op_wrap(lambda x, y: y / x, '__rdiv__')

    def __getitem__(self, key):
        """

        """
        try:
            return self._get_val_at(self.index.get_loc(key))

        except KeyError:
            if isinstance(key, (int, np.integer)):
                return self._get_val_at(key)
            raise Exception('Requested index not in this series!')

        except TypeError:
            # Could not hash item, must be array-like?
            pass

        # is there a case where this would NOT be an ndarray?
        # need to find an example, I took out the case for now

        dataSlice = self.values[key]
        new_index = Index(self.index.view(ndarray)[key])
        return self._constructor(dataSlice, index=new_index, name=self.name)

    def abs(self):
        """
        Return an object with absolute value taken. Only applicable to objects
        that are all numeric

        Returns
        -------
        abs: type of caller
        """
        res_sp_values = np.abs(self.sp_values)
        return SparseSeries(res_sp_values, index=self.index,
                            sparse_index=self.sp_index,
                            fill_value=self.fill_value)

    def get(self, label, default=None):
        """
        Returns value occupying requested label, default to specified
        missing value if not present. Analogous to dict.get

        Parameters
        ----------
        label : object
            Label value looking for
        default : object, optional
            Value to return if label not in index

        Returns
        -------
        y : scalar
        """
        if label in self.index:
            loc = self.index.get_loc(label)
            return self._get_val_at(loc)
        else:
            return default

    def get_value(self, label):
        """
        Retrieve single value at passed index label

        Parameters
        ----------
        index : label

        Returns
        -------
        value : scalar value
        """
        loc = self.index.get_loc(label)
        return self._get_val_at(loc)

    def set_value(self, label, value):
        """
        Quickly set single value at passed label. If label is not contained, a
        new object is created with the label placed at the end of the result
        index

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed
        value : object
            Scalar value

        Notes
        -----
        This method *always* returns a new object. It is not particularly
        efficient but is provided for API compatibility with Series

        Returns
        -------
        series : SparseSeries
        """
        dense = self.to_dense().set_value(label, value)
        return dense.to_sparse(kind=self.kind, fill_value=self.fill_value)

    def to_dense(self, sparse_only=False):
        """
        Convert SparseSeries to (dense) Series
        """
        if sparse_only:
            int_index = self.sp_index.to_int_index()
            index = self.index.take(int_index.indices)
            return Series(self.sp_values, index=index, name=self.name)
        else:
            return Series(self.values, index=self.index, name=self.name)

    @property
    def density(self):
        r = float(self.sp_index.npoints) / float(self.sp_index.length)
        return r

    def astype(self, dtype=None):
        """

        """
        if dtype is not None and dtype not in (np.float_, float):
            raise TypeError('Can only support floating point data')

        return self.copy()

    def copy(self, deep=True):
        """
        Make a copy of the SparseSeries. Only the actual sparse values need to
        be copied
        """
        if deep:
            values = self.sp_values.copy()
        else:
            values = self.sp_values
        return SparseSeries(values, index=self.index,
                            sparse_index=self.sp_index,
                            fill_value=self.fill_value, name=self.name)

    def reindex(self, index=None, method=None, copy=True, limit=None):
        """
        Conform SparseSeries to new Index

        See Series.reindex docstring for general behavior

        Returns
        -------
        reindexed : SparseSeries
        """
        new_index = _ensure_index(index)

        if self.index.equals(new_index):
            if copy:
                return self.copy()
            else:
                return self

        if len(self.index) == 0:
            # FIXME: inelegant / slow
            values = np.empty(len(new_index), dtype=np.float64)
            values.fill(nan)
            return SparseSeries(values, index=new_index,
                                fill_value=self.fill_value)

        new_index, fill_vec = self.index.reindex(index, method=method,
                                                 limit=limit)
        new_values = com.take_1d(self.values, fill_vec)
        return SparseSeries(new_values, index=new_index,
                            fill_value=self.fill_value, name=self.name)

    def sparse_reindex(self, new_index):
        """
        Conform sparse values to new SparseIndex

        Parameters
        ----------
        new_index : {BlockIndex, IntIndex}

        Returns
        -------
        reindexed : SparseSeries
        """
        if not (isinstance(new_index, splib.SparseIndex)):
            raise AssertionError()

        new_values = self.sp_index.to_int_index().reindex(self.sp_values,
                                                          self.fill_value,
                                                          new_index)
        return SparseSeries(new_values, index=self.index,
                            sparse_index=new_index,
                            fill_value=self.fill_value)

    @Appender(Series.fillna.__doc__)
    def fillna(self, value=None, method=None, inplace=False, limit=None):
        dense = self.to_dense()
        filled = dense.fillna(value=value, method=method, limit=limit)
        result = filled.to_sparse(kind=self.kind,
                                  fill_value=self.fill_value)

        if inplace:
            self.sp_values[:] = result.values
            return self
        else:
            return result

    def take(self, indices, axis=0, convert=True):
        """
        Sparse-compatible version of ndarray.take

        Returns
        -------
        taken : ndarray
        """
        new_values = SparseArray.take(self, indices)
        new_index = self.index.take(indices)
        return self._constructor(new_values, index=new_index)

    def cumsum(self, axis=0, dtype=None, out=None):
        """
        Cumulative sum of values. Preserves locations of NaN values

        Extra parameters are to preserve ndarray interface.

        Returns
        -------
        cumsum : Series or SparseSeries
        """
        result = SparseArray.cumsum(self)
        if isinstance(result, SparseArray):
            result = self._attach_meta(result)
        return result

    def _attach_meta(self, sparse_arr):
        sparse_series = sparse_arr.view(SparseSeries)
        sparse_series.index = self.index
        sparse_series.name = self.name
        return sparse_series

    def dropna(self):
        """
        Analogous to Series.dropna. If fill_value=NaN, returns a dense Series
        """
        # TODO: make more efficient
        dense_valid = self.to_dense().valid()
        if isnull(self.fill_value):
            return dense_valid
        else:
            return dense_valid.to_sparse(fill_value=self.fill_value)

    def shift(self, periods, freq=None, **kwds):
        """
        Analogous to Series.shift
        """
        from pandas.core.series import _resolve_offset

        offset = _resolve_offset(freq, kwds)

        # no special handling of fill values yet
        if not isnull(self.fill_value):
            dense_shifted = self.to_dense().shift(periods, freq=freq,
                                                  **kwds)
            return dense_shifted.to_sparse(fill_value=self.fill_value,
                                           kind=self.kind)

        if periods == 0:
            return self.copy()

        if offset is not None:
            return SparseSeries(self.sp_values,
                                sparse_index=self.sp_index,
                                index=self.index.shift(periods, offset),
                                fill_value=self.fill_value)

        int_index = self.sp_index.to_int_index()
        new_indices = int_index.indices + periods
        start, end = new_indices.searchsorted([0, int_index.length])

        new_indices = new_indices[start:end]

        new_sp_index = IntIndex(len(self), new_indices)
        if isinstance(self.sp_index, BlockIndex):
            new_sp_index = new_sp_index.to_block_index()

        return SparseSeries(self.sp_values[start:end].copy(),
                            index=self.index,
                            sparse_index=new_sp_index,
                            fill_value=self.fill_value)

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
        if isinstance(other, SparseSeries):
            other = other.to_dense()

        dense_combined = self.to_dense().combine_first(other)
        return dense_combined.to_sparse(fill_value=self.fill_value)


class SparseTimeSeries(SparseSeries, TimeSeries):
    """Data structure for labeled, sparse floating point data, with `TimeStamp`
    index labels

    Parameters
    ----------
    data : {array-like, Series, SparseSeries, dict}
    kind : {'block', 'integer'}
    fill_value : float
        Defaults to NaN (code for missing)
    sparse_index : {BlockIndex, IntIndex}, optional
        Only if you have one. Mainly used internally

    Notes
    -----
    SparseSeries objects are immutable via the typical Python means. If you
    must change values, convert to dense, make your changes, then convert back
    to sparse
    """

    pass
