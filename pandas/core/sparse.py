"""
Data structures for sparse float data. Life is made simpler by dealing only with
float64 data
"""

# pylint: disable=E1101,E1103,W0231

from numpy import nan, ndarray
import numpy as np

import operator

from pandas.core.common import (isnull, _pickle_array, _unpickle_array,
                                _mut_exclusive, _try_sort)
from pandas.core.index import Index, MultiIndex, NULL_INDEX, _ensure_index
from pandas.core.series import Series, TimeSeries, _maybe_match_name
from pandas.core.frame import (DataFrame, extract_index, _prep_ndarray,
                               _default_index)
from pandas.core.panel import Panel, LongPanel
import pandas.core.common as common
import pandas.core.datetools as datetools

from pandas.util import py3compat

from pandas._sparse import BlockIndex, IntIndex
import pandas._sparse as splib

def make_sparse(arr, kind='block', fill_value=nan):
    """
    Convert ndarray to sparse format

    Parameters
    ----------
    arr : ndarray
    kind : {'block', 'integer'}
    fill_value : NaN or another value

    Returns
    -------
    (sparse_values, index) : (ndarray, SparseIndex)
    """
    arr = np.asarray(arr)
    length = len(arr)

    if np.isnan(fill_value):
        mask = -np.isnan(arr)
    else:
        mask = arr != fill_value

    indices = np.arange(length, dtype=np.int32)[mask]

    if kind == 'block':
        locs, lens = splib.get_blocks(indices)
        index = BlockIndex(length, locs, lens)
    elif kind == 'integer':
        index = IntIndex(length, indices)
    else: # pragma: no cover
        raise ValueError('must be block or integer type')

    sparsified_values = arr[mask]
    return sparsified_values, index

#-------------------------------------------------------------------------------
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
            new_fill_value = op(np.float64(self.fill_value),
                                np.float64(other))

            return SparseSeries(op(self.sp_values, other),
                                index=self.index,
                                sparse_index=self.sp_index,
                                fill_value=new_fill_value,
                                name=self.name)
        else: # pragma: no cover
            raise TypeError('operation with %s not supported' % type(other))

    wrapper.__name__ = name
    return wrapper

def _sparse_series_op(left, right, op, name):
    if np.isnan(left.fill_value):
        sparse_op = lambda a, b: _sparse_nanop(a, b, name)
    else:
        sparse_op = lambda a, b: _sparse_fillop(a, b, name)

    new_index = left.index + right.index
    if not left.index.equals(new_index):
        left = left.reindex(new_index)

    if not right.index.equals(new_index):
        right = right.reindex(new_index)

    if left.sp_index.equals(right.sp_index):
        result = op(left.sp_values, right.sp_values)
        result_index = left.sp_index
    else:
        result, result_index = sparse_op(left, right)

    try:
        fill_value = op(left.fill_value, right.fill_value)
    except ZeroDivisionError:
        fill_value = nan

    new_name = _maybe_match_name(left, right)
    return SparseSeries(result, index=new_index,
                        sparse_index=result_index,
                        fill_value=fill_value, name=new_name)

def _sparse_nanop(this, other, name):
    sparse_op = getattr(splib, 'sparse_nan%s' % name)
    result, result_index = sparse_op(this.sp_values,
                                     this.sp_index,
                                     other.sp_values,
                                     other.sp_index)

    return result, result_index

def _sparse_fillop(this, other, name):
    sparse_op = getattr(splib, 'sparse_%s' % name)
    result, result_index = sparse_op(this.sp_values,
                                     this.sp_index,
                                     this.fill_value,
                                     other.sp_values,
                                     other.sp_index,
                                     other.fill_value)

    return result, result_index


class SparseSeries(Series):
    __array_priority__ = 15

    sp_index = None
    fill_value = None

    def __new__(cls, data, index=None, sparse_index=None, kind='block',
                fill_value=None, name=None, copy=False):

        is_sparse_series = isinstance(data, SparseSeries)
        if fill_value is None:
            if is_sparse_series:
                fill_value = data.fill_value
            else:
                fill_value = nan

        if is_sparse_series:
            if index is None:
                index = data.index
            else:
                assert(len(index) == len(data))

            sparse_index = data.sp_index
            values = np.asarray(data)
        elif isinstance(data, (Series, dict)):
            if index is None:
                index = data.index

            data = Series(data)
            values, sparse_index = make_sparse(data, kind=kind,
                                               fill_value=fill_value)
        elif np.isscalar(data): # pragma: no cover
            if index is None:
                raise Exception('must pass index!')

            values = np.empty(len(index))
            values.fill(data)

            # TODO: more efficient

            values, sparse_index = make_sparse(values, kind=kind,
                                               fill_value=fill_value)

        else:
            # array-like
            if sparse_index is None:
                values, sparse_index = make_sparse(data, kind=kind,
                                                   fill_value=fill_value)
            else:
                values = data
                assert(len(values) == sparse_index.npoints)

        if index is None:
            index = Index(np.arange(sparse_index.length))
        index = _ensure_index(index)

        # Create array, do *not* copy data by default
        if copy:
            subarr = np.array(values, dtype=np.float64, copy=True)
        else:
            subarr = np.asarray(values, dtype=np.float64)

        if index.is_all_dates():
            cls = SparseTimeSeries

        # Change the class of the array to be the subclass type.
        output = subarr.view(cls)
        output._sp_values = subarr
        output.sp_index = sparse_index
        output.fill_value = np.float64(fill_value)
        output.index = index
        output.name = name
        return output

    def __init__(self, data, index=None, sparse_index=None, kind='block',
                 fill_value=None, name=None, copy=False):
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

    def __repr__(self):
        series_rep = Series.__repr__(self)
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

    # Inplace operators
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __itruediv__ = __truediv__
    __ifloordiv__ = __floordiv__
    __ipow__ = __pow__

    # Python 2 division operators
    if not py3compat.PY3:
        __div__ = _sparse_op_wrap(operator.div, 'div')
        __rdiv__ = _sparse_op_wrap(lambda x, y: y / x, '__rdiv__')
        __idiv__ = __div__

    @property
    def values(self):
        output = np.empty(len(self), dtype=np.float64)
        int_index = self.sp_index.to_int_index()
        output.fill(self.fill_value)
        output.put(int_index.indices, self)
        return output

    @property
    def sp_values(self):
        try:
            return self._sp_values
        except AttributeError:
            self._sp_values = ret = np.asarray(self)
            return ret

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

    def _get_val_at(self, loc):
        n = len(self)
        if loc < 0:
            loc += n

        if loc >= len(self) or loc < 0:
            raise Exception('Out of bounds access')

        sp_loc = self.sp_index.lookup(loc)
        if sp_loc == -1:
            return self.fill_value
        else:
            return ndarray.__getitem__(self, sp_loc)

    def take(self, indices):
        """
        Sparse-compatible version of ndarray.take

        Returns
        -------
        taken : ndarray
        """
        indices = np.asarray(indices, dtype=int)

        n = len(self)
        if (indices < 0).any() or (indices >= n).any():
            raise Exception('out of bounds access')

        if self.sp_index.npoints > 0:
            locs = np.array([self.sp_index.lookup(loc) for loc in indices])
            result = self.sp_values.take(locs)
            result[locs == -1] = self.fill_value
        else:
            result = np.empty(len(indices))
            result.fill(self.fill_value)

        return result

    def __setitem__(self, key, value):
        raise Exception('SparseSeries objects are immutable')

    def __setslice__(self, i, j, value):
        raise Exception('SparseSeries objects are immutable')

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

    def astype(self, dtype=None):
        """

        """
        if dtype is not None and dtype not in (np.float_, float):
            raise Exception('Can only support floating point data')

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

    def reindex(self, index=None, method=None, copy=True):
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

        new_index, fill_vec = self.index.reindex(index, method=method)
        new_values = common.take_1d(self.values, fill_vec)
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
        assert(isinstance(new_index, splib.SparseIndex))

        new_values = self.sp_index.to_int_index().reindex(self.sp_values,
                                                          self.fill_value,
                                                          new_index)
        return SparseSeries(new_values, index=self.index,
                            sparse_index=new_index,
                            fill_value=self.fill_value)

    def count(self):
        """
        Compute sum of non-NA/null observations in SparseSeries. If the
        fill_value is not NaN, the "sparse" locations will be included in the
        observation count

        Returns
        -------
        nobs : int
        """
        sp_values = self.sp_values
        valid_spvals = np.isfinite(sp_values).sum()
        if self._null_fill_value:
            return valid_spvals
        else:
            return valid_spvals + (len(self) - len(sp_values))

    @property
    def _null_fill_value(self):
        return np.isnan(self.fill_value)

    @property
    def _valid_sp_values(self):
        sp_vals = self.sp_values
        mask = np.isfinite(sp_vals)
        return sp_vals[mask]

    def sum(self, axis=None, dtype=None, out=None):
        """
        Sum of non-NA/null values

        Returns
        -------
        sum : float
        """
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        if self._null_fill_value:
            return sp_sum
        else:
            nsparse = self.sp_index.npoints
            return sp_sum + self.fill_value * nsparse

    def cumsum(self, axis=0, dtype=None, out=None):
        """
        Cumulative sum of values. Preserves locations of NaN values

        Extra parameters are to preserve ndarray interface.

        Returns
        -------
        cumsum : Series
        """
        if not np.isnan(self.fill_value):
            return self.to_dense().cumsum()
        return SparseSeries(self.sp_values.cumsum(),
                            index=self.index,
                            sparse_index=self.sp_index,
                            fill_value=self.fill_value)

    def mean(self, axis=None, dtype=None, out=None):
        """
        Mean of non-NA/null values

        Returns
        -------
        mean : float
        """
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        ct = len(valid_vals)

        if self._null_fill_value:
            return sp_sum / ct
        else:
            nsparse = self.sp_index.npoints
            return (sp_sum + self.fill_value * nsparse) / (ct + nsparse)

    def valid(self):
        """
        Analogous to Series.valid
        """
        # TODO: make more efficient
        dense_valid = self.to_dense().valid()
        return dense_valid.to_sparse(fill_value=self.fill_value)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Analogous to Series.shift
        """
        # no special handling of fill values yet
        if not isnull(self.fill_value):
            dense_shifted = self.to_dense().shift(periods, offset=offset,
                                                  timeRule=timeRule)
            return dense_shifted.to_sparse(fill_value=self.fill_value,
                                           kind=self.kind)

        if periods == 0:
            return self.copy()

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

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
        dense_combined = self.to_dense().combine_first(other.to_dense())
        return dense_combined.to_sparse(fill_value=self.fill_value)

class SparseTimeSeries(SparseSeries, TimeSeries):
    pass

class SparseDataFrame(DataFrame):
    """
    DataFrame containing sparse floating point data in the form of SparseSeries
    objects

    Parameters
    ----------
    data : same types as can be passed to DataFrame
    index : array-like, optional
    column : array-like, optional
    default_kind : {'block', 'integer'}, default 'block'
        Default sparse kind for converting Series to SparseSeries. Will not
        override SparseSeries passed into constructor
    default_fill_value : float
        Default fill_value for converting Series to SparseSeries. Will not
        override SparseSeries passed in
    """
    _verbose_info = False
    _columns = None
    _series = None
    _is_mixed_type = False
    ndim = 2

    def __init__(self, data=None, index=None, columns=None,
                 default_kind='block', default_fill_value=None):
        if default_fill_value is None:
            default_fill_value = nan

        self.default_kind = default_kind
        self.default_fill_value = default_fill_value

        if isinstance(data, dict):
            sdict, columns, index = self._init_dict(data, index, columns)
        elif isinstance(data, (np.ndarray, list)):
            sdict, columns, index = self._init_matrix(data, index, columns)
        elif isinstance(data, DataFrame):
            sdict, columns, index = self._init_dict(data, data.index,
                                                    data.columns)
        elif data is None:
            sdict = {}

            if index is None:
                index = NULL_INDEX

            if columns is None:
                columns = NULL_INDEX
            else:
                for c in columns:
                    sdict[c] = Series(np.NaN, index=index)

        self._series = sdict
        self.columns = columns
        self.index = index

    def _get_numeric_columns(self):
        # everything is necessarily float64
        return self.columns

    def _consolidate_inplace(self):
        # do nothing when DataFrame calls this method
        pass

    @property
    def _constructor(self):
        def wrapper(data, index=None, columns=None):
            return SparseDataFrame(data, index=index, columns=columns,
                                   default_fill_value=self.default_fill_value,
                                   default_kind=self.default_kind)
        return wrapper

    def _init_dict(self, data, index, columns, dtype=None):
        # pre-filter out columns if we passed it
        if columns is not None:
            columns = _ensure_index(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in columns)
        else:
            columns = Index(_try_sort(data.keys()))

        if index is None:
            index = extract_index(data)

        sp_maker = lambda x: SparseSeries(x, index=index,
                                          kind=self.default_kind,
                                          fill_value=self.default_fill_value,
                                          copy=True)

        sdict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Force alignment, no copy necessary
                if not v.index.equals(index):
                    v = v.reindex(index)

                if not isinstance(v, SparseSeries):
                    v = sp_maker(v)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, nan) for i in index]

                v = sp_maker(v)
            sdict[k] = v

        # TODO: figure out how to handle this case, all nan's?
        # add in any other columns we want to have (completeness)
        nan_vec = np.empty(len(index))
        nan_vec.fill(nan)
        for c in columns:
            if c not in sdict:
                sdict[c] = sp_maker(nan_vec)

        return sdict, columns, index

    def _init_matrix(self, data, index, columns, dtype=None):
        data = _prep_ndarray(data, copy=False)
        N, K = data.shape
        if index is None:
            index = _default_index(N)
        if columns is None:
            columns = _default_index(K)

        if len(columns) != K:
            raise Exception('Column length mismatch: %d vs. %d' %
                            (len(columns), K))
        if len(index) != N:
            raise Exception('Index length mismatch: %d vs. %d' %
                            (len(index), N))

        data = dict([(idx, data[:, i]) for i, idx in enumerate(columns)])
        return self._init_dict(data, index, columns, dtype)

    def __array_wrap__(self, result):
        return SparseDataFrame(result, index=self.index, columns=self.columns,
                               default_kind=self.default_kind,
                               default_fill_value=self.default_fill_value)

    def __getstate__(self):
        series = dict((k, (v.sp_index, v.sp_values))
                      for k, v in self.iteritems())
        columns = _pickle_array(self.columns)
        index = _pickle_array(self.index)

        return (series, columns, index, self.default_fill_value,
                self.default_kind)

    def __setstate__(self, state):
        series, cols, idx, fv, kind = state
        columns = _unpickle_array(cols)
        index = _unpickle_array(idx)

        series_dict = {}
        for col, (sp_index, sp_values) in series.iteritems():
            series_dict[col] = SparseSeries(sp_values, sparse_index=sp_index,
                                            fill_value=fv)

        self._series = series_dict
        self.index = index
        self.columns = columns
        self.default_fill_value = fv
        self.default_kind = kind

    def to_dense(self):
        """
        Convert to dense DataFrame

        Returns
        -------
        df : DataFrame
        """
        data = dict((k, v.to_dense()) for k, v in self.iteritems())
        return DataFrame(data, index=self.index)

    def copy(self, deep=True):
        """
        Make a copy of this SparseDataFrame
        """
        series = self._series.copy()
        return SparseDataFrame(series, index=self.index, columns=self.columns,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

    @property
    def density(self):
        """
        Ratio of non-sparse points to total (dense) data points
        represented in the frame
        """
        tot_nonsparse = sum([ser.sp_index.npoints
                             for _, ser in self.iteritems()])
        tot = len(self.index) * len(self.columns)
        return tot_nonsparse / float(tot)

    #----------------------------------------------------------------------
    # Support different internal rep'n of SparseDataFrame

    def _set_item(self, key, value):
        sp_maker = lambda x: SparseSeries(x, index=self.index,
                                          fill_value=self.default_fill_value,
                                          kind=self.default_kind)
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                cleanSeries = value.reindex(self.index)
                if not isinstance(value, SparseSeries):
                    cleanSeries = sp_maker(cleanSeries)
            else:
                cleanSeries = sp_maker(value)

            self._series[key] = cleanSeries
        # Scalar
        else:
            self._series[key] = sp_maker(value)

        if key not in self.columns:
            self._insert_column(key)

    def _insert_column(self, key):
        self.columns = Index(np.concatenate((self.columns, [key])))

    def __delitem__(self, key):
        """
        Delete column from DataFrame
        """
        loc = self.columns.get_loc(key)
        del self._series[key]
        self._delete_column_index(loc)

    def _delete_column_index(self, loc):
        if loc == len(self.columns) - 1:
            new_columns = self.columns[:loc]
        else:
            new_columns = Index(np.concatenate((self.columns[:loc],
                                               self.columns[loc+1:])))
        self.columns = new_columns

    _index = None
    def _set_index(self, index):
        self._index = _ensure_index(index)
        for v in self._series.values():
            v.index = self._index

    def _get_index(self):
        return self._index

    def _get_columns(self):
        return self._columns

    def _set_columns(self, cols):
        if len(cols) != len(self._series):
            raise Exception('Columns length %d did not match data %d!' %
                            (len(cols), len(self._series)))
        self._columns = _ensure_index(cols)

    index = property(fget=_get_index, fset=_set_index)
    columns = property(fget=_get_columns, fset=_set_columns)

    def __getitem__(self, item):
        """
        Retrieve column or slice from DataFrame
        """
        try:
            # unsure about how kludgy this is
            s = self._series[item]
            s.name = item
            return s
        except (TypeError, KeyError):
            if isinstance(item, slice):
                dateRange = self.index[item]
                return self.reindex(dateRange)

            elif isinstance(item, np.ndarray):
                if len(item) != len(self.index):
                    raise Exception('Item wrong length %d instead of %d!' %
                                    (len(item), len(self.index)))
                newIndex = self.index[item]
                return self.reindex(newIndex)
            else: # pragma: no cover
                raise

    def _slice(self, slobj, axis=0):
        if axis == 0:
            new_index = self.index[slobj]
            new_columns = self.columns
        else:
            new_index = self.index
            new_columns = self.columns[slobj]

        return self.reindex(index=new_index, columns=new_columns)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        if columns is None:
            columns = self.columns

        if len(columns) == 0:
            return np.zeros((len(self.index), 0), dtype=float)

        return np.array([self[col].values for col in columns]).T

    values = property(as_matrix)

    def xs(self, key, axis=0, copy=False):
        """
        Returns a row (cross-section) from the SparseDataFrame as a Series
        object.

        Parameters
        ----------
        key : some index contained in the index

        Returns
        -------
        xs : Series
        """
        if axis == 1:
            data = self[key]
            return data

        i = self.index.get_loc(key)
        series = self._series
        values = [series[k][i] for k in self.columns]
        return Series(values, index=self.columns)

    #----------------------------------------------------------------------
    # Arithmetic-related methods

    def _combine_frame(self, other, func, fill_value=None):
        new_index = self.index.union(other.index)
        new_columns = self.columns.union(other.columns)

        if fill_value is not None:
            raise NotImplementedError

        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)
            other = other.reindex(new_index)

        if not self and not other:
            return SparseDataFrame(index=new_index)

        if not other:
            return self * nan

        if not self:
            return other * nan

        new_data = {}
        for col in new_columns:
            if col in this and col in other:
                new_data[col] = func(this[col], other[col])

        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns)

    def _combine_match_index(self, other, func, fill_value=None):
        new_data = {}

        if fill_value is not None:
            raise NotImplementedError

        new_index = self.index.union(other.index)
        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)

        if other.index is not new_index:
            other = other.reindex(new_index)

        for col, series in this.iteritems():
            new_data[col] = func(series.values, other.values)

        return self._constructor(new_data, index=new_index,
                                 columns=self.columns)

    def _combine_match_columns(self, other, func, fill_value):
        # patched version of DataFrame._combine_match_columns to account for
        # NumPy circumventing __rsub__ with float64 types, e.g.: 3.0 - series,
        # where 3.0 is numpy.float64 and series is a SparseSeries. Still
        # possible for this to happen, which is bothersome

        if fill_value is not None:
            raise NotImplementedError

        new_data = {}

        union = intersection = self.columns

        if not union.equals(other.index):
            union = other.index.union(self.columns)
            intersection = other.index.intersection(self.columns)

        for col in intersection:
            new_data[col] = func(self[col], float(other[col]))

        return self._constructor(new_data, index=self.index,
                                 columns=union)

    def _combine_const(self, other, func):
        new_data = {}
        for col, series in self.iteritems():
            new_data[col] = func(series, other)

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns)

    def _reindex_index(self, index, method, copy):
        if self.index.equals(index):
            if copy:
                return self.copy()
            else:
                return self

        if len(self.index) == 0:
            return SparseDataFrame(index=index, columns=self.columns)

        indexer = self.index.get_indexer(index, method)
        mask = indexer == -1
        need_mask = mask.any()

        new_series = {}
        for col, series in self.iteritems():
            values = series.values
            new = values.take(indexer)

            if need_mask:
                np.putmask(new, mask, nan)

            new_series[col] = new

        return SparseDataFrame(new_series, index=index, columns=self.columns,
                               default_fill_value=self.default_fill_value)

    def _reindex_columns(self, columns, copy):
        # TODO: fill value handling
        sdict = dict((k, v) for k, v in self.iteritems() if k in columns)
        return SparseDataFrame(sdict, index=self.index, columns=columns,
                               default_fill_value=self.default_fill_value)

    def _rename_index_inplace(self, mapper):
        self.index = [mapper(x) for x in self.index]

    def _rename_columns_inplace(self, mapper):
        new_series = {}
        new_columns = []

        for col in self.columns:
            new_col = mapper(col)
            if new_col in new_series: # pragma: no cover
                raise Exception('Non-unique mapping!')
            new_series[new_col] = self[col]
            new_columns.append(new_col)

        self.columns = new_columns
        self._series = new_series

    def _get_raw_column(self, col):
        return self._series[col].values

    def add_prefix(self, prefix):
        f = (('%s' % prefix) + '%s').__mod__
        return self.rename(columns=f)

    def add_suffix(self, suffix):
        f = ('%s' + ('%s' % suffix)).__mod__
        return self.rename(columns=f)

    def _join_on(self, other, on, how, lsuffix, rsuffix):
        # need to implement?
        raise NotImplementedError

    def _join_index(self, other, how, lsuffix, rsuffix):
        if isinstance(other, Series):
            assert(other.name is not None)
            other = SparseDataFrame({other.name : other},
                                    default_fill_value=self.default_fill_value)

        join_index = self.index.join(other.index, how=how)

        this = self.reindex(join_index)
        other = other.reindex(join_index)

        this, other = this._maybe_rename_join(other, lsuffix, rsuffix)

        result_series = this._series
        other_series = other._series
        result_series.update(other_series)

        return self._constructor(result_series, index=join_index)

    def _maybe_rename_join(self, other, lsuffix, rsuffix):
        intersection = self.columns.intersection(other.columns)

        if len(intersection) > 0:
            if not lsuffix and not rsuffix:
                raise Exception('columns overlap: %s' % intersection)

            def lrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, lsuffix)
                return x

            def rrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, rsuffix)
                return x

            this = self.rename(columns=lrenamer)
            other = other.rename(columns=rrenamer)
        else:
            this = self

        return this, other

    def transpose(self):
        """
        Returns a DataFrame with the rows/columns switched.
        """
        return SparseDataFrame(self.values.T, index=self.columns,
                               columns=self.index,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)
    T = property(transpose)

    def count(self, axis=0, **kwds):
        return self.apply(SparseSeries.count, axis=axis)
    count.__doc__ = DataFrame.count.__doc__

    def cumsum(self, axis=0):
        """
        Return SparseDataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : SparseDataFrame
        """
        return self.apply(SparseSeries.cumsum, axis=axis)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Analogous to DataFrame.shift
        """
        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        new_series = {}
        if offset is None:
            new_index = self.index
            for col, s in self.iteritems():
                new_series[col] = s.shift(periods)
        else:
            new_index = self.index.shift(periods, offset)
            for col, s in self.iteritems():
                new_series[col] = SparseSeries(s.sp_values, index=new_index,
                                               sparse_index=s.sp_index,
                                               fill_value=s.fill_value)

        return SparseDataFrame(new_series, index=new_index,
                               columns=self.columns,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

    def apply(self, func, axis=0, broadcast=False):
        """
        Analogous to DataFrame.apply, for SparseDataFrame

        Parameters
        ----------
        func : function
            Function to apply to each column
        axis : {0, 1}
        broadcast : bool, default False
            For aggregation functions, return object of same size with values
            propagated

        Returns
        -------
        applied : Series or SparseDataFrame
        """
        if not len(self.columns):
            return self

        if isinstance(func, np.ufunc):
            new_series = {}
            for k, v in self.iteritems():
                applied = func(v)
                applied.fill_value = func(applied.fill_value)
                new_series[k] = applied
            return SparseDataFrame(new_series, index=self.index,
                                   columns=self.columns,
                                   default_fill_value=self.default_fill_value,
                                   default_kind=self.default_kind)
        else:
            if not broadcast:
                return self._apply_standard(func, axis)
            else:
                return self._apply_broadcast(func, axis)

    def fillna(self, *args, **kwargs):
        raise NotImplementedError

def stack_sparse_frame(frame):
    """
    Only makes sense when fill_value is NaN
    """
    lengths = [s.sp_index.npoints for _, s in frame.iteritems()]
    nobs = sum(lengths)

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    inds_to_concat = []
    vals_to_concat = []
    for _, series in frame.iteritems():
        if not np.isnan(series.fill_value):
            raise Exception('This routine assumes NaN fill value')

        int_index = series.sp_index.to_int_index()
        inds_to_concat.append(int_index.indices)
        vals_to_concat.append(series.sp_values)

    major_labels = np.concatenate(inds_to_concat)
    stacked_values = np.concatenate(vals_to_concat)
    index = MultiIndex(levels=[frame.index, frame.columns],
                       labels=[major_labels, minor_labels])

    lp = LongPanel(stacked_values.reshape((nobs, 1)), index=index,
                   columns=['foo'])
    return lp.sortlevel(level=0)

def _stack_sparse_info(frame):
    lengths = [s.sp_index.npoints for _, s in frame.iteritems()]
    nobs = sum(lengths)

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    inds_to_concat = []
    vals_to_concat = []
    for col in frame.columns:
        series = frame[col]

        if not np.isnan(series.fill_value):
            raise Exception('This routine assumes NaN fill value')

        int_index = series.sp_index.to_int_index()
        inds_to_concat.append(int_index.indices)
        vals_to_concat.append(series.sp_values)

    major_labels = np.concatenate(inds_to_concat)
    sparse_values = np.concatenate(vals_to_concat)

    return sparse_values, major_labels, minor_labels


def homogenize(series_dict):
    """
    Conform a set of SparseSeries (with NaN fill_value) to a common SparseIndex
    corresponding to the locations where they all have data

    Parameters
    ----------
    series_dict : dict or DataFrame

    Notes
    -----
    Using the dumbest algorithm I could think of. Should put some more thought
    into this

    Returns
    -------
    homogenized : dict of SparseSeries
    """
    index = None

    need_reindex = False

    for _, series in series_dict.iteritems():
        if not np.isnan(series.fill_value):
            raise Exception('this method is only valid with NaN fill values')

        if index is None:
            index = series.sp_index
        elif not series.sp_index.equals(index):
            need_reindex = True
            index = index.intersect(series.sp_index)

    if need_reindex:
        output = {}
        for name, series in series_dict.iteritems():
            if not series.sp_index.equals(index):
                series = series.sparse_reindex(index)

            output[name] = series
    else:
        output = series_dict

    return output

class PanelAxis(object):

    def __init__(self, cache_field):
        self.cache_field = cache_field

    def __get__(self, obj, type=None):
        return getattr(obj, self.cache_field, None)

    def __set__(self, obj, value):
        value = _ensure_index(value)

        if isinstance(value, MultiIndex):
            raise NotImplementedError

        setattr(obj, self.cache_field, value)

class SparsePanel(Panel):
    """
    Sparse version of Panel

    Parameters
    ----------
    frames : dict of DataFrame objects
    items : array-like
    major_axis : array-like
    minor_axis : array-like
    default_kind : {'block', 'integer'}, default 'block'
        Default sparse kind for converting Series to SparseSeries. Will not
        override SparseSeries passed into constructor
    default_fill_value : float
        Default fill_value for converting Series to SparseSeries. Will not
        override SparseSeries passed in

    Notes
    -----
    """
    items = PanelAxis('_items')
    major_axis = PanelAxis('_major_axis')
    minor_axis = PanelAxis('_minor_axis')

    ndim = 3

    def __init__(self, frames, items=None, major_axis=None, minor_axis=None,
                 default_fill_value=nan, default_kind='block'):
        assert(isinstance(frames, dict))

        self.default_fill_value = fill_value = default_fill_value
        self.default_kind = kind = default_kind

        # pre-filter, if necessary
        if items is None:
            items = Index(sorted(frames.keys()))
        items = _ensure_index(items)

        (clean_frames,
         major_axis,
         minor_axis) = _convert_frames(frames, major_axis,
                                       minor_axis, kind=kind,
                                       fill_value=fill_value)

        self._frames = clean_frames

        # do we want to fill missing ones?
        for item in items:
            if item not in clean_frames:
                raise Exception('column %s not found in data' % item)

        self.items = items
        self.major_axis = major_axis
        self.minor_axis = minor_axis

    def _consolidate_inplace(self): # pragma: no cover
        # do nothing when DataFrame calls this method
        pass

    @classmethod
    def from_dict(cls, data):
        """
        Analogous to Panel.from_dict
        """
        return SparsePanel(data)

    def to_dense(self):
        """
        Convert SparsePanel to (dense) Panel

        Returns
        -------
        dense : Panel
        """
        return Panel(self.values, self.items, self.major_axis,
                     self.minor_axis)

    @property
    def values(self):
        # return dense values
        return np.array([self._frames[item].values
                         for item in self.items])

    def __getitem__(self, key):
        """
        """
        return self._frames[key]

    def __setitem__(self, key, value):
        if isinstance(value, DataFrame):
            value = value.reindex(index=self.major_axis,
                                  columns=self.minor_axis)
            if not isinstance(value, SparseDataFrame):
                value = value.to_sparse(fill_value=self.default_fill_value,
                                        kind=self.default_kind)
        else:
            raise ValueError('only DataFrame objects can be set currently')

        self._frames[key] = value

        if key not in self.items:
            self.items = Index(list(self.items) + [key])

    def __delitem__(self, key):
        loc = self.items.get_loc(key)
        indices = range(loc) + range(loc + 1, len(self.items))
        self.items = self.items[indices]
        del self._frames[key]

    def __getstate__(self):
        # pickling
        return (self._frames, _pickle_array(self.items),
                _pickle_array(self.major_axis), _pickle_array(self.minor_axis),
                self.default_fill_value, self.default_kind)

    def __setstate__(self, state):
        frames, items, major, minor, fv, kind = state

        self.default_fill_value = fv
        self.default_kind = kind
        self.items = _unpickle_array(items)
        self.major_axis = _unpickle_array(major)
        self.minor_axis = _unpickle_array(minor)
        self._frames = frames

    def copy(self):
        """
        Make a (shallow) copy of the sparse panel

        Returns
        -------
        copy : SparsePanel
        """
        return SparsePanel(self._frames.copy(), items=self.items,
                           major_axis=self.major_axis,
                           minor_axis=self.minor_axis,
                           default_fill_value=self.default_fill_value,
                           default_kind=self.default_kind)

    def to_long(self, filter_observations=True):
        """
        Convert SparsePanel to (dense) LongPanel

        Returns
        -------
        lp : LongPanel
        """
        if not filter_observations:
            raise Exception('filter_observations=False not supported for '
                            'SparsePanel.to_long')

        I, N, K = self.shape
        counts = np.zeros(N * K, dtype=int)

        d_values = {}
        d_indexer = {}

        for item in self.items:
            frame = self[item]

            values, major, minor = _stack_sparse_info(frame)

            # values are stacked column-major
            indexer = minor * N + major
            counts.put(indexer, counts.take(indexer) + 1) # cuteness

            d_values[item] = values
            d_indexer[item] = indexer

        # have full set of observations for each item
        mask = counts == I

        # for each item, take mask values at index locations for those sparse
        # values, and use that to select values
        values = np.column_stack([d_values[item][mask.take(d_indexer[item])]
                                  for item in self.items])

        inds, = mask.nonzero()

        # still column major
        major_labels = inds % N
        minor_labels = inds // N

        index = MultiIndex(levels=[self.major_axis, self.minor_axis],
                           labels=[major_labels, minor_labels])

        lp = LongPanel(values, index=index, columns=self.items)
        return lp.sortlevel(level=0)

    def reindex(self, major=None, items=None, minor=None, major_axis=None,
                minor_axis=None, copy=False):
        """
        Conform / reshape panel axis labels to new input labels

        Parameters
        ----------
        major : array-like, default None
        items : array-like, default None
        minor : array-like, default None
        copy : boolean, default False
            Copy underlying SparseDataFrame objects

        Returns
        -------
        reindexed : SparsePanel
        """
        major = _mut_exclusive(major, major_axis)
        minor = _mut_exclusive(minor, minor_axis)

        if None == major == items == minor:
            raise ValueError('Must specify at least one axis')

        major = self.major_axis if major is None else major
        minor = self.minor_axis if minor is None else minor

        if items is not None:
            new_frames = {}
            for item in items:
                if item in self._frames:
                    new_frames[item] = self._frames[item]
                else:
                    raise Exception('Reindexing with new items not yet '
                                    'supported')
        else:
            new_frames = self._frames

        if copy:
            new_frames = dict((k, v.copy()) for k, v in new_frames.iteritems())

        return SparsePanel(new_frames, items=items,
                           major_axis=major,
                           minor_axis=minor,
                           default_fill_value=self.default_fill_value,
                           default_kind=self.default_kind)

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            new_frames = dict((k, func(v, other))
                              for k, v in self.iterkv())
            return self._new_like(new_frames)

    def _combineFrame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._get_axis_number(axis)

        other = other.reindex(index=index, columns=columns)

        if axis == 0:
            new_values = func(self.values, other.values)
        elif axis == 1:
            new_values = func(self.values.swapaxes(0, 1), other.values.T)
            new_values = new_values.swapaxes(0, 1)
        elif axis == 2:
            new_values = func(self.values.swapaxes(0, 2), other.values)
            new_values = new_values.swapaxes(0, 2)

        # TODO: make faster!
        new_frames = {}
        for item, item_slice in zip(self.items, new_values):
            old_frame = self[item]
            ofv = old_frame.default_fill_value
            ok = old_frame.default_kind
            new_frames[item] = SparseDataFrame(item_slice,
                                               index=self.major_axis,
                                               columns=self.minor_axis,
                                               default_fill_value=ofv,
                                               default_kind=ok)

        return self._new_like(new_frames)

    def _new_like(self, new_frames):
        return SparsePanel(new_frames, self.items, self.major_axis,
                           self.minor_axis,
                           default_fill_value=self.default_fill_value,
                           default_kind=self.default_kind)

    def _combinePanel(self, other, func):
        # if isinstance(other, LongPanel):
        #     other = other.to_wide()
        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it

        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        new_frames = {}
        for item in items:
            new_frames[item] = func(this[item], other[item])

        # maybe unnecessary
        new_default_fill = func(self.default_fill_value,
                                other.default_fill_value)

        return SparsePanel(new_frames, items, major, minor,
                           default_fill_value=new_default_fill,
                           default_kind=self.default_kind)

    def major_xs(self, key):
        """
        Return slice of panel along major axis

        Parameters
        ----------
        key : object
            Major axis label

        Returns
        -------
        y : DataFrame
            index -> minor axis, columns -> items
        """
        slices = dict((k, v.xs(key)) for k, v in self.iterkv())
        return DataFrame(slices, index=self.minor_axis, columns=self.items)

    def minor_xs(self, key):
        """
        Return slice of panel along minor axis

        Parameters
        ----------
        key : object
            Minor axis label

        Returns
        -------
        y : SparseDataFrame
            index -> major axis, columns -> items
        """
        slices = dict((k, v[key]) for k, v in self.iterkv())
        return SparseDataFrame(slices, index=self.major_axis,
                               columns=self.items,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

SparseWidePanel = SparsePanel

def _convert_frames(frames, index, columns, fill_value=nan, kind='block'):
    from pandas.core.panel import _get_combined_index, _get_combined_columns
    output = {}
    for item, df in frames.iteritems():
        if not isinstance(df, SparseDataFrame):
            df = SparseDataFrame(df, default_kind=kind,
                                 default_fill_value=fill_value)

        output[item] = df

    if index is None:
        index = _get_combined_index(output)
    if columns is None:
        columns = _get_combined_columns(output)

    index = _ensure_index(index)
    columns = _ensure_index(columns)

    for item, df in output.iteritems():
        if not (df.index.equals(index) and df.columns.equals(columns)):
            output[item] = df.reindex(index=index, columns=columns)

    return output, index, columns

