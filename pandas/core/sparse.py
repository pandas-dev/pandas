"""
Data structures for sparse float data. Life is made simpler by dealing only with
float64 data
"""

from numpy import nan
import numpy as np

import operator

from pandas.core.index import Index, NULL_INDEX
from pandas.core.series import Series, TimeSeries, remove_na
from pandas.core.frame import DataFrame, extract_index, try_sort
import pandas.core.common as common

from pandas.lib.sparse import BlockIndex, IntIndex
import pandas.lib.sparse as splib
import pandas.lib.tseries as tseries

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
    vector : SparseVector
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
_MIRROR_OPS = {
    'add' : '__radd__',
    'sub' : '__rsub__',
    'div' : '__rdiv__',
    'truediv' : '__rdiv__',
    'mul' : '__rmul__',
}

def _sparse_op_wrap(name):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        py_op = getattr(operator, name)

        if isinstance(other, SparseSeries):
            if np.isnan(self.fill_value):
                sparse_op = lambda a, b: _sparse_nanop(a, b, name)
            else:
                sparse_op = lambda a, b: _sparse_fillop(a, b, name)

            new_index = self.index + other.index
            if self.index.equals(new_index):
                this = self
            else:
                this = self.reindex(new_index)

            if not other.index.equals(new_index):
                other = other.reindex(new_index)

            if self.sp_index.equals(other.sp_index):
                result = py_op(this.sp_values, other.sp_values)
                result_index = self.sp_index
            else:
                result, result_index = sparse_op(this, other)

            try:
                fill_value = py_op(this.fill_value, other.fill_value)
            except ZeroDivisionError:
                fill_value = nan

            return SparseSeries(result, index=new_index,
                                sparse_index=result_index,
                                fill_value=fill_value)

        elif isinstance(other, SparseDataFrame):
            reverse_op = _MIRROR_OPS.get(name)
            if reverse_op is None: # pragma: no cover
                raise Exception('Cannot do %s op, sorry!')
            return getattr(other, reverse_op)(self)
        elif np.isscalar(other):
            return SparseSeries(py_op(self.sp_values, other),
                                index=self.index,
                                sparse_index=self.sp_index,
                                fill_value=py_op(self.fill_value, other))

    wrapper.__name__ = name
    return wrapper

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

"""
Notes.

- Not sure if subclassing is the way to go, could just lead to trouble on down
  the road (e.g. if putting a SparseSeries in a regular DataFrame...). But on
  the other hand you do get many things for free...

- Will need to "disable" a number of methods?
"""

class SparseSeries(Series):
    """
    Data structure for labeled, sparse floating point data

    Parameters
    ----------
    data : {array-like, Series, SparseSeries, dict}
    kind : {'block', 'integer'}
    fill_value : float
        Defaults to NaN (code for missing)
    """
    sp_index = None
    fill_value = None

    def __new__(cls, data, index=None, sparse_index=None, kind='block',
                fill_value=None, copy=False):

        if isinstance(data, SparseSeries):
            if index is None:
                index = data.index

            if fill_value is None:
                fill_value = data.fill_value

            if index is not None:
                assert(len(index) == len(data))

            sparse_index = data.sp_index
            values = np.asarray(data)
        elif isinstance(data, (Series, dict)):
            if fill_value is None:
                fill_value = nan

            data = Series(data)
            if index is None:
                index = data.index
            values, sparse_index = make_sparse(data, kind=kind,
                                               fill_value=fill_value)
        elif np.isscalar(data): # pragma: no cover
            raise Exception('not supported yet')
        else:
            if fill_value is None:
                fill_value = nan

            # array-like
            if sparse_index is None:
                values, sparse_index = make_sparse(data, kind=kind,
                                                   fill_value=fill_value)
            else:
                values = data
                assert(len(values) == sparse_index.npoints)

        if index is None:
            index = Index(np.arange(sparse_index.length))
        elif not isinstance(index, Index):
            index = Index(index)

        # Create array, do *not* copy data by default
        subarr = np.array(values, dtype=np.float64, copy=copy)

        if index.is_all_dates():
            cls = SparseTimeSeries

        # Change the class of the array to be the subclass type.
        subarr = subarr.view(cls)
        subarr.sp_index = sparse_index
        subarr.fill_value = fill_value
        subarr.index = index
        return subarr

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)
        self.sp_index = getattr(obj, 'sp_index', None)
        self.fill_value = getattr(obj, 'fill_value', None)

    def __len__(self):
        return self.sp_index.length

    def __repr__(self):
        series_rep = Series.__repr__(self)
        rep = '%s\n%s' % (series_rep, repr(self.sp_index))
        return rep

    # Arithmetic operators

    __add__ = _sparse_op_wrap('add')
    __sub__ = _sparse_op_wrap('sub')
    __mul__ = _sparse_op_wrap('mul')
    __div__ = _sparse_op_wrap('div')
    __truediv__ = _sparse_op_wrap('div')
    __pow__ = _sparse_op_wrap('pow')

    # Inplace operators
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __idiv__ = __div__
    __ipow__ = __pow__

    @property
    def values(self):
        output = np.empty(self.sp_index.length, dtype=np.float64)
        int_index = self.sp_index.to_int_index()
        output.fill(self.fill_value)
        output.put(int_index.indices, self)
        return output

    @property
    def sp_values(self):
        return np.asarray(self)

    def to_dense(self, sparse_only=False):
        """
        Convert SparseSeries to (dense) Series
        """
        if sparse_only:
            int_index = self.sp_index.to_int_index()
            index = self.index.take(int_index.indices)
            return Series(self.sp_values, index=index)
        else:
            return Series(self.values, index=self.index)

    def astype(self, dtype):
        # HACK?
        return self.copy()

    def copy(self):
        values = self.sp_values.copy()
        return SparseSeries(values, index=self.index,
                            sparse_index=self.sp_index)

    def reindex(self, new_index, method=None):
        """
        Conform SparseSeries to new Index

        See Series.reindex docstring for general behavior

        Returns
        -------
        reindexed : SparseSeries
        """
        if not isinstance(new_index, Index):
            new_index = Index(new_index)

        if self.index.equals(new_index):
            return self.copy()

        if len(self.index) == 0:
            # FIXME: inelegant / slow
            values = np.empty(len(new_index), dtype=np.float64)
            values.fill(nan)
            return SparseSeries(values, index=new_index,
                                fill_value=self.fill_value)

        values = self.values
        indexer, mask = self.index.get_indexer(new_index, method=method)
        new_values = values.take(indexer)

        # TODO: always use NaN here?
        notmask = -mask
        if notmask.any():
            np.putmask(new_values, notmask, nan)

        return SparseSeries(new_values, index=new_index,
                            fill_value=self.fill_value)

    def take(self, indices):
        """
        Sparse-compatible version of ndarray.take

        Returns
        -------
        y : SparseSeries
        """
        pass

    def put(self, indices, values):
        """
        Sparse-compatible version of ndarray.put

        Returns
        -------
        y : SparseSeries
        """
        pass

    def count(self):
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
        Sum of non-null values
        """
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        if self._null_fill_value:
            return sp_sum
        else:
            nsparse = self.sp_index.npoints
            return sp_sum + self.fill_value * nsparse

    def mean(self, axis=None, dtype=None, out=None):
        """
        Mean of non-null values
        """
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        ct = len(valid_vals)

        if self._null_fill_value:
            return sp_sum / ct
        else:
            nsparse = self.sp_index.npoints
            return (sp_sum + self.fill_value * nsparse) / (ct + nsparse)

class SparseTimeSeries(SparseSeries, TimeSeries):
    pass

class SparseDataFrame(DataFrame):
    """
    DataFrame containing sparse floating point data in the form of SparseSeries
    objects
    """
    _columns = None

    def __init__(self, data=None, index=None, columns=None, kind='block',
                 default_fill_value=None):
        if default_fill_value is None:
            default_fill_value = nan

        self.default_kind = kind
        self.default_fill_value = default_fill_value

        DataFrame.__init__(self, data, index=index, columns=columns,
                           dtype=None)

    @property
    def _constructor(self):
        return SparseDataFrame

    def _init_dict(self, data, index, columns, dtype):
        # pre-filter out columns if we passed it
        if columns is not None:
            if not isinstance(columns, Index):
                columns = Index(columns)

            data = dict((k, v) for k, v in data.iteritems() if k in columns)
        else:
            columns = Index(try_sort(data.keys()))

        index = extract_index(data, index)

        sp_maker = lambda x: SparseSeries(x, index=index,
                                          kind=self.default_kind,
                                          fill_value=self.default_fill_value)

        sdict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Forces alignment and copies data
                v = v.reindex(index)

                if not isinstance(v, SparseSeries):
                    v = sp_maker(v)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, nan) for i in index]

                v = sp_maker(v).copy()
            sdict[k] = v

        # TODO: figure out how to handle this case, all nan's?
        # add in any other columns we want to have (completeness)
        for c in columns:
            if c not in sdict:
                sdict[c] = sp_maker([])

        return sdict, columns, index

    def _insert_item(self, key, value):
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                cleanSeries = value.reindex(self.index)
                if not isinstance(value, SparseSeries):
                    cleanSeries = SparseSeries(cleanSeries)
            else:
                cleanSeries = Series(value, index=self.index)

            self._series[key] = cleanSeries
        # Scalar
        else:
            self._series[key] = SparseSeries(value, index=self.index)

        if key not in self.columns:
            loc = self._get_insert_loc(key)
            self._insert_column_index(key, loc)

    def to_dense(self):
        """
        Convert to dense DataFrame

        Returns
        -------
        df : DataFrame
        """
        data = dict((k, v.to_dense()) for k, v in self.iteritems())
        return DataFrame(data, index=self.index)

    def _reindex_index(self, index, method):
        if self.index.equals(index):
            return self.copy()

        if not isinstance(index, Index):
            index = Index(index)

        if len(self.index) == 0:
            return SparseDataFrame(index=index, columns=self.columns)

        indexer, mask = common.get_indexer(self.index, index, method)
        notmask = -mask
        need_mask = notmask.any()

        new_series = {}
        for col, series in self.iteritems():
            values = series.values
            new = values.take(indexer)

            if need_mask:
                new[notmask] = nan

            new_series[col] = new

        return SparseDataFrame(new_series, index=index, columns=self.columns)

class SparsePanel(object):
    """

    """

    def __init__(self, frames):
        pass


