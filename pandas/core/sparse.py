"""
Data structures for sparse float data. Life is made simpler by dealing only with
float64 data
"""

from numpy import nan
import numpy as np

import operator

from pandas.core.index import Index, NULL_INDEX
from pandas.core.series import Series, TimeSeries
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

def to_sparse_series(series, kind='block', fill_value=nan):
    sp_values, sp_index = make_sparse(series, kind=kind, fill_value=fill_value)
    return SparseSeries(sp_values, index=series.index, sparse_index=sp_index,
                        fill_value=fill_value)

#-------------------------------------------------------------------------------
# Wrapper function for Series arithmetic methods
_MIRROR_OPS = {
    '__add__' : '__radd__',
    '__sub__' : '__rsub__',
    '__div__' : '__rdiv__',
    '__truediv__' : '__rdiv__',
    '__mul__' : '__rmul__',
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
            if reverse_op is None:
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

    def __new__(cls, data, index=None, sparse_index=None,
                kind='block', fill_value=None, copy=False):

        if isinstance(data, SparseSeries):
            if index is None:
                index = data.index

            if fill_value is None:
                sparse_index = data.fill_value

            if index is not None:
                assert(len(index) == data.length)

            values = np.asarray(data)
        elif isinstance(data, (Series, dict)):
            if fill_value is None:
                fill_value = nan

            data = Series(data)
            if index is None:
                index = data.index
            values, sparse_index = make_sparse(data, kind=kind,
                                               fill_value=fill_value)
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
        subarr = np.array(values, dtype=np.float64, copy=False)

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

    def to_dense(self):
        """
        Convert SparseSeries to (dense) Series
        """
        return Series(self.values, index=self.index)

    def astype(self, dtype):
        # HACK
        return self.copy()

    def copy(self):
        values = self.sp_values.copy()
        return SparseSeries(values, index=self.index,
                            sparse_index=self.sp_index)

    def reindex(self, new_index):
        return SparseSeries(self.to_dense().reindex(new_index))

        if self.index.equals(new_index):
            return self.copy()

        if not isinstance(new_index, Index):
            new_index = Index(new_index)

        if len(self.index) == 0:
            return Series(nan, index=new_index)

        indexer, mask = tseries.getFillVec(self.index, new_index,
                                           self.index.indexMap,
                                           new_index.indexMap)

    def take(self, indices):
        pass

    def put(self, indices, values):
        pass

class SparseTimeSeries(SparseSeries, TimeSeries):
    pass

class SparseDataFrame(DataFrame):
    """
    DataFrame containing sparse floating point data in the form of SparseSeries
    objects
    """
    _columns = None

    def __init__(self, data=None, index=None, columns=None, kind='block',
                 fill_value=None):
        self.kind = kind
        self.fill_value = fill_value
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

        sdict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Forces alignment and copies data
                v = v.reindex(index)
                if not isinstance(v, SparseSeries):
                    v = to_sparse_series(v, kind=self.kind,
                                         fill_value=self.fill_value)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, nan) for i in index]

                v = SparseSeries(v, index=index, kind=self.kind).copy()

            sdict[k] = v

        # TODO: figure out how to handle this case, all nan's?
        # add in any other columns we want to have (completeness)
        for c in columns:
            if c not in sdict:
                sdict[c] = SparseSeries([], index=index,
                                        fill_value=self.fill_value)

        return sdict, columns, index

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


