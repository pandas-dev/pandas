"""
Data structures for sparse float data. Life is made simpler by dealing only with
float64 data
"""

from numpy import nan, ndarray
import numpy as np

import operator

from pandas.core.index import Index, NULL_INDEX
from pandas.core.panel import _get_combined_index, _get_combined_columns
from pandas.core.series import Series, TimeSeries
from pandas.core.frame import DataFrame, extract_index, try_sort
from pandas.core.matrix import DataMatrix
import pandas.core.common as common

from pandas.lib.sparse import BlockIndex, IntIndex
import pandas.lib.sparse as splib

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
_MIRROR_OPS = {
    'add' : '__radd__',
    'sub' : '__rsub__',
    'div' : '__rdiv__',
    'truediv' : '__rdiv__',
    'mul' : '__rmul__',
}

def _sparse_op_wrap(op, name):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
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
                result = op(this.sp_values, other.sp_values)
                result_index = self.sp_index
            else:
                result, result_index = sparse_op(this, other)

            try:
                fill_value = op(this.fill_value, other.fill_value)
            except ZeroDivisionError:
                fill_value = nan

            return SparseSeries(result, index=new_index,
                                sparse_index=result_index,
                                fill_value=fill_value)

        elif isinstance(other, SparseDataFrame):
            reverse_op = _MIRROR_OPS.get(name)
            if reverse_op is None: # pragma: no cover
                raise Exception('Cannot do %s op, sorry!' % name)
            return getattr(other, reverse_op)(self)
        elif np.isscalar(other):
            new_fill_value = op(np.float64(self.fill_value),
                                np.float64(other))

            return SparseSeries(op(self.sp_values, other),
                                index=self.index,
                                sparse_index=self.sp_index,
                                fill_value=new_fill_value)
        else:
            raise Exception('operation with %s not supported' % type(other))

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

    Notes
    -----
    SparseSeries objects are immutable via the typical Python means. If you must
    change values, convert to dense, make your changes, then convert back to
    sparse
    """
    __array_priority__ = 15

    sp_index = None
    fill_value = None

    def __new__(cls, data, index=None, sparse_index=None, kind='block',
                fill_value=None, copy=False):

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

    @property
    def _constructor(self):
        def make_sp_series(data, index=None):
            return SparseSeries(data, index=index, fill_value=self.fill_value,
                                kind=self.kind)

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
        self.sp_index = getattr(obj, 'sp_index', None)
        self.fill_value = getattr(obj, 'fill_value', None)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))

        subclass_state = (self.index, self.fill_value, self.sp_index)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)

        index, fill_value, sp_index = own_state

        self.sp_index = sp_index
        self.fill_value = fill_value
        self.index = index

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
    __div__ = _sparse_op_wrap(operator.div, 'div')
    __truediv__ = _sparse_op_wrap(operator.truediv, 'truediv')
    __pow__ = _sparse_op_wrap(operator.pow, 'pow')

    # reverse operators
    __radd__ = _sparse_op_wrap(operator.add, '__radd__')
    __rmul__ = _sparse_op_wrap(operator.mul, '__rmul__')
    __rsub__ = _sparse_op_wrap(lambda x, y: y - x, '__rsub__')
    __rdiv__ = _sparse_op_wrap(lambda x, y: y / x, '__rdiv__')
    __rtruediv__ = _sparse_op_wrap(lambda x, y: y / x, '__rtruediv__')
    __rpow__ = _sparse_op_wrap(lambda x, y: y ** x, '__rpow__')

    # Inplace operators
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __idiv__ = __div__
    __ipow__ = __pow__

    @property
    def values(self):
        output = np.empty(len(self), dtype=np.float64)
        int_index = self.sp_index.to_int_index()
        output.fill(self.fill_value)
        output.put(int_index.indices, self)
        return output

    @property
    def sp_values(self):
        return np.asarray(self)

    def __getitem__(self, key):
        """

        """
        try:
            return self._get_val_at(self.index.indexMap[key])

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
        return self._constructor(dataSlice, index=new_index)

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
        y : ndarray
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

    def __getslice__(self, i, j):
        return self._constructor(self.values[i:j], index=self.index[i:j])

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
            return Series(self.sp_values, index=index)
        else:
            return Series(self.values, index=self.index)

    def astype(self, dtype):
        # HACK?
        return self.copy()

    def copy(self):
        values = self.sp_values.copy()
        return SparseSeries(values, index=self.index,
                            sparse_index=self.sp_index,
                            fill_value=self.fill_value)

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

        notmask = -mask
        if notmask.any():
            np.putmask(new_values, notmask, nan)

        return SparseSeries(new_values, index=new_index,
                            fill_value=self.fill_value)

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

        # indexer = self.sp_index.get_indexer(new_index)

        # new_values = self.sp_values.take(indexer)
        # new_values[indexer == -1] = self.fill_value

        return SparseSeries(new_values, index=self.index,
                            sparse_index=new_index,
                            fill_value=self.fill_value)

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

    def valid(self):
        """
        Analogous to Series.valid
        """
        # TODO: make more efficient
        dense_valid = self.to_dense().valid()
        return dense_valid.to_sparse(fill_value=self.fill_value)

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

    def __array_wrap__(self, result):
        return SparseDataFrame(result, index=self.index, columns=self.columns,
                               kind=self.default_kind,
                               default_fill_value=self.default_fill_value)

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
                                          fill_value=self.default_fill_value,
                                          copy=True)

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

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        from cStringIO import StringIO

        buf = StringIO()
        if len(self.index) < 500 and len(self.columns) < 10:
            self.toString(buf=buf)
        else:
            self.info(buf=buf, verbose=False)

        return buf.getvalue()

    def to_dense(self):
        """
        Convert to dense DataFrame

        Returns
        -------
        df : DataFrame
        """
        data = dict((k, v.to_dense()) for k, v in self.iteritems())
        return DataFrame(data, index=self.index)

    def copy(self):
        """
        Make a deep copy of this frame
        """
        return SparseDataFrame(self._series, index=self.index,
                               columns=self.columns,
                               default_fill_value=self.default_fill_value,
                               kind=self.default_kind)

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

    def _insert_item(self, key, value):
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
            loc = self._get_insert_loc(key)
            self._insert_column_index(key, loc)

    def _combine_match_columns(self, other, func):
        # patched version of DataFrame._combine_match_columns to account for
        # NumPy circumventing __rsub__ with float64 types, e.g.: 3.0 - series,
        # where 3.0 is numpy.float64 and series is a SparseSeries. Still
        # possible for this to happen, which is bothersome

        new_data = {}

        union = intersection = self.columns

        if not union.equals(other.index):
            union = other.index.union(self.columns)
            intersection = other.index.intersection(self.columns)

        for col in intersection:
            new_data[col] = func(self[col], float(other[col]))

        return self._constructor(new_data, index=self.index,
                                 columns=union)

    def _reindex_index(self, index, method):
        if self.index.equals(index):
            return self.copy()

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

        return SparseDataFrame(new_series, index=index, columns=self.columns,
                               default_fill_value=self.default_fill_value)

    def _reindex_columns(self, columns):
        # TODO: fill value handling
        sdict = dict((k, v) for k, v in self.iteritems() if k in columns)
        return SparseDataFrame(sdict, index=self.index, columns=columns,
                               default_fill_value=self.default_fill_value)

from pandas.core.panel import WidePanel

# from line_profiler import LineProfiler
# prof = LineProfiler()

def stack_sparse_frame(frame):
    """
    Only makes sense when fill_value is NaN
    """
    from pandas.core.panel import LongPanelIndex, LongPanel

    lengths = [s.sp_index.npoints for _, s in frame.iteritems()]
    nobs = sum(lengths)

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    # need to create
    major_labels = np.empty(nobs, dtype=int)
    stacked_values = np.empty(nobs, dtype=np.float64)

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
    index = LongPanelIndex(frame.index, frame.columns,
                           major_labels, minor_labels)

    lp = LongPanel(stacked_values.reshape((nobs, 1)), ['foo'], index)
    return lp.sort('major')

def homogenize(series_dict):
    """
    Conform a set of SparseSeries (with NaN fill_value) to a common SparseIndex
    corresponding to the locations where they all have data

    Notes
    -----
    Using the dumbest algorithm I could think of. Should put some more thought
    into this

    """
    index = None

    need_reindex = False

    for _, series in series_dict.iteritems():
        if index is None:
            index = series.sp_index
        elif not series.sp_index.equals(index):
            need_reindex = True
            index = index.intersect(series.sp_index)

    if need_reindex:
        output = {}
        for name, series in series_dict.iteritems():
            output[name] = series.sparse_reindex(index)
    else:
        output = series_dict

    return output

class SparseWidePanel(WidePanel):
    """



    """
    def __init__(self, frames, items=None, major_axis=None, minor_axis=None):
        assert(isinstance(frames, dict))

        self.frames = frames

        self.items = Index(sorted(frames))
        self.major_axis = foo

    @classmethod
    def from_dict(cls):
        pass

    @property
    def values(self):
        # return dense values
        pass

    def __getitem__(self, key):
        """
        """
        pass

    def __setitem__(self, key, value):
        pass

    def __delitem__(self, key):
        pass

    def pop(self, key):
        pass

    #----------------------------------------------------------------------
    # pickling

    def __getstate__(self):
        pass

    def __setstate__(self, state):
        pass

    @property
    def values(self):
        """

        """
        pass

    @classmethod
    def from_dict(cls, data, intersect=False):
        pass

    def copy(self):
        pass

    def to_long(self):
        pass

    def reindex(self, major=None, items=None, minor=None, method=None):
        result = self

        if major is not None:
            result = result._reindex_axis(major, method, 1)

        if minor is not None:
            result = result._reindex_axis(minor, method, 2)

        if items is not None:
            result = result._reindex_axis(items, method, 0)

        if result is self:
            raise ValueError('Must specify at least one axis')

        return result

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            newValues = func(self.values, other)

            return WidePanel(newValues, self.items, self.major_axis,
                             self.minor_axis)

    def _combineFrame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._get_axis_number(axis)

        other = other.reindex(index=index, columns=columns)

        if axis == 0:
            newValues = func(self.values, other.values)
        elif axis == 1:
            newValues = func(self.values.swapaxes(0, 1), other.values.T)
            newValues = newValues.swapaxes(0, 1)
        elif axis == 2:
            newValues = func(self.values.swapaxes(0, 2), other.values)
            newValues = newValues.swapaxes(0, 2)

        return WidePanel(newValues, self.items, self.major_axis,
                         self.minor_axis)

    def _combinePanel(self, other, func):
        if isinstance(other, LongPanel):
            other = other.to_wide()

        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it

        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        result_values = func(this.values, other.values)

        return WidePanel(result_values, items, major, minor)

    def major_xs(self, key):
        """
        Parameters
        ----------

        Returns
        -------
        y : DataMatrix
            index -> minor axis, columns -> items
        """
        try:
            loc = self.major_axis.indexMap[key]
        except KeyError:
            raise KeyError('%s not contained in major axis!' % key)

        mat = np.array(self.values[:, loc, :].T)
        return DataMatrix(mat, index=self.minor_axis, columns=self.items)

    def minor_xs(self, key):
        """
        Parameters
        ----------

        Returns
        -------
        y : DataMatrix
            index -> major axis, columns -> items
        """
        try:
            loc = self.minor_axis.indexMap[key]
        except KeyError:
            raise KeyError('%s not contained in minor axis!' % key)

        mat = np.array(self.values[:, :, loc].T)
        return DataMatrix(mat, index=self.major_axis, columns=self.items)

