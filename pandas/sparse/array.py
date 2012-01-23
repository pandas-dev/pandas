"""
SparseArray data structure
"""

# pylint: disable=E1101,E1103,W0231

from numpy import nan, ndarray
import numpy as np

import operator
import pandas.core.common as com

from pandas.util import py3compat

from pandas._sparse import BlockIndex, IntIndex
import pandas._sparse as splib
import pandas._engines as _gin


def _sparse_op_wrap(op, name):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        if isinstance(other, np.ndarray):
            assert(len(self) == len(other))
            if not isinstance(other, SparseArray):
                other = SparseArray(other, fill_value=self.fill_value)
            return _sparse_array_op(self, other, op, name)
        elif np.isscalar(other):
            new_fill_value = op(np.float64(self.fill_value),
                                np.float64(other))

            return SparseArray(op(self.sp_values, other),
                               sparse_index=self.sp_index,
                               fill_value=new_fill_value)
        else: # pragma: no cover
            raise TypeError('operation with %s not supported' % type(other))

    wrapper.__name__ = name
    return wrapper

def _sparse_array_op(left, right, op, name):
    if np.isnan(left.fill_value):
        sparse_op = lambda a, b: _sparse_nanop(a, b, name)
    else:
        sparse_op = lambda a, b: _sparse_fillop(a, b, name)

    if left.sp_index.equals(right.sp_index):
        result = op(left.sp_values, right.sp_values)
        result_index = left.sp_index
    else:
        result, result_index = sparse_op(left, right)

    try:
        fill_value = op(left.fill_value, right.fill_value)
    except ZeroDivisionError:
        fill_value = nan

    return SparseArray(result, sparse_index=result_index,
                       fill_value=fill_value)

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


class SparseArray(np.ndarray):
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

    def __new__(cls, data, sparse_index=None, kind='integer', fill_value=None,
                copy=False):

        is_sparse_array = isinstance(data, SparseArray)
        if fill_value is None:
            if is_sparse_array:
                fill_value = data.fill_value
            else:
                fill_value = nan

        if is_sparse_array:
            sparse_index = data.sp_index
            values = np.asarray(data)
        else:
            # array-like
            if sparse_index is None:
                values, sparse_index = make_sparse(data, kind=kind,
                                                   fill_value=fill_value)
            else:
                values = data
                assert(len(values) == sparse_index.npoints)

        # Create array, do *not* copy data by default
        if copy:
            subarr = np.array(values, dtype=np.float64, copy=True)
        else:
            subarr = np.asarray(values, dtype=np.float64)

        # Change the class of the array to be the subclass type.
        output = subarr.view(cls)
        output.sp_index = sparse_index
        output.fill_value = np.float64(fill_value)
        return output

    @property
    def _constructor(self):
        return lambda x: SparseArray(x, fill_value=self.fill_value,
                                     kind=self.kind)

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
        self.sp_index = getattr(obj, 'sp_index', None)
        self.fill_value = getattr(obj, 'fill_value', None)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))
        subclass_state = self.fill_value, self.sp_index
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)

        fill_value, sp_index = own_state[:2]
        self.sp_index = sp_index
        self.fill_value = fill_value

    def __len__(self):
        return self.sp_index.length

    def __repr__(self):
        return '%s\n%s' % (np.ndarray.__repr__(self),
                           repr(self.sp_index))

    # Arithmetic operators

    __add__ = _sparse_op_wrap(operator.add, 'add')
    __sub__ = _sparse_op_wrap(operator.sub, 'sub')
    __mul__ = _sparse_op_wrap(operator.mul, 'mul')
    __truediv__ = _sparse_op_wrap(operator.truediv, 'truediv')
    __floordiv__ = _sparse_op_wrap(operator.floordiv, 'floordiv')
    __pow__ = _sparse_op_wrap(operator.pow, 'pow')

    # reverse operators
    __radd__ = _sparse_op_wrap(operator.add, 'add')
    __rsub__ = _sparse_op_wrap(lambda x, y: y - x, 'rsub')
    __rmul__ = _sparse_op_wrap(operator.mul, 'mul')
    __rtruediv__ = _sparse_op_wrap(lambda x, y: y / x, 'rtruediv')
    __rfloordiv__ = _sparse_op_wrap(lambda x, y: y // x, 'rfloordiv')
    __rpow__ = _sparse_op_wrap(lambda x, y: y ** x, 'rpow')

    def disable(self, other):
        raise NotImplementedError('inplace binary ops not supported')
    # Inplace operators
    __iadd__ = disable
    __isub__ = disable
    __imul__ = disable
    __itruediv__ = disable
    __ifloordiv__ = disable
    __ipow__ = disable

    # Python 2 division operators
    if not py3compat.PY3:
        __div__ = _sparse_op_wrap(operator.div, 'div')
        __rdiv__ = _sparse_op_wrap(lambda x, y: y / x, '__rdiv__')
        __idiv__ = disable

    @property
    def values(self):
        """
        Dense values
        """
        output = np.empty(len(self), dtype=np.float64)
        int_index = self.sp_index.to_int_index()
        output.fill(self.fill_value)
        output.put(int_index.indices, self)
        return output

    @property
    def sp_values(self):
        # caching not an option, leaks memory
        return self.view(np.ndarray)

    def __getitem__(self, key):
        """

        """
        if com.is_integer(key):
            return self._get_val_at(key)
        else:
            data_slice = self.values[key]
            return self._constructor(data_slice)

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        if j < 0:
            j = 0
        slobj = slice(i, j)
        return self.__getitem__(slobj)

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
            return _gin.get_value_at(self, sp_loc)

    def take(self, indices, axis=0):
        """
        Sparse-compatible version of ndarray.take

        Returns
        -------
        taken : ndarray
        """
        assert(axis == 0)
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
        raise Exception('SparseArray objects are immutable')

    def __setslice__(self, i, j, value):
        raise Exception('SparseArray objects are immutable')

    def to_dense(self):
        """
        Convert SparseSeries to (dense) Series
        """
        return self.values

    def astype(self, dtype=None):
        """

        """
        dtype = np.dtype(dtype)
        if dtype is not None and dtype not in (np.float_, float):
            raise Exception('Can only support floating point data for now')
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
        return SparseArray(values, sparse_index=self.sp_index,
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
        if com.notnull(self.fill_value):
            return self.to_dense().cumsum()
        # TODO: what if sp_values contains NaN??
        return SparseArray(self.sp_values.cumsum(),
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
