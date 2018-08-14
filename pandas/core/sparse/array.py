"""
SparseArray data structure
"""
from __future__ import division
# pylint: disable=E1101,E1103,W0231

import operator
import numpy as np
import warnings

import pandas as pd
import collections
from pandas.core.base import PandasObject

from pandas import compat
from pandas.errors import PerformanceWarning
from pandas.compat.numpy import function as nv

from pandas.core.arrays.base import ExtensionArray, ExtensionOpsMixin
from pandas.core.common import is_bool_indexer
from pandas.core.dtypes.generic import (
    ABCSparseSeries, ABCSeries, ABCIndexClass
)
from pandas.core.dtypes.common import (
    is_integer,
    is_object_dtype,
    is_array_like,
    is_extension_array_dtype,
    pandas_dtype,
    is_bool_dtype,
    is_list_like,
    is_string_dtype,
    is_scalar, is_dtype_equal)
from pandas.core.dtypes.cast import (
    maybe_convert_platform,
    astype_nansafe, find_common_type, infer_dtype_from_scalar,
    construct_1d_arraylike_from_scalar)
from pandas.core.dtypes.missing import isna, notna, na_value_for_dtype
from pandas.core.missing import interpolate_2d

import pandas._libs.sparse as splib
from pandas._libs.sparse import BlockIndex, IntIndex
from pandas._libs import index as libindex
import pandas.core.algorithms as algos
import pandas.io.formats.printing as printing

from .dtype import SparseDtype


_sparray_doc_kwargs = dict(klass='SparseArray')


def _get_fill(arr):
    # type: (SparseArray) -> ndarray
    # coerce fill_value to arr dtype if possible
    # int64 SparseArray can have NaN as fill_value if there is no missing
    try:
        return np.asarray(arr.fill_value, dtype=arr.dtype.subdtype)
    except ValueError:
        return np.asarray(arr.fill_value)


def _sparse_array_op(left, right, op, name):
    # type: (SparseArray, SparseArray, Callable, str) -> Any
    if name.startswith('__'):
        # For lookups in _libs.sparse we need non-dunder op name
        name = name[2:-2]

    # dtype used to find corresponding sparse method
    ltype = left.dtype.subdtype
    rtype = right.dtype.subdtype

    if not is_dtype_equal(ltype, rtype):
        dtype = SparseDtype(find_common_type([ltype, rtype]))
        left = left.astype(dtype)
        right = right.astype(dtype)
        dtype = dtype.subdtype
    else:
        dtype = ltype

    # dtype the result must have
    result_dtype = None

    if left.sp_index.ngaps == 0 or right.sp_index.ngaps == 0:
        with np.errstate(all='ignore'):
            result = op(left.get_values(), right.get_values())
            fill = op(_get_fill(left), _get_fill(right))

        if left.sp_index.ngaps == 0:
            index = left.sp_index
        else:
            index = right.sp_index
    elif left.sp_index.equals(right.sp_index):
        with np.errstate(all='ignore'):
            result = op(left.sp_values, right.sp_values)
            fill = op(_get_fill(left), _get_fill(right))
        index = left.sp_index
    else:
        if name[0] == 'r':
            left, right = right, left
            name = name[1:]

        if name in ('and', 'or') and dtype == 'bool':
            opname = 'sparse_{name}_uint8'.format(name=name)
            # to make template simple, cast here
            left_sp_values = left.sp_values.view(np.uint8)
            right_sp_values = right.sp_values.view(np.uint8)
            result_dtype = np.bool
        else:
            opname = 'sparse_{name}_{dtype}'.format(name=name, dtype=dtype)
            left_sp_values = left.sp_values
            right_sp_values = right.sp_values

        sparse_op = getattr(splib, opname)
        with np.errstate(all='ignore'):
            result, index, fill = sparse_op(left_sp_values, left.sp_index,
                                            left.fill_value, right_sp_values,
                                            right.sp_index, right.fill_value)

    if result_dtype is None:
        result_dtype = result.dtype

    return _wrap_result(name, result, index, fill, dtype=result_dtype)


def _wrap_result(name, data, sparse_index, fill_value, dtype=None):
    """ wrap op result to have correct dtype """
    if name.startswith('__'):
        # e.g. __eq__ --> eq
        name = name[2:-2]

    if name in ('eq', 'ne', 'lt', 'gt', 'le', 'ge'):
        dtype = np.bool

    if not is_scalar(fill_value):
        fill_value = fill_value.item()

    if is_bool_dtype(dtype):
        # fill_value may be np.bool_
        fill_value = bool(fill_value)
    return SparseArray(data, sparse_index=sparse_index, fill_value=fill_value,
                       dtype=dtype)


class SparseArray(PandasObject, ExtensionArray, ExtensionOpsMixin):
    """
    An ExtensionArray for storing sparse data.

    Parameters
    ----------
    data : array-like
    sparse_index : SparseIndex, optional
    index : Index
    fill_value : scalar, optional
        The fill_value to use for this array. By default, this is depends
        on the dtype of data.

        ========== ==========
        data.dtype na_value
        ========== ==========
        float      ``np.nan``
        int        ``0``
        bool       False
        datetime64 ``pd.NaT``
        ========== ==========

        When ``data`` is already a ``SparseArray``, ``data.fill_value``
        is used unless specified, regardless of `data.dtype``.

    kind : {'integer', 'block'}, default 'integer'
        The type of storage for sparse locations.

        * 'block': Stores a `block` and `block_length` for each
          contiguous *span* of sparse values. This is best when
          sparse data tends to be clumped together, with large
          regsions of ``fill-value`` values between sparse values.
        * 'integer': uses an integer to store the location of
          each sparse value.

    dtype : np.dtype, optional
    copy : bool, default False
        Whether to explicitly copy the incoming `data` array.
    """

    __array_priority__ = 15
    _pandas_ftype = 'sparse'
    _subtyp = 'sparse_array'  # register ABCSparseArray

    def __init__(self, data, sparse_index=None, index=None, fill_value=None,
                 kind='integer', dtype=None, copy=False):
        from pandas.core.internals import SingleBlockManager

        if isinstance(data, SingleBlockManager):
            data = data.internal_values()

        if isinstance(data, (type(self), ABCSparseSeries)):
            # disable normal inference on dtype, sparse_index, & fill_value
            if sparse_index is None:
                sparse_index = data.sp_index
            if fill_value is None:
                fill_value = data.fill_value
            if dtype is None:
                dtype = data.dtype
            # TODO: make kind=None, and use data.kind?
            data = data.sp_values

        if isinstance(dtype, SparseDtype):
            dtype = dtype.subdtype

        # TODO: index feels strange... can we deprecate it?
        if index is not None:
            if data is None:
                data = np.nan
            if not is_scalar(data):
                raise Exception("must only pass scalars with an index ")
            dtype = infer_dtype_from_scalar(data)[0]
            data = construct_1d_arraylike_from_scalar(
                data, len(index), dtype)

        if dtype is not None:
            dtype = pandas_dtype(dtype)

        # TODO: disentangle the fill_value dtype inference from
        # dtype inference
        if data is None:
            # XXX: What should the empty dtype be? Object or float?
            data = np.array([], dtype=dtype)

        if not is_array_like(data):
            try:
                # probably shared code in sanitize_series
                from pandas.core.series import _sanitize_array
                data = _sanitize_array(data, index=None)
            except ValueError:
                # NumPy may raise a ValueError on data like [1, []]
                # we retry with object dtype here.
                if dtype is None:
                    dtype = object
                    data = np.atleast_1d(np.asarray(data, dtype=dtype))
                else:
                    raise

        if copy:
            # TODO: avoid double copy when dtype forces cast.
            data = data.copy()

        if fill_value is None:
            fill_value_dtype = dtype or data.dtype
            if fill_value_dtype is None:
                fill_value = np.nan
            else:
                fill_value = na_value_for_dtype(fill_value_dtype)

        if isinstance(data, type(self)) and sparse_index is None:
            sparse_index = data._sparse_index
            sparse_values = np.asarray(data.sp_values, dtype=dtype)
        elif sparse_index is None:
            sparse_values, sparse_index, fill_value = make_sparse(
                data, kind=kind, fill_value=fill_value, dtype=dtype
            )
        else:
            sparse_values = np.asarray(data, dtype=dtype)
            if len(sparse_values) != sparse_index.npoints:
                raise AssertionError("Non array-like type {type} must "
                                     "have the same length as the index"
                                     .format(type=type(sparse_values)))
        self._sparse_index = sparse_index
        self._sparse_values = sparse_values
        self._dtype = SparseDtype(sparse_values.dtype)
        self.fill_value = fill_value

    @classmethod
    def _simple_new(cls, sparse_array, sparse_index, fill_value=None):
        # type: (SparseArray, SparseIndex, Any) -> 'SparseArray'
        new = cls([])
        new._sparse_index = sparse_index
        new._sparse_values = sparse_array
        new._dtype = sparse_array.dtype

        if fill_value is None:
            fill_value = sparse_array.fill_value
        new.fill_value = fill_value
        return new

    def __array__(self, dtype=None, copy=True):
        if self.sp_index.ngaps == 0:
            # Compat for na dtype and int values.
            return self.sp_values
        if dtype is None:
            dtype = np.result_type(self.sp_values.dtype, self.fill_value)
        out = np.full(self.shape, self.fill_value, dtype=dtype)
        out[self.sp_index.to_int_index().indices] = self.sp_values
        return out

    def __setitem__(self, key, value):
        # I suppose we could allow setting of non-fill_value elements.
        msg = "SparseArray does not support item assignment via setitem"
        raise TypeError(msg)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls(scalars, dtype=dtype)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values)

    # ------------------------------------------------------------------------
    # Data
    # ------------------------------------------------------------------------
    @property
    def sp_index(self):
        return self._sparse_index

    @property
    def sp_values(self):
        return self._sparse_values

    @property
    def dtype(self):
        return self._dtype

    @property
    def fill_value(self):
        return self._fill_value

    @fill_value.setter
    def fill_value(self, value):
        if not is_scalar(value):
            raise ValueError('fill_value must be a scalar')
        # if the specified value triggers type promotion, raise ValueError
        # new_dtype, fill_value = maybe_promote(self.dtype.subdtype, value)
        # if is_dtype_equal(self.dtype, new_dtype):
        self._fill_value = value
        # else:
        #     msg = 'unable to set fill_value {fill} to {dtype} dtype'
        #     raise ValueError(msg.format(fill=value, dtype=self.dtype))

    @property
    def kind(self):
        """
        The kind of sparse index for this array. One of {'integer', 'block'}.
        """
        # TODO: make this an abstract attribute of SparseIndex
        if isinstance(self.sp_index, IntIndex):
            return 'integer'
        else:
            return 'block'

    @property
    def _valid_sp_values(self):
        sp_vals = self.sp_values
        mask = notna(sp_vals)
        return sp_vals[mask]

    def __len__(self):
        return self.sp_index.length

    @property
    def _null_fill_value(self):
        return isna(self.fill_value)

    def _fill_value_matches(self, fill_value):
        if self._null_fill_value:
            return pd.isna(fill_value)
        else:
            return self.fill_value == fill_value

    @property
    def nbytes(self):
        return self.sp_values.nbytes + self.sp_index.nbytes

    @property
    def values(self):
        """
        Dense values
        """
        return self.to_dense()

    def isna(self):
        if isna(self.fill_value):
            # Then just the sparse values
            mask = np.ones(len(self), dtype=bool)
            # TODO: avoid to_int_index
            mask[self.sp_index.to_int_index().indices] = False
        else:
            # This is inevitable expensive?
            mask = pd.isna(np.asarray(self))
        return mask

    def fillna(self, value=None, method=None, limit=None):
        # TODO: discussion on what the return type should be.
        # Does it make sense to always return a SparseArray?
        # We *could* have the return type depend on whether self.fill_value
        # is NA.
        # But I think that's probably a bad idea...
        if method is not None:
            warnings.warn("Converting to dense in fillna with 'method'",
                          PerformanceWarning)
            filled = interpolate_2d(np.asarray(self), method=method,
                                    limit=limit)
            return type(self)(filled, fill_value=self.fill_value)

        if issubclass(self.dtype.type, np.floating):
            value = float(value)

        new_values = np.where(isna(self.sp_values), value, self.sp_values)
        fill_value = value if self._null_fill_value else self.fill_value

        return type(self)(new_values, self.sp_index, fill_value=fill_value)

    def unique(self):
        # The EA API currently expects unique to return the same EA.
        # That doesn't really make sense for sparse.
        # Can we have it expect Union[EA, ndarray]?
        return type(self)(pd.unique(self.sp_values))

    def factorize(self, na_sentinel=-1):
        # hhhhhhhhhhhhhhhhhhhhhhhhhhhhmmmm
        # Ok. here's the plan...
        # We known that we'll share the same sparsity
        # so factorize our known values
        # and then rebuild using the same sparse index?
        if na_sentinel > 0:
            raise ValueError("na_sentinel must be less than 0. "
                             "Got {}".format(na_sentinel))

        known, uniques = pd.factorize(self.sp_values)
        new = SparseArray(known, sparse_index=self.sp_index,
                          fill_value=na_sentinel)
        # ah, but we have to go to sparse :/
        # so we're backwards in our sparsity her.
        return np.asarray(new), type(self)(uniques)

    def value_counts(self, dropna=True):
        """
        Returns a Series containing counts of unique values.

        Parameters
        ----------
        dropna : boolean, default True
            Don't include counts of NaN, even if NaN is in sp_values.

        Returns
        -------
        counts : Series
        """
        keys, counts = algos._value_counts_arraylike(self.sp_values,
                                                     dropna=dropna)
        fcounts = self.sp_index.ngaps
        if fcounts > 0:
            if self._null_fill_value and dropna:
                pass
            else:
                if self._null_fill_value:
                    mask = pd.isna(keys)
                else:
                    mask = keys == self.fill_value

                if mask.any():
                    counts[mask] += fcounts
                else:
                    keys = np.insert(keys, 0, self.fill_value)
                    counts = np.insert(counts, 0, fcounts)

        if not isinstance(keys, pd.Index):
            keys = pd.Index(keys)
        result = pd.Series(counts, index=keys)
        return result

    # --------
    # Indexing
    # --------

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) > 1:
                raise IndexError("too many indices for array.")
            key = key[0]

        if is_integer(key):
            return self._get_val_at(key)
        elif isinstance(key, tuple):
            data_slice = self.values[key]
        elif isinstance(key, slice):
            # special case to preserve dtypes
            if key == slice(None):
                return self.copy()
            # TODO: this logic is surely elsewhere
            # TODO: this could be more efficient
            indices = np.arange(len(self), dtype=np.int32)[key]
            return self.take(indices)
        else:
            if isinstance(key, SparseArray):
                if is_bool_dtype(key):
                    key = key.to_dense()
                else:
                    key = np.asarray(key)

            if hasattr(key, '__len__') and len(self) != len(key):
                return self.take(key)
            elif is_bool_indexer(key) and len(self) == len(key):
                return self.take(np.arange(len(key), dtype=np.int32)[key])
            else:
                # TODO: this densifies!
                data_slice = self.values[key]

        return type(self)(data_slice, kind=self.kind)

    def _get_val_at(self, loc):
        n = len(self)
        if loc < 0:
            loc += n

        if loc >= n or loc < 0:
            raise IndexError('Out of bounds access')

        sp_loc = self.sp_index.lookup(loc)
        if sp_loc == -1:
            return self.fill_value
        else:
            return libindex.get_value_at(self.sp_values, sp_loc)

    def _boolean_mask(self, key):
        # strategy:
        pass

    def take(self, indices, allow_fill=False, fill_value=None):
        indices = np.asarray(indices, dtype=np.int32)

        if indices.size == 0:
            result = []
        elif allow_fill:
            result = self._take_with_fill(indices, fill_value=fill_value)
        else:
            result = self._take_without_fill(indices)

        return type(self)(result, fill_value=self.fill_value, kind=self.kind)

    def _take_with_fill(self, indices, fill_value=None):
        if fill_value is None:
            fill_value = self.dtype.na_value

        if indices.min() < -1:
            raise ValueError("Invalid value in 'indices'. Must be between -1 "
                             "and the length of the array.")

        if indices.max() >= len(self):
            raise IndexError("out of bounds value in 'indices'.")

        if len(self) == 0:
            # Empty... Allow taking only if all empty
            if (indices == -1).all():
                taken = np.empty_like(indices, dtype=self.sp_values.dtype)
                taken.fill(fill_value)
                return taken
            else:
                raise IndexError('cannot do a non-empty take from an empty '
                                 'axes.')

        sp_indexer = self.sp_index.lookup_array(indices)

        if self.sp_index.npoints == 0:
            # Avoid taking from the empty self.sp_values
            taken = np.full(sp_indexer.shape, fill_value=fill_value)
        else:
            taken = self.sp_values.take(sp_indexer)

            # sp_indexer may be -1 for two reasons
            # 1.) we took for an index of -1 (new)
            # 2.) we took a value that was self.fill_value (old)
            new_fill_indices = indices == -1
            old_fill_indices = (sp_indexer == -1) & ~new_fill_indices

            # Fill in two steps.
            # Old fill values
            # New fill values
            # potentially coercing to a new dtype at each stage.

            m0 = sp_indexer[old_fill_indices] < 0
            m1 = sp_indexer[new_fill_indices] < 0

            result_type = taken.dtype

            if m0.any():
                result_type = np.result_type(result_type, self.fill_value)
                taken = taken.astype(result_type)
                taken[old_fill_indices] = self.fill_value

            if m1.any():
                result_type = np.result_type(result_type, fill_value)
                taken = taken.astype(result_type)
                taken[new_fill_indices] = fill_value

        return taken

    def _take_without_fill(self, indices):
        to_shift = indices < 0
        indices = indices.copy()

        n = len(self)

        if (indices.max() >= n) or (indices.min() < -n):
            if n == 0:
                raise IndexError("cannot do a non-empty take from an "
                                 "empty axes.")
            else:
                raise IndexError("out of bounds value in 'indices'.")

        if to_shift.any():
            indices[to_shift] += n

        if self.sp_index.npoints == 0:
            # edge case in take...
            # I think just return
            out = np.full(indices.shape, self.fill_value)
            arr, sp_index, fill_value = make_sparse(out,
                                                    fill_value=self.fill_value)
            return type(self)(arr, sparse_index=sp_index,
                              fill_value=fill_value)

        sp_indexer = self.sp_index.lookup_array(indices)
        taken = self.sp_values.take(sp_indexer)
        fillable = (sp_indexer < 0)

        if fillable.any():
            # TODO: may need to coerce array to fill value
            result_type = np.result_type(taken, self.fill_value)
            taken = taken.astype(result_type)
            taken[fillable] = self.fill_value

        return taken

    def copy(self, deep=False):
        if deep:
            values = self.sp_values.copy()
        else:
            values = self.sp_values

        return type(self)(values, sparse_index=self.sp_index, copy=False,
                          fill_value=self.fill_value)

    @classmethod
    def _concat_same_type(cls, to_concat):
        fill_values = list(x.fill_value for x in to_concat)

        fill_value = fill_values[0]

        if len(set(fill_values)) > 1:
            warnings.warn("Concatenating sparse arrays with multiple fill "
                          "values: '{}'. Picking the first and "
                          "converting the rest.".format(fill_values),
                          PerformanceWarning,
                          stacklevel=6)
            keep = to_concat[0]
            to_concat2 = [keep]

            for arr in to_concat[1:]:
                to_concat2.append(cls(np.asarray(arr), fill_value=fill_value))

            to_concat = to_concat2

        values = []
        length = 0

        if to_concat:
            sp_kind = to_concat[0].kind
        else:
            sp_kind = 'integer'

        if sp_kind == 'integer':
            indices = []

            for arr in to_concat:
                idx = arr.sp_index.to_int_index().indices.copy()
                idx += length  # TODO: wraparound
                length += arr.sp_index.length

                values.append(arr.sp_values)
                indices.append(idx)

            data = np.concatenate(values)
            indices = np.concatenate(indices)
            sp_index = IntIndex(length, indices)

        else:
            # when concatentating block indices, we don't claim that you'll
            # get an identical index as concating the values and then
            # creating a new index. We don't want to spend the time trying
            # to merge blocks across arrays in `to_concat`, so the resulting
            # BlockIndex may have more blocs.
            blengths = []
            blocs = []

            for arr in to_concat:
                idx = arr.sp_index.to_block_index()

                values.append(arr.sp_values)
                blocs.append(idx.blocs.copy() + length)
                blengths.append(idx.blengths)
                length += arr.sp_index.length

            data = np.concatenate(values)
            blocs = np.concatenate(blocs)
            blengths = np.concatenate(blengths)

            sp_index = BlockIndex(length, blocs, blengths)

        return cls(data, sparse_index=sp_index, fill_value=fill_value)

    def astype(self, dtype=None, copy=True):
        # TODO: Document API Change here: .astype(type) will densify
        # for non-sparse types
        dtype = pandas_dtype(dtype)

        if isinstance(dtype, SparseDtype):
            # Sparse -> Sparse
            sp_values = astype_nansafe(self.sp_values, dtype.subdtype,
                                       copy=copy)
            try:
                if is_bool_dtype(dtype):
                    # to avoid np.bool_ dtype
                    fill_value = bool(self.fill_value)
                else:
                    fill_value = dtype.type(self.fill_value)
            except ValueError:
                msg = ('unable to coerce current fill_value {fill} to '
                       '{dtype} dtype')
                raise ValueError(msg.format(fill=self.fill_value,
                                            dtype=dtype))
            return type(self)(sp_values, self.sp_index, fill_value=fill_value)
        elif is_extension_array_dtype(dtype):
            return dtype.construct_array_type()(self, copy=copy)
        else:
            return astype_nansafe(np.asarray(self), dtype=dtype)

    def map(self, mapper):
        # this is used in apply.
        # We get hit since we're an "is_extension_type" but regular extension
        # types are not hit...
        if isinstance(mapper, collections.Mapping):
            fill_value = mapper.get(self.fill_value, self.fill_value)
            sp_values = [mapper.get(x, None) for x in self.sp_values]
        else:
            fill_value = mapper(self.fill_value)
            sp_values = [mapper(x) for x in self.sp_values]

        # TODO: series?
        return type(self)(sp_values, sparse_index=self.sp_index,
                          fill_value=fill_value)

    def get_values(self, fill=None):
        """ return a dense representation """
        # TODO: deprecate for to_dense?
        return self.to_dense(fill=fill)

    def to_dense(self, fill=None):
        """
        Convert SparseArray to a NumPy array.

        Parameters
        ----------
        fill: float, default None
            .. deprecated:: 0.20.0
               This argument is not respected by this function.

        Returns
        -------
        arr : NumPy array
        """
        if fill is not None:
            warnings.warn(("The 'fill' parameter has been deprecated and "
                           "will be removed in a future version."),
                          FutureWarning, stacklevel=2)
        return np.asarray(self, dtype=self.sp_values.dtype)

    # ------------------------------------------------------------------------
    # IO
    # ------------------------------------------------------------------------
    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if isinstance(state, tuple):
            # Compat for pandas < 0.24.0
            nd_state, own_state = state
            sparse_values = np.array([])
            sparse_values.__setstate__(nd_state)

            self._sparse_values = sparse_values
            self.fill_value, self._sparse_index = own_state[:2]
            self._dtype = SparseDtype(sparse_values.dtype)
        else:
            self.__dict__.update(state)

    def nonzero(self):
        # TODO: Add to EA API? This is used by DataFrame.dropna
        if self.fill_value == 0:
            return self.sp_index.to_int_index().indices,
        else:
            return self.sp_index.to_int_index().indices[self.sp_values != 0],

    # ------------------------------------------------------------------------
    # Reductions
    # ------------------------------------------------------------------------

    def all(self, axis=None, *args, **kwargs):
        """
        Tests whether all elements evaluate True

        Returns
        -------
        all : bool

        See Also
        --------
        numpy.all
        """
        nv.validate_all(args, kwargs)

        values = self.sp_values

        if len(values) != len(self) and not np.all(self.fill_value):
            return False

        return values.all()

    def any(self, axis=0, *args, **kwargs):
        """
        Tests whether at least one of elements evaluate True

        Returns
        -------
        any : bool

        See Also
        --------
        numpy.any
        """
        nv.validate_any(args, kwargs)

        values = self.sp_values

        if len(values) != len(self) and np.any(self.fill_value):
            return True

        return values.any()

    def sum(self, axis=0, *args, **kwargs):
        """
        Sum of non-NA/null values

        Returns
        -------
        sum : float
        """
        nv.validate_sum(args, kwargs)
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        if self._null_fill_value:
            return sp_sum
        else:
            nsparse = self.sp_index.ngaps
            return sp_sum + self.fill_value * nsparse

    def cumsum(self, axis=0, *args, **kwargs):
        """
        Cumulative sum of non-NA/null values.

        When performing the cumulative summation, any non-NA/null values will
        be skipped. The resulting SparseArray will preserve the locations of
        NaN values, but the fill value will be `np.nan` regardless.

        Parameters
        ----------
        axis : int or None
            Axis over which to perform the cumulative summation. If None,
            perform cumulative summation over flattened array.

        Returns
        -------
        cumsum : SparseArray
        """
        nv.validate_cumsum(args, kwargs)

        if axis is not None and axis >= self.ndim:  # Mimic ndarray behaviour.
            raise ValueError("axis(={axis}) out of bounds".format(axis=axis))

        if not self._null_fill_value:
            return SparseArray(self.to_dense()).cumsum()

        return SparseArray(self.sp_values.cumsum(), sparse_index=self.sp_index,
                           fill_value=self.fill_value)

    def mean(self, axis=0, *args, **kwargs):
        """
        Mean of non-NA/null values

        Returns
        -------
        mean : float
        """
        nv.validate_mean(args, kwargs)
        valid_vals = self._valid_sp_values
        sp_sum = valid_vals.sum()
        ct = len(valid_vals)

        if self._null_fill_value:
            return sp_sum / ct
        else:
            nsparse = self.sp_index.ngaps
            return (sp_sum + self.fill_value * nsparse) / (ct + nsparse)

    def transpose(self, *axes):
        """Returns the SparseArray."""
        return self

    @property
    def T(self):
        """Returns the SparseArray."""
        return self

    # ------------------------------------------------------------------------
    # Ufuncs
    # ------------------------------------------------------------------------
    def __abs__(self):
        return np.abs(self)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        new_inputs = []
        new_fill_values = []

        special = {'add', 'sub', 'mul', 'pow', 'mod', 'floordiv', 'truediv',
                   'divmod', 'eq', 'ne', 'lt', 'gt', 'le', 'ge', 'remainder'}
        aliases = {
            'subtract': 'sub',
            'multiply': 'mul',
            'floor_divide': 'floordiv',
            'true_divide': 'truediv',
            'power': 'pow',
            'remainder': 'mod',
        }
        op_name = ufunc.__name__
        op_name = aliases.get(op_name, op_name)

        if op_name in special:
            if isinstance(inputs[0], type(self)):
                # this is surely incorrect...
                return getattr(self, '__{}__'.format(op_name))(inputs[1])
            else:
                return getattr(self, '__r{}__'.format(op_name))(inputs[0])

        for input in inputs:
            if isinstance(input, type(self)):
                new_inputs.append(self.sp_values)
                new_fill_values.append(self.fill_value)
            else:
                new_inputs.append(input)
                new_fill_values.append(input)

        new_values = ufunc(*new_inputs, **kwargs)
        new_fill = ufunc(*new_fill_values, **kwargs)
        # TODO:
        # call ufunc on fill_value?
        # What about a new sparse index?
        return type(self)(new_values, sparse_index=self.sp_index,
                          fill_value=new_fill)

    # ------------------------------------------------------------------------
    # Ops
    # ------------------------------------------------------------------------

    @classmethod
    def _create_arithmetic_method(cls, op):
        def sparse_arithmetic_method(self, other):
            op_name = op.__name__

            if isinstance(other, (ABCSeries, ABCIndexClass)):
                other = getattr(other, 'values', other)

            if isinstance(other, SparseArray):
                return _sparse_array_op(self, other, op, op_name)

            elif is_scalar(other):
                with np.errstate(all='ignore'):
                    fill = op(_get_fill(self), np.asarray(other))
                    result = op(self.sp_values, other)
                return _wrap_result(op_name, result, self.sp_index, fill)

            else:
                with np.errstate(all='ignore'):
                    # TODO: delete sparse stuff in core/ops.py
                    # TODO: look into _wrap_result
                    if len(self) != len(other):
                        raise AssertionError(
                            ("length mismatch: {self} vs. {other}".format(
                                self=len(self), other=len(other))))
                    if not isinstance(other, SparseArray):
                        dtype = getattr(other, 'dtype', None)
                        other = SparseArray(other, fill_value=self.fill_value,
                                            dtype=dtype)
                    return _sparse_array_op(self, other, op, op_name)
                    # fill_value = op(self.fill_value, other)
                    # result = op(self.sp_values, other)

                # TODO: is self.sp_index right? An op could change what's
                # sparse...
                # return type(self)(result, sparse_index=self.sp_index,
                #                   fill_value=fill_value)

        name = '__{name}__'.format(name=op.__name__)
        return compat.set_function_name(sparse_arithmetic_method, name, cls)

    @classmethod
    def _create_comparison_method(cls, op):
        def cmp_method(self, other):
            op_name = op.__name__

            if op_name in {'and_', 'or_'}:
                op_name = op_name[:-1]

            if isinstance(other, (ABCSeries, ABCIndexClass)):
                other = getattr(other, 'values', other)

            if isinstance(other, np.ndarray):
                # TODO: make this more flexible than just ndarray...
                if len(self) != len(other):
                    raise AssertionError("length mismatch: {self} vs. {other}"
                                         .format(self=len(self),
                                                 other=len(other)))
                other = SparseArray(other, fill_value=self.fill_value)

            if isinstance(other, SparseArray):
                return _sparse_array_op(self, other, op, op_name)
            else:
                with np.errstate(all='ignore'):
                    fill_value = op(self.fill_value, other)
                    result = op(self.sp_values, other)

                return type(self)(result,
                                  sparse_index=self.sp_index,
                                  fill_value=fill_value,
                                  dtype=np.bool_)

        name = '__{name}__'.format(name=op.__name__)
        return compat.set_function_name(cmp_method, name, cls)

    # ----------
    # Formatting
    # -----------
    def __unicode__(self):
        return '{self}\nFill: {fill}\n{index}'.format(
            self=printing.pprint_thing(self),
            fill=printing.pprint_thing(self.fill_value),
            index=printing.pprint_thing(self.sp_index))


SparseArray._add_arithmetic_ops()
SparseArray._add_comparison_ops()
SparseArray.__and__ = SparseArray._create_comparison_method(operator.and_)
SparseArray.__or__ = SparseArray._create_comparison_method(operator.or_)


def _maybe_to_dense(obj):
    """ try to convert to dense """
    if hasattr(obj, 'to_dense'):
        return obj.to_dense()
    return obj


def _maybe_to_sparse(array):
    """ array must be SparseSeries or SparseArray """
    if isinstance(array, ABCSparseSeries):
        array = array.values.copy()
    return array


def _sanitize_values(arr):
    """
    return an ndarray for our input,
    in a platform independent manner
    """

    if hasattr(arr, 'values'):
        arr = arr.values
    else:

        # scalar
        if is_scalar(arr):
            arr = [arr]

        # ndarray
        if isinstance(arr, np.ndarray):
            pass

        elif is_list_like(arr) and len(arr) > 0:
            arr = maybe_convert_platform(arr)

        else:
            arr = np.asarray(arr)

    return arr


def make_sparse(arr, kind='block', fill_value=None, dtype=None, copy=False):
    """
    Convert ndarray to sparse format

    Parameters
    ----------
    arr : ndarray
    kind : {'block', 'integer'}
    fill_value : NaN or another value
    dtype : np.dtype, optional
    copy : bool, default False

    Returns
    -------
    (sparse_values, index) : (ndarray, SparseIndex)
    """

    arr = _sanitize_values(arr)

    if arr.ndim > 1:
        raise TypeError("expected dimension <= 1 data")

    if fill_value is None:
        fill_value = na_value_for_dtype(arr.dtype)

    if isna(fill_value):
        mask = notna(arr)
    else:
        # For str arrays in NumPy 1.12.0, operator!= below isn't
        # element-wise but just returns False if fill_value is not str,
        # so cast to object comparison to be safe
        if is_string_dtype(arr):
            arr = arr.astype(object)

        if is_object_dtype(arr.dtype):
            # element-wise equality check method in numpy doesn't treat
            # each element type, eg. 0, 0.0, and False are treated as
            # same. So we have to check the both of its type and value.
            mask = splib.make_mask_object_ndarray(arr, fill_value)
        else:
            mask = arr != fill_value

    length = len(arr)
    if length != mask.size:
        # the arr is a SparseArray
        indices = mask.sp_index.indices
    else:
        indices = mask.nonzero()[0].astype(np.int32)

    index = _make_index(length, indices, kind)
    sparsified_values = arr[mask]

    sparsified_values = np.asarray(sparsified_values, dtype=dtype)
    # TODO: copy
    return sparsified_values, index, fill_value


def _make_index(length, indices, kind):

    if kind == 'block' or isinstance(kind, BlockIndex):
        locs, lens = splib.get_blocks(indices)
        index = BlockIndex(length, locs, lens)
    elif kind == 'integer' or isinstance(kind, IntIndex):
        index = IntIndex(length, indices)
    else:  # pragma: no cover
        raise ValueError('must be block or integer type')
    return index
