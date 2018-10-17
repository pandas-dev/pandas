"""
SparseArray data structure
"""
from __future__ import division
# pylint: disable=E1101,E1103,W0231

import re
import operator
import numbers
import numpy as np
import warnings

import pandas as pd
from pandas.core.base import PandasObject

from pandas import compat
from pandas.errors import PerformanceWarning
from pandas.compat.numpy import function as nv

from pandas.core.arrays import ExtensionArray, ExtensionOpsMixin
import pandas.core.common as com
from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import register_extension_dtype
from pandas.core.dtypes.generic import (
    ABCSparseSeries, ABCSeries, ABCIndexClass
)
from pandas.core.dtypes.common import (
    is_datetime64_any_dtype,
    is_integer,
    is_object_dtype,
    is_array_like,
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
from pandas._libs import lib
import pandas.core.algorithms as algos
import pandas.io.formats.printing as printing


# ----------------------------------------------------------------------------
# Dtype

@register_extension_dtype
class SparseDtype(ExtensionDtype):
    """
    Dtype for data stored in :class:`SparseArray`.

    This dtype implements the pandas ExtensionDtype interface.

    .. versionadded:: 0.24.0

    Parameters
    ----------
    dtype : str, ExtensionDtype, numpy.dtype, type, default numpy.float64
        The dtype of the underlying array storing the non-fill value values.
    fill_value : scalar, optional.
        The scalar value not stored in the SparseArray. By default, this
        depends on `dtype`.

        ========== ==========
        dtype      na_value
        ========== ==========
        float      ``np.nan``
        int        ``0``
        bool       ``False``
        datetime64 ``pd.NaT``
        timedelta64 ``pd.NaT``
        ========== ==========

        The default value may be overridden by specifying a `fill_value`.
    """
    # We include `_is_na_fill_value` in the metadata to avoid hash collisions
    # between SparseDtype(float, 0.0) and SparseDtype(float, nan).
    # Without is_na_fill_value in the comparison, those would be equal since
    # hash(nan) is (sometimes?) 0.
    _metadata = ('_dtype', '_fill_value', '_is_na_fill_value')

    def __init__(self, dtype=np.float64, fill_value=None):
        # type: (Union[str, np.dtype, 'ExtensionDtype', type], Any) -> None
        from pandas.core.dtypes.missing import na_value_for_dtype
        from pandas.core.dtypes.common import (
            pandas_dtype, is_string_dtype, is_scalar
        )

        if isinstance(dtype, type(self)):
            if fill_value is None:
                fill_value = dtype.fill_value
            dtype = dtype.subtype

        dtype = pandas_dtype(dtype)
        if is_string_dtype(dtype):
            dtype = np.dtype('object')

        if fill_value is None:
            fill_value = na_value_for_dtype(dtype)

        if not is_scalar(fill_value):
            raise ValueError("fill_value must be a scalar. Got {} "
                             "instead".format(fill_value))
        self._dtype = dtype
        self._fill_value = fill_value

    def __hash__(self):
        # Python3 doesn't inherit __hash__ when a base class overrides
        # __eq__, so we explicitly do it here.
        return super(SparseDtype, self).__hash__()

    def __eq__(self, other):
        # We have to override __eq__ to handle NA values in _metadata.
        # The base class does simple == checks, which fail for NA.
        if isinstance(other, compat.string_types):
            try:
                other = self.construct_from_string(other)
            except TypeError:
                return False

        if isinstance(other, type(self)):
            subtype = self.subtype == other.subtype
            if self._is_na_fill_value:
                # this case is complicated by two things:
                # SparseDtype(float, float(nan)) == SparseDtype(float, np.nan)
                # SparseDtype(float, np.nan)     != SparseDtype(float, pd.NaT)
                # i.e. we want to treat any floating-point NaN as equal, but
                # not a floating-point NaN and a datetime NaT.
                fill_value = (
                    other._is_na_fill_value and
                    isinstance(self.fill_value, type(other.fill_value)) or
                    isinstance(other.fill_value, type(self.fill_value))
                )
            else:
                fill_value = self.fill_value == other.fill_value

            return subtype and fill_value
        return False

    @property
    def fill_value(self):
        """
        The fill value of the array.

        Converting the SparseArray to a dense ndarray will fill the
        array with this value.

        .. warning::

           It's possible to end up with a SparseArray that has ``fill_value``
           values in ``sp_values``. This can occur, for example, when setting
           ``SparseArray.fill_value`` directly.
        """
        return self._fill_value

    @property
    def _is_na_fill_value(self):
        from pandas.core.dtypes.missing import isna
        return isna(self.fill_value)

    @property
    def _is_numeric(self):
        from pandas.core.dtypes.common import is_object_dtype
        return not is_object_dtype(self.subtype)

    @property
    def _is_boolean(self):
        from pandas.core.dtypes.common import is_bool_dtype
        return is_bool_dtype(self.subtype)

    @property
    def kind(self):
        return self.subtype.kind

    @property
    def type(self):
        return self.subtype.type

    @property
    def subtype(self):
        return self._dtype

    @property
    def name(self):
        return 'Sparse[{}, {}]'.format(self.subtype.name, self.fill_value)

    def __repr__(self):
        return self.name

    @classmethod
    def construct_array_type(cls):
        return SparseArray

    @classmethod
    def construct_from_string(cls, string):
        """
        Construct a SparseDtype from a string form.

        Parameters
        ----------
        string : str
            Can take the following forms.

            string           dtype
            ================ ============================
            'int'            SparseDtype[np.int64, 0]
            'Sparse'         SparseDtype[np.float64, nan]
            'Sparse[int]'    SparseDtype[np.int64, 0]
            'Sparse[int, 0]' SparseDtype[np.int64, 0]
            ================ ============================

            It is not possible to specify non-default fill values
            with a string. An argument like ``'Sparse[int, 1]'``
            will raise a ``TypeError`` because the default fill value
            for integers is 0.

        Returns
        -------
        SparseDtype
        """
        msg = "Could not construct SparseDtype from '{}'".format(string)
        if string.startswith("Sparse"):
            try:
                sub_type, has_fill_value = cls._parse_subtype(string)
                result = SparseDtype(sub_type)
            except Exception:
                raise TypeError(msg)
            else:
                msg = ("Could not construct SparseDtype from '{}'.\n\nIt "
                       "looks like the fill_value in the string is not "
                       "the default for the dtype. Non-default fill_values "
                       "are not supported. Use the 'SparseDtype()' "
                       "constructor instead.")
                if has_fill_value and str(result) != string:
                    raise TypeError(msg.format(string))
                return result
        else:
            raise TypeError(msg)

    @staticmethod
    def _parse_subtype(dtype):
        """
        Parse a string to get the subtype

        Parameters
        ----------
        dtype : str
            A string like

            * Sparse[subtype]
            * Sparse[subtype, fill_value]

        Returns
        -------
        subtype : str

        Raises
        ------
        ValueError
            When the subtype cannot be extracted.
        """
        xpr = re.compile(
            r"Sparse\[(?P<subtype>[^,]*)(, )?(?P<fill_value>.*?)?\]$"
        )
        m = xpr.match(dtype)
        has_fill_value = False
        if m:
            subtype = m.groupdict()['subtype']
            has_fill_value = m.groupdict()['fill_value'] or has_fill_value
        elif dtype == "Sparse":
            subtype = 'float64'
        else:
            raise ValueError("Cannot parse {}".format(dtype))
        return subtype, has_fill_value

    @classmethod
    def is_dtype(cls, dtype):
        dtype = getattr(dtype, 'dtype', dtype)
        if (isinstance(dtype, compat.string_types) and
                dtype.startswith("Sparse")):
            sub_type, _ = cls._parse_subtype(dtype)
            dtype = np.dtype(sub_type)
        elif isinstance(dtype, cls):
            return True
        return isinstance(dtype, np.dtype) or dtype == 'Sparse'

# ----------------------------------------------------------------------------
# Array


_sparray_doc_kwargs = dict(klass='SparseArray')


def _get_fill(arr):
    # type: (SparseArray) -> ndarray
    """
    Create a 0-dim ndarray containing the fill value

    Parameters
    ----------
    arr : SparseArray

    Returns
    -------
    fill_value : ndarray
        0-dim ndarray with just the fill value.

    Notes
    -----
    coerce fill_value to arr dtype if possible
    int64 SparseArray can have NaN as fill_value if there is no missing
    """
    try:
        return np.asarray(arr.fill_value, dtype=arr.dtype.subtype)
    except ValueError:
        return np.asarray(arr.fill_value)


def _sparse_array_op(left, right, op, name):
    """
    Perform a binary operation between two arrays.

    Parameters
    ----------
    left : Union[SparseArray, ndarray]
    right : Union[SparseArray, ndarray]
    op : Callable
        The binary operation to perform
    name str
        Name of the callable.

    Returns
    -------
    SparseArray
    """
    # type: (SparseArray, SparseArray, Callable, str) -> Any
    if name.startswith('__'):
        # For lookups in _libs.sparse we need non-dunder op name
        name = name[2:-2]

    # dtype used to find corresponding sparse method
    ltype = left.dtype.subtype
    rtype = right.dtype.subtype

    if not is_dtype_equal(ltype, rtype):
        subtype = find_common_type([ltype, rtype])
        ltype = SparseDtype(subtype, left.fill_value)
        rtype = SparseDtype(subtype, right.fill_value)

        # TODO(GH-23092): pass copy=False. Need to fix astype_nansafe
        left = left.astype(ltype)
        right = right.astype(rtype)
        dtype = ltype.subtype
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
            result, index, fill = sparse_op(
                left_sp_values, left.sp_index, left.fill_value,
                right_sp_values, right.sp_index, right.fill_value)

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

    fill_value = lib.item_from_zerodim(fill_value)

    if is_bool_dtype(dtype):
        # fill_value may be np.bool_
        fill_value = bool(fill_value)
    return SparseArray(data,
                       sparse_index=sparse_index,
                       fill_value=fill_value,
                       dtype=dtype)


class SparseArray(PandasObject, ExtensionArray, ExtensionOpsMixin):
    """
    An ExtensionArray for storing sparse data.

    .. versionchanged:: 0.24.0

       Implements the ExtensionArray interface.

    Parameters
    ----------
    data : array-like
        A dense array of values to store in the SparseArray. This may contain
        `fill_value`.
    sparse_index : SparseIndex, optional
    index : Index
    fill_value : scalar, optional
        Elements in `data` that are `fill_value` are not stored in the
        SparseArray. For memory savings, this should be the most common value
        in `data`. By default, `fill_value` depends on the dtype of `data`:

        =========== ==========
        data.dtype  na_value
        =========== ==========
        float       ``np.nan``
        int         ``0``
        bool        False
        datetime64  ``pd.NaT``
        timedelta64 ``pd.NaT``
        =========== ==========

        The fill value is potentiall specified in three ways. In order of
        precedence, these are

        1. The `fill_value` argument
        2. ``dtype.fill_value`` if `fill_value` is None and `dtype` is
           a ``SparseDtype``
        3. ``data.dtype.fill_value`` if `fill_value` is None and `dtype`
           is not a ``SparseDtype`` and `data` is a ``SparseArray``.


    kind : {'integer', 'block'}, default 'integer'
        The type of storage for sparse locations.

        * 'block': Stores a `block` and `block_length` for each
          contiguous *span* of sparse values. This is best when
          sparse data tends to be clumped together, with large
          regsions of ``fill-value`` values between sparse values.
        * 'integer': uses an integer to store the location of
          each sparse value.

    dtype : np.dtype or SparseDtype, optional
        The dtype to use for the SparseArray. For numpy dtypes, this
        determines the dtype of ``self.sp_values``. For SparseDtype,
        this determines ``self.sp_values`` and ``self.fill_value``.
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

        if fill_value is None and isinstance(dtype, SparseDtype):
            fill_value = dtype.fill_value

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

        # Handle use-provided dtype
        if isinstance(dtype, compat.string_types):
            # Two options: dtype='int', regular numpy dtype
            # or dtype='Sparse[int]', a sparse dtype
            try:
                dtype = SparseDtype.construct_from_string(dtype)
            except TypeError:
                dtype = pandas_dtype(dtype)

        if isinstance(dtype, SparseDtype):
            if fill_value is None:
                fill_value = dtype.fill_value
            dtype = dtype.subtype

        if index is not None and not is_scalar(data):
            raise Exception("must only pass scalars with an index ")

        if is_scalar(data):
            if index is not None:
                if data is None:
                    data = np.nan

            if index is not None:
                npoints = len(index)
            elif sparse_index is None:
                npoints = 1
            else:
                npoints = sparse_index.length

            dtype = infer_dtype_from_scalar(data)[0]
            data = construct_1d_arraylike_from_scalar(
                data, npoints, dtype
            )

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
            fill_value_dtype = data.dtype if dtype is None else dtype
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
        self._dtype = SparseDtype(sparse_values.dtype, fill_value)

    @classmethod
    def _simple_new(cls, sparse_array, sparse_index, dtype):
        # type: (np.ndarray, SparseIndex, SparseDtype) -> 'SparseArray'
        new = cls([])
        new._sparse_index = sparse_index
        new._sparse_values = sparse_array
        new._dtype = dtype
        return new

    def __array__(self, dtype=None, copy=True):
        fill_value = self.fill_value

        if self.sp_index.ngaps == 0:
            # Compat for na dtype and int values.
            return self.sp_values
        if dtype is None:
            # Can NumPy represent this type?
            # If not, `np.result_type` will raise. We catch that
            # and return object.
            if is_datetime64_any_dtype(self.sp_values.dtype):
                # However, we *do* special-case the common case of
                # a datetime64 with pandas NaT.
                if fill_value is pd.NaT:
                    # Can't put pd.NaT in a datetime64[ns]
                    fill_value = np.datetime64('NaT')
            try:
                dtype = np.result_type(self.sp_values.dtype, fill_value)
            except TypeError:
                dtype = object

        out = np.full(self.shape, fill_value, dtype=dtype)
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
        return cls(values, dtype=original.dtype)

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
        """
        Elements in `data` that are `fill_value` are not stored.

        For memory savings, this should be the most common value in the array.
        """
        return self.dtype.fill_value

    @fill_value.setter
    def fill_value(self, value):
        self._dtype = SparseDtype(self.dtype.subtype, value)

    @property
    def kind(self):
        """
        The kind of sparse index for this array. One of {'integer', 'block'}.
        """
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
        return self._dtype._is_na_fill_value

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
        from pandas import isna
        # If null fill value, we want SparseDtype[bool, true]
        # to preserve the same memory usage.
        dtype = SparseDtype(bool, self._null_fill_value)
        return type(self)._simple_new(isna(self.sp_values),
                                      self.sp_index, dtype)

    def fillna(self, value=None, method=None, limit=None):
        """
        Fill missing values with `value`.

        Parameters
        ----------
        value : scalar, optional
        method : str, optional

            .. warning::

               Using 'method' will result in high memory use,
               as all `fill_value` methods will be converted to
               an in-memory ndarray

        limit : int, optional

        Returns
        -------
        SparseArray

        Notes
        -----
        When `value` is specified, the result's ``fill_value`` depends on
        ``self.fill_value``. The goal is to maintain low-memory use.

        If ``self.fill_value`` is NA, the result dtype will be
        ``SparseDtype(self.dtype, fill_value=value)``. This will preserve
        amount of memory used before and after filling.

        When ``self.fill_value`` is not NA, the result dtype will be
        ``self.dtype``. Again, this preserves the amount of memory used.
        """
        if ((method is None and value is None) or
                (method is not None and value is not None)):
            raise ValueError("Must specify one of 'method' or 'value'.")

        elif method is not None:
            msg = "fillna with 'method' requires high memory usage."
            warnings.warn(msg, PerformanceWarning)
            filled = interpolate_2d(np.asarray(self), method=method,
                                    limit=limit)
            return type(self)(filled, fill_value=self.fill_value)

        else:
            new_values = np.where(isna(self.sp_values), value, self.sp_values)

            if self._null_fill_value:
                # This is essentially just updating the dtype.
                new_dtype = SparseDtype(self.dtype.subtype, fill_value=value)
            else:
                new_dtype = self.dtype

        return self._simple_new(new_values, self._sparse_index, new_dtype)

    def shift(self, periods=1):

        if periods == 0:
            return self.copy()

        subtype = np.result_type(np.nan, self.dtype.subtype)

        if subtype != self.dtype.subtype:
            # just coerce up front
            arr = self.astype(SparseDtype(subtype, self.fill_value))
        else:
            arr = self

        empty = self._from_sequence([self.dtype.na_value] * abs(periods),
                                    dtype=arr.dtype)
        if periods > 0:
            a = empty
            b = arr[:-periods]
        else:
            a = arr[abs(periods):]
            b = empty
        return arr._concat_same_type([a, b])

    def _first_fill_value_loc(self):
        """
        Get the location of the first missing value.

        Returns
        -------
        int
        """
        if len(self) == 0 or self.sp_index.npoints == len(self):
            return -1

        indices = self.sp_index.to_int_index().indices
        if indices[0] > 0:
            return 0

        diff = indices[1:] - indices[:-1]
        return np.searchsorted(diff, 2) + 1

    def unique(self):
        uniques = list(pd.unique(self.sp_values))
        fill_loc = self._first_fill_value_loc()
        if fill_loc >= 0:
            uniques.insert(fill_loc, self.fill_value)
        return type(self)._from_sequence(uniques, dtype=self.dtype)

    def _values_for_factorize(self):
        # Still override this for hash_pandas_object
        return np.asarray(self), self.fill_value

    def factorize(self, na_sentinel=-1):
        # Currently, ExtensionArray.factorize -> Tuple[ndarray, EA]
        # The sparsity on this is backwards from what Sparse would want. Want
        # ExtensionArray.factorize -> Tuple[EA, EA]
        # Given that we have to return a dense array of labels, why bother
        # implementing an efficient factorize?
        labels, uniques = pd.factorize(np.asarray(self),
                                       na_sentinel=na_sentinel)
        uniques = SparseArray(uniques, dtype=self.dtype)
        return labels, uniques

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
            # TODO: I think we can avoid densifying when masking a
            # boolean SparseArray with another. Need to look at the
            # key's fill_value for True / False, and then do an intersection
            # on the indicies of the sp_values.
            if isinstance(key, SparseArray):
                if is_bool_dtype(key):
                    key = key.to_dense()
                else:
                    key = np.asarray(key)

            if com.is_bool_indexer(key) and len(self) == len(key):
                # TODO(numpy 1.11): Remove this asarray.
                # Old NumPy didn't treat array-like as boolean masks.
                key = np.asarray(key)
                return self.take(np.arange(len(key), dtype=np.int32)[key])
            elif hasattr(key, '__len__'):
                return self.take(key)
            else:
                raise ValueError("Cannot slice with '{}'".format(key))

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

    def take(self, indices, allow_fill=False, fill_value=None):
        if is_scalar(indices):
            raise ValueError("'indices' must be an array, not a "
                             "scalar '{}'.".format(indices))
        indices = np.asarray(indices, dtype=np.int32)

        if indices.size == 0:
            result = []
            kwargs = {'dtype': self.dtype}
        elif allow_fill:
            result = self._take_with_fill(indices, fill_value=fill_value)
            kwargs = {}
        else:
            result = self._take_without_fill(indices)
            kwargs = {'dtype': self.dtype}

        return type(self)(result, fill_value=self.fill_value, kind=self.kind,
                          **kwargs)

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
                dtype = np.result_type(self.sp_values, fill_value)
                taken = np.empty_like(indices, dtype=dtype)
                taken.fill(fill_value)
                return taken
            else:
                raise IndexError('cannot do a non-empty take from an empty '
                                 'axes.')

        sp_indexer = self.sp_index.lookup_array(indices)

        if self.sp_index.npoints == 0:
            # Avoid taking from the empty self.sp_values
            taken = np.full(sp_indexer.shape, fill_value=fill_value,
                            dtype=np.result_type(fill_value))
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
            out = np.full(indices.shape, self.fill_value,
                          dtype=np.result_type(self.fill_value))
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

        return self._simple_new(values, self.sp_index, self.dtype)

    @classmethod
    def _concat_same_type(cls, to_concat):
        fill_values = [x.fill_value for x in to_concat]

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
        """
        Change the dtype of a SparseArray.

        The output will always be a SparseArray. To convert to a dense
        ndarray with a certain dtype, use :meth:`numpy.asarray`.

        Parameters
        ----------
        dtype : np.dtype or ExtensionDtype
            For SparseDtype, this changes the dtype of
            ``self.sp_values`` and the ``self.fill_value``.

            For other dtypes, this only changes the dtype of
            ``self.sp_values``.

        copy : bool, default True
            Whether to ensure a copy is made, even if not necessary.

        Returns
        -------
        SparseArray

        Examples
        --------
        >>> arr = SparseArray([0, 0, 1, 2])
        >>> arr
        [0, 0, 1, 2]
        Fill: 0
        IntIndex
        Indices: array([2, 3], dtype=int32)

        >>> arr.astype(np.dtype('int32'))
        [0, 0, 1, 2]
        Fill: 0
        IntIndex
        Indices: array([2, 3], dtype=int32)

        Using a NumPy dtype with a different kind (e.g. float) will coerce
        just ``self.sp_values``.

        >>> arr.astype(np.dtype('float64'))
        ... # doctest: +NORMALIZE_WHITESPACE
        [0, 0, 1.0, 2.0]
        Fill: 0
        IntIndex
        Indices: array([2, 3], dtype=int32)

        Use a SparseDtype if you wish to be change the fill value as well.

        >>> arr.astype(SparseDtype("float64", fill_value=np.nan))
        ... # doctest: +NORMALIZE_WHITESPACE
        [nan, nan, 1.0, 2.0]
        Fill: nan
        IntIndex
        Indices: array([2, 3], dtype=int32)
        """
        dtype = pandas_dtype(dtype)

        if not isinstance(dtype, SparseDtype):
            dtype = SparseDtype(dtype, fill_value=self.fill_value)

        sp_values = astype_nansafe(self.sp_values,
                                   dtype.subtype,
                                   copy=copy)
        if sp_values is self.sp_values and copy:
            sp_values = sp_values.copy()

        return self._simple_new(sp_values,
                                self.sp_index,
                                dtype)

    def map(self, mapper):
        """
        Map categories using input correspondence (dict, Series, or function).

        Parameters
        ----------
        mapper : dict, Series, callable
            The correspondence from old values to new.

        Returns
        -------
        SparseArray
            The output array will have the same density as the input.
            The output fill value will be the result of applying the
            mapping to ``self.fill_value``

        Examples
        --------
        >>> arr = pd.SparseArray([0, 1, 2])
        >>> arr.apply(lambda x: x + 10)
        [10, 11, 12]
        Fill: 10
        IntIndex
        Indices: array([1, 2], dtype=int32)

        >>> arr.apply({0: 10, 1: 11, 2: 12})
        [10, 11, 12]
        Fill: 10
        IntIndex
        Indices: array([1, 2], dtype=int32)

        >>> arr.apply(pd.Series([10, 11, 12], index=[0, 1, 2]))
        [10, 11, 12]
        Fill: 10
        IntIndex
        Indices: array([1, 2], dtype=int32)
        """
        # this is used in apply.
        # We get hit since we're an "is_extension_type" but regular extension
        # types are not hit. This may be worth adding to the interface.
        if isinstance(mapper, ABCSeries):
            mapper = mapper.to_dict()

        if isinstance(mapper, compat.Mapping):
            fill_value = mapper.get(self.fill_value, self.fill_value)
            sp_values = [mapper.get(x, None) for x in self.sp_values]
        else:
            fill_value = mapper(self.fill_value)
            sp_values = [mapper(x) for x in self.sp_values]

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
            nd_state, (fill_value, sp_index) = state
            sparse_values = np.array([])
            sparse_values.__setstate__(nd_state)

            self._sparse_values = sparse_values
            self._sparse_index = sp_index
            self._dtype = SparseDtype(sparse_values.dtype, fill_value)
        else:
            self.__dict__.update(state)

    def nonzero(self):
        if self.fill_value == 0:
            return self.sp_index.to_int_index().indices,
        else:
            return self.sp_index.to_int_index().indices[self.sp_values != 0],

    # ------------------------------------------------------------------------
    # Reductions
    # ------------------------------------------------------------------------

    def _reduce(self, name, skipna=True, **kwargs):
        method = getattr(self, name, None)

        if method is None:
            raise TypeError("cannot perform {name} with type {dtype}".format(
                name=name, dtype=self.dtype))

        if skipna:
            arr = self
        else:
            arr = self.dropna()

        # we don't support these kwargs.
        # They should only be present when called via pandas, so do it here.
        # instead of in `any` / `all` (which will raise if they're present,
        # thanks to nv.validate
        kwargs.pop('filter_type', None)
        kwargs.pop('numeric_only', None)
        kwargs.pop('op', None)
        return getattr(arr, name)(**kwargs)

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

        return values.any().item()

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

    def __array_wrap__(self, array, context=None):
        from pandas.core.dtypes.generic import ABCSparseSeries

        ufunc, inputs, _ = context
        inputs = tuple(x.values if isinstance(x, ABCSparseSeries) else x
                       for x in inputs)
        return self.__array_ufunc__(ufunc, '__call__', *inputs)

    _HANDLED_TYPES = (np.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.get('out', ())

        for x in inputs + out:
            if not isinstance(x, self._HANDLED_TYPES + (SparseArray,)):
                return NotImplemented

        special = {'add', 'sub', 'mul', 'pow', 'mod', 'floordiv', 'truediv',
                   'divmod', 'eq', 'ne', 'lt', 'gt', 'le', 'ge', 'remainder'}
        if compat.PY2:
            special.add('div')
        aliases = {
            'subtract': 'sub',
            'multiply': 'mul',
            'floor_divide': 'floordiv',
            'true_divide': 'truediv',
            'power': 'pow',
            'remainder': 'mod',
            'divide': 'div',
            'equal': 'eq',
            'not_equal': 'ne',
            'less': 'lt',
            'less_equal': 'le',
            'greater': 'gt',
            'greater_equal': 'ge',
        }

        flipped = {
            'lt': '__gt__',
            'le': '__ge__',
            'gt': '__lt__',
            'ge': '__le__',
            'eq': '__eq__',
            'ne': '__ne__',
        }

        op_name = ufunc.__name__
        op_name = aliases.get(op_name, op_name)

        if op_name in special and kwargs.get('out') is None:
            if isinstance(inputs[0], type(self)):
                return getattr(self, '__{}__'.format(op_name))(inputs[1])
            else:
                name = flipped.get(op_name, '__r{}__'.format(op_name))
                return getattr(self, name)(inputs[0])

        if len(inputs) == 1:
            # No alignment necessary.
            sp_values = getattr(ufunc, method)(self.sp_values, **kwargs)
            fill_value = getattr(ufunc, method)(self.fill_value, **kwargs)
            return self._simple_new(sp_values,
                                    self.sp_index,
                                    SparseDtype(sp_values.dtype, fill_value))

        result = getattr(ufunc, method)(*[np.asarray(x) for x in inputs],
                                        **kwargs)
        if out:
            if len(out) == 1:
                out = out[0]
            return out

        if type(result) is tuple:
            return tuple(type(self)(x) for x in result)
        elif method == 'at':
            # no return value
            return None
        else:
            return type(self)(result)

    def __abs__(self):
        return np.abs(self)

    # ------------------------------------------------------------------------
    # Ops
    # ------------------------------------------------------------------------

    @classmethod
    def _create_unary_method(cls, op):
        def sparse_unary_method(self):
            fill_value = op(np.array(self.fill_value)).item()
            values = op(self.sp_values)
            dtype = SparseDtype(values.dtype, fill_value)
            return cls._simple_new(values, self.sp_index, dtype)

        name = '__{name}__'.format(name=op.__name__)
        return compat.set_function_name(sparse_unary_method, name, cls)

    @classmethod
    def _create_arithmetic_method(cls, op):
        def sparse_arithmetic_method(self, other):
            op_name = op.__name__

            if isinstance(other, (ABCSeries, ABCIndexClass)):
                # Rely on pandas to dispatch to us.
                return NotImplemented

            if isinstance(other, SparseArray):
                return _sparse_array_op(self, other, op, op_name)

            elif is_scalar(other):
                with np.errstate(all='ignore'):
                    fill = op(_get_fill(self), np.asarray(other))
                    result = op(self.sp_values, other)

                if op_name == 'divmod':
                    left, right = result
                    lfill, rfill = fill
                    return (_wrap_result(op_name, left, self.sp_index, lfill),
                            _wrap_result(op_name, right, self.sp_index, rfill))

                return _wrap_result(op_name, result, self.sp_index, fill)

            else:
                other = np.asarray(other)
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

        name = '__{name}__'.format(name=op.__name__)
        return compat.set_function_name(sparse_arithmetic_method, name, cls)

    @classmethod
    def _create_comparison_method(cls, op):
        def cmp_method(self, other):
            op_name = op.__name__

            if op_name in {'and_', 'or_'}:
                op_name = op_name[:-1]

            if isinstance(other, (ABCSeries, ABCIndexClass)):
                # Rely on pandas to unbox and dispatch to us.
                return NotImplemented

            if not is_scalar(other) and not isinstance(other, type(self)):
                # convert list-like to ndarray
                other = np.asarray(other)

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

    @classmethod
    def _add_unary_ops(cls):
        cls.__pos__ = cls._create_unary_method(operator.pos)
        cls.__neg__ = cls._create_unary_method(operator.neg)
        cls.__invert__ = cls._create_unary_method(operator.invert)

    @classmethod
    def _add_comparison_ops(cls):
        cls.__and__ = cls._create_comparison_method(operator.and_)
        cls.__or__ = cls._create_comparison_method(operator.or_)
        super(SparseArray, cls)._add_comparison_ops()

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
SparseArray._add_unary_ops()


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
    (sparse_values, index, fill_value) : (ndarray, SparseIndex, Scalar)
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
    if length != len(mask):
        # the arr is a SparseArray
        indices = mask.sp_index.indices
    else:
        indices = mask.nonzero()[0].astype(np.int32)

    index = _make_index(length, indices, kind)
    sparsified_values = arr[mask]
    if dtype is not None:
        sparsified_values = astype_nansafe(sparsified_values, dtype=dtype)
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
