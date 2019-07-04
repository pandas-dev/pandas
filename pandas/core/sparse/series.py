"""
Data structures for sparse float data. Life is made simpler by dealing only
with float64 data
"""
from collections import abc
import warnings

import numpy as np

import pandas._libs.index as libindex
import pandas._libs.sparse as splib
from pandas._libs.sparse import BlockIndex, IntIndex
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution

from pandas.core.dtypes.common import is_integer, is_scalar
from pandas.core.dtypes.generic import ABCSeries, ABCSparseSeries
from pandas.core.dtypes.missing import isna, notna

from pandas.core import generic
from pandas.core.arrays import SparseArray
from pandas.core.arrays.sparse import SparseAccessor
from pandas.core.index import Index
from pandas.core.internals import SingleBlockManager
import pandas.core.ops as ops
from pandas.core.series import Series
from pandas.core.sparse.scipy_sparse import _coo_to_sparse_series, _sparse_series_to_coo

_shared_doc_kwargs = dict(
    axes="index",
    klass="SparseSeries",
    axes_single_arg="{0, 'index'}",
    optional_labels="",
    optional_axis="",
)


depr_msg = """\
SparseSeries is deprecated and will be removed in a future version.
Use a Series with sparse values instead.

    >>> series = pd.Series(pd.SparseArray(...))

See http://pandas.pydata.org/pandas-docs/stable/\
user_guide/sparse.html#migrating for more.
"""


class SparseSeries(Series):
    """Data structure for labeled, sparse floating point data

    .. deprecated:: 0.25.0

       Use a Series with sparse values instead.

    Parameters
    ----------
    data : {array-like, Series, SparseSeries, dict}
        .. versionchanged :: 0.23.0
           If data is a dict, argument order is maintained for Python 3.6
           and later.

    kind : {'block', 'integer'}
    fill_value : float
        Code for missing value. Defaults depends on dtype.
        0 for int dtype, False for bool dtype, and NaN for other dtypes
    sparse_index : {BlockIndex, IntIndex}, optional
        Only if you have one. Mainly used internally

    Notes
    -----
    SparseSeries objects are immutable via the typical Python means. If you
    must change values, convert to dense, make your changes, then convert back
    to sparse
    """

    _subtyp = "sparse_series"

    def __init__(
        self,
        data=None,
        index=None,
        sparse_index=None,
        kind="block",
        fill_value=None,
        name=None,
        dtype=None,
        copy=False,
        fastpath=False,
    ):
        warnings.warn(depr_msg, FutureWarning, stacklevel=2)
        # TODO: Most of this should be refactored and shared with Series
        # 1. BlockManager -> array
        # 2. Series.index, Series.name, index, name reconciliation
        # 3. Implicit reindexing
        # 4. Implicit broadcasting
        # 5. Dict construction
        if data is None:
            data = []
        elif isinstance(data, SingleBlockManager):
            index = data.index
            data = data.blocks[0].values
        elif isinstance(data, (ABCSeries, ABCSparseSeries)):
            index = data.index if index is None else index
            dtype = data.dtype if dtype is None else dtype
            name = data.name if name is None else name

            if index is not None:
                data = data.reindex(index)

        elif isinstance(data, abc.Mapping):
            data, index = Series()._init_dict(data, index=index)

        elif is_scalar(data) and index is not None:
            data = np.full(len(index), fill_value=data)

        super().__init__(
            SparseArray(
                data,
                sparse_index=sparse_index,
                kind=kind,
                dtype=dtype,
                fill_value=fill_value,
                copy=copy,
            ),
            index=index,
            name=name,
            copy=False,
            fastpath=fastpath,
        )

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # avoid infinite recursion for other SparseSeries inputs
        inputs = tuple(x.values if isinstance(x, type(self)) else x for x in inputs)
        result = self.values.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        return self._constructor(
            result,
            index=self.index,
            sparse_index=self.sp_index,
            fill_value=result.fill_value,
            copy=False,
        ).__finalize__(self)

    # unary ops
    # TODO: See if this can be shared
    def __pos__(self):
        result = self.values.__pos__()
        return self._constructor(
            result,
            index=self.index,
            sparse_index=self.sp_index,
            fill_value=result.fill_value,
            copy=False,
        ).__finalize__(self)

    def __neg__(self):
        result = self.values.__neg__()
        return self._constructor(
            result,
            index=self.index,
            sparse_index=self.sp_index,
            fill_value=result.fill_value,
            copy=False,
        ).__finalize__(self)

    def __invert__(self):
        result = self.values.__invert__()
        return self._constructor(
            result,
            index=self.index,
            sparse_index=self.sp_index,
            fill_value=result.fill_value,
            copy=False,
        ).__finalize__(self)

    @property
    def block(self):
        warnings.warn("SparseSeries.block is deprecated.", FutureWarning, stacklevel=2)
        return self._data._block

    @property
    def fill_value(self):
        return self.values.fill_value

    @fill_value.setter
    def fill_value(self, v):
        self.values.fill_value = v

    @property
    def sp_index(self):
        return self.values.sp_index

    @property
    def sp_values(self):
        return self.values.sp_values

    @property
    def npoints(self):
        return self.values.npoints

    @classmethod
    def from_array(
        cls, arr, index=None, name=None, copy=False, fill_value=None, fastpath=False
    ):
        """Construct SparseSeries from array.

        .. deprecated:: 0.23.0
            Use the pd.SparseSeries(..) constructor instead.
        """
        warnings.warn(
            "'from_array' is deprecated and will be removed in a "
            "future version. Please use the pd.SparseSeries(..) "
            "constructor instead.",
            FutureWarning,
            stacklevel=2,
        )
        return cls(
            arr,
            index=index,
            name=name,
            copy=copy,
            fill_value=fill_value,
            fastpath=fastpath,
        )

    @property
    def _constructor(self):
        return SparseSeries

    @property
    def _constructor_expanddim(self):
        from pandas.core.sparse.api import SparseDataFrame

        return SparseDataFrame

    @property
    def kind(self):
        if isinstance(self.sp_index, BlockIndex):
            return "block"
        elif isinstance(self.sp_index, IntIndex):
            return "integer"

    def as_sparse_array(self, kind=None, fill_value=None, copy=False):
        """ return my self as a sparse array, do not copy by default """

        if fill_value is None:
            fill_value = self.fill_value
        if kind is None:
            kind = self.kind
        return SparseArray(
            self.values,
            sparse_index=self.sp_index,
            fill_value=fill_value,
            kind=kind,
            copy=copy,
        )

    def __repr__(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Sparse")
            series_rep = Series.__repr__(self)
            rep = "{series}\n{index!r}".format(series=series_rep, index=self.sp_index)
            return rep

    def _reduce(
        self, op, name, axis=0, skipna=True, numeric_only=None, filter_type=None, **kwds
    ):
        """ perform a reduction operation """
        return op(self.array.to_dense(), skipna=skipna, **kwds)

    def __getstate__(self):
        # pickling
        return dict(
            _typ=self._typ,
            _subtyp=self._subtyp,
            _data=self._data,
            fill_value=self.fill_value,
            name=self.name,
        )

    def _unpickle_series_compat(self, state):

        nd_state, own_state = state

        # recreate the ndarray
        data = np.empty(nd_state[1], dtype=nd_state[2])
        np.ndarray.__setstate__(data, nd_state)

        index, fill_value, sp_index = own_state[:3]
        name = None
        if len(own_state) > 3:
            name = own_state[3]

        # create a sparse array
        if not isinstance(data, SparseArray):
            data = SparseArray(
                data, sparse_index=sp_index, fill_value=fill_value, copy=False
            )

        # recreate
        data = SingleBlockManager(data, index, fastpath=True)
        generic.NDFrame.__init__(self, data)

        self._set_axis(0, index)
        self.name = name

    def _set_subtyp(self, is_all_dates):
        if is_all_dates:
            object.__setattr__(self, "_subtyp", "sparse_time_series")
        else:
            object.__setattr__(self, "_subtyp", "sparse_series")

    def _ixs(self, i, axis=0):
        """
        Return the i-th value or values in the SparseSeries by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Returns
        -------
        value : scalar (int) or Series (slice, sequence)
        """
        label = self.index[i]
        if isinstance(label, Index):
            return self.take(i, axis=axis)
        else:
            return self._get_val_at(i)

    def _get_val_at(self, loc):
        """ forward to the array """
        return self.values._get_val_at(loc)

    def __getitem__(self, key):
        # TODO: Document difference from Series.__getitem__, deprecate,
        # and remove!
        if is_integer(key) and key not in self.index:
            return self._get_val_at(key)
        else:
            return super().__getitem__(key)

    def _get_values(self, indexer):
        try:
            return self._constructor(
                self._data.get_slice(indexer), fastpath=True
            ).__finalize__(self)
        except Exception:
            return self[indexer]

    def _set_with_engine(self, key, value):
        return self._set_value(key, value)

    def abs(self):
        """
        Return an object with absolute value taken. Only applicable to objects
        that are all numeric

        Returns
        -------
        abs: same type as caller
        """
        return self._constructor(np.abs(self.values), index=self.index).__finalize__(
            self
        )

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

    def get_value(self, label, takeable=False):
        """
        Retrieve single value at passed index label

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        index : label
        takeable : interpret the index as indexers, default False

        Returns
        -------
        value : scalar value
        """
        warnings.warn(
            "get_value is deprecated and will be removed "
            "in a future release. Please use "
            ".at[] or .iat[] accessors instead",
            FutureWarning,
            stacklevel=2,
        )

        return self._get_value(label, takeable=takeable)

    def _get_value(self, label, takeable=False):
        loc = label if takeable is True else self.index.get_loc(label)
        return self._get_val_at(loc)

    _get_value.__doc__ = get_value.__doc__

    def set_value(self, label, value, takeable=False):
        """
        Quickly set single value at passed label. If label is not contained, a
        new object is created with the label placed at the end of the result
        index

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed
        value : object
            Scalar value
        takeable : interpret the index as indexers, default False

        Notes
        -----
        This method *always* returns a new object. It is not particularly
        efficient but is provided for API compatibility with Series

        Returns
        -------
        series : SparseSeries
        """
        warnings.warn(
            "set_value is deprecated and will be removed "
            "in a future release. Please use "
            ".at[] or .iat[] accessors instead",
            FutureWarning,
            stacklevel=2,
        )
        return self._set_value(label, value, takeable=takeable)

    def _set_value(self, label, value, takeable=False):
        values = self.to_dense()

        # if the label doesn't exist, we will create a new object here
        # and possibly change the index
        new_values = values._set_value(label, value, takeable=takeable)
        if new_values is not None:
            values = new_values
        new_index = values.index
        values = SparseArray(values, fill_value=self.fill_value, kind=self.kind)
        self._data = SingleBlockManager(values, new_index)
        self._index = new_index

    _set_value.__doc__ = set_value.__doc__

    def _set_values(self, key, value):

        # this might be inefficient as we have to recreate the sparse array
        # rather than setting individual elements, but have to convert
        # the passed slice/boolean that's in dense space into a sparse indexer
        # not sure how to do that!
        if isinstance(key, Series):
            key = key.values

        values = self.values.to_dense()
        values[key] = libindex.convert_scalar(values, value)
        values = SparseArray(values, fill_value=self.fill_value, kind=self.kind)
        self._data = SingleBlockManager(values, self.index)

    def to_dense(self):
        """
        Convert SparseSeries to a Series.

        Returns
        -------
        s : Series
        """
        return Series(self.values.to_dense(), index=self.index, name=self.name)

    @property
    def density(self):
        return self.values.density

    def copy(self, deep=True):
        """
        Make a copy of the SparseSeries. Only the actual sparse values need to
        be copied
        """
        # TODO: https://github.com/pandas-dev/pandas/issues/22314
        # We skip the block manager till that is resolved.
        new_data = self.values
        if deep:
            new_data = new_data.copy()
        return self._constructor(
            new_data,
            sparse_index=self.sp_index,
            fill_value=self.fill_value,
            index=self.index.copy(),
            name=self.name,
        ).__finalize__(self)

    @Substitution(**_shared_doc_kwargs)
    @Appender(generic.NDFrame.reindex.__doc__)
    def reindex(self, index=None, method=None, copy=True, limit=None, **kwargs):
        # TODO: remove?
        return super().reindex(
            index=index, method=method, copy=copy, limit=limit, **kwargs
        )

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
        if not isinstance(new_index, splib.SparseIndex):
            raise TypeError("new index must be a SparseIndex")
        values = self.values
        values = values.sp_index.to_int_index().reindex(
            values.sp_values.astype("float64"), values.fill_value, new_index
        )
        values = SparseArray(
            values, sparse_index=new_index, fill_value=self.values.fill_value
        )
        return self._constructor(values, index=self.index).__finalize__(self)

    def cumsum(self, axis=0, *args, **kwargs):
        """
        Cumulative sum of non-NA/null values.

        When performing the cumulative summation, any non-NA/null values will
        be skipped. The resulting SparseSeries will preserve the locations of
        NaN values, but the fill value will be `np.nan` regardless.

        Parameters
        ----------
        axis : {0}

        Returns
        -------
        cumsum : SparseSeries
        """
        nv.validate_cumsum(args, kwargs)
        # Validate axis
        if axis is not None:
            self._get_axis_number(axis)

        new_array = self.values.cumsum()

        return self._constructor(
            new_array, index=self.index, sparse_index=new_array.sp_index
        ).__finalize__(self)

    # TODO: SparseSeries.isna is Sparse, while Series.isna is dense
    @Appender(generic._shared_docs["isna"] % _shared_doc_kwargs)
    def isna(self):
        arr = SparseArray(
            isna(self.values.sp_values),
            sparse_index=self.values.sp_index,
            fill_value=isna(self.fill_value),
        )
        return self._constructor(arr, index=self.index).__finalize__(self)

    isnull = isna

    @Appender(generic._shared_docs["notna"] % _shared_doc_kwargs)
    def notna(self):
        arr = SparseArray(
            notna(self.values.sp_values),
            sparse_index=self.values.sp_index,
            fill_value=notna(self.fill_value),
        )
        return self._constructor(arr, index=self.index).__finalize__(self)

    notnull = notna

    def dropna(self, axis=0, inplace=False, **kwargs):
        """
        Analogous to Series.dropna. If fill_value=NaN, returns a dense Series
        """
        # TODO: make more efficient
        # Validate axis
        self._get_axis_number(axis or 0)
        dense_valid = self.to_dense().dropna()
        if inplace:
            raise NotImplementedError(
                "Cannot perform inplace dropna" " operations on a SparseSeries"
            )
        if isna(self.fill_value):
            return dense_valid
        else:
            dense_valid = dense_valid[dense_valid != self.fill_value]
            return dense_valid.to_sparse(fill_value=self.fill_value)

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

    @Appender(SparseAccessor.to_coo.__doc__)
    def to_coo(self, row_levels=(0,), column_levels=(1,), sort_labels=False):
        A, rows, columns = _sparse_series_to_coo(
            self, row_levels, column_levels, sort_labels=sort_labels
        )
        return A, rows, columns

    @classmethod
    @Appender(SparseAccessor.from_coo.__doc__)
    def from_coo(cls, A, dense_index=False):
        return _coo_to_sparse_series(A, dense_index=dense_index)


# overwrite series methods with unaccelerated Sparse-specific versions
ops.add_flex_arithmetic_methods(SparseSeries)
ops.add_special_arithmetic_methods(SparseSeries)
