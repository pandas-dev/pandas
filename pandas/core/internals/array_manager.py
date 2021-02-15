"""
Experimental manager based on storing a collection of 1D arrays
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, List, Optional, Tuple, TypeVar, Union

import numpy as np

from pandas._libs import algos as libalgos, lib
from pandas._typing import ArrayLike, DtypeObj, Hashable
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import find_common_type, infer_dtype_from_scalar
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_numeric_dtype,
)
from pandas.core.dtypes.dtypes import ExtensionDtype, PandasDtype
from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries
from pandas.core.dtypes.missing import array_equals, isna

import pandas.core.algorithms as algos
from pandas.core.arrays import ExtensionArray
from pandas.core.arrays.sparse import SparseDtype
from pandas.core.construction import extract_array
from pandas.core.indexers import maybe_convert_indices
from pandas.core.indexes.api import Index, ensure_index
from pandas.core.internals.base import DataManager
from pandas.core.internals.blocks import make_block

if TYPE_CHECKING:
    from pandas.core.internals.managers import SingleBlockManager


T = TypeVar("T", bound="ArrayManager")


class ArrayManager(DataManager):
    """
    Core internal data structure to implement DataFrame and Series.

    Alternative to the BlockManager, storing a list of 1D arrays instead of
    Blocks.

    This is *not* a public API class

    Parameters
    ----------
    arrays : Sequence of arrays
    axes : Sequence of Index
    do_integrity_check : bool, default True

    """

    __slots__ = [
        "_axes",  # private attribute, because 'axes' has different order, see below
        "arrays",
    ]

    arrays: List[Union[np.ndarray, ExtensionArray]]
    _axes: List[Index]

    def __init__(
        self,
        arrays: List[Union[np.ndarray, ExtensionArray]],
        axes: List[Index],
        do_integrity_check: bool = True,
    ):
        # Note: we are storing the axes in "_axes" in the (row, columns) order
        # which contrasts the order how it is stored in BlockManager
        self._axes = axes
        self.arrays = arrays

        if do_integrity_check:
            self._axes = [ensure_index(ax) for ax in axes]
            self._verify_integrity()

    def make_empty(self: T, axes=None) -> T:
        """Return an empty ArrayManager with the items axis of len 0 (no columns)"""
        if axes is None:
            axes = [self.axes[1:], Index([])]

        arrays: List[Union[np.ndarray, ExtensionArray]] = []
        return type(self)(arrays, axes)

    @property
    def items(self) -> Index:
        return self._axes[1]

    @property
    def axes(self) -> List[Index]:  # type: ignore[override]
        # mypy doesn't work to override attribute with property
        # see https://github.com/python/mypy/issues/4125
        """Axes is BlockManager-compatible order (columns, rows)"""
        return [self._axes[1], self._axes[0]]

    @property
    def shape(self) -> Tuple[int, ...]:
        # this still gives the BlockManager-compatible transposed shape
        return tuple(len(ax) for ax in self.axes)

    @property
    def shape_proper(self) -> Tuple[int, ...]:
        # this returns (n_rows, n_columns)
        return tuple(len(ax) for ax in self._axes)

    @staticmethod
    def _normalize_axis(axis):
        # switch axis
        axis = 1 if axis == 0 else 0
        return axis

    # TODO can be shared
    def set_axis(self, axis: int, new_labels: Index) -> None:
        # Caller is responsible for ensuring we have an Index object.
        axis = self._normalize_axis(axis)
        old_len = len(self._axes[axis])
        new_len = len(new_labels)

        if new_len != old_len:
            raise ValueError(
                f"Length mismatch: Expected axis has {old_len} elements, new "
                f"values have {new_len} elements"
            )

        self._axes[axis] = new_labels

    def consolidate(self) -> ArrayManager:
        return self

    def is_consolidated(self) -> bool:
        return True

    def _consolidate_inplace(self) -> None:
        pass

    def get_dtypes(self):
        return np.array([arr.dtype for arr in self.arrays], dtype="object")

    # TODO setstate getstate

    def __repr__(self) -> str:
        output = type(self).__name__
        output += f"\nIndex: {self._axes[0]}"
        output += f"\nColumns: {self._axes[1]}"
        output += f"\n{len(self.arrays)} arrays:"
        for arr in self.arrays:
            output += f"\n{arr.dtype}"
        return output

    def _verify_integrity(self) -> None:
        n_rows, n_columns = self.shape_proper
        if not len(self.arrays) == n_columns:
            raise ValueError(
                "Number of passed arrays must equal the size of the column Index: "
                f"{len(self.arrays)} arrays vs {n_columns} columns."
            )
        for arr in self.arrays:
            if not len(arr) == n_rows:
                raise ValueError(
                    "Passed arrays should have the same length as the rows Index: "
                    f"{len(arr)} vs {n_rows} rows"
                )
            if not isinstance(arr, (np.ndarray, ExtensionArray)):
                raise ValueError(
                    "Passed arrays should be np.ndarray or ExtensionArray instances, "
                    f"got {type(arr)} instead"
                )

    def reduce(
        self: T, func: Callable, ignore_failures: bool = False
    ) -> Tuple[T, np.ndarray]:
        # TODO this still fails because `func` assumes to work on 2D arrays
        # TODO implement ignore_failures
        assert self.ndim == 2

        res_arrays = []
        for arr in self.arrays:
            res = func(arr, axis=0)
            res_arrays.append(np.array([res]))

        index = Index([None])  # placeholder
        new_mgr = type(self)(res_arrays, [index, self.items])
        indexer = np.arange(self.shape[0])
        return new_mgr, indexer

    def operate_blockwise(self, other: ArrayManager, array_op) -> ArrayManager:
        """
        Apply array_op blockwise with another (aligned) BlockManager.
        """
        # TODO what if `other` is BlockManager ?
        left_arrays = self.arrays
        right_arrays = other.arrays
        result_arrays = [
            array_op(left, right) for left, right in zip(left_arrays, right_arrays)
        ]
        return type(self)(result_arrays, self._axes)

    def apply(
        self: T,
        f,
        align_keys: Optional[List[str]] = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        """
        Iterate over the arrays, collect and create a new ArrayManager.

        Parameters
        ----------
        f : str or callable
            Name of the Array method to apply.
        align_keys: List[str] or None, default None
        ignore_failures: bool, default False
        **kwargs
            Keywords to pass to `f`

        Returns
        -------
        ArrayManager
        """
        assert "filter" not in kwargs

        align_keys = align_keys or []
        result_arrays: List[np.ndarray] = []
        result_indices: List[int] = []
        # fillna: Series/DataFrame is responsible for making sure value is aligned

        aligned_args = {k: kwargs[k] for k in align_keys}

        if f == "apply":
            f = kwargs.pop("func")

        for i, arr in enumerate(self.arrays):

            if aligned_args:

                for k, obj in aligned_args.items():
                    if isinstance(obj, (ABCSeries, ABCDataFrame)):
                        # The caller is responsible for ensuring that
                        #  obj.axes[-1].equals(self.items)
                        if obj.ndim == 1:
                            kwargs[k] = obj.iloc[i]
                        else:
                            kwargs[k] = obj.iloc[:, i]._values
                    else:
                        # otherwise we have an array-like
                        kwargs[k] = obj[i]

            try:
                if callable(f):
                    applied = f(arr, **kwargs)
                else:
                    applied = getattr(arr, f)(**kwargs)
            except (TypeError, NotImplementedError):
                if not ignore_failures:
                    raise
                continue
            # if not isinstance(applied, ExtensionArray):
            #     # TODO not all EA operations return new EAs (eg astype)
            #     applied = array(applied)
            result_arrays.append(applied)
            result_indices.append(i)

        new_axes: List[Index]
        if ignore_failures:
            # TODO copy?
            new_axes = [self._axes[0], self._axes[1][result_indices]]
        else:
            new_axes = self._axes

        if len(result_arrays) == 0:
            return self.make_empty(new_axes)

        return type(self)(result_arrays, new_axes)

    def apply_with_block(self: T, f, align_keys=None, **kwargs) -> T:

        align_keys = align_keys or []
        aligned_args = {k: kwargs[k] for k in align_keys}

        result_arrays = []

        for i, arr in enumerate(self.arrays):

            if aligned_args:
                for k, obj in aligned_args.items():
                    if isinstance(obj, (ABCSeries, ABCDataFrame)):
                        # The caller is responsible for ensuring that
                        #  obj.axes[-1].equals(self.items)
                        if obj.ndim == 1:
                            kwargs[k] = obj.iloc[[i]]
                        else:
                            kwargs[k] = obj.iloc[:, [i]]._values
                    else:
                        # otherwise we have an ndarray
                        kwargs[k] = obj[[i]]

            if hasattr(arr, "tz") and arr.tz is None:  # type: ignore[union-attr]
                # DatetimeArray needs to be converted to ndarray for DatetimeBlock
                arr = arr._data  # type: ignore[union-attr]
            elif arr.dtype.kind == "m":
                # TimedeltaArray needs to be converted to ndarray for TimedeltaBlock
                arr = arr._data  # type: ignore[union-attr]
            if isinstance(arr, np.ndarray):
                arr = np.atleast_2d(arr)
            block = make_block(arr, placement=slice(0, 1, 1), ndim=2)
            applied = getattr(block, f)(**kwargs)
            if isinstance(applied, list):
                applied = applied[0]
            arr = applied.values
            if isinstance(arr, np.ndarray):
                arr = arr[0, :]
            result_arrays.append(arr)

        return type(self)(result_arrays, self._axes)

    # TODO quantile

    def isna(self, func) -> ArrayManager:
        return self.apply("apply", func=func)

    def where(self, other, cond, align: bool, errors: str, axis: int) -> ArrayManager:
        if align:
            align_keys = ["other", "cond"]
        else:
            align_keys = ["cond"]
            other = extract_array(other, extract_numpy=True)

        return self.apply_with_block(
            "where",
            align_keys=align_keys,
            other=other,
            cond=cond,
            errors=errors,
            axis=axis,
        )

    # TODO what is this used for?
    # def setitem(self, indexer, value) -> ArrayManager:
    #     return self.apply_with_block("setitem", indexer=indexer, value=value)

    def putmask(self, mask, new, align: bool = True):

        if align:
            align_keys = ["new", "mask"]
        else:
            align_keys = ["mask"]
            new = extract_array(new, extract_numpy=True)

        return self.apply_with_block(
            "putmask",
            align_keys=align_keys,
            mask=mask,
            new=new,
        )

    def diff(self, n: int, axis: int) -> ArrayManager:
        return self.apply_with_block("diff", n=n, axis=axis)

    def interpolate(self, **kwargs) -> ArrayManager:
        return self.apply_with_block("interpolate", **kwargs)

    def shift(self, periods: int, axis: int, fill_value) -> ArrayManager:
        if fill_value is lib.no_default:
            fill_value = None

        if axis == 0 and self.ndim == 2:
            # TODO column-wise shift
            raise NotImplementedError

        return self.apply_with_block(
            "shift", periods=periods, axis=axis, fill_value=fill_value
        )

    def fillna(self, value, limit, inplace: bool, downcast) -> ArrayManager:
        # TODO implement downcast
        inplace = validate_bool_kwarg(inplace, "inplace")

        def array_fillna(array, value, limit, inplace):

            mask = isna(array)
            if limit is not None:
                limit = libalgos.validate_limit(None, limit=limit)
                mask[mask.cumsum() > limit] = False

            # TODO could optimize for arrays that cannot hold NAs
            # (like _can_hold_na on Blocks)
            if not inplace:
                array = array.copy()

            # np.putmask(array, mask, value)
            if np.any(mask):
                # TODO allow invalid value if there is nothing to fill?
                array[mask] = value
            return array

        return self.apply(array_fillna, value=value, limit=limit, inplace=inplace)

    def downcast(self) -> ArrayManager:
        return self.apply_with_block("downcast")

    def astype(self, dtype, copy: bool = False, errors: str = "raise") -> ArrayManager:
        return self.apply("astype", dtype=dtype, copy=copy)  # , errors=errors)

    def convert(
        self,
        copy: bool = True,
        datetime: bool = True,
        numeric: bool = True,
        timedelta: bool = True,
    ) -> ArrayManager:
        return self.apply_with_block(
            "convert",
            copy=copy,
            datetime=datetime,
            numeric=numeric,
            timedelta=timedelta,
        )

    def replace(self, value, **kwargs) -> ArrayManager:
        assert np.ndim(value) == 0, value
        # TODO "replace" is right now implemented on the blocks, we should move
        # it to general array algos so it can be reused here
        return self.apply_with_block("replace", value=value, **kwargs)

    def replace_list(
        self: T,
        src_list: List[Any],
        dest_list: List[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> T:
        """ do a list replace """
        inplace = validate_bool_kwarg(inplace, "inplace")

        return self.apply_with_block(
            "_replace_list",
            src_list=src_list,
            dest_list=dest_list,
            inplace=inplace,
            regex=regex,
        )

    def to_native_types(self, **kwargs):
        return self.apply_with_block("to_native_types", **kwargs)

    @property
    def is_mixed_type(self) -> bool:
        return True

    @property
    def is_numeric_mixed_type(self) -> bool:
        return False

    @property
    def any_extension_types(self) -> bool:
        """Whether any of the blocks in this manager are extension blocks"""
        return False  # any(block.is_extension for block in self.blocks)

    @property
    def is_view(self) -> bool:
        """ return a boolean if we are a single block and are a view """
        # TODO what is this used for?
        return False

    @property
    def is_single_block(self) -> bool:
        return False

    def get_bool_data(self, copy: bool = False) -> ArrayManager:
        """
        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """
        mask = np.array([is_bool_dtype(t) for t in self.get_dtypes()], dtype="object")
        arrays = [self.arrays[i] for i in np.nonzero(mask)[0]]
        # TODO copy?
        new_axes = [self._axes[0], self._axes[1][mask]]
        return type(self)(arrays, new_axes)

    def get_numeric_data(self, copy: bool = False) -> ArrayManager:
        """
        Parameters
        ----------
        copy : bool, default False
            Whether to copy the blocks
        """
        mask = np.array([is_numeric_dtype(t) for t in self.get_dtypes()])
        arrays = [self.arrays[i] for i in np.nonzero(mask)[0]]
        # TODO copy?
        new_axes = [self._axes[0], self._axes[1][mask]]
        return type(self)(arrays, new_axes)

    def copy(self: T, deep=True) -> T:
        """
        Make deep or shallow copy of ArrayManager

        Parameters
        ----------
        deep : bool or string, default True
            If False, return shallow copy (do not copy data)
            If 'all', copy data and a deep copy of the index

        Returns
        -------
        BlockManager
        """
        # this preserves the notion of view copying of axes
        if deep:
            # hit in e.g. tests.io.json.test_pandas

            def copy_func(ax):
                return ax.copy(deep=True) if deep == "all" else ax.view()

            new_axes = [copy_func(ax) for ax in self._axes]
        else:
            new_axes = list(self._axes)

        if deep:
            new_arrays = [arr.copy() for arr in self.arrays]
        else:
            new_arrays = self.arrays
        return type(self)(new_arrays, new_axes)

    def as_array(
        self,
        transpose: bool = False,
        dtype=None,
        copy: bool = False,
        na_value=lib.no_default,
    ) -> np.ndarray:
        """
        Convert the blockmanager data into an numpy array.

        Parameters
        ----------
        transpose : bool, default False
            If True, transpose the return array.
        dtype : object, default None
            Data type of the return array.
        copy : bool, default False
            If True then guarantee that a copy is returned. A value of
            False does not guarantee that the underlying data is not
            copied.
        na_value : object, default lib.no_default
            Value to be used as the missing value sentinel.

        Returns
        -------
        arr : ndarray
        """
        if len(self.arrays) == 0:
            arr = np.empty(self.shape, dtype=float)
            return arr.transpose() if transpose else arr

        # We want to copy when na_value is provided to avoid
        # mutating the original object
        copy = copy or na_value is not lib.no_default

        if not dtype:
            dtype = _interleaved_dtype(self.arrays)

        if isinstance(dtype, SparseDtype):
            dtype = dtype.subtype
        elif isinstance(dtype, PandasDtype):
            dtype = dtype.numpy_dtype
        elif is_extension_array_dtype(dtype):
            dtype = "object"
        elif is_dtype_equal(dtype, str):
            dtype = "object"

        result = np.empty(self.shape_proper, dtype=dtype)

        for i, arr in enumerate(self.arrays):
            arr = arr.astype(dtype, copy=copy)
            result[:, i] = arr

        if na_value is not lib.no_default:
            result[isna(result)] = na_value

        return result
        # return arr.transpose() if transpose else arr

    def get_slice(self, slobj: slice, axis: int = 0) -> ArrayManager:
        axis = self._normalize_axis(axis)

        if axis == 0:
            arrays = [arr[slobj] for arr in self.arrays]
        elif axis == 1:
            arrays = self.arrays[slobj]

        new_axes = list(self._axes)
        new_axes[axis] = new_axes[axis][slobj]

        return type(self)(arrays, new_axes, do_integrity_check=False)

    def fast_xs(self, loc: int) -> ArrayLike:
        """
        Return the array corresponding to `frame.iloc[loc]`.

        Parameters
        ----------
        loc : int

        Returns
        -------
        np.ndarray or ExtensionArray
        """
        dtype = _interleaved_dtype(self.arrays)

        if isinstance(dtype, SparseDtype):
            temp_dtype = dtype.subtype
        elif isinstance(dtype, PandasDtype):
            temp_dtype = dtype.numpy_dtype
        elif is_extension_array_dtype(dtype):
            temp_dtype = "object"
        elif is_dtype_equal(dtype, str):
            temp_dtype = "object"
        else:
            temp_dtype = dtype

        result = np.array([arr[loc] for arr in self.arrays], dtype=temp_dtype)
        if isinstance(dtype, ExtensionDtype):
            result = dtype.construct_array_type()._from_sequence(result, dtype=dtype)
        return result

    def iget(self, i: int) -> SingleBlockManager:
        """
        Return the data as a SingleBlockManager.
        """
        from pandas.core.internals.managers import SingleBlockManager

        values = self.arrays[i]
        block = make_block(values, placement=slice(0, len(values)), ndim=1)

        return SingleBlockManager(block, self._axes[0])

    def iget_values(self, i: int) -> ArrayLike:
        """
        Return the data for column i as the values (ndarray or ExtensionArray).
        """
        return self.arrays[i]

    def idelete(self, indexer):
        """
        Delete selected locations in-place (new block and array, same BlockManager)
        """
        to_keep = np.ones(self.shape[0], dtype=np.bool_)
        to_keep[indexer] = False

        self.arrays = [self.arrays[i] for i in np.nonzero(to_keep)[0]]
        self._axes = [self._axes[0], self._axes[1][to_keep]]

    def iset(self, loc: Union[int, slice, np.ndarray], value):
        """
        Set new column(s).

        This changes the ArrayManager in-place, but replaces (an) existing
        column(s), not changing column values in-place).

        Parameters
        ----------
        loc : integer, slice or boolean mask
            Positional location (already bounds checked)
        value : array-like
        """
        # single column -> single integer index
        if lib.is_integer(loc):
            # TODO the extract array should in theory not be needed?
            value = extract_array(value, extract_numpy=True)

            # TODO can we avoid needing to unpack this here? That means converting
            # DataFrame into 1D array when loc is an integer
            if isinstance(value, np.ndarray) and value.ndim == 2:
                assert value.shape[1] == 1
                value = value[0, :]

            assert isinstance(value, (np.ndarray, ExtensionArray))
            assert value.ndim == 1
            assert len(value) == len(self._axes[0])
            self.arrays[loc] = value
            return

        # multiple columns -> convert slice or array to integer indices
        elif isinstance(loc, slice):
            indices = range(
                loc.start if loc.start is not None else 0,
                loc.stop if loc.stop is not None else self.shape_proper[1],
                loc.step if loc.step is not None else 1,
            )
        else:
            assert isinstance(loc, np.ndarray)
            assert loc.dtype == "bool"
            indices = np.nonzero(loc)[0]

        assert value.ndim == 2
        assert value.shape[0] == len(self._axes[0])

        for value_idx, mgr_idx in enumerate(indices):
            value_arr = value[:, value_idx]
            self.arrays[mgr_idx] = value_arr
        return

    def insert(self, loc: int, item: Hashable, value, allow_duplicates: bool = False):
        """
        Insert item at selected position.

        Parameters
        ----------
        loc : int
        item : hashable
        value : array_like
        allow_duplicates: bool
            If False, trying to insert non-unique item will raise

        """
        if not allow_duplicates and item in self.items:
            # Should this be a different kind of error??
            raise ValueError(f"cannot insert {item}, already exists")

        if not isinstance(loc, int):
            raise TypeError("loc must be int")

        # insert to the axis; this could possibly raise a TypeError
        new_axis = self.items.insert(loc, item)

        value = extract_array(value, extract_numpy=True)
        if value.ndim == 2:
            value = value[0, :]
        # TODO self.arrays can be empty
        # assert len(value) == len(self.arrays[0])

        # TODO is this copy needed?
        arrays = self.arrays.copy()
        arrays.insert(loc, value)

        self.arrays = arrays
        self._axes[1] = new_axis

    def reindex_indexer(
        self: T,
        new_axis,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
        # ignored keywords
        consolidate: bool = True,
        only_slice: bool = False,
    ) -> T:
        axis = self._normalize_axis(axis)
        return self._reindex_indexer(
            new_axis, indexer, axis, fill_value, allow_dups, copy
        )

    def _reindex_indexer(
        self: T,
        new_axis,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
    ) -> T:
        """
        Parameters
        ----------
        new_axis : Index
        indexer : ndarray of int64 or None
        axis : int
        fill_value : object, default None
        allow_dups : bool, default False
        copy : bool, default True


        pandas-indexer with -1's only.
        """
        if indexer is None:
            if new_axis is self._axes[axis] and not copy:
                return self

            result = self.copy(deep=copy)
            result._axes = list(self._axes)
            result._axes[axis] = new_axis
            return result

        # some axes don't allow reindexing with dups
        if not allow_dups:
            self._axes[axis]._can_reindex(indexer)

        # if axis >= self.ndim:
        #     raise IndexError("Requested axis not found in manager")

        if axis == 1:
            new_arrays = []
            for i in indexer:
                if i == -1:
                    arr = self._make_na_array(fill_value=fill_value)
                else:
                    arr = self.arrays[i]
                new_arrays.append(arr)

        else:
            new_arrays = [
                algos.take(
                    arr,
                    indexer,
                    allow_fill=True,
                    fill_value=fill_value,
                    # if fill_value is not None else blk.fill_value
                )
                for arr in self.arrays
            ]

        new_axes = list(self._axes)
        new_axes[axis] = new_axis

        return type(self)(new_arrays, new_axes)

    def take(self, indexer, axis: int = 1, verify: bool = True, convert: bool = True):
        """
        Take items along any axis.
        """
        axis = self._normalize_axis(axis)

        indexer = (
            np.arange(indexer.start, indexer.stop, indexer.step, dtype="int64")
            if isinstance(indexer, slice)
            else np.asanyarray(indexer, dtype="int64")
        )

        n = self.shape_proper[axis]
        if convert:
            indexer = maybe_convert_indices(indexer, n)

        if verify:
            if ((indexer == -1) | (indexer >= n)).any():
                raise Exception("Indices must be nonzero and less than the axis length")

        new_labels = self._axes[axis].take(indexer)
        return self._reindex_indexer(
            new_axis=new_labels, indexer=indexer, axis=axis, allow_dups=True
        )

    def _make_na_array(self, fill_value=None):
        if fill_value is None:
            fill_value = np.nan

        dtype, fill_value = infer_dtype_from_scalar(fill_value)
        values = np.empty(self.shape_proper[0], dtype=dtype)
        values.fill(fill_value)
        return values

    def _equal_values(self, other) -> bool:
        """
        Used in .equals defined in base class. Only check the column values
        assuming shape and indexes have already been checked.
        """
        for left, right in zip(self.arrays, other.arrays):
            if not array_equals(left, right):
                return False
        else:
            return True

    def unstack(self, unstacker, fill_value) -> ArrayManager:
        """
        Return a BlockManager with all blocks unstacked..

        Parameters
        ----------
        unstacker : reshape._Unstacker
        fill_value : Any
            fill_value for newly introduced missing values.

        Returns
        -------
        unstacked : BlockManager
        """
        indexer, _ = unstacker._indexer_and_to_sort
        new_indexer = np.full(unstacker.mask.shape, -1)
        new_indexer[unstacker.mask] = indexer
        new_indexer2D = new_indexer.reshape(*unstacker.full_shape)

        new_arrays = []
        for arr in self.arrays:
            for i in range(unstacker.full_shape[1]):
                new_arr = algos.take(
                    arr, new_indexer2D[:, i], allow_fill=True, fill_value=fill_value
                )
                new_arrays.append(new_arr)

        new_index = unstacker.new_index
        new_columns = unstacker.get_new_columns(self._axes[1])
        new_axes = [new_index, new_columns]

        return type(self)(new_arrays, new_axes, do_integrity_check=False)

    # TODO
    # equals
    # to_dict
    # quantile


def _interleaved_dtype(blocks) -> Optional[DtypeObj]:
    """
    Find the common dtype for `blocks`.

    Parameters
    ----------
    blocks : List[Block]

    Returns
    -------
    dtype : np.dtype, ExtensionDtype, or None
        None is returned when `blocks` is empty.
    """
    if not len(blocks):
        return None

    return find_common_type([b.dtype for b in blocks])
