from __future__ import annotations

import functools
from itertools import product
from numbers import Integral

import numpy as np
from toolz import pluck

from dask._task_spec import Alias, DataNode, Task, TaskRef
from dask.array._array_expr._expr import ArrayExpr
from dask.array.chunk import getitem
from dask.array.slicing import (
    _slice_1d,
    check_index,
    expander,
    new_blockdim,
    normalize_slice,
    posify_index,
    replace_ellipsis,
    sanitize_index,
)
from dask.array.utils import meta_from_array
from dask.layers import ArrayBlockwiseDep
from dask.tokenize import tokenize
from dask.utils import cached_cumsum, is_arraylike


def slice_with_int_dask_array(x, index):
    """Slice x with at most one 1D dask arrays of ints.

    This is a helper function of :meth:`Array.__getitem__`.

    Parameters
    ----------
    x: Array
    index: tuple with as many elements as x.ndim, among which there are
           one or more Array's with dtype=int

    Returns
    -------
    tuple of (sliced x, new index)

    where the new index is the same as the input, but with slice(None)
    replaced to the original slicer where a 1D filter has been applied and
    one less element where a zero-dimensional filter has been applied.
    """
    from dask.array._array_expr._collection import Array

    assert len(index) == x.ndim
    fancy_indexes = [
        isinstance(idx, (tuple, list))
        or (isinstance(idx, (np.ndarray, Array)) and idx.ndim > 0)
        for idx in index
    ]
    if sum(fancy_indexes) > 1:
        raise NotImplementedError("Don't yet support nd fancy indexing")

    out_index = []
    dropped_axis_cnt = 0
    for in_axis, idx in enumerate(index):
        out_axis = in_axis - dropped_axis_cnt
        if isinstance(idx, Array) and idx.dtype.kind in "iu":
            if idx.ndim == 0:
                idx = idx[np.newaxis]
                x = slice_with_int_dask_array_on_axis(x, idx, out_axis, in_axis)
                x = x[tuple(0 if i == out_axis else slice(None) for i in range(x.ndim))]
                dropped_axis_cnt += 1
            elif idx.ndim == 1:
                x = slice_with_int_dask_array_on_axis(x, idx, out_axis, in_axis)
                out_index.append(slice(None))
            else:
                raise NotImplementedError(
                    "Slicing with dask.array of ints only permitted when "
                    "the indexer has zero or one dimensions"
                )
        else:
            out_index.append(idx)
    return x, tuple(out_index)


def normalize_index(idx, shape):
    """Normalize slicing indexes

    1.  Replaces ellipses with many full slices
    2.  Adds full slices to end of index
    3.  Checks bounding conditions
    4.  Replace multidimensional numpy arrays with dask arrays
    5.  Replaces numpy arrays with lists
    6.  Posify's integers and lists
    7.  Normalizes slices to canonical form

    Examples
    --------
    >>> normalize_index(1, (10,))
    (1,)
    >>> normalize_index(-1, (10,))
    (9,)
    >>> normalize_index([-1], (10,))
    (array([9]),)
    >>> normalize_index(slice(-3, 10, 1), (10,))
    (slice(7, None, None),)
    >>> normalize_index((Ellipsis, None), (10,))
    (slice(None, None, None), None)
    >>> normalize_index(np.array([[True, False], [False, True], [True, True]]), (3, 2))
    (dask.array<array, shape=(3, 2), dtype=bool, chunksize=(3, 2), chunktype=numpy.ndarray>,)
    """
    from dask.array._array_expr._collection import Array, from_array

    if not isinstance(idx, tuple):
        idx = (idx,)

    # if a > 1D numpy.array is provided, cast it to a dask array
    if len(idx) > 0 and len(shape) > 1:
        i = idx[0]
        if is_arraylike(i) and not isinstance(i, Array) and i.shape == shape:
            idx = (from_array(i), *idx[1:])

    idx = replace_ellipsis(len(shape), idx)
    n_sliced_dims = 0
    for i in idx:
        if hasattr(i, "ndim") and i.ndim >= 1:
            n_sliced_dims += i.ndim
        elif i is None:
            continue
        else:
            n_sliced_dims += 1

    idx = idx + (slice(None),) * (len(shape) - n_sliced_dims)
    if len([i for i in idx if i is not None]) > len(shape):
        raise IndexError("Too many indices for array")

    none_shape = []
    i = 0
    for ind in idx:
        if ind is not None:
            none_shape.append(shape[i])
            i += 1
        else:
            none_shape.append(None)

    for axis, (i, d) in enumerate(zip(idx, none_shape)):
        if d is not None:
            check_index(axis, i, d)
    idx = tuple(map(sanitize_index, idx))
    idx = tuple(map(normalize_slice, idx, none_shape))
    idx = posify_index(none_shape, idx)
    return idx


def slice_with_int_dask_array_on_axis(x, idx, axis, in_axis):
    """Slice a ND dask array with a 1D dask arrays of ints along the given
    axis.

    This is a helper function of :func:`slice_with_int_dask_array`.
    """
    from dask.array import chunk
    from dask.array._array_expr._collection import blockwise

    assert 0 <= axis < x.ndim

    if np.isnan(x.chunks[axis]).any():
        raise NotImplementedError(
            "Slicing an array with unknown chunks with "
            "a dask.array of ints is not supported"
        )
    x_axes = tuple(range(x.ndim))
    idx_axes = (x.ndim,)  # arbitrary index not already in x_axes

    # Calculate the offset at which each chunk starts along axis
    # e.g. chunks=(..., (5, 3, 4), ...) -> offset=[0, 5, 8]
    offset = np.roll(np.cumsum(np.asarray(x.chunks[axis], like=x._meta)), 1)
    offset[0] = 0
    offset = ArrayOffsetDep(x.chunks, offset, in_axis)
    # Define axis labels for blockwise

    p_axes = x_axes[: axis + 1] + idx_axes + x_axes[axis + 1 :]
    y_axes = x_axes[:axis] + idx_axes + x_axes[axis + 1 :]

    # Calculate the cartesian product of every chunk of x vs every chunk of idx
    p = blockwise(
        chunk.slice_with_int_dask_array,
        p_axes,
        x,
        x_axes,
        idx,
        idx_axes,
        offset,
        p_axes,
        x_size=x.shape[axis],
        axis=axis,
        dtype=x.dtype,
        meta=x._meta,
    )

    # Aggregate on the chunks of x along axis
    y = blockwise(
        chunk.slice_with_int_dask_array_aggregate,
        y_axes,
        idx,
        idx_axes,
        p,
        p_axes,
        concatenate=True,
        x_chunks=x.chunks[axis],
        axis=axis,
        dtype=x.dtype,
        meta=x._meta,
    )
    return y


class ArrayOffsetDep(ArrayBlockwiseDep):
    def __init__(
        self, chunks: tuple[tuple[int, ...], ...], values: np.ndarray | dict, axis: int
    ):
        super().__init__(chunks)
        self.values = values
        self.axis = axis

    def __getitem__(self, idx: tuple):
        return self.values[idx[self.axis]]


def slice_array(x, index):
    """
    slice_with_newaxis : handle None/newaxis case
    slice_wrap_lists : handle fancy indexing with lists
    slice_slices_and_integers : handle everything else
    """
    if all(
        isinstance(index, slice) and index == slice(None, None, None) for index in index
    ):
        # all none slices
        return x.expr

    # Add in missing colons at the end as needed.  x[5] -> x[5, :, :]
    not_none_count = sum(i is not None for i in index)
    missing = len(x.chunks) - not_none_count
    index += (slice(None, None, None),) * missing
    return slice_with_newaxes(x, index)


def slice_with_newaxes(x, index):
    """
    Handle indexing with Nones

    Strips out Nones then hands off to slice_wrap_lists
    """
    # Strip Nones from index
    index2 = tuple(ind for ind in index if ind is not None)
    where_none = [i for i, ind in enumerate(index) if ind is None]
    for i, xx in enumerate(where_none):
        n = sum(isinstance(ind, Integral) for ind in index[:xx])
        if n:
            where_none[i] -= n

    # Pass down and do work
    x = slice_wrap_lists(x, index2, not where_none)

    if where_none:
        return SlicesWrapNone(
            x.array, x.index, x.allow_getitem_optimization, where_none
        )

    else:
        return x


def slice_wrap_lists(x, index, allow_getitem_optimization):
    """
    Fancy indexing along blocked array dasks

    Handles index of type list.  Calls slice_slices_and_integers for the rest

    See Also
    --------

    take : handle slicing with lists ("fancy" indexing)
    slice_slices_and_integers : handle slicing with slices and integers
    """
    assert all(isinstance(i, (slice, list, Integral)) or is_arraylike(i) for i in index)
    if not len(x.chunks) == len(index):
        raise IndexError("Too many indices for array")

    # Do we have more than one list in the index?
    where_list = [
        i for i, ind in enumerate(index) if is_arraylike(ind) and ind.ndim > 0
    ]
    if len(where_list) > 1:
        raise NotImplementedError("Don't yet support nd fancy indexing")
    # Is the single list an empty list? In this case just treat it as a zero
    # length slice
    if where_list and not index[where_list[0]].size:
        index = list(index)
        index[where_list.pop()] = slice(0, 0, 1)
        index = tuple(index)

    # No lists, hooray! just use slice_slices_and_integers
    if not where_list:
        return slice_slices_and_integers(x, index, allow_getitem_optimization)

    # Replace all lists with full slices  [3, 1, 0] -> slice(None, None, None)
    index_without_list = tuple(
        slice(None, None, None) if is_arraylike(i) else i for i in index
    )

    # lists and full slices.  Just use take
    if all(is_arraylike(i) or i == slice(None, None, None) for i in index):
        axis = where_list[0]
        x = take(x, index[where_list[0]], axis=axis)
    # Mixed case. Both slices/integers and lists. slice/integer then take
    else:
        x = slice_slices_and_integers(
            x,
            index_without_list,
            allow_getitem_optimization=False,
        )
        axis = where_list[0]
        axis2 = axis - sum(
            1 for i, ind in enumerate(index) if i < axis and isinstance(ind, Integral)
        )
        x = take(x, index[axis], axis=axis2)

    return x


def slice_slices_and_integers(x, index, allow_getitem_optimization=False):
    from dask.array.core import unknown_chunk_message

    shape = tuple(cached_cumsum(dim, initial_zero=True)[-1] for dim in x.chunks)

    for dim, ind in zip(shape, index):
        if np.isnan(dim) and ind != slice(None, None, None):
            raise ValueError(
                f"Arrays chunk sizes are unknown: {shape}{unknown_chunk_message}"
            )
    assert all(isinstance(ind, (slice, Integral)) for ind in index)
    return SliceSlicesIntegers(x, index, allow_getitem_optimization)


def take(x, index, axis=0):
    if not np.isnan(x.chunks[axis]).any():
        from dask.array._array_expr._shuffle import _shuffle
        from dask.array.utils import arange_safe, asarray_safe

        arange = arange_safe(np.sum(x.chunks[axis]), like=index)
        if len(index) == len(arange) and np.abs(index - arange).sum() == 0:
            # no-op
            return x

        average_chunk_size = int(sum(x.chunks[axis]) / len(x.chunks[axis]))

        indexer = []
        index = asarray_safe(index, like=index)
        for i in range(0, len(index), average_chunk_size):
            indexer.append(index[i : i + average_chunk_size].tolist())
        return _shuffle(x, indexer, axis, "getitem-")
    elif len(x.chunks[axis]) == 1:
        return TakeUnknownOneChunk(x, index, axis)
    else:
        from dask.array.core import unknown_chunk_message

        raise ValueError(
            f"Array chunk size or shape is unknown. {unknown_chunk_message}"
        )


class Slice(ArrayExpr):
    @functools.cached_property
    def _name(self):
        return f"getitem-{self.deterministic_token}"

    @functools.cached_property
    def _meta(self):
        meta = meta_from_array(self.array._meta, ndim=len(self.chunks))
        if np.isscalar(meta):
            meta = np.array(meta)
        return meta


class SliceSlicesIntegers(Slice):
    _parameters = ["array", "index", "allow_getitem_optimization"]

    @functools.cached_property
    def chunks(self):
        new_blockdims = [
            new_blockdim(d, db, i)
            for d, i, db in zip(self.array.shape, self.index, self.array.chunks)
            if not isinstance(i, Integral)
        ]
        return tuple(map(tuple, new_blockdims))

    def _layer(self) -> dict:
        # Get a list (for each dimension) of dicts{blocknum: slice()}
        block_slices = list(
            map(_slice_1d, self.array.shape, self.array.chunks, self.index)
        )
        sorted_block_slices = [sorted(i.items()) for i in block_slices]

        # (in_name, 1, 1, 2), (in_name, 1, 1, 4), (in_name, 2, 1, 2), ...
        in_names = list(
            product([self.array._name], *[pluck(0, s) for s in sorted_block_slices])
        )

        # (out_name, 0, 0, 0), (out_name, 0, 0, 1), (out_name, 0, 1, 0), ...
        out_names = list(
            product(
                [self._name],
                *[
                    range(len(d))[::-1] if i.step and i.step < 0 else range(len(d))
                    for d, i in zip(block_slices, self.index)
                    if not isinstance(i, Integral)
                ],
            )
        )

        all_slices = list(product(*[pluck(1, s) for s in sorted_block_slices]))

        dsk_out = {
            out_name: (
                Task(out_name, getitem, TaskRef(in_name), slices)
                if not self.allow_getitem_optimization
                or not all(sl == slice(None, None, None) for sl in slices)
                else Alias(out_name, in_name)
            )
            for out_name, in_name, slices in zip(out_names, in_names, all_slices)
        }
        return dsk_out


class SlicesWrapNone(SliceSlicesIntegers):
    _parameters = ["array", "index", "allow_getitem_optimization", "where_none"]

    @functools.cached_property
    def chunks(self):
        return self.expand(super().chunks, (1,))

    @functools.cached_property
    def expand(self):
        return expander(self.where_none)

    def _layer(self) -> dict:
        dsk = super()._layer()

        where_none_orig = list(self.where_none)
        expand_orig = expander(where_none_orig)

        # Insert ",0" into the key:  ('x', 2, 3) -> ('x', 0, 2, 0, 3)
        dsk2: dict = {}
        for k, v in dsk.items():
            if k[0] == self._name:
                k2 = (self._name,) + self.expand(k[1:], 0)
                if isinstance(v.args[1], Alias):
                    # positional indexing with newaxis
                    indexer = expand_orig(dsk[v.args[1].key].value[1], None)
                    tok = "shuffle-taker-" + tokenize(indexer)
                    dsk2[tok] = DataNode(tok, (1, indexer))
                    arg = TaskRef(tok)
                else:
                    arg = expand_orig(v.args[1], None)
                # raise NotImplementedError
                dsk2[k2] = Task(k2, v.func, v.args[0], arg)
            else:
                dsk2[k] = v
        return dsk2


class TakeUnknownOneChunk(Slice):
    _parameters = ["array", "index", "axis"]

    @functools.cached_property
    def chunks(self):
        return self.array.chunks

    def _layer(self) -> dict:
        slices = [slice(None)] * len(self.array.chunks)
        slices[self.axis] = list(self.index)
        sl = tuple(slices)
        chunk_tuples = list(
            product(*(range(len(c)) for i, c in enumerate(self.array.chunks)))
        )
        dsk = {
            (self._name,)
            + ct: Task(
                (self._name,) + ct, getitem, TaskRef((self.array.name,) + ct), sl
            )
            for ct in chunk_tuples
        }
        return dsk
