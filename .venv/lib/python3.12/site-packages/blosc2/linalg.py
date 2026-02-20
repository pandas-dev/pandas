#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import builtins
import math
import warnings
from itertools import product
from typing import TYPE_CHECKING, Any

import numpy as np

import blosc2

from .utils import get_intersecting_chunks, nptranspose, npvecdot, slice_to_chunktuple

if TYPE_CHECKING:
    from collections.abc import Sequence


def matmul(x1: blosc2.Array, x2: blosc2.NDArray, **kwargs: Any) -> blosc2.NDArray:
    """
    Computes the matrix product between two Blosc2 NDArrays.

    Parameters
    ----------
    x1: :ref:`NDArray` | np.ndarray
        The first input array.
    x2: :ref:`NDArray` | np.ndarray
        The second input array.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        The matrix product of the inputs. This is a scalar only when both x1,
        x2 are 1-d vectors.

    Raises
    ------
    ValueError
        If the last dimension of ``x1`` is not the same size as
        the second-to-last dimension of ``x2``.

        If a scalar value is passed in.

    References
    ----------
    `numpy.matmul <https://numpy.org/doc/stable/reference/generated/numpy.matmul.html>`_

    Examples
    --------
    For 2-D arrays it is the matrix product:

    >>> import numpy as np
    >>> import blosc2
    >>> a = np.array([[1, 2],
    ...               [3, 4]])
    >>> nd_a = blosc2.asarray(a)
    >>> b = np.array([[2, 3],
    ...               [2, 1]])
    >>> nd_b = blosc2.asarray(b)
    >>> blosc2.matmul(nd_a, nd_b)
    array([[ 6,  5],
           [14, 13]])

    For 2-D mixed with 1-D, the result is the usual.

    >>> a = np.array([[1, 3],
    ...               [0, 1]])
    >>> nd_a = blosc2.asarray(a)
    >>> v = np.array([1, 2])
    >>> nd_v = blosc2.asarray(v)
    >>> blosc2.matmul(nd_a, nd_v)
    array([7, 2])
    >>> blosc2.matmul(nd_v, nd_a)
    array([1, 5])

    """
    # Validate arguments are not scalars
    if np.isscalar(x1) or np.isscalar(x2):
        raise ValueError("Arguments can't be scalars.")

    # Makes a SimpleProxy if inputs are not blosc2 arrays
    x1, x2 = blosc2.as_simpleproxy(x1, x2)

    # Validate matrix multiplication compatibility
    if x1.shape[builtins.max(-1, -len(x2.shape))] != x2.shape[builtins.max(-2, -len(x2.shape))]:
        raise ValueError("Shapes are not aligned for matrix multiplication.")

    # Promote 1D arrays to 2D if necessary
    x1_is_vector = False
    x2_is_vector = False
    if x1.ndim == 1:
        x1 = blosc2.expand_dims(x1, axis=0)  # (N,) -> (1, N)
        x1_is_vector = True
    if x2.ndim == 1:
        x2 = blosc2.expand_dims(x2, axis=1)  # (M,) -> (M, 1)
        x2_is_vector = True

    n, k = x1.shape[-2:]
    m = x2.shape[-1]
    result_shape = np.broadcast_shapes(x1.shape[:-2], x2.shape[:-2]) + (n, m)
    # For matmul, we don't want to reduce the chunksize, as experiments show that
    # the larger, the better (as long as some limits are not exceeded).
    kwargs["_chunksize_reduc_factor"] = 1
    result = blosc2.zeros(result_shape, dtype=blosc2.result_type(x1, x2), **kwargs)

    if 0 not in result.shape + x1.shape + x2.shape:  # if any array is empty, return array of 0s
        p, q = result.chunks[-2:]
        r = x2.chunks[-1]

        intersecting_chunks = get_intersecting_chunks((), result.shape[:-2], result.chunks[:-2])
        for chunk in intersecting_chunks:
            chunk = chunk.raw
            for row in range(0, n, p):
                row_end = builtins.min(row + p, n)
                for col in range(0, m, q):
                    col_end = builtins.min(col + q, m)
                    for aux in range(0, k, r):
                        aux_end = builtins.min(aux + r, k)
                        bx1 = (
                            x1[chunk[-x1.ndim + 2 :] + (slice(row, row_end), slice(aux, aux_end))]
                            if x1.ndim > 2
                            else x1[row:row_end, aux:aux_end]
                        )
                        bx2 = (
                            x2[chunk[-x2.ndim + 2 :] + (slice(aux, aux_end), slice(col, col_end))]
                            if x2.ndim > 2
                            else x2[aux:aux_end, col:col_end]
                        )
                        result[chunk + (slice(row, row_end), slice(col, col_end))] += np.matmul(bx1, bx2)

    if x1_is_vector:
        result = result.squeeze(axis=-2)
    if x2_is_vector:
        result = result.squeeze(axis=-1)

    return result


def tensordot(
    x1: blosc2.NDArray,
    x2: blosc2.NDArray,
    axes: int | tuple[Sequence[int], Sequence[int]] = 2,
    **kwargs: Any,
) -> blosc2.NDArray:
    """
    Returns a tensor contraction of x1 and x2 over specific axes. The tensordot function corresponds to the
    generalized matrix product. Note: Neither argument is complex-conjugated or transposed. If conjugation and/or transposition is desired, these operations should be explicitly
    performed prior to computing the generalized matrix product.

    Parameters
    ----------
    x1: blosc2.NDArray
        First input array. Should have a numeric data type.

    x2: blosc2.NDArray
        Second input array. Should have a numeric data type. Corresponding contracted axes of x1 and x2
        must be equal.

    axes: int | tuple[Sequence[int], Sequence[int]]
        Number of axes (dimensions) to contract or explicit sequences of axis (dimension) indices for x1 and x2,
        respectively.

        * If axes is an int equal to N, then contraction is performed over the last N axes of x1 and the first N axes of x2 in order. The size of each corresponding axis (dimension) must match. Must be nonnegative.

        * If N equals 0, the result is the tensor (outer) product.

        * If N equals 1, the result is the tensor dot product.

        * If N equals 2, the result is the tensor double contraction (default).

        * If axes is a tuple of two sequences (x1_axes, x2_axes), the first sequence applies to x1 and the second sequence to x2.
        Both sequences must have the same length. Each axis (dimension) x1_axes[i] for x1 must have the same size as the respective
        axis (dimension) x2_axes[i] for x2. Each index referred to in a sequence must be unique. If x1 has rank (i.e, number of dimensions) N,
        a valid x1 axis must reside on the half-open interval [-N, N). If x2 has rank M, a valid x2 axis must reside on the half-open interval [-M, M).

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: blosc2.NDArray
        An array containing the tensor contraction whose shape consists of the non-contracted axes (dimensions) of the first array x1, followed by
        the non-contracted axes (dimensions) of the second array x2.
    """
    fast_path = kwargs.pop("fast_path", None)  # for testing purposes
    # TODO: add fast path for when don't need to change chunkshapes

    # Makes a SimpleProxy if inputs are not blosc2 arrays
    x1, x2 = blosc2.as_simpleproxy(x1, x2)

    if isinstance(axes, tuple):
        a_axes, b_axes = axes
        a_axes = list(a_axes)
        b_axes = list(b_axes)
        if len(a_axes) != len(b_axes):
            raise ValueError("Lengths of reduction axes for x1 and x2 must be equal!")
        # need to track order of b_axes; later we cycle through a_axes sorted for op_chunk
        # a_sorted[inv_sort][b_sort] matches b_sorted since b_axes matches a_axes
        inv_sort = np.argsort(np.argsort(a_axes))
        b_sort = np.argsort(b_axes)
        order = inv_sort[b_sort]
        a_keep, b_keep = [True] * x1.ndim, [True] * x2.ndim
        for i, j in zip(a_axes, b_axes, strict=False):
            i = x1.ndim + i if i < 0 else i
            j = x2.ndim + j if j < 0 else j
            a_keep[i] = False
            b_keep[j] = False
        a_axes = [] if a_axes == () else a_axes  # handle no reduction
        b_axes = [] if b_axes == () else b_axes  # handle no reduction
    elif isinstance(axes, int):
        if axes < 0:
            raise ValueError("Integer axes argument must be nonnegative!")
        order = np.arange(axes, dtype=int)  # no reordering required
        a_axes = list(range(x1.ndim - axes, x1.ndim))
        b_axes = list(range(0, axes))
        a_keep = [i + axes < x1.ndim for i in range(x1.ndim)]
        b_keep = [i >= axes for i in range(x2.ndim)]
    else:
        raise ValueError("Axes argument must be two element tuple of sequences or an integer.")
    x1shape = np.array(x1.shape)
    x2shape = np.array(x2.shape)
    a_chunks_red = tuple(c for i, c in enumerate(x1.chunks) if not a_keep[i])
    a_shape_red = tuple(c for i, c in enumerate(x1.shape) if not a_keep[i])

    if np.any(x1shape[a_axes] != x2shape[b_axes]):
        raise ValueError("x1 and x2 must have same shapes along reduction dimensions")

    result_shape = tuple(x1shape[a_keep]) + tuple(x2shape[b_keep])
    result = blosc2.zeros(result_shape, dtype=blosc2.result_type(x1, x2), **kwargs)

    op_chunks = [
        slice_to_chunktuple(slice(0, s, 1), c) for s, c in zip(x1shape[a_axes], a_chunks_red, strict=True)
    ]
    res_chunks = [
        slice_to_chunktuple(s, c)
        for s, c in zip([slice(0, r, 1) for r in result.shape], result.chunks, strict=True)
    ]
    a_selection = (slice(None, None, 1),) * x1.ndim
    b_selection = (slice(None, None, 1),) * x2.ndim

    chunk_memory = np.prod(result.chunks) * (
        np.prod(x1shape[a_axes]) * x1.dtype.itemsize + np.prod(x2shape[b_axes]) * x2.dtype.itemsize
    )
    if chunk_memory < blosc2.MAX_FAST_PATH_SIZE:
        fast_path = True if fast_path is None else fast_path
    fast_path = False if fast_path is None else fast_path  # fast_path set via kwargs for testing

    # adapted from numpy.tensordot
    a_keep_axes = [i for i, k in enumerate(a_keep) if k]
    b_keep_axes = [i for i, k in enumerate(b_keep) if k]
    newaxes_a = a_keep_axes + a_axes
    newaxes_b = b_axes + b_keep_axes

    for rchunk in product(*res_chunks):
        res_chunk = tuple(
            slice(rc * rcs, builtins.min((rc + 1) * rcs, rshape), 1)
            for rc, rcs, rshape in zip(rchunk, result.chunks, result.shape, strict=True)
        )
        rchunk_iter = iter(res_chunk)
        a_selection = tuple(next(rchunk_iter) if a else slice(None, None, 1) for a in a_keep)
        b_selection = tuple(next(rchunk_iter) if b else slice(None, None, 1) for b in b_keep)
        res_chunks = tuple(s.stop - s.start for s in res_chunk)
        for ochunk in product(*op_chunks):
            if not fast_path:  # operands too big, have to go chunk-by-chunk
                op_chunk = tuple(
                    slice(rc * rcs, builtins.min((rc + 1) * rcs, x1s), 1)
                    for rc, rcs, x1s in zip(ochunk, a_chunks_red, a_shape_red, strict=True)
                )  # use x1 chunk shape to iterate over reduction axes
                ochunk_iter = iter(op_chunk)
                a_selection = tuple(
                    next(ochunk_iter) if not a else as_ for as_, a in zip(a_selection, a_keep, strict=True)
                )
                # have to permute to match order of a_axes
                order_iter = iter(order)
                b_selection = tuple(
                    op_chunk[next(order_iter)] if not b else bs_
                    for bs_, b in zip(b_selection, b_keep, strict=True)
                )
            bx1 = x1[a_selection]
            bx2 = x2[b_selection]
            # adapted from numpy tensordot
            newshape_a = (
                math.prod([bx1.shape[i] for i in a_keep_axes]),
                math.prod([bx1.shape[a] for a in a_axes]),
            )
            newshape_b = (
                math.prod([bx2.shape[b] for b in b_axes]),
                math.prod([bx2.shape[i] for i in b_keep_axes]),
            )
            at = nptranspose(bx1, newaxes_a).reshape(newshape_a)
            bt = nptranspose(bx2, newaxes_b).reshape(newshape_b)
            res = np.dot(at, bt)
            result[res_chunk] += res.reshape(res_chunks)
            if fast_path:  # already done everything
                break
    return result


def vecdot(x1: blosc2.NDArray, x2: blosc2.NDArray, axis: int = -1, **kwargs) -> blosc2.NDArray:
    """
    Computes the (vector) dot product of two arrays. Complex conjugates x1.

    Parameters
    ----------
    x1: blosc2.NDArray
        First input array. Must have floating-point data type.

    x2: blosc2.NDArray
        Second input array. Must be compatible with x1 for all non-contracted axes (via broadcasting).
        The size of the axis over which to compute the dot product must be the same size as the respective axis in x1.
        Must have a floating-point data type.

    axis: int
        The axis (dimension) of x1 and x2 containing the vectors for which to compute the dot product.
        Should be an integer on the interval [-N, -1], where N is min(x1.ndim, x2.ndim). Default: -1.

    Returns
    -------
    out: blosc2.NDArray
        If x1 and x2 are both one-dimensional arrays, a zero-dimensional containing the dot product;
        otherwise, a non-zero-dimensional array containing the dot products and having rank N-1,
        where N is the rank (number of dimensions) of the shape determined according to broadcasting
        along the non-contracted axes.
    """
    fast_path = kwargs.pop("fast_path", None)  # for testing purposes
    # Added this to pass array-api tests (which use internal getitem to check results)
    if isinstance(x1, np.ndarray) and isinstance(x2, np.ndarray):
        return npvecdot(x1, x2, axis=axis)

    # Makes a SimpleProxy if inputs are not blosc2 arrays
    x1, x2 = blosc2.as_simpleproxy(x1, x2)

    N = builtins.min(x1.ndim, x2.ndim)
    if axis < -N or axis > -1:
        raise ValueError("axis must be on interval [-N,-1].")
    a_axes = axis + x1.ndim
    b_axes = axis + x2.ndim
    a_keep = [True] * x1.ndim
    a_keep[a_axes] = False
    b_keep = [True] * x2.ndim
    b_keep[b_axes] = False

    x1shape = np.array(x1.shape)
    x2shape = np.array(x2.shape)
    a_chunks_red = x1.chunks[a_axes]
    a_shape_red = x1.shape[a_axes]

    if np.any(x1shape[a_axes] != x2shape[b_axes]):
        raise ValueError("x1 and x2 must have same shapes along reduction dimensions")

    result_shape = np.broadcast_shapes(x1shape[a_keep], x2shape[b_keep])
    result = blosc2.zeros(result_shape, dtype=blosc2.result_type(x1, x2), **kwargs)

    res_chunks = [
        slice_to_chunktuple(s, c)
        for s, c in zip([slice(0, r, 1) for r in result.shape], result.chunks, strict=True)
    ]
    a_selection = (slice(None, None, 1),) * x1.ndim
    b_selection = (slice(None, None, 1),) * x2.ndim

    chunk_memory = np.prod(result.chunks) * (
        x1shape[a_axes] * x1.dtype.itemsize + x2shape[b_axes] * x2.dtype.itemsize
    )
    if chunk_memory < blosc2.MAX_FAST_PATH_SIZE:
        fast_path = True if fast_path is None else fast_path
    fast_path = False if fast_path is None else fast_path  # fast_path set via kwargs for testing

    for rchunk in product(*res_chunks):
        res_chunk = tuple(
            slice(rc * rcs, builtins.min((rc + 1) * rcs, rshape), 1)
            for rc, rcs, rshape in zip(rchunk, result.chunks, result.shape, strict=True)
        )
        # handle broadcasting - if x1, x2 different ndim, could have to prepend 1s
        rchunk_iter = (
            slice(0, 1, 1) if s == 1 else r
            for r, s in zip(res_chunk[-x1.ndim + 1 :], x1shape[a_keep], strict=True)
        )
        a_selection = tuple(next(rchunk_iter) if a else slice(None, None, 1) for a in a_keep)
        rchunk_iter = (
            slice(0, 1, 1) if s == 1 else r
            for r, s in zip(res_chunk[-x2.ndim + 1 :], x2shape[b_keep], strict=True)
        )
        b_selection = tuple(next(rchunk_iter) if b else slice(None, None, 1) for b in b_keep)

        for ochunk in range(0, a_shape_red, a_chunks_red):
            if not fast_path:  # operands too big, go chunk-by-chunk
                op_chunk = (slice(ochunk, builtins.min(ochunk + a_chunks_red, x1.shape[a_axes]), 1),)
                a_selection = a_selection[:a_axes] + op_chunk + a_selection[a_axes + 1 :]
                b_selection = b_selection[:b_axes] + op_chunk + b_selection[b_axes + 1 :]
            bx1 = x1[a_selection]
            bx2 = x2[b_selection]
            res = npvecdot(bx1, bx2, axis=axis)  # handles conjugation of bx1
            result[res_chunk] += res
            if fast_path:  # already done everything
                break
    return result


def permute_dims(
    arr: blosc2.Array, axes: tuple[int] | list[int] | None = None, **kwargs: Any
) -> blosc2.NDArray:
    """
    Permutes the axes (dimensions) of an array.

    Parameters
    ----------
    arr: :ref:`blosc2.NDArray` | np.ndarray
        The input array.
    axes: tuple[int], list[int], optional
        The desired permutation of axes. If None, the axes are reversed by default.
        If specified, axes must be a tuple or list representing a permutation of
        ``[0, 1, ..., N-1]``, where ``N`` is the number of dimensions of the input array.
        Negative indices are also supported. The *i*-th axis of the result will correspond
        to the axis numbered ``axes[i]`` of the input.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`blosc2.NDArray`
        A Blosc2 :ref:`blosc2.NDArray` with axes transposed.

    Raises
    ------
    ValueError
        If ``axes`` is not a valid permutation of the dimensions of ``arr``.

    References
    ----------
    `numpy.transpose <https://numpy.org/doc/2.2/reference/generated/numpy.transpose.html>`_

    `permute_dims <https://data-apis.org/array-api/latest/API_specification/generated/array_api.permute_dims.html#permute-dims>`_

    Examples
    --------
    For 2-D arrays it is the matrix transposition as usual:

    >>> import blosc2
    >>> a = blosc2.arange(1, 10).reshape((3, 3))
    >>> a[:]
    array([[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]])
    >>> at = blosc2.permute_dims(a)
    >>> at[:]
    array([[1, 4, 7],
           [2, 5, 8],
           [3, 6, 9]])

    For 3-D arrays:

    >>> import blosc2
    >>> a = blosc2.arange(1, 25).reshape((2, 3, 4))
    >>> a[:]
    array([[[ 1,  2,  3,  4],
            [ 5,  6,  7,  8],
            [ 9, 10, 11, 12]],
           [[13, 14, 15, 16],
            [17, 18, 19, 20],
            [21, 22, 23, 24]]])

    >>> at = blosc2.permute_dims(a, axes=(1, 0, 2))
    >>> at[:]
    array([[[ 1,  2,  3,  4],
            [13, 14, 15, 16]],
           [[ 5,  6,  7,  8],
            [17, 18, 19, 20]],
           [[ 9, 10, 11, 12],
            [21, 22, 23, 24]]])
    """
    if np.isscalar(arr) or arr.ndim < 2:
        return arr

    # Makes a SimpleProxy if input is not blosc2 array
    arr = blosc2.as_simpleproxy(arr)

    ndim = arr.ndim

    if axes is None:
        axes = tuple(range(ndim))[::-1]
    else:
        axes = tuple(axis if axis >= 0 else ndim + axis for axis in axes)
        if sorted(axes) != list(range(ndim)):
            raise ValueError(f"axes {axes} is not a valid permutation of {ndim} dimensions")

    new_shape = tuple(arr.shape[axis] for axis in axes)
    if "chunks" not in kwargs or kwargs["chunks"] is None:
        kwargs["chunks"] = tuple(arr.chunks[axis] for axis in axes)

    result = blosc2.empty(shape=new_shape, dtype=arr.dtype, **kwargs)

    chunks = arr.chunks
    shape = arr.shape
    # handle SimpleProxy which doesn't have iterchunks_info
    if hasattr(arr, "iterchunks_info"):
        my_it = arr.iterchunks_info()
        _get_el = lambda x: x.coords  # noqa: E731
    else:
        my_it = get_intersecting_chunks((), shape, chunks)
        _get_el = lambda x: x.raw  # noqa: E731
    for info in my_it:
        coords = _get_el(info)
        start_stop = [
            (coord * chunk, builtins.min(chunk * (coord + 1), dim))
            for coord, chunk, dim in zip(coords, chunks, shape, strict=False)
        ]

        src_slice = tuple(slice(start, stop) for start, stop in start_stop)
        dst_slice = tuple(slice(start_stop[ax][0], start_stop[ax][1]) for ax in axes)

        transposed = nptranspose(arr[src_slice], axes=axes)
        result[dst_slice] = np.ascontiguousarray(transposed)

    return result


def transpose(x, **kwargs: Any) -> blosc2.NDArray:
    """
    Returns a Blosc2 blosc2.NDArray with axes transposed.

    Only 2D arrays are supported for now.  Other dimensions raise an error.

    Parameters
    ----------
    x: :ref:`blosc2.NDArray`
        The input array.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`blosc2.NDArray`
        The Blosc2 blosc2.NDArray with axes transposed.

    References
    ----------
    `numpy.transpose <https://numpy.org/doc/2.2/reference/generated/numpy.transpose.html>`_
    """
    warnings.warn(
        "transpose is deprecated and will be removed in a future version. "
        "Use matrix_transpose or permute_dims instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    # If arguments are dimension < 2, they are returned
    if np.isscalar(x) or x.ndim < 2:
        return x
    # Makes a SimpleProxy if input is not blosc2 array
    x = blosc2.as_simpleproxy(x)
    # Validate arguments are dimension 2
    if x.ndim > 2:
        raise ValueError("Transposing arrays with dimension greater than 2 is not supported yet.")
    return permute_dims(x, **kwargs)


def matrix_transpose(arr: blosc2.Array, **kwargs: Any) -> blosc2.NDArray:
    """
    Transposes a matrix (or a stack of matrices).

    Parameters
    ----------
    arr: :ref:`blosc2.NDArray` | np.ndarray
        The input blosc2.NDArray having shape ``(..., M, N)`` and whose innermost two dimensions form
        ``MxN`` matrices.

    Returns
    -------
    out: :ref:`blosc2.NDArray`
        A new :ref:`blosc2.NDArray` containing the transpose for each matrix and having shape
        ``(..., N, M)``.
    """
    axes = None
    # Makes a SimpleProxy if input is not blosc2 array
    arr = blosc2.as_simpleproxy(arr)
    if not np.isscalar(arr) and arr.ndim > 2:
        axes = list(range(arr.ndim))
        axes[-2], axes[-1] = axes[-1], axes[-2]
    return permute_dims(arr, axes, **kwargs)


def diagonal(x: blosc2.blosc2.NDArray, offset: int = 0) -> blosc2.blosc2.NDArray:
    """
    Returns the specified diagonals of a matrix (or a stack of matrices) x.

    Parameters
    ----------
    x: blosc2.NDArray
        Input array having shape (..., M, N) and whose innermost two dimensions form MxN matrices.

    offset: int
        Offset specifying the off-diagonal relative to the main diagonal.

        * offset = 0: the main diagonal.
        * offset > 0: off-diagonal above the main diagonal.
        * offset < 0: off-diagonal below the main diagonal.

        Default: 0.

    Returns
    -------
    out: blosc2.NDArray
        An array containing the diagonals and whose shape is determined by
        removing the last two dimensions and appending a dimension equal to the size of the
        resulting diagonals.

    Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.diag.html#diag
    """
    # Makes a SimpleProxy if input is not blosc2 array
    x = blosc2.as_simpleproxy(x)
    n_rows, n_cols = x.shape[-2:]
    min_idx = builtins.min(n_rows, n_cols)
    if offset < 0:
        start = -offset
        rows = np.arange(start, builtins.min(start + n_cols, n_rows))
        cols = np.arange(len(rows))
    elif offset > 0:
        cols = np.arange(offset, builtins.min(offset + n_rows, n_cols))
        rows = np.arange(len(cols))
    else:
        rows = cols = np.arange(min_idx)
    key = tuple(slice(None, None, 1) for i in range(x.ndim - 2)) + (rows, cols)
    # TODO: change to use slice to give optimised compressing
    return blosc2.asarray(x[key])


def outer(x1: blosc2.blosc2.NDArray, x2: blosc2.blosc2.NDArray, **kwargs: Any) -> blosc2.blosc2.NDArray:
    """
    Returns the outer product of two vectors x1 and x2.

    Parameters
    ----------
    x1: blosc2.NDArray
        First one-dimensional input array of size N. Must have a numeric data type.

    x2: blosc2.NDArray
        Second one-dimensional input array of size M. Must have a numeric data type.

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: blosc2.NDArray
        A two-dimensional array containing the outer product and whose shape is (N, M).
    """
    x1, x2 = blosc2.as_simpleproxy(x1, x2)
    if (x1.ndim != 1) or (x2.ndim != 1):
        raise ValueError("outer only valid for 1D inputs.")
    return tensordot(x1, x2, ((), ()), **kwargs)  # for testing purposes


def cholesky(x: blosc2.blosc2.NDArray, upper: bool = False) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.cholesky.html#cholesky
    # """
    raise NotImplementedError


def cross(x1: blosc2.blosc2.NDArray, x2: blosc2.blosc2.NDArray, axis: int = -1) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.cross.html#cross
    # """
    raise NotImplementedError


def det(x: blosc2.blosc2.NDArray) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.det.html#det
    # """
    raise NotImplementedError


def eigh(x: blosc2.blosc2.NDArray) -> tuple[blosc2.blosc2.NDArray, blosc2.blosc2.NDArray]:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.eigh.html#eigh
    # """
    raise NotImplementedError


def eigvalsh(x: blosc2.blosc2.NDArray) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.eigvalsh.html#eigvalsh
    # """
    raise NotImplementedError


def inv(x: blosc2.blosc2.NDArray) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.inv.html#inv
    # """
    raise NotImplementedError


def matrix_norm(
    x: blosc2.blosc2.NDArray, keepdims: bool = False, ord: int | float | str | None = "fro"
) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented but could be doable. ord may take values:
    #     * 'fro' - Frobenius norm
    #     * 'nuc' - nuclear norm
    #     * 1 - max(sum(abs(x), axis=-2))
    #     * 2 - largest singular value (sum(x**2, axis=[-1,-2]))
    #     * inf - max(sum(abs(x), axis=-1))
    #     * -1 - min(sum(abs(x), axis=-2))
    #     * -2 - smallest singular value
    #     * -inf - min(sum(abs(x), axis=-1))
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.matrix_norm.html#matrix_norm
    # """
    raise NotImplementedError


def matrix_power(x: blosc2.blosc2.NDArray, n: int) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.matrix_power.html#matrix_power
    # """
    raise NotImplementedError


def matrix_rank(
    x: blosc2.blosc2.NDArray, rtol: float | blosc2.blosc2.NDArray | None = None
) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.matrix_rank.html#matrix_rank
    # """
    raise NotImplementedError


def pinv(
    x: blosc2.blosc2.NDArray, rtol: float | blosc2.blosc2.NDArray | None = None
) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.pinv.html#pinv
    # """
    raise NotImplementedError


def qr(
    x: blosc2.blosc2.NDArray, mode: str = "reduced"
) -> tuple[blosc2.blosc2.NDArray, blosc2.blosc2.NDArray]:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.qr.html#qr
    # """
    raise NotImplementedError


def slogdet(x: blosc2.blosc2.NDArray) -> tuple[blosc2.blosc2.NDArray, blosc2.blosc2.NDArray]:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.slogdet.html#slogdet
    # """
    raise NotImplementedError


def solve(x1: blosc2.blosc2.NDArray, x2: blosc2.blosc2.NDArray) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.solve.html#solve
    # """
    raise NotImplementedError


def svd(
    x: blosc2.blosc2.NDArray, full_matrices: bool = True
) -> tuple[blosc2.blosc2.NDArray, blosc2.blosc2.NDArray, blosc2.blosc2.NDArray]:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.svd.html#svd
    # """
    raise NotImplementedError


def svdvals(x: blosc2.blosc2.NDArray) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.svdvals.html#svdvals
    # """
    raise NotImplementedError


def trace(x: blosc2.blosc2.NDArray, offset: int = 0, dtype: np.dtype | None = None) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented
    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.trace.html#trace
    # """
    raise NotImplementedError


def vector_norm(
    x: blosc2.blosc2.NDArray,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    ord: int | float = 2,
) -> blosc2.blosc2.NDArray:
    # """
    # Not Implemented but could be doable. ord may take values:
    #     * p: int - p-norm
    #     * inf - max(x)
    #     * -inf - min(abs(x))

    # Reference: https://data-apis.org/array-api/latest/extensions/generated/array_api.linalg.vector_norm.html#vector_norm
    # """
    raise NotImplementedError
