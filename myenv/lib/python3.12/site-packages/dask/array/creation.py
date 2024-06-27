from __future__ import annotations

import itertools
from collections.abc import Sequence
from functools import partial
from itertools import product
from numbers import Integral, Number

import numpy as np
from tlz import sliding_window

from dask.array import chunk
from dask.array.backends import array_creation_dispatch
from dask.array.core import (
    Array,
    asarray,
    block,
    blockwise,
    broadcast_arrays,
    broadcast_to,
    concatenate,
    normalize_chunks,
    stack,
)
from dask.array.numpy_compat import NUMPY_GE_200, AxisError
from dask.array.ufunc import greater_equal, rint
from dask.array.utils import meta_from_array
from dask.array.wrap import empty, full, ones, zeros
from dask.base import tokenize
from dask.highlevelgraph import HighLevelGraph
from dask.utils import cached_cumsum, derived_from, is_cupy_type


def to_backend(x: Array, backend: str | None = None, **kwargs):
    """Move an Array collection to a new backend

    Parameters
    ----------
    x : Array
        The input Array collection.
    backend : str, Optional
        The name of the new backend to move to. The default
        is the current "array.backend" configuration.

    Returns
    -------
    dask.Array
        A new Array collection with the backend specified
        by ``backend``.
    """
    # Get desired backend
    backend = backend or array_creation_dispatch.backend
    # Check that "backend" has a registered entrypoint
    backend_entrypoint = array_creation_dispatch.dispatch(backend)
    # Call `ArrayBackendEntrypoint.to_backend`
    return backend_entrypoint.to_backend(x, **kwargs)


def empty_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return a new array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of the
        returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of uninitialized (arbitrary) data with the same
        shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    zeros_like : Return an array of zeros with shape and type of input.
    empty : Return a new uninitialized array.
    ones : Return a new array setting values to one.
    zeros : Return a new array setting values to zero.

    Notes
    -----
    This function does *not* initialize the returned array; to do that use
    `zeros_like` or `ones_like` instead.  It may be marginally faster than
    the functions that do set the array values.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular empty function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.empty_like, dtype=(dtype or a.dtype)))

    return empty(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def ones_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return an array of ones with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of ones with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular ones function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.ones_like, dtype=(dtype or a.dtype)))

    return ones(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def zeros_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return an array of zeros with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of zeros with the same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular zeros function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.zeros_like, dtype=(dtype or a.dtype)))

    return zeros(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def full_like(a, fill_value, order="C", dtype=None, chunks=None, name=None, shape=None):
    """
    Return a full array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    fill_value : scalar
        Fill value.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of `fill_value` with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    full : Fill a new array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular full function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.full_like, dtype=(dtype or a.dtype)), fill_value)

    return full(
        shape,
        fill_value,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def _get_like_function_shapes_chunks(a, chunks, shape):
    """
    Helper function for finding shapes and chunks for *_like()
    array creation functions.
    """
    if shape is None:
        shape = a.shape
        if chunks is None:
            chunks = a.chunks
    elif chunks is None:
        chunks = "auto"
    return shape, chunks


def linspace(
    start, stop, num=50, endpoint=True, retstep=False, chunks="auto", dtype=None
):
    """
    Return `num` evenly spaced values over the closed interval [`start`,
    `stop`].

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The last value of the sequence.
    num : int, optional
        Number of samples to include in the returned dask array, including the
        endpoints. Default is 50.
    endpoint : bool, optional
        If True, ``stop`` is the last sample. Otherwise, it is not included.
        Default is True.
    retstep : bool, optional
        If True, return (samples, step), where step is the spacing between
        samples. Default is False.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if `num % blocksize != 0`
    dtype : dtype, optional
        The type of the output array.

    Returns
    -------
    samples : dask array
    step : float, optional
        Only returned if ``retstep`` is True. Size of spacing between samples.


    See Also
    --------
    dask.array.arange
    """
    num = int(num)

    if dtype is None:
        dtype = np.linspace(0, 1, 1).dtype

    chunks = normalize_chunks(chunks, (num,), dtype=dtype)

    range_ = stop - start

    div = (num - 1) if endpoint else num
    if div == 0:
        div = 1

    step = float(range_) / div

    name = "linspace-" + tokenize((start, stop, num, endpoint, chunks, dtype))

    dsk = {}
    blockstart = start

    for i, bs in enumerate(chunks[0]):
        bs_space = bs - 1 if endpoint else bs
        blockstop = blockstart + (bs_space * step)
        task = (
            partial(chunk.linspace, endpoint=endpoint, dtype=dtype),
            blockstart,
            blockstop,
            bs,
        )
        blockstart = blockstart + (step * bs)
        dsk[(name, i)] = task

    if retstep:
        return Array(dsk, name, chunks, dtype=dtype), step
    else:
        return Array(dsk, name, chunks, dtype=dtype)


@array_creation_dispatch.register_inplace("numpy")
def arange(*args, chunks="auto", like=None, dtype=None, **kwargs):
    """
    Return evenly spaced values from `start` to `stop` with step size `step`.

    The values are half-open [start, stop), so including start and excluding
    stop. This is basically the same as python's range function but for dask
    arrays.

    When using a non-integer step, such as 0.1, the results will often not be
    consistent. It is better to use linspace for these cases.

    Parameters
    ----------
    start : int, optional
        The starting value of the sequence. The default is 0.
    stop : int
        The end of the interval, this value is excluded from the interval.
    step : int, optional
        The spacing between the values. The default is 1 when not specified.
        The last value of the sequence.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
        Defaults to "auto" which will automatically determine chunk sizes.
    dtype : numpy.dtype
        Output dtype. Omit to infer it from start, stop, step
        Defaults to ``None``.
    like : array type or ``None``
        Array to extract meta from. Defaults to ``None``.

    Returns
    -------
    samples : dask array

    See Also
    --------
    dask.array.linspace
    """
    if len(args) == 1:
        start = 0
        stop = args[0]
        step = 1
    elif len(args) == 2:
        start = args[0]
        stop = args[1]
        step = 1
    elif len(args) == 3:
        start, stop, step = args
    else:
        raise TypeError(
            """
        arange takes 3 positional arguments: arange([start], stop, [step])
        """
        )

    num = int(max(np.ceil((stop - start) / step), 0))

    meta = meta_from_array(like) if like is not None else None

    if dtype is None:
        dtype = np.arange(start, stop, step * num if num else step).dtype

    chunks = normalize_chunks(chunks, (num,), dtype=dtype)

    if kwargs:
        raise TypeError("Unexpected keyword argument(s): %s" % ",".join(kwargs.keys()))

    name = "arange-" + tokenize((start, stop, step, chunks, dtype))
    dsk = {}
    elem_count = 0

    for i, bs in enumerate(chunks[0]):
        blockstart = start + (elem_count * step)
        blockstop = start + ((elem_count + bs) * step)
        task = (
            partial(chunk.arange, like=like),
            blockstart,
            blockstop,
            step,
            bs,
            dtype,
        )
        dsk[(name, i)] = task
        elem_count += bs

    return Array(dsk, name, chunks, dtype=dtype, meta=meta)


@derived_from(np)
def meshgrid(*xi, sparse=False, indexing="xy", **kwargs):
    sparse = bool(sparse)

    if "copy" in kwargs:
        raise NotImplementedError("`copy` not supported")

    if kwargs:
        raise TypeError("unsupported keyword argument(s) provided")

    if indexing not in ("ij", "xy"):
        raise ValueError("`indexing` must be `'ij'` or `'xy'`")

    xi = [asarray(e) for e in xi]
    xi = [e.flatten() for e in xi]

    if indexing == "xy" and len(xi) > 1:
        xi[0], xi[1] = xi[1], xi[0]

    grid = []
    for i in range(len(xi)):
        s = len(xi) * [None]
        s[i] = slice(None)
        s = tuple(s)

        r = xi[i][s]

        grid.append(r)

    if not sparse:
        grid = broadcast_arrays(*grid)

    if indexing == "xy" and len(xi) > 1:
        grid = (grid[1], grid[0], *grid[2:])

    out_type = tuple if NUMPY_GE_200 else list
    return out_type(grid)


def indices(dimensions, dtype=int, chunks="auto"):
    """
    Implements NumPy's ``indices`` for Dask Arrays.

    Generates a grid of indices covering the dimensions provided.

    The final array has the shape ``(len(dimensions), *dimensions)``. The
    chunks are used to specify the chunking for axis 1 up to
    ``len(dimensions)``. The 0th axis always has chunks of length 1.

    Parameters
    ----------
    dimensions : sequence of ints
        The shape of the index grid.
    dtype : dtype, optional
        Type to use for the array. Default is ``int``.
    chunks : sequence of ints, str
        The size of each block.  Must be one of the following forms:

        - A blocksize like (500, 1000)
        - A size in bytes, like "100 MiB" which will choose a uniform
          block-like shape
        - The word "auto" which acts like the above, but uses a configuration
          value ``array.chunk-size`` for the chunk size

        Note that the last block will have fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    grid : dask array
    """
    dimensions = tuple(dimensions)
    dtype = np.dtype(dtype)
    chunks = normalize_chunks(chunks, shape=dimensions, dtype=dtype)

    if len(dimensions) != len(chunks):
        raise ValueError("Need same number of chunks as dimensions.")

    xi = []
    for i in range(len(dimensions)):
        xi.append(arange(dimensions[i], dtype=dtype, chunks=(chunks[i],)))

    grid = []
    if all(dimensions):
        grid = meshgrid(*xi, indexing="ij")

    if grid:
        grid = stack(grid)
    else:
        grid = empty((len(dimensions),) + dimensions, dtype=dtype, chunks=(1,) + chunks)

    return grid


def eye(N, chunks="auto", M=None, k=0, dtype=float):
    """
    Return a 2-D Array with ones on the diagonal and zeros elsewhere.

    Parameters
    ----------
    N : int
      Number of rows in the output.
    chunks : int, str
        How to chunk the array. Must be one of the following forms:

        -   A blocksize like 1000.
        -   A size in bytes, like "100 MiB" which will choose a uniform
            block-like shape
        -   The word "auto" which acts like the above, but uses a configuration
            value ``array.chunk-size`` for the chunk size
    M : int, optional
      Number of columns in the output. If None, defaults to `N`.
    k : int, optional
      Index of the diagonal: 0 (the default) refers to the main diagonal,
      a positive value refers to an upper diagonal, and a negative value
      to a lower diagonal.
    dtype : data-type, optional
      Data-type of the returned array.

    Returns
    -------
    I : Array of shape (N,M)
      An array where all elements are equal to zero, except for the `k`-th
      diagonal, whose values are equal to one.
    """
    eye = {}
    if M is None:
        M = N
    if dtype is None:
        dtype = float

    if not isinstance(chunks, (int, str)):
        raise ValueError("chunks must be an int or string")

    vchunks, hchunks = normalize_chunks(chunks, shape=(N, M), dtype=dtype)
    chunks = vchunks[0]

    token = tokenize(N, chunks, M, k, dtype)
    name_eye = "eye-" + token

    for i, vchunk in enumerate(vchunks):
        for j, hchunk in enumerate(hchunks):
            if (j - i - 1) * chunks <= k <= (j - i + 1) * chunks:
                eye[name_eye, i, j] = (
                    np.eye,
                    vchunk,
                    hchunk,
                    k - (j - i) * chunks,
                    dtype,
                )
            else:
                eye[name_eye, i, j] = (np.zeros, (vchunk, hchunk), dtype)
    return Array(eye, name_eye, shape=(N, M), chunks=(chunks, chunks), dtype=dtype)


@derived_from(np)
def diag(v, k=0):
    if not isinstance(v, np.ndarray) and not isinstance(v, Array):
        raise TypeError(f"v must be a dask array or numpy array, got {type(v)}")

    name = "diag-" + tokenize(v, k)

    meta = meta_from_array(v, 2 if v.ndim == 1 else 1)

    if isinstance(v, np.ndarray) or (
        hasattr(v, "__array_function__") and not isinstance(v, Array)
    ):
        if v.ndim == 1:
            m = abs(k)
            chunks = ((v.shape[0] + m,), (v.shape[0] + m,))
            dsk = {(name, 0, 0): (np.diag, v, k)}
        elif v.ndim == 2:
            kdiag_row_start = max(0, -k)
            kdiag_row_stop = min(v.shape[0], v.shape[1] - k)
            len_kdiag = kdiag_row_stop - kdiag_row_start
            chunks = ((0,),) if len_kdiag <= 0 else ((len_kdiag,),)
            dsk = {(name, 0): (np.diag, v, k)}
        else:
            raise ValueError("Array must be 1d or 2d only")
        return Array(dsk, name, chunks, meta=meta)

    if v.ndim != 1:
        if v.ndim != 2:
            raise ValueError("Array must be 1d or 2d only")
        if k == 0 and v.chunks[0] == v.chunks[1]:
            dsk = {
                (name, i): (np.diag, row[i]) for i, row in enumerate(v.__dask_keys__())
            }
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[v])
            return Array(graph, name, (v.chunks[0],), meta=meta)
        else:
            return diagonal(v, k)

    if k == 0:
        chunks_1d = v.chunks[0]
        blocks = v.__dask_keys__()
        dsk = {}
        for i, m in enumerate(chunks_1d):
            for j, n in enumerate(chunks_1d):
                key = (name, i, j)
                if i == j:
                    dsk[key] = (np.diag, blocks[i])
                else:
                    dsk[key] = (np.zeros, (m, n))
                    dsk[key] = (partial(np.zeros_like, shape=(m, n)), meta)

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[v])
        return Array(graph, name, (chunks_1d, chunks_1d), meta=meta)

    elif k > 0:
        return pad(diag(v), [[0, k], [k, 0]], mode="constant")
    elif k < 0:
        return pad(diag(v), [[-k, 0], [0, -k]], mode="constant")


@derived_from(np)
def diagonal(a, offset=0, axis1=0, axis2=1):
    name = "diagonal-" + tokenize(a, offset, axis1, axis2)

    if a.ndim < 2:
        # NumPy uses `diag` as we do here.
        raise ValueError("diag requires an array of at least two dimensions")

    def _axis_fmt(axis, name, ndim):
        if axis < 0:
            t = ndim + axis
            if t < 0:
                msg = "{}: axis {} is out of bounds for array of dimension {}"
                raise AxisError(msg.format(name, axis, ndim))
            axis = t
        return axis

    def pop_axes(chunks, axis1, axis2):
        chunks = list(chunks)
        chunks.pop(axis2)
        chunks.pop(axis1)
        return tuple(chunks)

    axis1 = _axis_fmt(axis1, "axis1", a.ndim)
    axis2 = _axis_fmt(axis2, "axis2", a.ndim)

    if axis1 == axis2:
        raise ValueError("axis1 and axis2 cannot be the same")

    a = asarray(a)
    k = offset
    if axis1 > axis2:
        axis1, axis2 = axis2, axis1
        k = -offset

    free_axes = set(range(a.ndim)) - {axis1, axis2}
    free_indices = list(product(*(range(a.numblocks[i]) for i in free_axes)))
    ndims_free = len(free_axes)

    # equation of diagonal: i = j - k
    kdiag_row_start = max(0, -k)
    kdiag_col_start = max(0, k)
    kdiag_row_stop = min(a.shape[axis1], a.shape[axis2] - k)
    len_kdiag = kdiag_row_stop - kdiag_row_start

    if len_kdiag <= 0:
        xp = np

        if is_cupy_type(a._meta):
            import cupy

            xp = cupy

        out_chunks = pop_axes(a.chunks, axis1, axis2) + ((0,),)
        dsk = dict()
        for free_idx in free_indices:
            shape = tuple(
                out_chunks[axis][free_idx[axis]] for axis in range(ndims_free)
            )
            dsk[(name,) + free_idx + (0,)] = (
                partial(xp.empty, dtype=a.dtype),
                shape + (0,),
            )

        meta = meta_from_array(a, ndims_free + 1)
        return Array(dsk, name, out_chunks, meta=meta)

    # compute row index ranges for chunks along axis1:
    row_stops_ = np.cumsum(a.chunks[axis1])
    row_starts = np.roll(row_stops_, 1)
    row_starts[0] = 0

    # compute column index ranges for chunks along axis2:
    col_stops_ = np.cumsum(a.chunks[axis2])
    col_starts = np.roll(col_stops_, 1)
    col_starts[0] = 0

    # locate first chunk containing diagonal:
    row_blockid = np.arange(a.numblocks[axis1])
    col_blockid = np.arange(a.numblocks[axis2])

    row_filter = (row_starts <= kdiag_row_start) & (kdiag_row_start < row_stops_)
    col_filter = (col_starts <= kdiag_col_start) & (kdiag_col_start < col_stops_)
    (I,) = row_blockid[row_filter]
    (J,) = col_blockid[col_filter]

    # follow k-diagonal through chunks while constructing dask graph:
    dsk = dict()
    i = 0
    kdiag_chunks = ()
    while kdiag_row_start < a.shape[axis1] and kdiag_col_start < a.shape[axis2]:
        # localize block info:
        nrows, ncols = a.chunks[axis1][I], a.chunks[axis2][J]
        kdiag_row_start -= row_starts[I]
        kdiag_col_start -= col_starts[J]
        k = -kdiag_row_start if kdiag_row_start > 0 else kdiag_col_start
        kdiag_row_end = min(nrows, ncols - k)
        kdiag_len = kdiag_row_end - kdiag_row_start

        # increment dask graph:
        for free_idx in free_indices:
            input_idx = (
                free_idx[:axis1]
                + (I,)
                + free_idx[axis1 : axis2 - 1]
                + (J,)
                + free_idx[axis2 - 1 :]
            )
            output_idx = free_idx + (i,)
            dsk[(name,) + output_idx] = (
                np.diagonal,
                (a.name,) + input_idx,
                k,
                axis1,
                axis2,
            )

        kdiag_chunks += (kdiag_len,)
        # prepare for next iteration:
        i += 1
        kdiag_row_start = kdiag_row_end + row_starts[I]
        kdiag_col_start = min(ncols, nrows + k) + col_starts[J]
        I = I + 1 if kdiag_row_start == row_stops_[I] else I
        J = J + 1 if kdiag_col_start == col_stops_[J] else J

    out_chunks = pop_axes(a.chunks, axis1, axis2) + (kdiag_chunks,)
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[a])
    meta = meta_from_array(a, ndims_free + 1)
    return Array(graph, name, out_chunks, meta=meta)


@derived_from(np)
def tri(N, M=None, k=0, dtype=float, chunks="auto", *, like=None):
    if M is None:
        M = N

    chunks = normalize_chunks(chunks, shape=(N, M), dtype=dtype)

    m = greater_equal(
        arange(N, chunks=chunks[0][0], like=like).reshape(1, N).T,
        arange(-k, M - k, chunks=chunks[1][0], like=like),
    )

    # Avoid making a copy if the requested type is already bool
    m = m.astype(dtype, copy=False)

    return m


@derived_from(np)
def fromfunction(func, chunks="auto", shape=None, dtype=None, **kwargs):
    dtype = dtype or float
    chunks = normalize_chunks(chunks, shape, dtype=dtype)

    inds = tuple(range(len(shape)))

    arrs = [arange(s, dtype=dtype, chunks=c) for s, c in zip(shape, chunks)]
    arrs = meshgrid(*arrs, indexing="ij")

    args = sum(zip(arrs, itertools.repeat(inds)), ())

    res = blockwise(func, inds, *args, token="fromfunction", **kwargs)

    return res


@derived_from(np)
def repeat(a, repeats, axis=None):
    if axis is None:
        if a.ndim == 1:
            axis = 0
        else:
            raise NotImplementedError("Must supply an integer axis value")

    if not isinstance(repeats, Integral):
        raise NotImplementedError("Only integer valued repeats supported")

    if -a.ndim <= axis < 0:
        axis += a.ndim
    elif not 0 <= axis <= a.ndim - 1:
        raise ValueError("axis(=%d) out of bounds" % axis)

    if repeats == 0:
        return a[tuple(slice(None) if d != axis else slice(0) for d in range(a.ndim))]
    elif repeats == 1:
        return a

    cchunks = cached_cumsum(a.chunks[axis], initial_zero=True)
    slices = []
    for c_start, c_stop in sliding_window(2, cchunks):
        ls = np.linspace(c_start, c_stop, repeats).round(0)
        for ls_start, ls_stop in sliding_window(2, ls):
            if ls_start != ls_stop:
                slices.append(slice(ls_start, ls_stop))

    all_slice = slice(None, None, None)
    slices = [
        (all_slice,) * axis + (s,) + (all_slice,) * (a.ndim - axis - 1) for s in slices
    ]

    slabs = [a[slc] for slc in slices]

    out = []
    for slab in slabs:
        chunks = list(slab.chunks)
        assert len(chunks[axis]) == 1
        chunks[axis] = (chunks[axis][0] * repeats,)
        chunks = tuple(chunks)
        result = slab.map_blocks(
            np.repeat, repeats, axis=axis, chunks=chunks, dtype=slab.dtype
        )
        out.append(result)

    return concatenate(out, axis=axis)


@derived_from(np)
def tile(A, reps):
    try:
        tup = tuple(reps)
    except TypeError:
        tup = (reps,)
    if any(i < 0 for i in tup):
        raise ValueError("Negative `reps` are not allowed.")
    c = asarray(A)

    if all(tup):
        for nrep in tup[::-1]:
            c = nrep * [c]
        return block(c)

    d = len(tup)
    if d < c.ndim:
        tup = (1,) * (c.ndim - d) + tup
    if c.ndim < d:
        shape = (1,) * (d - c.ndim) + c.shape
    else:
        shape = c.shape
    shape_out = tuple(s * t for s, t in zip(shape, tup))
    return empty(shape=shape_out, dtype=c.dtype)


def expand_pad_value(array, pad_value):
    if isinstance(pad_value, Number) or getattr(pad_value, "ndim", None) == 0:
        pad_value = array.ndim * ((pad_value, pad_value),)
    elif (
        isinstance(pad_value, Sequence)
        and all(isinstance(pw, Number) for pw in pad_value)
        and len(pad_value) == 1
    ):
        pad_value = array.ndim * ((pad_value[0], pad_value[0]),)
    elif (
        isinstance(pad_value, Sequence)
        and len(pad_value) == 2
        and all(isinstance(pw, Number) for pw in pad_value)
    ):
        pad_value = array.ndim * (tuple(pad_value),)
    elif (
        isinstance(pad_value, Sequence)
        and len(pad_value) == array.ndim
        and all(isinstance(pw, Sequence) for pw in pad_value)
        and all((len(pw) == 2) for pw in pad_value)
        and all(all(isinstance(w, Number) for w in pw) for pw in pad_value)
    ):
        pad_value = tuple(tuple(pw) for pw in pad_value)
    elif (
        isinstance(pad_value, Sequence)
        and len(pad_value) == 1
        and isinstance(pad_value[0], Sequence)
        and len(pad_value[0]) == 2
        and all(isinstance(pw, Number) for pw in pad_value[0])
    ):
        pad_value = array.ndim * (tuple(pad_value[0]),)
    else:
        raise TypeError("`pad_value` must be composed of integral typed values.")

    return pad_value


def get_pad_shapes_chunks(array, pad_width, axes):
    """
    Helper function for finding shapes and chunks of end pads.
    """

    pad_shapes = [list(array.shape), list(array.shape)]
    pad_chunks = [list(array.chunks), list(array.chunks)]

    for d in axes:
        for i in range(2):
            pad_shapes[i][d] = pad_width[d][i]
            pad_chunks[i][d] = (pad_width[d][i],)

    pad_shapes = [tuple(s) for s in pad_shapes]
    pad_chunks = [tuple(c) for c in pad_chunks]

    return pad_shapes, pad_chunks


def linear_ramp_chunk(start, stop, num, dim, step):
    """
    Helper function to find the linear ramp for a chunk.
    """
    num1 = num + 1

    shape = list(start.shape)
    shape[dim] = num
    shape = tuple(shape)

    dtype = np.dtype(start.dtype)

    result = np.empty_like(start, shape=shape, dtype=dtype)
    for i in np.ndindex(start.shape):
        j = list(i)
        j[dim] = slice(None)
        j = tuple(j)

        result[j] = np.linspace(start[i], stop, num1, dtype=dtype)[1:][::step]

    return result


def pad_edge(array, pad_width, mode, **kwargs):
    """
    Helper function for padding edges.

    Handles the cases where the only the values on the edge are needed.
    """

    kwargs = {k: expand_pad_value(array, v) for k, v in kwargs.items()}

    result = array
    for d in range(array.ndim):
        pad_shapes, pad_chunks = get_pad_shapes_chunks(result, pad_width, (d,))
        pad_arrays = [result, result]

        if mode == "constant":
            from dask.array.utils import asarray_safe

            constant_values = kwargs["constant_values"][d]
            constant_values = [
                asarray_safe(c, like=meta_from_array(array), dtype=result.dtype)
                for c in constant_values
            ]

            pad_arrays = [
                broadcast_to(v, s, c)
                for v, s, c in zip(constant_values, pad_shapes, pad_chunks)
            ]
        elif mode in ["edge", "linear_ramp"]:
            pad_slices = [result.ndim * [slice(None)], result.ndim * [slice(None)]]
            pad_slices[0][d] = slice(None, 1, None)
            pad_slices[1][d] = slice(-1, None, None)
            pad_slices = [tuple(sl) for sl in pad_slices]

            pad_arrays = [result[sl] for sl in pad_slices]

            if mode == "edge":
                pad_arrays = [
                    broadcast_to(a, s, c)
                    for a, s, c in zip(pad_arrays, pad_shapes, pad_chunks)
                ]
            elif mode == "linear_ramp":
                end_values = kwargs["end_values"][d]

                pad_arrays = [
                    a.map_blocks(
                        linear_ramp_chunk,
                        ev,
                        pw,
                        chunks=c,
                        dtype=result.dtype,
                        dim=d,
                        step=(2 * i - 1),
                    )
                    for i, (a, ev, pw, c) in enumerate(
                        zip(pad_arrays, end_values, pad_width[d], pad_chunks)
                    )
                ]
        elif mode == "empty":
            pad_arrays = [
                empty_like(array, shape=s, dtype=array.dtype, chunks=c)
                for s, c in zip(pad_shapes, pad_chunks)
            ]

        result = concatenate([pad_arrays[0], result, pad_arrays[1]], axis=d)

    return result


def pad_reuse(array, pad_width, mode, **kwargs):
    """
    Helper function for padding boundaries with values in the array.

    Handles the cases where the padding is constructed from values in
    the array. Namely by reflecting them or tiling them to create periodic
    boundary constraints.
    """

    if mode in {"reflect", "symmetric"}:
        reflect_type = kwargs.get("reflect", "even")
        if reflect_type == "odd":
            raise NotImplementedError("`pad` does not support `reflect_type` of `odd`.")
        if reflect_type != "even":
            raise ValueError(
                "unsupported value for reflect_type, must be one of (`even`, `odd`)"
            )

    result = np.empty(array.ndim * (3,), dtype=object)
    for idx in np.ndindex(result.shape):
        select = []
        orient = []
        for i, s, pw in zip(idx, array.shape, pad_width):
            if mode == "wrap":
                pw = pw[::-1]

            if i < 1:
                if mode == "reflect":
                    select.append(slice(1, pw[0] + 1, None))
                else:
                    select.append(slice(None, pw[0], None))
            elif i > 1:
                if mode == "reflect":
                    select.append(slice(s - pw[1] - 1, s - 1, None))
                else:
                    select.append(slice(s - pw[1], None, None))
            else:
                select.append(slice(None))

            if i != 1 and mode in ["reflect", "symmetric"]:
                orient.append(slice(None, None, -1))
            else:
                orient.append(slice(None))

        select = tuple(select)
        orient = tuple(orient)

        if mode == "wrap":
            idx = tuple(2 - i for i in idx)

        result[idx] = array[select][orient]

    result = block(result.tolist())

    return result


def pad_stats(array, pad_width, mode, stat_length):
    """
    Helper function for padding boundaries with statistics from the array.

    In cases where the padding requires computations of statistics from part
    or all of the array, this function helps compute those statistics as
    requested and then adds those statistics onto the boundaries of the array.
    """

    if mode == "median":
        raise NotImplementedError("`pad` does not support `mode` of `median`.")

    stat_length = expand_pad_value(array, stat_length)

    result = np.empty(array.ndim * (3,), dtype=object)
    for idx in np.ndindex(result.shape):
        axes = []
        select = []
        pad_shape = []
        pad_chunks = []
        for d, (i, s, c, w, l) in enumerate(
            zip(idx, array.shape, array.chunks, pad_width, stat_length)
        ):
            if i < 1:
                axes.append(d)
                select.append(slice(None, l[0], None))
                pad_shape.append(w[0])
                pad_chunks.append(w[0])
            elif i > 1:
                axes.append(d)
                select.append(slice(s - l[1], None, None))
                pad_shape.append(w[1])
                pad_chunks.append(w[1])
            else:
                select.append(slice(None))
                pad_shape.append(s)
                pad_chunks.append(c)

        axes = tuple(axes)
        select = tuple(select)
        pad_shape = tuple(pad_shape)
        pad_chunks = tuple(pad_chunks)

        result_idx = array[select]
        if axes:
            if mode == "maximum":
                result_idx = result_idx.max(axis=axes, keepdims=True)
            elif mode == "mean":
                result_idx = result_idx.mean(axis=axes, keepdims=True)
            elif mode == "minimum":
                result_idx = result_idx.min(axis=axes, keepdims=True)

            result_idx = broadcast_to(result_idx, pad_shape, chunks=pad_chunks)

            if mode == "mean":
                if np.issubdtype(array.dtype, np.integer):
                    result_idx = rint(result_idx)
                result_idx = result_idx.astype(array.dtype)

        result[idx] = result_idx

    result = block(result.tolist())

    return result


def wrapped_pad_func(array, pad_func, iaxis_pad_width, iaxis, pad_func_kwargs):
    result = np.empty_like(array)
    for i in np.ndindex(array.shape[:iaxis] + array.shape[iaxis + 1 :]):
        i = i[:iaxis] + (slice(None),) + i[iaxis:]
        result[i] = pad_func(array[i], iaxis_pad_width, iaxis, pad_func_kwargs)

    return result


def pad_udf(array, pad_width, mode, **kwargs):
    """
    Helper function for padding boundaries with a user defined function.

    In cases where the padding requires a custom user defined function be
    applied to the array, this function assists in the prepping and
    application of this function to the Dask Array to construct the desired
    boundaries.
    """

    result = pad_edge(array, pad_width, "constant", constant_values=0)

    chunks = result.chunks
    for d in range(result.ndim):
        result = result.rechunk(
            chunks[:d] + (result.shape[d : d + 1],) + chunks[d + 1 :]
        )

        result = result.map_blocks(
            wrapped_pad_func,
            name="pad",
            dtype=result.dtype,
            pad_func=mode,
            iaxis_pad_width=pad_width[d],
            iaxis=d,
            pad_func_kwargs=kwargs,
        )

        result = result.rechunk(chunks)

    return result


@derived_from(np)
def pad(array, pad_width, mode="constant", **kwargs):
    array = asarray(array)

    pad_width = expand_pad_value(array, pad_width)

    if callable(mode):
        return pad_udf(array, pad_width, mode, **kwargs)

    # Make sure that no unsupported keywords were passed for the current mode
    allowed_kwargs = {
        "empty": [],
        "edge": [],
        "wrap": [],
        "constant": ["constant_values"],
        "linear_ramp": ["end_values"],
        "maximum": ["stat_length"],
        "mean": ["stat_length"],
        "median": ["stat_length"],
        "minimum": ["stat_length"],
        "reflect": ["reflect_type"],
        "symmetric": ["reflect_type"],
    }
    try:
        unsupported_kwargs = set(kwargs) - set(allowed_kwargs[mode])
    except KeyError as e:
        raise ValueError(f"mode '{mode}' is not supported") from e
    if unsupported_kwargs:
        raise ValueError(
            "unsupported keyword arguments for mode '{}': {}".format(
                mode, unsupported_kwargs
            )
        )

    if mode in {"maximum", "mean", "median", "minimum"}:
        stat_length = kwargs.get("stat_length", tuple((n, n) for n in array.shape))
        return pad_stats(array, pad_width, mode, stat_length)
    elif mode == "constant":
        kwargs.setdefault("constant_values", 0)
        return pad_edge(array, pad_width, mode, **kwargs)
    elif mode == "linear_ramp":
        kwargs.setdefault("end_values", 0)
        return pad_edge(array, pad_width, mode, **kwargs)
    elif mode in {"edge", "empty"}:
        return pad_edge(array, pad_width, mode)
    elif mode in ["reflect", "symmetric", "wrap"]:
        return pad_reuse(array, pad_width, mode, **kwargs)

    raise RuntimeError("unreachable")
