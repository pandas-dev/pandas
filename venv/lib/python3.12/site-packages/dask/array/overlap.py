from __future__ import annotations

import warnings
from functools import reduce
from numbers import Integral, Number
from operator import mul

import numpy as np
from tlz import concat, get, partial
from tlz.curried import map

from dask._compatibility import import_optional_dependency
from dask.array import chunk
from dask.array._shuffle import _calculate_new_chunksizes
from dask.array.core import Array, broadcast_to, concatenate, map_blocks, unify_chunks
from dask.array.creation import arange, empty_like, full_like, repeat
from dask.array.numpy_compat import normalize_axis_tuple
from dask.array.reductions import cumreduction
from dask.array.routines import notnull, where
from dask.base import tokenize
from dask.highlevelgraph import HighLevelGraph
from dask.layers import ArrayOverlapLayer
from dask.utils import derived_from


def _overlap_internal_chunks(original_chunks, axes):
    """Get new chunks for array with overlap."""
    chunks = []
    for i, bds in enumerate(original_chunks):
        depth = axes.get(i, 0)
        if isinstance(depth, tuple):
            left_depth = depth[0]
            right_depth = depth[1]
        else:
            left_depth = depth
            right_depth = depth

        if len(bds) == 1:
            chunks.append(bds)
        else:
            left = [bds[0] + right_depth]
            right = [bds[-1] + left_depth]
            mid = []
            for bd in bds[1:-1]:
                mid.append(bd + left_depth + right_depth)
            chunks.append(left + mid + right)
    return chunks


def overlap_internal(x, axes):
    """Share boundaries between neighboring blocks

    Parameters
    ----------

    x: da.Array
        A dask array
    axes: dict
        The size of the shared boundary per axis

    The axes input informs how many cells to overlap between neighboring blocks
    {0: 2, 2: 5} means share two cells in 0 axis, 5 cells in 2 axis
    """
    token = tokenize(x, axes)
    name = "overlap-" + token

    graph = ArrayOverlapLayer(
        name=x.name,
        axes=axes,
        chunks=x.chunks,
        numblocks=x.numblocks,
        token=token,
    )
    graph = HighLevelGraph.from_collections(name, graph, dependencies=[x])
    chunks = _overlap_internal_chunks(x.chunks, axes)

    return Array(graph, name, chunks, meta=x)


def trim_overlap(x, depth, boundary=None):
    """Trim sides from each block.

    This couples well with the ``map_overlap`` operation which may leave
    excess data on each block.

    See also
    --------
    dask.array.overlap.map_overlap

    """

    # parameter to be passed to trim_internal
    axes = coerce_depth(x.ndim, depth)
    return trim_internal(x, axes=axes, boundary=boundary)


def trim_internal(x, axes, boundary=None):
    """Trim sides from each block

    This couples well with the overlap operation, which may leave excess data on
    each block

    See also
    --------
    dask.array.chunk.trim
    dask.array.map_blocks
    """
    boundary = coerce_boundary(x.ndim, boundary)

    olist = []
    for i, bd in enumerate(x.chunks):
        bdy = boundary.get(i, "none")
        overlap = axes.get(i, 0)
        ilist = []
        for j, d in enumerate(bd):
            if bdy != "none":
                if isinstance(overlap, tuple):
                    d = d - sum(overlap)
                else:
                    d = d - overlap * 2

            else:
                if isinstance(overlap, tuple):
                    d = d - overlap[0] if j != 0 else d
                    d = d - overlap[1] if j != len(bd) - 1 else d
                else:
                    d = d - overlap if j != 0 else d
                    d = d - overlap if j != len(bd) - 1 else d

            ilist.append(d)
        olist.append(tuple(ilist))
    chunks = tuple(olist)

    return map_blocks(
        partial(_trim, axes=axes, boundary=boundary),
        x,
        chunks=chunks,
        dtype=x.dtype,
        meta=x._meta,
    )


def _trim(x, axes, boundary, _overlap_trim_info):
    """Similar to dask.array.chunk.trim but requires one to specify the
    boundary condition.

    ``axes``, and ``boundary`` are assumed to have been coerced.

    """
    chunk_location = _overlap_trim_info[0]
    num_chunks = _overlap_trim_info[1]
    axes = [axes.get(i, 0) for i in range(x.ndim)]
    axes_front = (ax[0] if isinstance(ax, tuple) else ax for ax in axes)
    axes_back = (
        (
            -ax[1]
            if isinstance(ax, tuple) and ax[1]
            else -ax if isinstance(ax, Integral) and ax else None
        )
        for ax in axes
    )

    trim_front = (
        0 if (chunk_location == 0 and boundary.get(i, "none") == "none") else ax
        for i, (chunk_location, ax) in enumerate(zip(chunk_location, axes_front))
    )
    trim_back = (
        (
            None
            if (chunk_location == chunks - 1 and boundary.get(i, "none") == "none")
            else ax
        )
        for i, (chunks, chunk_location, ax) in enumerate(
            zip(num_chunks, chunk_location, axes_back)
        )
    )
    ind = tuple(slice(front, back) for front, back in zip(trim_front, trim_back))
    return x[ind]


def periodic(x, axis, depth):
    """Copy a slice of an array around to its other side

    Useful to create periodic boundary conditions for overlap
    """

    left = (
        (slice(None, None, None),) * axis
        + (slice(0, depth),)
        + (slice(None, None, None),) * (x.ndim - axis - 1)
    )
    right = (
        (slice(None, None, None),) * axis
        + (slice(-depth, None),)
        + (slice(None, None, None),) * (x.ndim - axis - 1)
    )
    l = x[left]
    r = x[right]

    l, r = _remove_overlap_boundaries(l, r, axis, depth)

    return concatenate([r, x, l], axis=axis)


def reflect(x, axis, depth):
    """Reflect boundaries of array on the same side

    This is the converse of ``periodic``
    """
    if depth == 1:
        left = (
            (slice(None, None, None),) * axis
            + (slice(0, 1),)
            + (slice(None, None, None),) * (x.ndim - axis - 1)
        )
    else:
        left = (
            (slice(None, None, None),) * axis
            + (slice(depth - 1, None, -1),)
            + (slice(None, None, None),) * (x.ndim - axis - 1)
        )
    right = (
        (slice(None, None, None),) * axis
        + (slice(-1, -depth - 1, -1),)
        + (slice(None, None, None),) * (x.ndim - axis - 1)
    )
    l = x[left]
    r = x[right]

    l, r = _remove_overlap_boundaries(l, r, axis, depth)

    return concatenate([l, x, r], axis=axis)


def nearest(x, axis, depth):
    """Each reflect each boundary value outwards

    This mimics what the skimage.filters.gaussian_filter(... mode="nearest")
    does.
    """
    left = (
        (slice(None, None, None),) * axis
        + (slice(0, 1),)
        + (slice(None, None, None),) * (x.ndim - axis - 1)
    )
    right = (
        (slice(None, None, None),) * axis
        + (slice(-1, -2, -1),)
        + (slice(None, None, None),) * (x.ndim - axis - 1)
    )

    l = repeat(x[left], depth, axis=axis)
    r = repeat(x[right], depth, axis=axis)

    l, r = _remove_overlap_boundaries(l, r, axis, depth)

    return concatenate([l, x, r], axis=axis)


def constant(x, axis, depth, value):
    """Add constant slice to either side of array"""
    chunks = list(x.chunks)
    chunks[axis] = (depth,)

    c = full_like(
        x,
        value,
        shape=tuple(map(sum, chunks)),
        chunks=tuple(chunks),
        dtype=x.dtype,
    )

    return concatenate([c, x, c], axis=axis)


def _remove_overlap_boundaries(l, r, axis, depth):
    lchunks = list(l.chunks)
    lchunks[axis] = (depth,)
    rchunks = list(r.chunks)
    rchunks[axis] = (depth,)

    l = l.rechunk(tuple(lchunks))
    r = r.rechunk(tuple(rchunks))
    return l, r


def boundaries(x, depth=None, kind=None):
    """Add boundary conditions to an array before overlapping

    See Also
    --------
    periodic
    constant
    """
    if not isinstance(kind, dict):
        kind = dict.fromkeys(range(x.ndim), kind)
    if not isinstance(depth, dict):
        depth = dict.fromkeys(range(x.ndim), depth)

    for i in range(x.ndim):
        d = depth.get(i, 0)
        if d == 0:
            continue

        this_kind = kind.get(i, "none")
        if this_kind == "none":
            continue
        elif this_kind == "periodic":
            x = periodic(x, i, d)
        elif this_kind == "reflect":
            x = reflect(x, i, d)
        elif this_kind == "nearest":
            x = nearest(x, i, d)
        elif i in kind:
            x = constant(x, i, d, kind[i])

    return x


def ensure_minimum_chunksize(size, chunks):
    """Determine new chunks to ensure that every chunk >= size

    Parameters
    ----------
    size: int
        The maximum size of any chunk.
    chunks: tuple
        Chunks along one axis, e.g. ``(3, 3, 2)``

    Examples
    --------
    >>> ensure_minimum_chunksize(10, (20, 20, 1))
    (20, 11, 10)
    >>> ensure_minimum_chunksize(3, (1, 1, 3))
    (5,)

    See Also
    --------
    overlap
    """
    if size <= min(chunks):
        return chunks

    # add too-small chunks to chunks before them
    output = []
    new = 0
    for c in chunks:
        if c < size:
            if new > size + (size - c):
                output.append(new - (size - c))
                new = size
            else:
                new += c
        if new >= size:
            output.append(new)
            new = 0
        if c >= size:
            new += c
    if new >= size:
        output.append(new)
    elif len(output) >= 1:
        output[-1] += new
    else:
        raise ValueError(
            f"The overlapping depth {size} is larger than your array {sum(chunks)}."
        )

    return tuple(output)


def _get_overlap_rechunked_chunks(x, depth2):
    depths = [max(d) if isinstance(d, tuple) else d for d in depth2.values()]
    # rechunk if new chunks are needed to fit depth in every chunk
    return tuple(ensure_minimum_chunksize(size, c) for size, c in zip(depths, x.chunks))


def overlap(x, depth, boundary, *, allow_rechunk=True):
    """Share boundaries between neighboring blocks

    Parameters
    ----------

    x: da.Array
        A dask array
    depth: dict
        The size of the shared boundary per axis
    boundary: dict
        The boundary condition on each axis. Options are 'reflect', 'periodic',
        'nearest', 'none', or an array value.  Such a value will fill the
        boundary with that value.
    allow_rechunk: bool, keyword only
        Allows rechunking, otherwise chunk sizes need to match and core
        dimensions are to consist only of one chunk.

    The depth input informs how many cells to overlap between neighboring
    blocks ``{0: 2, 2: 5}`` means share two cells in 0 axis, 5 cells in 2 axis.
    Axes missing from this input will not be overlapped.

    Any axis containing chunks smaller than depth will be rechunked if
    possible, provided the keyword ``allow_rechunk`` is True (recommended).

    Examples
    --------
    >>> import numpy as np
    >>> import dask.array as da

    >>> x = np.arange(64).reshape((8, 8))
    >>> d = da.from_array(x, chunks=(4, 4))
    >>> d.chunks
    ((4, 4), (4, 4))

    >>> g = da.overlap.overlap(d, depth={0: 2, 1: 1},
    ...                       boundary={0: 100, 1: 'reflect'})
    >>> g.chunks
    ((8, 8), (6, 6))

    >>> np.array(g)
    array([[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [  0,   0,   1,   2,   3,   4,   3,   4,   5,   6,   7,   7],
           [  8,   8,   9,  10,  11,  12,  11,  12,  13,  14,  15,  15],
           [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
           [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
           [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
           [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
           [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
           [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
           [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
           [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
           [ 48,  48,  49,  50,  51,  52,  51,  52,  53,  54,  55,  55],
           [ 56,  56,  57,  58,  59,  60,  59,  60,  61,  62,  63,  63],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]])
    """
    depth2 = coerce_depth(x.ndim, depth)
    boundary2 = coerce_boundary(x.ndim, boundary)

    depths = [max(d) if isinstance(d, tuple) else d for d in depth2.values()]
    if allow_rechunk:
        # rechunk if new chunks are needed to fit depth in every chunk
        x1 = x.rechunk(
            _get_overlap_rechunked_chunks(x, depth2)
        )  # this is a no-op if x.chunks == new_chunks

    else:
        original_chunks_too_small = any(min(c) < d for d, c in zip(depths, x.chunks))
        if original_chunks_too_small:
            raise ValueError(
                "Overlap depth is larger than smallest chunksize.\n"
                "Please set allow_rechunk=True to rechunk automatically.\n"
                f"Overlap depths required: {depths}\n"
                f"Input chunks: {x.chunks}\n"
            )
        x1 = x

    x2 = boundaries(x1, depth2, boundary2)
    x3 = overlap_internal(x2, depth2)
    trim = {
        k: v * 2 if boundary2.get(k, "none") != "none" else 0 for k, v in depth2.items()
    }
    x4 = chunk.trim(x3, trim)
    return x4


def add_dummy_padding(x, depth, boundary):
    """
    Pads an array which has 'none' as the boundary type.
    Used to simplify trimming arrays which use 'none'.

    >>> import dask.array as da
    >>> x = da.arange(6, chunks=3)
    >>> add_dummy_padding(x, {0: 1}, {0: 'none'}).compute()  # doctest: +NORMALIZE_WHITESPACE
    array([..., 0, 1, 2, 3, 4, 5, ...])
    """
    for k, v in boundary.items():
        d = depth.get(k, 0)
        if v == "none" and d > 0:
            empty_shape = list(x.shape)
            empty_shape[k] = d

            empty_chunks = list(x.chunks)
            empty_chunks[k] = (d,)

            empty = empty_like(
                getattr(x, "_meta", x),
                shape=empty_shape,
                chunks=empty_chunks,
                dtype=x.dtype,
            )

            out_chunks = list(x.chunks)
            ax_chunks = list(out_chunks[k])
            ax_chunks[0] += d
            ax_chunks[-1] += d
            out_chunks[k] = tuple(ax_chunks)

            x = concatenate([empty, x, empty], axis=k)
            x = x.rechunk(out_chunks)
    return x


def map_overlap(
    func,
    *args,
    depth=None,
    boundary=None,
    trim=True,
    align_arrays=True,
    allow_rechunk=True,
    **kwargs,
):
    """Map a function over blocks of arrays with some overlap

    We share neighboring zones between blocks of the array, map a
    function, and then trim away the neighboring strips. If depth is
    larger than any chunk along a particular axis, then the array is
    rechunked.

    Note that this function will attempt to automatically determine the output
    array type before computing it, please refer to the ``meta`` keyword argument
    in ``map_blocks`` if you expect that the function will not succeed when
    operating on 0-d arrays.

    Parameters
    ----------
    func: function
        The function to apply to each extended block.
        If multiple arrays are provided, then the function should expect to
        receive chunks of each array in the same order.
    args : dask arrays
    depth: int, tuple, dict or list, keyword only
        The number of elements that each block should share with its neighbors
        If a tuple or dict then this can be different per axis.
        If a list then each element of that list must be an int, tuple or dict
        defining depth for the corresponding array in `args`.
        Asymmetric depths may be specified using a dict value of (-/+) tuples.
        Note that asymmetric depths are currently only supported when
        ``boundary`` is 'none'.
        The default value is 0.
    boundary: str, tuple, dict or list, keyword only
        How to handle the boundaries.
        Values include 'reflect', 'periodic', 'nearest', 'none',
        or any constant value like 0 or np.nan.
        If a list then each element must be a str, tuple or dict defining the
        boundary for the corresponding array in `args`.
        The default value is 'reflect'.
    trim: bool, keyword only
        Whether or not to trim ``depth`` elements from each block after
        calling the map function.
        Set this to False if your mapping function already does this for you
    align_arrays: bool, keyword only
        Whether or not to align chunks along equally sized dimensions when
        multiple arrays are provided.  This allows for larger chunks in some
        arrays to be broken into smaller ones that match chunk sizes in other
        arrays such that they are compatible for block function mapping. If
        this is false, then an error will be thrown if arrays do not already
        have the same number of blocks in each dimension.
    allow_rechunk: bool, keyword only
        Allows rechunking, otherwise chunk sizes need to match and core
        dimensions are to consist only of one chunk.
    **kwargs:
        Other keyword arguments valid in ``map_blocks``

    Examples
    --------
    >>> import numpy as np
    >>> import dask.array as da

    >>> x = np.array([1, 1, 2, 3, 3, 3, 2, 1, 1])
    >>> x = da.from_array(x, chunks=5)
    >>> def derivative(x):
    ...     return x - np.roll(x, 1)

    >>> y = x.map_overlap(derivative, depth=1, boundary=0)
    >>> y.compute()
    array([ 1,  0,  1,  1,  0,  0, -1, -1,  0])

    >>> x = np.arange(16).reshape((4, 4))
    >>> d = da.from_array(x, chunks=(2, 2))
    >>> d.map_overlap(lambda x: x + x.size, depth=1, boundary='reflect').compute()
    array([[16, 17, 18, 19],
           [20, 21, 22, 23],
           [24, 25, 26, 27],
           [28, 29, 30, 31]])

    >>> func = lambda x: x + x.size
    >>> depth = {0: 1, 1: 1}
    >>> boundary = {0: 'reflect', 1: 'none'}
    >>> d.map_overlap(func, depth, boundary).compute()  # doctest: +NORMALIZE_WHITESPACE
    array([[12,  13,  14,  15],
           [16,  17,  18,  19],
           [20,  21,  22,  23],
           [24,  25,  26,  27]])

    The ``da.map_overlap`` function can also accept multiple arrays.

    >>> func = lambda x, y: x + y
    >>> x = da.arange(8).reshape(2, 4).rechunk((1, 2))
    >>> y = da.arange(4).rechunk(2)
    >>> da.map_overlap(func, x, y, depth=1, boundary='reflect').compute() # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  2,  4,  6],
           [ 4,  6,  8,  10]])

    When multiple arrays are given, they do not need to have the
    same number of dimensions but they must broadcast together.
    Arrays are aligned block by block (just as in ``da.map_blocks``)
    so the blocks must have a common chunk size.  This common chunking
    is determined automatically as long as ``align_arrays`` is True.

    >>> x = da.arange(8, chunks=4)
    >>> y = da.arange(8, chunks=2)
    >>> r = da.map_overlap(func, x, y, depth=1, boundary='reflect', align_arrays=True)
    >>> len(r.to_delayed())
    4

    >>> da.map_overlap(func, x, y, depth=1, boundary='reflect', align_arrays=False).compute()
    Traceback (most recent call last):
        ...
    ValueError: Shapes do not align {'.0': {2, 4}}

    Note also that this function is equivalent to ``map_blocks``
    by default.  A non-zero ``depth`` must be defined for any
    overlap to appear in the arrays provided to ``func``.

    >>> func = lambda x: x.sum()
    >>> x = da.ones(10, dtype='int')
    >>> block_args = dict(chunks=(), drop_axis=0)
    >>> da.map_blocks(func, x, **block_args).compute()
    np.int64(10)
    >>> da.map_overlap(func, x, **block_args, boundary='reflect').compute()
    np.int64(10)
    >>> da.map_overlap(func, x, **block_args, depth=1, boundary='reflect').compute()
    np.int64(12)

    For functions that may not handle 0-d arrays, it's also possible to specify
    ``meta`` with an empty array matching the type of the expected result. In
    the example below, ``func`` will result in an ``IndexError`` when computing
    ``meta``:

    >>> x = np.arange(16).reshape((4, 4))
    >>> d = da.from_array(x, chunks=(2, 2))
    >>> y = d.map_overlap(lambda x: x + x[2], depth=1, boundary='reflect', meta=np.array(()))
    >>> y
    dask.array<_trim, shape=(4, 4), dtype=float64, chunksize=(2, 2), chunktype=numpy.ndarray>
    >>> y.compute()
    array([[ 4,  6,  8, 10],
           [ 8, 10, 12, 14],
           [20, 22, 24, 26],
           [24, 26, 28, 30]])

    Similarly, it's possible to specify a non-NumPy array to ``meta``:

    >>> import cupy  # doctest: +SKIP
    >>> x = cupy.arange(16).reshape((4, 4))  # doctest: +SKIP
    >>> d = da.from_array(x, chunks=(2, 2))  # doctest: +SKIP
    >>> y = d.map_overlap(lambda x: x + x[2], depth=1, boundary='reflect', meta=cupy.array(()))  # doctest: +SKIP
    >>> y  # doctest: +SKIP
    dask.array<_trim, shape=(4, 4), dtype=float64, chunksize=(2, 2), chunktype=cupy.ndarray>
    >>> y.compute()  # doctest: +SKIP
    array([[ 4,  6,  8, 10],
           [ 8, 10, 12, 14],
           [20, 22, 24, 26],
           [24, 26, 28, 30]])
    """
    # Look for invocation using deprecated single-array signature
    # map_overlap(x, func, depth, boundary=None, trim=True, **kwargs)
    if isinstance(func, Array) and callable(args[0]):
        warnings.warn(
            "The use of map_overlap(array, func, **kwargs) is deprecated since dask 2.17.0 "
            "and will be an error in a future release. To silence this warning, use the syntax "
            "map_overlap(func, array0,[ array1, ...,] **kwargs) instead.",
            FutureWarning,
        )
        sig = ["func", "depth", "boundary", "trim"]
        depth = get(sig.index("depth"), args, depth)
        boundary = get(sig.index("boundary"), args, boundary)
        trim = get(sig.index("trim"), args, trim)
        func, args = args[0], [func]

    if not callable(func):
        raise TypeError(
            "First argument must be callable function, not {}\n"
            "Usage:   da.map_overlap(function, x)\n"
            "   or:   da.map_overlap(function, x, y, z)".format(type(func).__name__)
        )
    if not all(isinstance(x, Array) for x in args):
        raise TypeError(
            "All variadic arguments must be arrays, not {}\n"
            "Usage:   da.map_overlap(function, x)\n"
            "   or:   da.map_overlap(function, x, y, z)".format(
                [type(x).__name__ for x in args]
            )
        )

    # Coerce depth and boundary arguments to lists of individual
    # specifications for each array argument
    def coerce(xs, arg, fn):
        if not isinstance(arg, list):
            arg = [arg] * len(xs)
        return [fn(x.ndim, a) for x, a in zip(xs, arg)]

    depth = coerce(args, depth, coerce_depth)
    boundary = coerce(args, boundary, coerce_boundary)

    # Align chunks in each array to a common size
    if align_arrays:
        # Reverse unification order to allow block broadcasting
        inds = [list(reversed(range(x.ndim))) for x in args]
        _, args = unify_chunks(*list(concat(zip(args, inds))), warn=False)

    # Escape to map_blocks if depth is zero (a more efficient computation)
    if all(all(depth_val == 0 for depth_val in d.values()) for d in depth):
        return map_blocks(func, *args, **kwargs)

    for i, x in enumerate(args):
        for j in range(x.ndim):
            if isinstance(depth[i][j], tuple) and boundary[i][j] != "none":
                raise NotImplementedError(
                    "Asymmetric overlap is currently only implemented "
                    "for boundary='none', however boundary for dimension "
                    "{} in array argument {} is {}".format(j, i, boundary[i][j])
                )

    def assert_int_chunksize(xs):
        assert all(type(c) is int for x in xs for cc in x.chunks for c in cc)

    assert_int_chunksize(args)
    if not trim and "chunks" not in kwargs:
        if allow_rechunk:
            # Adjust chunks based on the rechunking result
            kwargs["chunks"] = _get_overlap_rechunked_chunks(
                args[0], coerce_depth(args[0].ndim, depth[0])
            )
        else:
            kwargs["chunks"] = args[0].chunks
    args = [
        overlap(x, depth=d, boundary=b, allow_rechunk=allow_rechunk)
        for x, d, b in zip(args, depth, boundary)
    ]
    assert_int_chunksize(args)
    x = map_blocks(func, *args, **kwargs)
    assert_int_chunksize([x])
    if trim:
        # Find index of array argument with maximum rank and break ties by choosing first provided
        i = sorted(enumerate(args), key=lambda v: (v[1].ndim, -v[0]))[-1][0]
        # Trim using depth/boundary setting for array of highest rank
        depth = depth[i]
        boundary = boundary[i]
        # remove any dropped axes from depth and boundary variables
        drop_axis = kwargs.pop("drop_axis", None)
        if drop_axis is not None:
            if isinstance(drop_axis, Number):
                drop_axis = [drop_axis]

            # convert negative drop_axis to equivalent positive value
            ndim_out = max(a.ndim for a in args if isinstance(a, Array))
            drop_axis = [d % ndim_out for d in drop_axis]

            kept_axes = tuple(ax for ax in range(args[i].ndim) if ax not in drop_axis)
            # note that keys are relabeled to match values in range(x.ndim)
            depth = {n: depth[ax] for n, ax in enumerate(kept_axes)}
            boundary = {n: boundary[ax] for n, ax in enumerate(kept_axes)}

        # add any new axes to depth and boundary variables
        new_axis = kwargs.pop("new_axis", None)
        if new_axis is not None:
            if isinstance(new_axis, Number):
                new_axis = [new_axis]

            # convert negative new_axis to equivalent positive value
            ndim_out = max(a.ndim for a in args if isinstance(a, Array))
            new_axis = [d % ndim_out for d in new_axis]

            for axis in new_axis:
                for existing_axis in list(depth.keys()):
                    if existing_axis >= axis:
                        # Shuffle existing axis forward to give room to insert new_axis
                        depth[existing_axis + 1] = depth[existing_axis]
                        boundary[existing_axis + 1] = boundary[existing_axis]

                depth[axis] = 0
                boundary[axis] = "none"

        return trim_internal(x, depth, boundary)
    else:
        return x


def coerce_depth(ndim, depth):
    default = 0
    if depth is None:
        depth = default
    if isinstance(depth, Integral):
        depth = (depth,) * ndim
    if isinstance(depth, tuple):
        depth = dict(zip(range(ndim), depth))
    if isinstance(depth, dict):
        depth = {ax: depth.get(ax, default) for ax in range(ndim)}
    return coerce_depth_type(ndim, depth)


def coerce_depth_type(ndim, depth):
    for i in range(ndim):
        if isinstance(depth[i], tuple):
            depth[i] = tuple(int(d) for d in depth[i])
        else:
            depth[i] = int(depth[i])
    return depth


def coerce_boundary(ndim, boundary):
    default = "none"
    if boundary is None:
        boundary = default
    if not isinstance(boundary, (tuple, dict)):
        boundary = (boundary,) * ndim
    if isinstance(boundary, tuple):
        boundary = dict(zip(range(ndim), boundary))
    if isinstance(boundary, dict):
        boundary = {ax: boundary.get(ax, default) for ax in range(ndim)}
    return boundary


@derived_from(np.lib.stride_tricks)
def sliding_window_view(x, window_shape, axis=None, automatic_rechunk=True):
    window_shape = tuple(window_shape) if np.iterable(window_shape) else (window_shape,)

    window_shape_array = np.array(window_shape)
    if np.any(window_shape_array <= 0):
        raise ValueError("`window_shape` must contain values > 0")

    if axis is None:
        axis = tuple(range(x.ndim))
        if len(window_shape) != len(axis):
            raise ValueError(
                f"Since axis is `None`, must provide "
                f"window_shape for all dimensions of `x`; "
                f"got {len(window_shape)} window_shape elements "
                f"and `x.ndim` is {x.ndim}."
            )
    else:
        axis = normalize_axis_tuple(axis, x.ndim, allow_duplicate=True)
        if len(window_shape) != len(axis):
            raise ValueError(
                f"Must provide matching length window_shape and "
                f"axis; got {len(window_shape)} window_shape "
                f"elements and {len(axis)} axes elements."
            )

    depths = [0] * x.ndim
    for ax, window in zip(axis, window_shape):
        depths[ax] += window - 1

    # Ensure that each chunk is big enough to leave at least a size-1 chunk
    # after windowing (this is only really necessary for the last chunk).
    safe_chunks = list(
        ensure_minimum_chunksize(d + 1, c) for d, c in zip(depths, x.chunks)
    )
    if automatic_rechunk:
        safe_chunks = [
            s if d != 0 else c for d, c, s in zip(depths, x.chunks, safe_chunks)
        ]
        # safe chunks is our output chunks, so add the new dimensions
        safe_chunks.extend([(w,) for w in window_shape])
        max_chunk = reduce(mul, map(max, x.chunks))
        new_chunks = _calculate_new_chunksizes(
            x.chunks,
            safe_chunks.copy(),
            {i for i, d in enumerate(depths) if d == 0},
            max_chunk,
        )
        x = x.rechunk(tuple(new_chunks))
    else:
        x = x.rechunk(tuple(safe_chunks))

    # result.shape = x_shape_trimmed + window_shape,
    # where x_shape_trimmed is x.shape with every entry
    # reduced by one less than the corresponding window size.
    # trim chunks to match x_shape_trimmed
    newchunks = tuple(c[:-1] + (c[-1] - d,) for d, c in zip(depths, x.chunks)) + tuple(
        (window,) for window in window_shape
    )

    return map_overlap(
        np.lib.stride_tricks.sliding_window_view,
        x,
        depth=tuple((0, d) for d in depths),  # Overlap on +ve side only
        boundary="none",
        meta=x._meta,
        new_axis=range(x.ndim, x.ndim + len(axis)),
        chunks=newchunks,
        trim=False,
        align_arrays=False,
        window_shape=window_shape,
        axis=axis,
    )


def push(array, n, axis):
    """
    Dask-version of bottleneck.push

    .. note::

        Requires bottleneck to be installed.
    """
    import_optional_dependency("bottleneck", min_version="1.3.7")

    if n is not None and 0 < n < array.shape[axis] - 1:
        arr = broadcast_to(
            arange(
                array.shape[axis], chunks=array.chunks[axis], dtype=array.dtype
            ).reshape(
                tuple(size if i == axis else 1 for i, size in enumerate(array.shape))
            ),
            array.shape,
            array.chunks,
        )
        valid_arange = where(notnull(array), arr, np.nan)
        valid_limits = (arr - push(valid_arange, None, axis)) <= n
        # omit the forward fill that violate the limit
        return where(valid_limits, push(array, None, axis), np.nan)

    # The method parameter makes that the tests for python 3.7 fails.
    return cumreduction(
        func=_push,
        binop=_fill_with_last_one,
        ident=np.nan,
        x=array,
        axis=axis,
        dtype=array.dtype,
    )


def _fill_with_last_one(a, b):
    # cumreduction apply the push func over all the blocks first so, the only missing part is filling
    # the missing values using the last data of the previous chunk
    return np.where(~np.isnan(b), b, a)


def _push(array, n: int | None = None, axis: int = -1):
    # work around for bottleneck 178
    limit = n if n is not None else array.shape[axis]

    import bottleneck as bn

    return bn.push(array, limit, axis)
