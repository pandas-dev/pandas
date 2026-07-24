from __future__ import annotations

import warnings
from itertools import product
from numbers import Number

from toolz import concat

from dask.array._array_expr._collection import Array, blockwise
from dask.array.core import _pass_extra_kwargs, apply_and_enforce, apply_infer_dtype
from dask.array.utils import compute_meta
from dask.layers import ArrayBlockIdDep, ArrayValuesDep
from dask.utils import cached_cumsum, funcname, has_keyword


def map_blocks(
    func,
    *args,
    name=None,
    token=None,
    dtype=None,
    chunks=None,
    drop_axis=None,
    new_axis=None,
    enforce_ndim=False,
    meta=None,
    **kwargs,
):
    """Map a function across all blocks of a dask array.

    Note that ``map_blocks`` will attempt to automatically determine the output
    array type by calling ``func`` on 0-d versions of the inputs. Please refer to
    the ``meta`` keyword argument below if you expect that the function will not
    succeed when operating on 0-d arrays.

    Parameters
    ----------
    func : callable
        Function to apply to every block in the array.
        If ``func`` accepts ``block_info=`` or ``block_id=``
        as keyword arguments, these will be passed dictionaries
        containing information about input and output chunks/arrays
        during computation. See examples for details.
    args : dask arrays or other objects
    dtype : np.dtype, optional
        The ``dtype`` of the output array. It is recommended to provide this.
        If not provided, will be inferred by applying the function to a small
        set of fake data.
    chunks : tuple, optional
        Chunk shape of resulting blocks if the function does not preserve
        shape. If not provided, the resulting array is assumed to have the same
        block structure as the first input array.
    drop_axis : number or iterable, optional
        Dimensions lost by the function.
    new_axis : number or iterable, optional
        New dimensions created by the function. Note that these are applied
        after ``drop_axis`` (if present). The size of each chunk along this
        dimension will be set to 1. Please specify ``chunks`` if the individual
        chunks have a different size.
    enforce_ndim : bool, default False
        Whether to enforce at runtime that the dimensionality of the array
        produced by ``func`` actually matches that of the array returned by
        ``map_blocks``.
        If True, this will raise an error when there is a mismatch.
    token : string, optional
        The key prefix to use for the output array. If not provided, will be
        determined from the function name.
    name : string, optional
        The key name to use for the output array. Note that this fully
        specifies the output key name, and must be unique. If not provided,
        will be determined by a hash of the arguments.
    meta : array-like, optional
        The ``meta`` of the output array, when specified is expected to be an
        array of the same type and dtype of that returned when calling ``.compute()``
        on the array returned by this function. When not provided, ``meta`` will be
        inferred by applying the function to a small set of fake data, usually a
        0-d array. It's important to ensure that ``func`` can successfully complete
        computation without raising exceptions when 0-d is passed to it, providing
        ``meta`` will be required otherwise. If the output type is known beforehand
        (e.g., ``np.ndarray``, ``cupy.ndarray``), an empty array of such type dtype
        can be passed, for example: ``meta=np.array((), dtype=np.int32)``.
    **kwargs :
        Other keyword arguments to pass to function. Values must be constants
        (not dask.arrays)

    See Also
    --------
    dask.array.map_overlap : Generalized operation with overlap between neighbors.
    dask.array.blockwise : Generalized operation with control over block alignment.

    Examples
    --------
    >>> import dask.array as da
    >>> x = da.arange(6, chunks=3)

    >>> x.map_blocks(lambda x: x * 2).compute()
    array([ 0,  2,  4,  6,  8, 10])

    The ``da.map_blocks`` function can also accept multiple arrays.

    >>> d = da.arange(5, chunks=2)
    >>> e = da.arange(5, chunks=2)

    >>> f = da.map_blocks(lambda a, b: a + b**2, d, e)
    >>> f.compute()
    array([ 0,  2,  6, 12, 20])

    If the function changes shape of the blocks then you must provide chunks
    explicitly.

    >>> y = x.map_blocks(lambda x: x[::2], chunks=((2, 2),))

    You have a bit of freedom in specifying chunks.  If all of the output chunk
    sizes are the same, you can provide just that chunk size as a single tuple.

    >>> a = da.arange(18, chunks=(6,))
    >>> b = a.map_blocks(lambda x: x[:3], chunks=(3,))

    If the function changes the dimension of the blocks you must specify the
    created or destroyed dimensions.

    >>> b = a.map_blocks(lambda x: x[None, :, None], chunks=(1, 6, 1),
    ...                  new_axis=[0, 2])

    If ``chunks`` is specified but ``new_axis`` is not, then it is inferred to
    add the necessary number of axes on the left.

    Note that ``map_blocks()`` will concatenate chunks along axes specified by
    the keyword parameter ``drop_axis`` prior to applying the function.
    This is illustrated in the figure below:

    .. image:: /images/map_blocks_drop_axis.png

    Due to memory-size-constraints, it is often not advisable to use ``drop_axis``
    on an axis that is chunked.  In that case, it is better not to use
    ``map_blocks`` but rather
    ``dask.array.reduction(..., axis=dropped_axes, concatenate=False)`` which
    maintains a leaner memory footprint while it drops any axis.

    Map_blocks aligns blocks by block positions without regard to shape. In the
    following example we have two arrays with the same number of blocks but
    with different shape and chunk sizes.

    >>> x = da.arange(1000, chunks=(100,))
    >>> y = da.arange(100, chunks=(10,))

    The relevant attribute to match is numblocks.

    >>> x.numblocks
    (10,)
    >>> y.numblocks
    (10,)

    If these match (up to broadcasting rules) then we can map arbitrary
    functions across blocks

    >>> def func(a, b):
    ...     return np.array([a.max(), b.max()])

    >>> da.map_blocks(func, x, y, chunks=(2,), dtype='i8')
    dask.array<func, shape=(20,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>

    >>> _.compute()
    array([ 99,   9, 199,  19, 299,  29, 399,  39, 499,  49, 599,  59, 699,
            69, 799,  79, 899,  89, 999,  99])

    Your block function can get information about where it is in the array by
    accepting a special ``block_info`` or ``block_id`` keyword argument.
    During computation, they will contain information about each of the input
    and output chunks (and dask arrays) relevant to each call of ``func``.

    >>> def func(block_info=None):
    ...     pass

    This will receive the following information:

    >>> block_info  # doctest: +SKIP
    {0: {'shape': (1000,),
         'num-chunks': (10,),
         'chunk-location': (4,),
         'array-location': [(400, 500)]},
     None: {'shape': (1000,),
            'num-chunks': (10,),
            'chunk-location': (4,),
            'array-location': [(400, 500)],
            'chunk-shape': (100,),
            'dtype': dtype('float64')}}

    The keys to the ``block_info`` dictionary indicate which is the input and
    output Dask array:

    - **Input Dask array(s):** ``block_info[0]`` refers to the first input Dask array.
      The dictionary key is ``0`` because that is the argument index corresponding
      to the first input Dask array.
      In cases where multiple Dask arrays have been passed as input to the function,
      you can access them with the number corresponding to the input argument,
      eg: ``block_info[1]``, ``block_info[2]``, etc.
      (Note that if you pass multiple Dask arrays as input to map_blocks,
      the arrays must match each other by having matching numbers of chunks,
      along corresponding dimensions up to broadcasting rules.)
    - **Output Dask array:** ``block_info[None]`` refers to the output Dask array,
      and contains information about the output chunks.
      The output chunk shape and dtype may may be different than the input chunks.

    For each dask array, ``block_info`` describes:

    - ``shape``: the shape of the full Dask array,
    - ``num-chunks``: the number of chunks of the full array in each dimension,
    - ``chunk-location``: the chunk location (for example the fourth chunk over
      in the first dimension), and
    - ``array-location``: the array location within the full Dask array
      (for example the slice corresponding to ``40:50``).

    In addition to these, there are two extra parameters described by
    ``block_info`` for the output array (in ``block_info[None]``):

    - ``chunk-shape``: the output chunk shape, and
    - ``dtype``: the output dtype.

    These features can be combined to synthesize an array from scratch, for
    example:

    >>> def func(block_info=None):
    ...     loc = block_info[None]['array-location'][0]
    ...     return np.arange(loc[0], loc[1])

    >>> da.map_blocks(func, chunks=((4, 4),), dtype=np.float64)
    dask.array<func, shape=(8,), dtype=float64, chunksize=(4,), chunktype=numpy.ndarray>

    >>> _.compute()
    array([0, 1, 2, 3, 4, 5, 6, 7])

    ``block_id`` is similar to ``block_info`` but contains only the ``chunk_location``:

    >>> def func(block_id=None):
    ...     pass

    This will receive the following information:

    >>> block_id  # doctest: +SKIP
    (4, 3)

    You may specify the key name prefix of the resulting task in the graph with
    the optional ``token`` keyword argument.

    >>> x.map_blocks(lambda x: x + 1, name='increment')
    dask.array<increment, shape=(1000,), dtype=int64, chunksize=(100,), chunktype=numpy.ndarray>

    For functions that may not handle 0-d arrays, it's also possible to specify
    ``meta`` with an empty array matching the type of the expected result. In
    the example below, ``func`` will result in an ``IndexError`` when computing
    ``meta``:

    >>> rng = da.random.default_rng()
    >>> da.map_blocks(lambda x: x[2], rng.random(5), meta=np.array(()))
    dask.array<lambda, shape=(5,), dtype=float64, chunksize=(5,), chunktype=numpy.ndarray>

    Similarly, it's possible to specify a non-NumPy array to ``meta``, and provide
    a ``dtype``:

    >>> import cupy  # doctest: +SKIP
    >>> rng = da.random.default_rng(cupy.random.default_rng())  # doctest: +SKIP
    >>> dt = np.float32
    >>> da.map_blocks(lambda x: x[2], rng.random(5, dtype=dt), meta=cupy.array((), dtype=dt))  # doctest: +SKIP
    dask.array<lambda, shape=(5,), dtype=float32, chunksize=(5,), chunktype=cupy.ndarray>
    """
    if drop_axis is None:
        drop_axis = []

    if not callable(func):
        msg = (
            "First argument must be callable function, not %s\n"
            "Usage:   da.map_blocks(function, x)\n"
            "   or:   da.map_blocks(function, x, y, z)"
        )
        raise TypeError(msg % type(func).__name__)
    if token:
        warnings.warn(
            "The `token=` keyword to `map_blocks` has been moved to `name=`. "
            "Please use `name=` instead as the `token=` keyword will be removed "
            "in a future release.",
            category=FutureWarning,
        )
        name = token

    name = f"{name or funcname(func)}"
    new_axes = {}

    if isinstance(drop_axis, Number):
        drop_axis = [drop_axis]
    if isinstance(new_axis, Number):
        new_axis = [new_axis]  # TODO: handle new_axis

    arrs = [a for a in args if isinstance(a, Array)]

    argpairs = [
        (a, tuple(range(a.ndim))[::-1]) if isinstance(a, Array) else (a, None)
        for a in args
    ]
    if arrs:
        out_ind = tuple(range(max(a.ndim for a in arrs)))[::-1]
    else:
        out_ind = ()

    original_kwargs = kwargs

    if dtype is None and meta is None:
        try:
            meta = compute_meta(func, dtype, *args, **kwargs)
        except Exception:
            pass

        dtype = apply_infer_dtype(func, args, original_kwargs, "map_blocks")

    if drop_axis:
        ndim_out = len(out_ind)
        if any(i < -ndim_out or i >= ndim_out for i in drop_axis):
            raise ValueError(
                f"drop_axis out of range (drop_axis={drop_axis}, "
                f"but output is {ndim_out}d)."
            )
        drop_axis = [i % ndim_out for i in drop_axis]
        out_ind = tuple(x for i, x in enumerate(out_ind) if i not in drop_axis)
    if new_axis is None and chunks is not None and len(out_ind) < len(chunks):
        new_axis = range(len(chunks) - len(out_ind))
    if new_axis:
        # new_axis = [x + len(drop_axis) for x in new_axis]
        out_ind = list(out_ind)
        for ax in sorted(new_axis):
            n = len(out_ind) + len(drop_axis)
            out_ind.insert(ax, n)
            if chunks is not None:
                new_axes[n] = chunks[ax]
            else:
                new_axes[n] = 1
        out_ind = tuple(out_ind)
        if max(new_axis) > max(out_ind):
            raise ValueError("New_axis values do not fill in all dimensions")

    if chunks is not None:
        if len(chunks) != len(out_ind):
            raise ValueError(
                f"Provided chunks have {len(chunks)} dims; expected {len(out_ind)} dims"
            )
        adjust_chunks = dict(zip(out_ind, chunks))
    else:
        adjust_chunks = None

    if enforce_ndim:
        out = blockwise(
            apply_and_enforce,
            out_ind,
            *concat(argpairs),
            expected_ndim=len(out_ind),
            _func=func,
            token=name,
            new_axes=new_axes,
            dtype=dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=adjust_chunks,
            meta=meta,
            **kwargs,
        )
    else:
        out = blockwise(
            func,
            out_ind,
            *concat(argpairs),
            token=name,
            new_axes=new_axes,
            dtype=dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=adjust_chunks,
            meta=meta,
            **kwargs,
        )

    extra_argpairs = []
    extra_names = []

    # If func has block_id as an argument, construct an object to inject it.
    if has_keyword(func, "block_id"):
        extra_argpairs.append((ArrayBlockIdDep(out.chunks), out_ind))
        extra_names.append("block_id")

    if has_keyword(func, "_overlap_trim_info"):
        # Internal for map overlap to reduce size of graph
        num_chunks = out.numblocks
        block_id_dict = {
            block_id: (block_id, num_chunks)
            for block_id in product(*(range(len(c)) for c in out.chunks))
        }
        extra_argpairs.append((ArrayValuesDep(out.chunks, block_id_dict), out_ind))
        extra_names.append("_overlap_trim_info")

    # If func has block_info as an argument, construct a dict of block info
    # objects and prepare to inject it.
    if has_keyword(func, "block_info"):
        starts = {}
        num_chunks = {}
        shapes = {}

        for i, (arg, in_ind) in enumerate(argpairs):
            if in_ind is not None:
                shapes[i] = arg.shape
                if drop_axis:
                    # We concatenate along dropped axes, so we need to treat them
                    # as if there is only a single chunk.
                    starts[i] = [
                        (
                            cached_cumsum(arg.chunks[j], initial_zero=True)
                            if ind in out_ind
                            else [0, arg.shape[j]]
                        )
                        for j, ind in enumerate(in_ind)
                    ]
                    num_chunks[i] = tuple(len(s) - 1 for s in starts[i])
                else:
                    starts[i] = [
                        cached_cumsum(c, initial_zero=True) for c in arg.chunks
                    ]
                    num_chunks[i] = arg.numblocks
        out_starts = [cached_cumsum(c, initial_zero=True) for c in out.chunks]

        block_info_dict = {}
        for block_id in product(*(range(len(c)) for c in out.chunks)):
            # Get position of chunk, indexed by axis labels
            location = {out_ind[i]: loc for i, loc in enumerate(block_id)}
            info = {}
            for i, shape in shapes.items():
                # Compute chunk key in the array, taking broadcasting into
                # account. We don't directly know which dimensions are
                # broadcast, but any dimension with only one chunk can be
                # treated as broadcast.
                arr_k = tuple(
                    location.get(ind, 0) if num_chunks[i][j] > 1 else 0
                    for j, ind in enumerate(argpairs[i][1])
                )
                info[i] = {
                    "shape": shape,
                    "num-chunks": num_chunks[i],
                    "array-location": [
                        (starts[i][ij][j], starts[i][ij][j + 1])
                        for ij, j in enumerate(arr_k)
                    ],
                    "chunk-location": arr_k,
                }

            info[None] = {
                "shape": out.shape,
                "num-chunks": out.numblocks,
                "array-location": [
                    (out_starts[ij][j], out_starts[ij][j + 1])
                    for ij, j in enumerate(block_id)
                ],
                "chunk-location": block_id,
                "chunk-shape": tuple(
                    out.chunks[ij][j] for ij, j in enumerate(block_id)
                ),
                "dtype": dtype,
            }
            block_info_dict[block_id] = info

        extra_argpairs.append((ArrayValuesDep(out.chunks, block_info_dict), out_ind))
        extra_names.append("block_info")

    if extra_argpairs:
        # Rewrite the Blockwise layer. It would be nice to find a way to
        # avoid doing it twice, but it's currently needed to determine
        # out.chunks from the first pass. Since it constructs a Blockwise
        # rather than an expanded graph, it shouldn't be too expensive.
        out = blockwise(
            _pass_extra_kwargs,
            out_ind,
            func,
            None,
            tuple(extra_names),
            None,
            *concat(extra_argpairs),
            *concat(argpairs),
            dtype=out.dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=dict(zip(out_ind, out.chunks)),
            meta=meta,
            **kwargs,
        )

    return out
