from __future__ import annotations

import numbers
import warnings
from math import isnan

import numpy as np
import tlz as toolz

from dask import base, utils
from dask.blockwise import _blockwise_unpack_collections_task_spec
from dask.blockwise import blockwise as core_blockwise
from dask.highlevelgraph import HighLevelGraph
from dask.layers import ArrayBlockwiseDep


def blockwise(
    func,
    out_ind,
    *args,
    name=None,
    token=None,
    dtype=None,
    adjust_chunks=None,
    new_axes=None,
    align_arrays=True,
    concatenate=None,
    meta=None,
    **kwargs,
):
    """Tensor operation: Generalized inner and outer products

    A broad class of blocked algorithms and patterns can be specified with a
    concise multi-index notation.  The ``blockwise`` function applies an in-memory
    function across multiple blocks of multiple inputs in a variety of ways.
    Many dask.array operations are special cases of blockwise including
    elementwise, broadcasting, reductions, tensordot, and transpose.

    Parameters
    ----------
    func : callable
        Function to apply to individual tuples of blocks
    out_ind : iterable
        Block pattern of the output, something like 'ijk' or (1, 2, 3)
    *args : sequence of Array, index pairs
        You may also pass literal arguments, accompanied by None index
        e.g. (x, 'ij', y, 'jk', z, 'i', some_literal, None)
    **kwargs : dict
        Extra keyword arguments to pass to function
    dtype : np.dtype
        Datatype of resulting array.
    concatenate : bool, keyword only
        If true concatenate arrays along dummy indices, else provide lists
    adjust_chunks : dict
        Dictionary mapping index to function to be applied to chunk sizes
    new_axes : dict, keyword only
        New indexes and their dimension lengths
    align_arrays: bool
        Whether or not to align chunks along equally sized dimensions when
        multiple arrays are provided.  This allows for larger chunks in some
        arrays to be broken into smaller ones that match chunk sizes in other
        arrays such that they are compatible for block function mapping. If
        this is false, then an error will be thrown if arrays do not already
        have the same number of blocks in each dimension.

    Examples
    --------
    2D embarrassingly parallel operation from two arrays, x, and y.

    >>> import operator, numpy as np, dask.array as da
    >>> x = da.from_array([[1, 2],
    ...                    [3, 4]], chunks=(1, 2))
    >>> y = da.from_array([[10, 20],
    ...                    [0, 0]])
    >>> z = blockwise(operator.add, 'ij', x, 'ij', y, 'ij', dtype='f8')
    >>> z.compute()
    array([[11, 22],
           [ 3,  4]])

    Outer product multiplying a by b, two 1-d vectors

    >>> a = da.from_array([0, 1, 2], chunks=1)
    >>> b = da.from_array([10, 50, 100], chunks=1)
    >>> z = blockwise(np.outer, 'ij', a, 'i', b, 'j', dtype='f8')
    >>> z.compute()
    array([[  0,   0,   0],
           [ 10,  50, 100],
           [ 20, 100, 200]])

    z = x.T

    >>> z = blockwise(np.transpose, 'ji', x, 'ij', dtype=x.dtype)
    >>> z.compute()
    array([[1, 3],
           [2, 4]])

    The transpose case above is illustrative because it does transposition
    both on each in-memory block by calling ``np.transpose`` and on the order
    of the blocks themselves, by switching the order of the index ``ij -> ji``.

    We can compose these same patterns with more variables and more complex
    in-memory functions

    z = X + Y.T

    >>> z = blockwise(lambda x, y: x + y.T, 'ij', x, 'ij', y, 'ji', dtype='f8')
    >>> z.compute()
    array([[11,  2],
           [23,  4]])

    Any index, like ``i`` missing from the output index is interpreted as a
    contraction (note that this differs from Einstein convention; repeated
    indices do not imply contraction.)  In the case of a contraction the passed
    function should expect an iterable of blocks on any array that holds that
    index.  To receive arrays concatenated along contracted dimensions instead
    pass ``concatenate=True``.

    Inner product multiplying a by b, two 1-d vectors

    >>> def sequence_dot(a_blocks, b_blocks):
    ...     result = 0
    ...     for a, b in zip(a_blocks, b_blocks):
    ...         result += a.dot(b)
    ...     return result

    >>> z = blockwise(sequence_dot, '', a, 'i', b, 'i', dtype='f8')
    >>> z.compute()
    np.int64(250)

    Add new single-chunk dimensions with the ``new_axes=`` keyword, including
    the length of the new dimension.  New dimensions will always be in a single
    chunk.

    >>> def f(a):
    ...     return a[:, None] * np.ones((1, 5))

    >>> z = blockwise(f, 'az', a, 'a', new_axes={'z': 5}, dtype=a.dtype)

    New dimensions can also be multi-chunk by specifying a tuple of chunk
    sizes.  This has limited utility as is (because the chunks are all the
    same), but the resulting graph can be modified to achieve more useful
    results (see ``da.map_blocks``).

    >>> z = blockwise(f, 'az', a, 'a', new_axes={'z': (5, 5)}, dtype=x.dtype)
    >>> z.chunks
    ((1, 1, 1), (5, 5))

    If the applied function changes the size of each chunk you can specify this
    with a ``adjust_chunks={...}`` dictionary holding a function for each index
    that modifies the dimension size in that index.

    >>> def double(x):
    ...     return np.concatenate([x, x])

    >>> y = blockwise(double, 'ij', x, 'ij',
    ...               adjust_chunks={'i': lambda n: 2 * n}, dtype=x.dtype)
    >>> y.chunks
    ((2, 2), (2,))

    Include literals by indexing with None

    >>> z = blockwise(operator.add, 'ij', x, 'ij', 1234, None, dtype=x.dtype)
    >>> z.compute()
    array([[1235, 1236],
           [1237, 1238]])
    """
    out = name
    new_axes = new_axes or {}

    # Input Validation
    if len(set(out_ind)) != len(out_ind):
        raise ValueError(
            "Repeated elements not allowed in output index",
            [k for k, v in toolz.frequencies(out_ind).items() if v > 1],
        )
    new = (
        set(out_ind)
        - {a for arg in args[1::2] if arg is not None for a in arg}
        - set(new_axes or ())
    )
    if new:
        raise ValueError("Unknown dimension", new)

    from dask.array.core import normalize_arg, unify_chunks

    if align_arrays:
        chunkss, arrays = unify_chunks(*args)
    else:
        arginds = [(a, i) for (a, i) in toolz.partition(2, args) if i is not None]
        chunkss = {}
        # For each dimension, use the input chunking that has the most blocks;
        # this will ensure that broadcasting works as expected, and in
        # particular the number of blocks should be correct if the inputs are
        # consistent.
        for arg, ind in arginds:
            for c, i in zip(arg.chunks, ind):
                if i not in chunkss or len(c) > len(chunkss[i]):
                    chunkss[i] = c
        arrays = args[::2]

    for k, v in new_axes.items():
        if not isinstance(v, tuple):
            v = (v,)
        chunkss[k] = v

    arginds = list(zip(arrays, args[1::2]))
    numblocks = {}

    dependencies = []
    arrays = []

    # Normalize arguments
    argindsstr = []

    for arg, ind in arginds:
        if ind is None:
            arg = normalize_arg(arg)
            arg, collections = _blockwise_unpack_collections_task_spec(arg)

            dependencies.extend(collections)
        else:
            if (
                hasattr(arg, "ndim")
                and hasattr(ind, "__len__")
                and arg.ndim != len(ind)
            ):
                raise ValueError(
                    "Index string %s does not match array dimension %d"
                    % (ind, arg.ndim)
                )
            if not isinstance(arg, ArrayBlockwiseDep):
                numblocks[arg.name] = arg.numblocks
                arrays.append(arg)
                arg = arg.name
        argindsstr.extend((arg, ind))

    # Normalize keyword arguments
    kwargs2 = {}
    for k, v in kwargs.items():
        v = normalize_arg(v)
        v, collections = _blockwise_unpack_collections_task_spec(v)
        dependencies.extend(collections)
        kwargs2[k] = v

    # Finish up the name
    if not out:
        out = "{}-{}".format(
            token or utils.funcname(func).strip("_"),
            base.tokenize(
                func,
                out_ind,
                argindsstr,
                adjust_chunks,
                new_axes,
                align_arrays,
                concatenate,
                meta,
                dtype,
                **kwargs,
            ),
        )

    graph = core_blockwise(
        func,
        out,
        out_ind,
        *argindsstr,
        numblocks=numblocks,
        dependencies=dependencies,
        new_axes=new_axes,
        concatenate=concatenate,
        **kwargs2,
    )
    graph = HighLevelGraph.from_collections(
        out, graph, dependencies=arrays + dependencies
    )

    chunks = [chunkss[i] for i in out_ind]
    if adjust_chunks:
        for i, ind in enumerate(out_ind):
            if ind in adjust_chunks:
                if callable(adjust_chunks[ind]):
                    chunks[i] = tuple(map(adjust_chunks[ind], chunks[i]))
                elif isinstance(adjust_chunks[ind], numbers.Integral) or (
                    np.isscalar(adjust_chunks[ind]) and isnan(adjust_chunks[ind])
                ):
                    chunks[i] = tuple(adjust_chunks[ind] for _ in chunks[i])
                elif isinstance(adjust_chunks[ind], (tuple, list)):
                    if len(adjust_chunks[ind]) != len(chunks[i]):
                        raise ValueError(
                            f"Dimension {i} has {len(chunks[i])} blocks, adjust_chunks "
                            f"specified with {len(adjust_chunks[ind])} blocks"
                        )
                    chunks[i] = tuple(adjust_chunks[ind])
                else:
                    raise NotImplementedError(
                        "adjust_chunks values must be callable, int, or tuple"
                    )
    chunks = tuple(chunks)

    if meta is None:
        from dask.array.utils import compute_meta

        meta = compute_meta(func, dtype, *args[::2], **kwargs)
    return new_da_object(graph, out, chunks, meta=meta, dtype=dtype)


def atop(*args, **kwargs):
    warnings.warn("The da.atop function has moved to da.blockwise")
    return blockwise(*args, **kwargs)


from dask.array.core import new_da_object
