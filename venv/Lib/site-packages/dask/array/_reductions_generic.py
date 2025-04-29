from __future__ import annotations

import builtins
import contextlib
import math
import warnings
from functools import partial
from itertools import product
from numbers import Integral

import numpy as np
from tlz import compose, get, partition_all

from dask import config
from dask.array.blockwise import blockwise
from dask.array.core import Array, _concatenate2, asanyarray, broadcast_to, handle_out
from dask.array.numpy_compat import ComplexWarning
from dask.array.utils import compute_meta, is_arraylike, validate_axis
from dask.base import tokenize
from dask.blockwise import lol_tuples
from dask.highlevelgraph import HighLevelGraph
from dask.utils import funcname, getargspec, is_series_like


def reduction(
    x,
    chunk,
    aggregate,
    axis=None,
    keepdims=False,
    dtype=None,
    split_every=None,
    combine=None,
    name=None,
    out=None,
    concatenate=True,
    output_size=1,
    meta=None,
    weights=None,
):
    """General version of reductions

    Parameters
    ----------
    x: Array
        Data being reduced along one or more axes
    chunk: callable(x_chunk, [weights_chunk=None], axis, keepdims)
        First function to be executed when resolving the dask graph.
        This function is applied in parallel to all original chunks of x.
        See below for function parameters.
    combine: callable(x_chunk, axis, keepdims), optional
        Function used for intermediate recursive aggregation (see
        split_every below). If omitted, it defaults to aggregate.
        If the reduction can be performed in less than 3 steps, it will not
        be invoked at all.
    aggregate: callable(x_chunk, axis, keepdims)
        Last function to be executed when resolving the dask graph,
        producing the final output. It is always invoked, even when the reduced
        Array counts a single chunk along the reduced axes.
    axis: int or sequence of ints, optional
        Axis or axes to aggregate upon. If omitted, aggregate along all axes.
    keepdims: boolean, optional
        Whether the reduction function should preserve the reduced axes,
        leaving them at size ``output_size``, or remove them.
    dtype: np.dtype
        data type of output. This argument was previously optional, but
        leaving as ``None`` will now raise an exception.
    split_every: int >= 2 or dict(axis: int), optional
        Determines the depth of the recursive aggregation. If set to or more
        than the number of input chunks, the aggregation will be performed in
        two steps, one ``chunk`` function per input chunk and a single
        ``aggregate`` function at the end. If set to less than that, an
        intermediate ``combine`` function will be used, so that any one
        ``combine`` or ``aggregate`` function has no more than ``split_every``
        inputs. The depth of the aggregation graph will be
        :math:`log_{split_every}(input chunks along reduced axes)`. Setting to
        a low value can reduce cache size and network transfers, at the cost of
        more CPU and a larger dask graph.

        Omit to let dask heuristically decide a good default. A default can
        also be set globally with the ``split_every`` key in
        :mod:`dask.config`.
    name: str, optional
        Prefix of the keys of the intermediate and output nodes. If omitted it
        defaults to the function names.
    out: Array, optional
        Another dask array whose contents will be replaced. Omit to create a
        new one. Note that, unlike in numpy, this setting gives no performance
        benefits whatsoever, but can still be useful  if one needs to preserve
        the references to a previously existing Array.
    concatenate: bool, optional
        If True (the default), the outputs of the ``chunk``/``combine``
        functions are concatenated into a single np.array before being passed
        to the ``combine``/``aggregate`` functions. If False, the input of
        ``combine`` and ``aggregate`` will be either a list of the raw outputs
        of the previous step or a single output, and the function will have to
        concatenate it itself. It can be useful to set this to False if the
        chunk and/or combine steps do not produce np.arrays.
    output_size: int >= 1, optional
        Size of the output of the ``aggregate`` function along the reduced
        axes. Ignored if keepdims is False.
    weights : array_like, optional
        Weights to be used in the reduction of `x`. Will be
        automatically broadcast to the shape of `x`, and so must have
        a compatible shape. For instance, if `x` has shape ``(3, 4)``
        then acceptable shapes for `weights` are ``(3, 4)``, ``(4,)``,
        ``(3, 1)``, ``(1, 1)``, ``(1)``, and ``()``.

    Returns
    -------
    dask array

    **Function Parameters**

    x_chunk: numpy.ndarray
        Individual input chunk. For ``chunk`` functions, it is one of the
        original chunks of x. For ``combine`` and ``aggregate`` functions, it's
        the concatenation of the outputs produced by the previous ``chunk`` or
        ``combine`` functions. If concatenate=False, it's a list of the raw
        outputs from the previous functions.
    weights_chunk: numpy.ndarray, optional
        Only applicable to the ``chunk`` function. Weights, with the
        same shape as `x_chunk`, to be applied during the reduction of
        the individual input chunk. If ``weights`` have not been
        provided then the function may omit this parameter. When
        `weights_chunk` is included then it must occur immediately
        after the `x_chunk` parameter, and must also have a default
        value for cases when ``weights`` are not provided.
    axis: tuple
        Normalized list of axes to reduce upon, e.g. ``(0, )``
        Scalar, negative, and None axes have been normalized away.
        Note that some numpy reduction functions cannot reduce along multiple
        axes at once and strictly require an int in input. Such functions have
        to be wrapped to cope.
    keepdims: bool
        Whether the reduction function should preserve the reduced axes or
        remove them.

    """
    if axis is None:
        axis = tuple(range(x.ndim))
    if isinstance(axis, Integral):
        axis = (axis,)
    axis = validate_axis(axis, x.ndim)

    if dtype is None:
        raise ValueError("Must specify dtype")
    if "dtype" in getargspec(chunk).args:
        chunk = partial(chunk, dtype=dtype)
    if "dtype" in getargspec(aggregate).args:
        aggregate = partial(aggregate, dtype=dtype)
    if is_series_like(x):
        x = x.values

    # Map chunk across all blocks
    inds = tuple(range(x.ndim))

    args = (x, inds)
    if weights is not None:
        # Broadcast weights to x and add to args
        wgt = asanyarray(weights)
        try:
            wgt = broadcast_to(wgt, x.shape)
        except ValueError:
            raise ValueError(
                f"Weights with shape {wgt.shape} are not broadcastable "
                f"to x with shape {x.shape}"
            )

        args += (wgt, inds)

    # The dtype of `tmp` doesn't actually matter, and may be incorrect.
    tmp = blockwise(
        chunk, inds, *args, axis=axis, keepdims=True, token=name, dtype=dtype or float
    )
    tmp._chunks = tuple(
        (output_size,) * len(c) if i in axis else c for i, c in enumerate(tmp.chunks)
    )

    if meta is None and hasattr(x, "_meta"):
        try:
            reduced_meta = compute_meta(
                chunk, x.dtype, x._meta, axis=axis, keepdims=True, computing_meta=True
            )
        except TypeError:
            reduced_meta = compute_meta(
                chunk, x.dtype, x._meta, axis=axis, keepdims=True
            )
        except ValueError:
            pass
    else:
        reduced_meta = None

    result = _tree_reduce(
        tmp,
        aggregate,
        axis,
        keepdims,
        dtype,
        split_every,
        combine,
        name=name,
        concatenate=concatenate,
        reduced_meta=reduced_meta,
    )
    if keepdims and output_size != 1:
        result._chunks = tuple(
            (output_size,) if i in axis else c for i, c in enumerate(tmp.chunks)
        )
    if meta is not None:
        result._meta = meta
    return handle_out(out, result)


def _tree_reduce(
    x,
    aggregate,
    axis,
    keepdims,
    dtype,
    split_every=None,
    combine=None,
    name=None,
    concatenate=True,
    reduced_meta=None,
):
    """Perform the tree reduction step of a reduction.

    Lower level, users should use ``reduction`` or ``arg_reduction`` directly.
    """
    # Normalize split_every
    split_every = split_every or config.get("split_every", 4)
    if isinstance(split_every, dict):
        split_every = {k: split_every.get(k, 2) for k in axis}
    elif isinstance(split_every, Integral):
        n = builtins.max(int(split_every ** (1 / (len(axis) or 1))), 2)
        split_every = dict.fromkeys(axis, n)
    else:
        raise ValueError("split_every must be a int or a dict")

    # Reduce across intermediates
    depth = 1
    for i, n in enumerate(x.numblocks):
        if i in split_every and split_every[i] != 1:
            depth = int(builtins.max(depth, math.ceil(math.log(n, split_every[i]))))
    func = partial(combine or aggregate, axis=axis, keepdims=True)
    if concatenate:
        func = compose(func, partial(_concatenate2, axes=sorted(axis)))
    for _ in range(depth - 1):
        x = partial_reduce(
            func,
            x,
            split_every,
            True,
            dtype=dtype,
            name=(name or funcname(combine or aggregate)) + "-partial",
            reduced_meta=reduced_meta,
        )
    func = partial(aggregate, axis=axis, keepdims=keepdims)
    if concatenate:
        func = compose(func, partial(_concatenate2, axes=sorted(axis)))
    return partial_reduce(
        func,
        x,
        split_every,
        keepdims=keepdims,
        dtype=dtype,
        name=(name or funcname(aggregate)) + "-aggregate",
        reduced_meta=reduced_meta,
    )


def partial_reduce(
    func, x, split_every, keepdims=False, dtype=None, name=None, reduced_meta=None
):
    """Partial reduction across multiple axes.

    Parameters
    ----------
    func : function
    x : Array
    split_every : dict
        Maximum reduction block sizes in each dimension.

    Examples
    --------
    Reduce across axis 0 and 2, merging a maximum of 1 block in the 0th
    dimension, and 3 blocks in the 2nd dimension:

    >>> partial_reduce(np.min, x, {0: 1, 2: 3})    # doctest: +SKIP
    """
    name = (
        (name or funcname(func)) + "-" + tokenize(func, x, split_every, keepdims, dtype)
    )
    parts = [
        list(partition_all(split_every.get(i, 1), range(n)))
        for (i, n) in enumerate(x.numblocks)
    ]
    keys = product(*map(range, map(len, parts)))
    out_chunks = [
        tuple(1 for p in partition_all(split_every[i], c)) if i in split_every else c
        for (i, c) in enumerate(x.chunks)
    ]
    if not keepdims:
        out_axis = [i for i in range(x.ndim) if i not in split_every]
        getter = lambda k: get(out_axis, k)
        keys = map(getter, keys)
        out_chunks = list(getter(out_chunks))
    dsk = {}
    for k, p in zip(keys, product(*parts)):
        free = {
            i: j[0] for (i, j) in enumerate(p) if len(j) == 1 and i not in split_every
        }
        dummy = dict(i for i in enumerate(p) if i[0] in split_every)
        g = lol_tuples((x.name,), range(x.ndim), free, dummy)
        dsk[(name,) + k] = (func, g)
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[x])

    meta = x._meta
    if reduced_meta is not None:
        try:
            meta = func(reduced_meta, computing_meta=True)
        # no meta keyword argument exists for func, and it isn't required
        except TypeError:
            try:
                meta = func(reduced_meta)
            except ValueError as e:
                # min/max functions have no identity, don't apply function to meta
                if "zero-size array to reduction operation" in str(e):
                    meta = reduced_meta
        # when no work can be computed on the empty array (e.g., func is a ufunc)
        except ValueError:
            pass

    # some functions can't compute empty arrays (those for which reduced_meta
    # fall into the ValueError exception) and we have to rely on reshaping
    # the array according to len(out_chunks)
    if is_arraylike(meta) and meta.ndim != len(out_chunks):
        if len(out_chunks) == 0:
            meta = meta.sum()
        else:
            meta = meta.reshape((0,) * len(out_chunks))

    if np.isscalar(meta):
        return Array(graph, name, out_chunks, dtype=dtype)
    else:
        with contextlib.suppress(AttributeError), warnings.catch_warnings():
            if name.startswith("var") or name.startswith("moment"):
                warnings.simplefilter("ignore", ComplexWarning)
            meta = meta.astype(dtype)
        return Array(graph, name, out_chunks, meta=meta)
