from __future__ import annotations

import builtins
import contextlib
import math
import operator
import warnings
from collections.abc import Iterable
from functools import partial
from itertools import product, repeat
from numbers import Integral, Number

import numpy as np
from tlz import accumulate, compose, drop, get, partition_all, pluck

from dask import config
from dask.array import chunk
from dask.array.blockwise import blockwise
from dask.array.core import (
    Array,
    _concatenate2,
    asanyarray,
    broadcast_to,
    handle_out,
    implements,
    unknown_chunk_message,
)
from dask.array.creation import arange, diagonal
from dask.array.dispatch import divide_lookup, nannumel_lookup, numel_lookup
from dask.array.numpy_compat import ComplexWarning
from dask.array.utils import (
    array_safe,
    asarray_safe,
    compute_meta,
    is_arraylike,
    meta_from_array,
    validate_axis,
)
from dask.array.wrap import ones, zeros
from dask.base import tokenize
from dask.blockwise import lol_tuples
from dask.highlevelgraph import HighLevelGraph
from dask.utils import (
    apply,
    deepmap,
    derived_from,
    funcname,
    getargspec,
    is_series_like,
)


def divide(a, b, dtype=None):
    key = lambda x: getattr(x, "__array_priority__", float("-inf"))
    f = divide_lookup.dispatch(type(builtins.max(a, b, key=key)))
    return f(a, b, dtype=dtype)


def numel(x, **kwargs):
    return numel_lookup(x, **kwargs)


def nannumel(x, **kwargs):
    return nannumel_lookup(x, **kwargs)


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


@derived_from(np)
def sum(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is None:
        dtype = getattr(np.zeros(1, dtype=a.dtype).sum(), "dtype", object)
    result = reduction(
        a,
        chunk.sum,
        chunk.sum,
        axis=axis,
        keepdims=keepdims,
        dtype=dtype,
        split_every=split_every,
        out=out,
    )
    return result


@derived_from(np)
def prod(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(np.ones((1,), dtype=a.dtype).prod(), "dtype", object)
    return reduction(
        a,
        chunk.prod,
        chunk.prod,
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        out=out,
    )


@implements(np.min, np.amin)
@derived_from(np)
def min(a, axis=None, keepdims=False, split_every=None, out=None):
    return reduction(
        a,
        chunk_min,
        chunk.min,
        combine=chunk_min,
        axis=axis,
        keepdims=keepdims,
        dtype=a.dtype,
        split_every=split_every,
        out=out,
    )


def chunk_min(x, axis=None, keepdims=None):
    """Version of np.min which ignores size 0 arrays"""
    if x.size == 0:
        return array_safe([], x, ndmin=x.ndim, dtype=x.dtype)
    else:
        return np.min(x, axis=axis, keepdims=keepdims)


@implements(np.max, np.amax)
@derived_from(np)
def max(a, axis=None, keepdims=False, split_every=None, out=None):
    return reduction(
        a,
        chunk_max,
        chunk.max,
        combine=chunk_max,
        axis=axis,
        keepdims=keepdims,
        dtype=a.dtype,
        split_every=split_every,
        out=out,
    )


def chunk_max(x, axis=None, keepdims=None):
    """Version of np.max which ignores size 0 arrays"""
    if x.size == 0:
        return array_safe([], x, ndmin=x.ndim, dtype=x.dtype)
    else:
        return np.max(x, axis=axis, keepdims=keepdims)


@derived_from(np)
def any(a, axis=None, keepdims=False, split_every=None, out=None):
    return reduction(
        a,
        chunk.any,
        chunk.any,
        axis=axis,
        keepdims=keepdims,
        dtype="bool",
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def all(a, axis=None, keepdims=False, split_every=None, out=None):
    return reduction(
        a,
        chunk.all,
        chunk.all,
        axis=axis,
        keepdims=keepdims,
        dtype="bool",
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def nansum(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(chunk.nansum(np.ones((1,), dtype=a.dtype)), "dtype", object)
    return reduction(
        a,
        chunk.nansum,
        chunk.sum,
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def nanprod(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(chunk.nansum(np.ones((1,), dtype=a.dtype)), "dtype", object)
    return reduction(
        a,
        chunk.nanprod,
        chunk.prod,
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def nancumsum(x, axis, dtype=None, out=None, *, method="sequential"):
    """Dask added an additional keyword-only argument ``method``.

    method : {'sequential', 'blelloch'}, optional
        Choose which method to use to perform the cumsum.  Default is 'sequential'.

        * 'sequential' performs the cumsum of each prior block before the current block.
        * 'blelloch' is a work-efficient parallel cumsum.  It exposes parallelism by
            first taking the sum of each block and combines the sums via a binary tree.
            This method may be faster or more memory efficient depending on workload,
            scheduler, and hardware.  More benchmarking is necessary.
    """
    return cumreduction(
        chunk.nancumsum,
        operator.add,
        0,
        x,
        axis,
        dtype,
        out=out,
        method=method,
        preop=np.nansum,
    )


@derived_from(np)
def nancumprod(x, axis, dtype=None, out=None, *, method="sequential"):
    """Dask added an additional keyword-only argument ``method``.

    method : {'sequential', 'blelloch'}, optional
        Choose which method to use to perform the cumprod.  Default is 'sequential'.

        * 'sequential' performs the cumprod of each prior block before the current block.
        * 'blelloch' is a work-efficient parallel cumprod.  It exposes parallelism by first
            taking the product of each block and combines the products via a binary tree.
            This method may be faster or more memory efficient depending on workload,
            scheduler, and hardware.  More benchmarking is necessary.
    """
    return cumreduction(
        chunk.nancumprod,
        operator.mul,
        1,
        x,
        axis,
        dtype,
        out=out,
        method=method,
        preop=np.nanprod,
    )


@derived_from(np)
def nanmin(a, axis=None, keepdims=False, split_every=None, out=None):
    if np.isnan(a.size):
        raise ValueError(f"Arrays chunk sizes are unknown. {unknown_chunk_message}")
    if a.size == 0:
        raise ValueError(
            "zero-size array to reduction operation fmin which has no identity"
        )
    return reduction(
        a,
        _nanmin_skip,
        _nanmin_skip,
        axis=axis,
        keepdims=keepdims,
        dtype=a.dtype,
        split_every=split_every,
        out=out,
    )


def _nanmin_skip(x_chunk, axis, keepdims):
    if x_chunk.size > 0:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "All-NaN slice encountered", RuntimeWarning
            )
            return np.nanmin(x_chunk, axis=axis, keepdims=keepdims)
    else:
        return asarray_safe(
            np.array([], dtype=x_chunk.dtype), like=meta_from_array(x_chunk)
        )


@derived_from(np)
def nanmax(a, axis=None, keepdims=False, split_every=None, out=None):
    if np.isnan(a.size):
        raise ValueError(f"Arrays chunk sizes are unknown. {unknown_chunk_message}")
    if a.size == 0:
        raise ValueError(
            "zero-size array to reduction operation fmax which has no identity"
        )
    return reduction(
        a,
        _nanmax_skip,
        _nanmax_skip,
        axis=axis,
        keepdims=keepdims,
        dtype=a.dtype,
        split_every=split_every,
        out=out,
    )


def _nanmax_skip(x_chunk, axis, keepdims):
    if x_chunk.size > 0:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "All-NaN slice encountered", RuntimeWarning
            )
            return np.nanmax(x_chunk, axis=axis, keepdims=keepdims)
    else:
        return asarray_safe(
            np.array([], dtype=x_chunk.dtype), like=meta_from_array(x_chunk)
        )


def mean_chunk(
    x, sum=chunk.sum, numel=numel, dtype="f8", computing_meta=False, **kwargs
):
    if computing_meta:
        return x
    n = numel(x, dtype=dtype, **kwargs)

    total = sum(x, dtype=dtype, **kwargs)

    return {"n": n, "total": total}


def mean_combine(
    pairs,
    sum=chunk.sum,
    numel=numel,
    dtype="f8",
    axis=None,
    computing_meta=False,
    **kwargs,
):
    if not isinstance(pairs, list):
        pairs = [pairs]

    ns = deepmap(lambda pair: pair["n"], pairs) if not computing_meta else pairs
    n = _concatenate2(ns, axes=axis).sum(axis=axis, **kwargs)

    if computing_meta:
        return n

    totals = deepmap(lambda pair: pair["total"], pairs)
    total = _concatenate2(totals, axes=axis).sum(axis=axis, **kwargs)

    return {"n": n, "total": total}


def mean_agg(pairs, dtype="f8", axis=None, computing_meta=False, **kwargs):
    ns = deepmap(lambda pair: pair["n"], pairs) if not computing_meta else pairs
    n = _concatenate2(ns, axes=axis)
    n = np.sum(n, axis=axis, dtype=dtype, **kwargs)

    if computing_meta:
        return n

    totals = deepmap(lambda pair: pair["total"], pairs)
    total = _concatenate2(totals, axes=axis).sum(axis=axis, dtype=dtype, **kwargs)

    with np.errstate(divide="ignore", invalid="ignore"):
        return divide(total, n, dtype=dtype)


@derived_from(np)
def mean(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    elif a.dtype == object:
        dt = object
    else:
        dt = getattr(np.mean(np.zeros(shape=(1,), dtype=a.dtype)), "dtype", object)
    return reduction(
        a,
        mean_chunk,
        mean_agg,
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        combine=mean_combine,
        out=out,
        concatenate=False,
    )


@derived_from(np)
def nanmean(a, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(np.mean(np.ones(shape=(1,), dtype=a.dtype)), "dtype", object)
    return reduction(
        a,
        partial(mean_chunk, sum=chunk.nansum, numel=nannumel),
        mean_agg,
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        out=out,
        concatenate=False,
        combine=partial(mean_combine, sum=chunk.nansum, numel=nannumel),
    )


def moment_chunk(
    A,
    order=2,
    sum=chunk.sum,
    numel=numel,
    dtype="f8",
    computing_meta=False,
    implicit_complex_dtype=False,
    **kwargs,
):
    if computing_meta:
        return A
    n = numel(A, **kwargs)

    n = n.astype(np.int64)
    if implicit_complex_dtype:
        total = sum(A, **kwargs)
    else:
        total = sum(A, dtype=dtype, **kwargs)

    with np.errstate(divide="ignore", invalid="ignore"):
        u = total / n
    d = A - u
    if np.issubdtype(A.dtype, np.complexfloating):
        d = np.abs(d)
    xs = [sum(d**i, dtype=dtype, **kwargs) for i in range(2, order + 1)]
    M = np.stack(xs, axis=-1)
    return {"total": total, "n": n, "M": M}


def _moment_helper(Ms, ns, inner_term, order, sum, axis, kwargs):
    M = Ms[..., order - 2].sum(axis=axis, **kwargs) + sum(
        ns * inner_term**order, axis=axis, **kwargs
    )
    for k in range(1, order - 1):
        coeff = math.factorial(order) / (math.factorial(k) * math.factorial(order - k))
        M += coeff * sum(Ms[..., order - k - 2] * inner_term**k, axis=axis, **kwargs)
    return M


def moment_combine(
    pairs,
    order=2,
    ddof=0,
    dtype="f8",
    sum=np.sum,
    axis=None,
    computing_meta=False,
    **kwargs,
):
    if not isinstance(pairs, list):
        pairs = [pairs]

    kwargs["dtype"] = None
    kwargs["keepdims"] = True

    ns = deepmap(lambda pair: pair["n"], pairs) if not computing_meta else pairs
    ns = _concatenate2(ns, axes=axis)
    n = ns.sum(axis=axis, **kwargs)

    if computing_meta:
        return n

    totals = _concatenate2(deepmap(lambda pair: pair["total"], pairs), axes=axis)
    Ms = _concatenate2(deepmap(lambda pair: pair["M"], pairs), axes=axis)

    total = totals.sum(axis=axis, **kwargs)

    with np.errstate(divide="ignore", invalid="ignore"):
        if np.issubdtype(total.dtype, np.complexfloating):
            mu = divide(total, n)
            inner_term = np.abs(divide(totals, ns) - mu)
        else:
            mu = divide(total, n, dtype=dtype)
            inner_term = divide(totals, ns, dtype=dtype) - mu

    xs = [
        _moment_helper(Ms, ns, inner_term, o, sum, axis, kwargs)
        for o in range(2, order + 1)
    ]
    M = np.stack(xs, axis=-1)
    return {"total": total, "n": n, "M": M}


def moment_agg(
    pairs,
    order=2,
    ddof=0,
    dtype="f8",
    sum=np.sum,
    axis=None,
    computing_meta=False,
    **kwargs,
):
    if not isinstance(pairs, list):
        pairs = [pairs]

    kwargs["dtype"] = dtype
    # To properly handle ndarrays, the original dimensions need to be kept for
    # part of the calculation.
    keepdim_kw = kwargs.copy()
    keepdim_kw["keepdims"] = True
    keepdim_kw["dtype"] = None

    ns = deepmap(lambda pair: pair["n"], pairs) if not computing_meta else pairs
    ns = _concatenate2(ns, axes=axis)
    n = ns.sum(axis=axis, **keepdim_kw)

    if computing_meta:
        return n

    totals = _concatenate2(deepmap(lambda pair: pair["total"], pairs), axes=axis)
    Ms = _concatenate2(deepmap(lambda pair: pair["M"], pairs), axes=axis)

    mu = divide(totals.sum(axis=axis, **keepdim_kw), n)

    with np.errstate(divide="ignore", invalid="ignore"):
        if np.issubdtype(totals.dtype, np.complexfloating):
            inner_term = np.abs(divide(totals, ns) - mu)
        else:
            inner_term = divide(totals, ns, dtype=dtype) - mu

    M = _moment_helper(Ms, ns, inner_term, order, sum, axis, kwargs)

    denominator = n.sum(axis=axis, **kwargs) - ddof

    # taking care of the edge case with empty or all-nans array with ddof > 0
    if isinstance(denominator, Number):
        if denominator < 0:
            denominator = np.nan
    elif denominator is not np.ma.masked:
        denominator[denominator < 0] = np.nan

    return divide(M, denominator, dtype=dtype)


def moment(
    a, order, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None
):
    """Calculate the nth centralized moment.

    Parameters
    ----------
    a : Array
        Data over which to compute moment
    order : int
        Order of the moment that is returned, must be >= 2.
    axis : int, optional
        Axis along which the central moment is computed. The default is to
        compute the moment of the flattened array.
    dtype : data-type, optional
        Type to use in computing the moment. For arrays of integer type the
        default is float64; for arrays of float types it is the same as the
        array type.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original array.
    ddof : int, optional
        "Delta Degrees of Freedom": the divisor used in the calculation is
        N - ddof, where N represents the number of elements. By default
        ddof is zero.

    Returns
    -------
    moment : Array

    References
    ----------
    .. [1] Pebay, Philippe (2008), "Formulas for Robust, One-Pass Parallel
        Computation of Covariances and Arbitrary-Order Statistical Moments",
        Technical Report SAND2008-6212, Sandia National Laboratories.

    """
    if not isinstance(order, Integral) or order < 0:
        raise ValueError("Order must be an integer >= 0")

    if order < 2:
        reduced = a.sum(axis=axis)  # get reduced shape and chunks
        if order == 0:
            # When order equals 0, the result is 1, by definition.
            return ones(
                reduced.shape, chunks=reduced.chunks, dtype="f8", meta=reduced._meta
            )
        # By definition the first order about the mean is 0.
        return zeros(
            reduced.shape, chunks=reduced.chunks, dtype="f8", meta=reduced._meta
        )

    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(np.var(np.ones(shape=(1,), dtype=a.dtype)), "dtype", object)

    implicit_complex_dtype = dtype is None and np.iscomplexobj(a)

    return reduction(
        a,
        partial(
            moment_chunk, order=order, implicit_complex_dtype=implicit_complex_dtype
        ),
        partial(moment_agg, order=order, ddof=ddof),
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        out=out,
        concatenate=False,
        combine=partial(moment_combine, order=order),
    )


@derived_from(np)
def var(a, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(np.var(np.ones(shape=(1,), dtype=a.dtype)), "dtype", object)

    implicit_complex_dtype = dtype is None and np.iscomplexobj(a)

    return reduction(
        a,
        partial(moment_chunk, implicit_complex_dtype=implicit_complex_dtype),
        partial(moment_agg, ddof=ddof),
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        combine=moment_combine,
        name="var",
        out=out,
        concatenate=False,
    )


@derived_from(np)
def nanvar(
    a, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None
):
    if dtype is not None:
        dt = dtype
    else:
        dt = getattr(np.var(np.ones(shape=(1,), dtype=a.dtype)), "dtype", object)

    implicit_complex_dtype = dtype is None and np.iscomplexobj(a)

    return reduction(
        a,
        partial(
            moment_chunk,
            sum=chunk.nansum,
            numel=nannumel,
            implicit_complex_dtype=implicit_complex_dtype,
        ),
        partial(moment_agg, sum=np.nansum, ddof=ddof),
        axis=axis,
        keepdims=keepdims,
        dtype=dt,
        split_every=split_every,
        combine=partial(moment_combine, sum=np.nansum),
        out=out,
        concatenate=False,
    )


def _sqrt(a):
    if isinstance(a, np.ma.masked_array) and not a.shape and a.mask.all():
        return np.ma.masked
    return np.sqrt(a)


def safe_sqrt(a):
    """A version of sqrt that properly handles scalar masked arrays.

    To mimic ``np.ma`` reductions, we need to convert scalar masked arrays that
    have an active mask to the ``np.ma.masked`` singleton. This is properly
    handled automatically for reduction code, but not for ufuncs. We implement
    a simple version here, since calling `np.ma.sqrt` everywhere is
    significantly more expensive.
    """
    if hasattr(a, "_elemwise"):
        return a._elemwise(_sqrt, a)
    return _sqrt(a)


@derived_from(np)
def std(a, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None):
    result = safe_sqrt(
        var(
            a,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            ddof=ddof,
            split_every=split_every,
            out=out,
        )
    )
    if dtype and dtype != result.dtype:
        result = result.astype(dtype)
    return result


@derived_from(np)
def nanstd(
    a, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None
):
    result = safe_sqrt(
        nanvar(
            a,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            ddof=ddof,
            split_every=split_every,
            out=out,
        )
    )
    if dtype and dtype != result.dtype:
        result = result.astype(dtype)
    return result


def _arg_combine(data, axis, argfunc, keepdims=False):
    """Merge intermediate results from ``arg_*`` functions"""
    if isinstance(data, dict):
        # Array type doesn't support structured arrays (e.g., CuPy),
        # therefore `data` is stored in a `dict`.
        assert data["vals"].ndim == data["arg"].ndim
        axis = (
            None
            if len(axis) == data["vals"].ndim or data["vals"].ndim == 1
            else axis[0]
        )
    else:
        axis = None if len(axis) == data.ndim or data.ndim == 1 else axis[0]

    vals = data["vals"]
    arg = data["arg"]
    if axis is None:
        local_args = argfunc(vals, axis=axis, keepdims=keepdims)
        vals = vals.ravel()[local_args]
        arg = arg.ravel()[local_args]
    else:
        local_args = argfunc(vals, axis=axis)
        inds = list(np.ogrid[tuple(map(slice, local_args.shape))])
        inds.insert(axis, local_args)
        inds = tuple(inds)
        vals = vals[inds]
        arg = arg[inds]
        if keepdims:
            vals = np.expand_dims(vals, axis)
            arg = np.expand_dims(arg, axis)
    return arg, vals


def arg_chunk(func, argfunc, x, axis, offset_info):
    arg_axis = None if len(axis) == x.ndim or x.ndim == 1 else axis[0]
    vals = func(x, axis=arg_axis, keepdims=True)
    arg = argfunc(x, axis=arg_axis, keepdims=True)
    if x.ndim > 0:
        if arg_axis is None:
            offset, total_shape = offset_info
            ind = np.unravel_index(arg.ravel()[0], x.shape)
            total_ind = tuple(o + i for (o, i) in zip(offset, ind))
            arg[:] = np.ravel_multi_index(total_ind, total_shape)
        else:
            arg += offset_info

    if isinstance(vals, np.ma.masked_array):
        if "min" in argfunc.__name__:
            fill_value = np.ma.minimum_fill_value(vals)
        else:
            fill_value = np.ma.maximum_fill_value(vals)
        vals = np.ma.filled(vals, fill_value)

    try:
        result = np.empty_like(
            vals, shape=vals.shape, dtype=[("vals", vals.dtype), ("arg", arg.dtype)]
        )
    except TypeError:
        # Array type doesn't support structured arrays (e.g., CuPy)
        result = dict()

    result["vals"] = vals
    result["arg"] = arg
    return result


def arg_combine(argfunc, data, axis=None, **kwargs):
    arg, vals = _arg_combine(data, axis, argfunc, keepdims=True)

    try:
        result = np.empty_like(
            vals, shape=vals.shape, dtype=[("vals", vals.dtype), ("arg", arg.dtype)]
        )
    except TypeError:
        # Array type doesn't support structured arrays (e.g., CuPy).
        result = dict()

    result["vals"] = vals
    result["arg"] = arg
    return result


def arg_agg(argfunc, data, axis=None, keepdims=False, **kwargs):
    return _arg_combine(data, axis, argfunc, keepdims=keepdims)[0]


def nanarg_agg(argfunc, data, axis=None, keepdims=False, **kwargs):
    arg, vals = _arg_combine(data, axis, argfunc, keepdims=keepdims)
    if np.any(np.isnan(vals)):
        raise ValueError("All NaN slice encountered")
    return arg


def arg_reduction(
    x, chunk, combine, agg, axis=None, keepdims=False, split_every=None, out=None
):
    """Generic function for argreduction.

    Parameters
    ----------
    x : Array
    chunk : callable
        Partialed ``arg_chunk``.
    combine : callable
        Partialed ``arg_combine``.
    agg : callable
        Partialed ``arg_agg``.
    axis : int, optional
    split_every : int or dict, optional
    """
    if axis is None:
        axis = tuple(range(x.ndim))
        ravel = True
    elif isinstance(axis, Integral):
        axis = validate_axis(axis, x.ndim)
        axis = (axis,)
        ravel = x.ndim == 1
    else:
        raise TypeError(f"axis must be either `None` or int, got '{axis}'")

    for ax in axis:
        chunks = x.chunks[ax]
        if len(chunks) > 1 and np.isnan(chunks).any():
            raise ValueError(
                "Arg-reductions do not work with arrays that have "
                "unknown chunksizes. At some point in your computation "
                "this array lost chunking information.\n\n"
                "A possible solution is with \n"
                "  x.compute_chunk_sizes()"
            )

    # Map chunk across all blocks
    name = f"arg-reduce-{tokenize(axis, x, chunk, combine, split_every)}"
    old = x.name
    keys = list(product(*map(range, x.numblocks)))
    offsets = list(product(*(accumulate(operator.add, bd[:-1], 0) for bd in x.chunks)))
    if ravel:
        offset_info = zip(offsets, repeat(x.shape))
    else:
        offset_info = pluck(axis[0], offsets)

    chunks = tuple((1,) * len(c) if i in axis else c for (i, c) in enumerate(x.chunks))
    dsk = {
        (name,) + k: (chunk, (old,) + k, axis, off)
        for (k, off) in zip(keys, offset_info)
    }

    dtype = np.argmin(asarray_safe([1], like=meta_from_array(x)))
    meta = None
    if is_arraylike(dtype):
        # This case occurs on non-NumPy types (e.g., CuPy), where the returned
        # value is an ndarray rather than a scalar.
        meta = dtype
        dtype = meta.dtype

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[x])
    tmp = Array(graph, name, chunks, dtype=dtype, meta=meta)

    result = _tree_reduce(
        tmp,
        agg,
        axis,
        keepdims=keepdims,
        dtype=dtype,
        split_every=split_every,
        combine=combine,
    )
    return handle_out(out, result)


def _nanargmin(x, axis, **kwargs):
    try:
        return chunk.nanargmin(x, axis, **kwargs)
    except ValueError:
        return chunk.nanargmin(np.where(np.isnan(x), np.inf, x), axis, **kwargs)


def _nanargmax(x, axis, **kwargs):
    try:
        return chunk.nanargmax(x, axis, **kwargs)
    except ValueError:
        return chunk.nanargmax(np.where(np.isnan(x), -np.inf, x), axis, **kwargs)


@derived_from(np)
def argmax(a, axis=None, keepdims=False, split_every=None, out=None):
    return arg_reduction(
        a,
        partial(arg_chunk, chunk.max, chunk.argmax),
        partial(arg_combine, chunk.argmax),
        partial(arg_agg, chunk.argmax),
        axis=axis,
        keepdims=keepdims,
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def argmin(a, axis=None, keepdims=False, split_every=None, out=None):
    return arg_reduction(
        a,
        partial(arg_chunk, chunk.min, chunk.argmin),
        partial(arg_combine, chunk.argmin),
        partial(arg_agg, chunk.argmin),
        axis=axis,
        keepdims=keepdims,
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def nanargmax(a, axis=None, keepdims=False, split_every=None, out=None):
    return arg_reduction(
        a,
        partial(arg_chunk, chunk.nanmax, _nanargmax),
        partial(arg_combine, _nanargmax),
        partial(nanarg_agg, _nanargmax),
        axis=axis,
        keepdims=keepdims,
        split_every=split_every,
        out=out,
    )


@derived_from(np)
def nanargmin(a, axis=None, keepdims=False, split_every=None, out=None):
    return arg_reduction(
        a,
        partial(arg_chunk, chunk.nanmin, _nanargmin),
        partial(arg_combine, _nanargmin),
        partial(nanarg_agg, _nanargmin),
        axis=axis,
        keepdims=keepdims,
        split_every=split_every,
        out=out,
    )


def _prefixscan_combine(func, binop, pre, x, axis, dtype):
    """Combine results of a parallel prefix scan such as cumsum

    Parameters
    ----------
    func : callable
        Cumulative function (e.g. ``np.cumsum``)
    binop : callable
        Associative function (e.g. ``add``)
    pre : np.array
        The value calculated in parallel from ``preop``.
        For example, the sum of all the previous blocks.
    x : np.array
        Current block
    axis : int
    dtype : dtype

    Returns
    -------
    np.array
    """
    # We could compute this in two tasks.
    # This would allow us to do useful work (i.e., func), while waiting on `pre`.
    # Using one task may guide the scheduler to do better and reduce scheduling overhead.
    return binop(pre, func(x, axis=axis, dtype=dtype))


def _prefixscan_first(func, x, axis, dtype):
    """Compute the prefix scan (e.g., cumsum) on the first block

    Parameters
    ----------
    func : callable
        Cumulative function (e.g. ``np.cumsum``)
    x : np.array
        Current block
    axis : int
    dtype : dtype

    Returns
    -------
    np.array
    """
    return func(x, axis=axis, dtype=dtype)


def prefixscan_blelloch(func, preop, binop, x, axis=None, dtype=None, out=None):
    """Generic function to perform parallel cumulative scan (a.k.a prefix scan)

    The Blelloch prefix scan is work-efficient and exposes parallelism.
    A parallel cumsum works by first taking the sum of each block, then do a binary tree
    merge followed by a fan-out (i.e., the Brent-Kung pattern).  We then take the cumsum
    of each block and add the sum of the previous blocks.

    When performing a cumsum across N chunks, this method has 2 * lg(N) levels of dependencies.
    In contrast, the sequential method has N levels of dependencies.

    Floating point operations should be more accurate with this method compared to sequential.

    Parameters
    ----------
    func : callable
        Cumulative function (e.g. ``np.cumsum``)
    preop : callable
        Function to get the final value of a cumulative function (e.g., ``np.sum``)
    binop : callable
        Associative function (e.g. ``add``)
    x : dask array
    axis : int
    dtype : dtype

    Returns
    -------
    dask array
    """
    if axis is None:
        x = x.flatten().rechunk(chunks=x.npartitions)
        axis = 0
    if dtype is None:
        dtype = getattr(func(np.ones((0,), dtype=x.dtype)), "dtype", object)
    assert isinstance(axis, Integral)
    axis = validate_axis(axis, x.ndim)
    name = f"{func.__name__}-{tokenize(func, axis, preop, binop, x, dtype)}"
    base_key = (name,)

    # Right now, the metadata for batches is incorrect, but this should be okay
    batches = x.map_blocks(preop, axis=axis, keepdims=True, dtype=dtype)
    # We don't need the last index until the end
    *indices, last_index = full_indices = [
        list(
            product(
                *[range(nb) if j != axis else [i] for j, nb in enumerate(x.numblocks)]
            )
        )
        for i in range(x.numblocks[axis])
    ]
    prefix_vals = [[(batches.name,) + index for index in vals] for vals in indices]
    dsk = {}
    n_vals = len(prefix_vals)
    level = 0
    if n_vals >= 2:
        # Upsweep
        stride = 1
        stride2 = 2
        while stride2 <= n_vals:
            for i in range(stride2 - 1, n_vals, stride2):
                new_vals = []
                for index, left_val, right_val in zip(
                    indices[i], prefix_vals[i - stride], prefix_vals[i]
                ):
                    key = base_key + index + (level, i)
                    dsk[key] = (binop, left_val, right_val)
                    new_vals.append(key)
                prefix_vals[i] = new_vals
            stride = stride2
            stride2 *= 2
            level += 1

        # Downsweep
        # With `n_vals == 3`, we would have `stride = 1` and `stride = 0`, but we need
        # to do a downsweep iteration, so make sure stride2 is at least 2.
        stride2 = builtins.max(2, 2 ** math.ceil(math.log2(n_vals // 2)))
        stride = stride2 // 2
        while stride > 0:
            for i in range(stride2 + stride - 1, n_vals, stride2):
                new_vals = []
                for index, left_val, right_val in zip(
                    indices[i], prefix_vals[i - stride], prefix_vals[i]
                ):
                    key = base_key + index + (level, i)
                    dsk[key] = (binop, left_val, right_val)
                    new_vals.append(key)
                prefix_vals[i] = new_vals
            stride2 = stride
            stride //= 2
            level += 1

    if full_indices:
        for index in full_indices[0]:
            dsk[base_key + index] = (
                _prefixscan_first,
                func,
                (x.name,) + index,
                axis,
                dtype,
            )
        for indexes, vals in zip(drop(1, full_indices), prefix_vals):
            for index, val in zip(indexes, vals):
                dsk[base_key + index] = (
                    _prefixscan_combine,
                    func,
                    binop,
                    val,
                    (x.name,) + index,
                    axis,
                    dtype,
                )
    if len(full_indices) < 2:
        deps = [x]
    else:
        deps = [x, batches]
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=deps)
    result = Array(graph, name, x.chunks, batches.dtype)
    return handle_out(out, result)


def cumreduction(
    func,
    binop,
    ident,
    x,
    axis=None,
    dtype=None,
    out=None,
    method="sequential",
    preop=None,
):
    """Generic function for cumulative reduction

    Parameters
    ----------
    func: callable
        Cumulative function like np.cumsum or np.cumprod
    binop: callable
        Associated binary operator like ``np.cumsum->add`` or ``np.cumprod->mul``
    ident: Number
        Associated identity like ``np.cumsum->0`` or ``np.cumprod->1``
    x: dask Array
    axis: int
    dtype: dtype
    method : {'sequential', 'blelloch'}, optional
        Choose which method to use to perform the cumsum.  Default is 'sequential'.

        * 'sequential' performs the scan of each prior block before the current block.
        * 'blelloch' is a work-efficient parallel scan.  It exposes parallelism by first
          calling ``preop`` on each block and combines the values via a binary tree.
          This method may be faster or more memory efficient depending on workload,
          scheduler, and hardware.  More benchmarking is necessary.
    preop: callable, optional
        Function used by 'blelloch' method,
        like ``np.cumsum->np.sum`` or ``np.cumprod->np.prod``

    Returns
    -------
    dask array

    See also
    --------
    cumsum
    cumprod
    """
    if method == "blelloch":
        if preop is None:
            raise TypeError(
                'cumreduction with "blelloch" method required `preop=` argument'
            )
        return prefixscan_blelloch(func, preop, binop, x, axis, dtype, out=out)
    elif method != "sequential":
        raise ValueError(
            f'Invalid method for cumreduction.  Expected "sequential" or "blelloch".  Got: {method!r}'
        )

    if axis is None:
        x = x.flatten().rechunk(chunks=x.npartitions)
        axis = 0
    if dtype is None:
        dtype = getattr(func(np.ones((0,), dtype=x.dtype)), "dtype", object)
    assert isinstance(axis, Integral)
    axis = validate_axis(axis, x.ndim)

    m = x.map_blocks(func, axis=axis, dtype=dtype)

    name = f"{func.__name__}-{tokenize(func, axis, binop, ident, x, dtype)}"
    n = x.numblocks[axis]
    full = slice(None, None, None)
    slc = (full,) * axis + (slice(-1, None),) + (full,) * (x.ndim - axis - 1)

    indices = list(
        product(*[range(nb) if i != axis else [0] for i, nb in enumerate(x.numblocks)])
    )
    dsk = dict()
    for ind in indices:
        shape = tuple(x.chunks[i][ii] if i != axis else 1 for i, ii in enumerate(ind))
        dsk[(name, "extra") + ind] = (
            apply,
            np.full_like,
            (x._meta, ident, m.dtype),
            {"shape": shape},
        )
        dsk[(name,) + ind] = (m.name,) + ind

    for i in range(1, n):
        last_indices = indices
        indices = list(
            product(
                *[range(nb) if ii != axis else [i] for ii, nb in enumerate(x.numblocks)]
            )
        )
        for old, ind in zip(last_indices, indices):
            this_slice = (name, "extra") + ind
            dsk[this_slice] = (
                binop,
                (name, "extra") + old,
                (operator.getitem, (m.name,) + old, slc),
            )
            dsk[(name,) + ind] = (binop, this_slice, (m.name,) + ind)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[m])
    result = Array(graph, name, x.chunks, m.dtype, meta=x._meta)
    return handle_out(out, result)


def _cumsum_merge(a, b):
    if isinstance(a, np.ma.masked_array) or isinstance(b, np.ma.masked_array):
        values = np.ma.getdata(a) + np.ma.getdata(b)
        return np.ma.masked_array(values, mask=np.ma.getmaskarray(b))
    return a + b


def _cumprod_merge(a, b):
    if isinstance(a, np.ma.masked_array) or isinstance(b, np.ma.masked_array):
        values = np.ma.getdata(a) * np.ma.getdata(b)
        return np.ma.masked_array(values, mask=np.ma.getmaskarray(b))
    return a * b


@derived_from(np)
def cumsum(x, axis=None, dtype=None, out=None, method="sequential"):
    """Dask added an additional keyword-only argument ``method``.

    method : {'sequential', 'blelloch'}, optional
        Choose which method to use to perform the cumsum.  Default is 'sequential'.

        * 'sequential' performs the cumsum of each prior block before the current block.
        * 'blelloch' is a work-efficient parallel cumsum.  It exposes parallelism by
          first taking the sum of each block and combines the sums via a binary tree.
          This method may be faster or more memory efficient depending on workload,
          scheduler, and hardware.  More benchmarking is necessary.
    """
    return cumreduction(
        np.cumsum,
        _cumsum_merge,
        0,
        x,
        axis,
        dtype,
        out=out,
        method=method,
        preop=np.sum,
    )


@derived_from(np)
def cumprod(x, axis=None, dtype=None, out=None, method="sequential"):
    """Dask added an additional keyword-only argument ``method``.

    method : {'sequential', 'blelloch'}, optional
        Choose which method to use to perform the cumprod.  Default is 'sequential'.

        * 'sequential' performs the cumprod of each prior block before the current block.
        * 'blelloch' is a work-efficient parallel cumprod.  It exposes parallelism by first
          taking the product of each block and combines the products via a binary tree.
          This method may be faster or more memory efficient depending on workload,
          scheduler, and hardware.  More benchmarking is necessary.
    """
    return cumreduction(
        np.cumprod,
        _cumprod_merge,
        1,
        x,
        axis,
        dtype,
        out=out,
        method=method,
        preop=np.prod,
    )


def topk(a, k, axis=-1, split_every=None):
    """Extract the k largest elements from a on the given axis,
    and return them sorted from largest to smallest.
    If k is negative, extract the -k smallest elements instead,
    and return them sorted from smallest to largest.

    This performs best when ``k`` is much smaller than the chunk size. All
    results will be returned in a single chunk along the given axis.

    Parameters
    ----------
    x: Array
        Data being sorted
    k: int
    axis: int, optional
    split_every: int >=2, optional
        See :func:`reduce`. This parameter becomes very important when k is
        on the same order of magnitude of the chunk size or more, as it
        prevents getting the whole or a significant portion of the input array
        in memory all at once, with a negative impact on network transfer
        too when running on distributed.

    Returns
    -------
    Selection of x with size abs(k) along the given axis.

    Examples
    --------
    >>> import dask.array as da
    >>> x = np.array([5, 1, 3, 6])
    >>> d = da.from_array(x, chunks=2)
    >>> d.topk(2).compute()
    array([6, 5])
    >>> d.topk(-2).compute()
    array([1, 3])
    """
    axis = validate_axis(axis, a.ndim)

    # chunk and combine steps of the reduction, which recursively invoke
    # np.partition to pick the top/bottom k elements from the previous step.
    # The selection is not sorted internally.
    chunk_combine = partial(chunk.topk, k=k)
    # aggregate step of the reduction. Internally invokes the chunk/combine
    # function, then sorts the results internally.
    aggregate = partial(chunk.topk_aggregate, k=k)

    return reduction(
        a,
        chunk=chunk_combine,
        combine=chunk_combine,
        aggregate=aggregate,
        axis=axis,
        keepdims=True,
        dtype=a.dtype,
        split_every=split_every,
        output_size=abs(k),
    )


def argtopk(a, k, axis=-1, split_every=None):
    """Extract the indices of the k largest elements from a on the given axis,
    and return them sorted from largest to smallest. If k is negative, extract
    the indices of the -k smallest elements instead, and return them sorted
    from smallest to largest.

    This performs best when ``k`` is much smaller than the chunk size. All
    results will be returned in a single chunk along the given axis.

    Parameters
    ----------
    x: Array
        Data being sorted
    k: int
    axis: int, optional
    split_every: int >=2, optional
        See :func:`topk`. The performance considerations for topk also apply
        here.

    Returns
    -------
    Selection of np.intp indices of x with size abs(k) along the given axis.

    Examples
    --------
    >>> import dask.array as da
    >>> x = np.array([5, 1, 3, 6])
    >>> d = da.from_array(x, chunks=2)
    >>> d.argtopk(2).compute()
    array([3, 0])
    >>> d.argtopk(-2).compute()
    array([1, 2])
    """
    axis = validate_axis(axis, a.ndim)

    # Generate nodes where every chunk is a tuple of (a, original index of a)
    idx = arange(a.shape[axis], chunks=(a.chunks[axis],), dtype=np.intp)
    idx = idx[tuple(slice(None) if i == axis else np.newaxis for i in range(a.ndim))]
    a_plus_idx = a.map_blocks(chunk.argtopk_preprocess, idx, dtype=object)

    # chunk and combine steps of the reduction. They acquire in input a tuple
    # of (a, original indices of a) and return another tuple containing the top
    # k elements of a and the matching original indices. The selection is not
    # sorted internally, as in np.argpartition.
    chunk_combine = partial(chunk.argtopk, k=k)
    # aggregate step of the reduction. Internally invokes the chunk/combine
    # function, then sorts the results internally, drops a and returns the
    # index only.
    aggregate = partial(chunk.argtopk_aggregate, k=k)

    if isinstance(axis, Number):
        naxis = 1
    else:
        naxis = len(axis)

    meta = a._meta.astype(np.intp).reshape((0,) * (a.ndim - naxis + 1))

    return reduction(
        a_plus_idx,
        chunk=chunk_combine,
        combine=chunk_combine,
        aggregate=aggregate,
        axis=axis,
        keepdims=True,
        dtype=np.intp,
        split_every=split_every,
        concatenate=False,
        output_size=abs(k),
        meta=meta,
    )


@derived_from(np)
def trace(a, offset=0, axis1=0, axis2=1, dtype=None):
    return diagonal(a, offset=offset, axis1=axis1, axis2=axis2).sum(-1, dtype=dtype)


@derived_from(np)
def median(a, axis=None, keepdims=False, out=None):
    """
    This works by automatically chunking the reduced axes to a single chunk if necessary
    and then calling ``numpy.median`` function across the remaining dimensions
    """
    if axis is None:
        raise NotImplementedError(
            "The da.median function only works along an axis.  "
            "The full algorithm is difficult to do in parallel"
        )

    if not isinstance(axis, Iterable):
        axis = (axis,)

    axis = [ax + a.ndim if ax < 0 else ax for ax in axis]

    # rechunk if reduced axes are not contained in a single chunk
    if builtins.any(a.numblocks[ax] > 1 for ax in axis):
        a = a.rechunk({ax: -1 if ax in axis else "auto" for ax in range(a.ndim)})

    result = a.map_blocks(
        np.median,
        axis=axis,
        keepdims=keepdims,
        drop_axis=axis if not keepdims else None,
        chunks=[1 if ax in axis else c for ax, c in enumerate(a.chunks)]
        if keepdims
        else None,
    )

    result = handle_out(out, result)
    return result


@derived_from(np)
def nanmedian(a, axis=None, keepdims=False, out=None):
    """
    This works by automatically chunking the reduced axes to a single chunk
    and then calling ``numpy.nanmedian`` function across the remaining dimensions
    """
    if axis is None:
        raise NotImplementedError(
            "The da.nanmedian function only works along an axis or a subset of axes.  "
            "The full algorithm is difficult to do in parallel"
        )

    if not isinstance(axis, Iterable):
        axis = (axis,)

    axis = [ax + a.ndim if ax < 0 else ax for ax in axis]

    # rechunk if reduced axes are not contained in a single chunk
    if builtins.any(a.numblocks[ax] > 1 for ax in axis):
        a = a.rechunk({ax: -1 if ax in axis else "auto" for ax in range(a.ndim)})

    result = a.map_blocks(
        np.nanmedian,
        axis=axis,
        keepdims=keepdims,
        drop_axis=axis if not keepdims else None,
        chunks=[1 if ax in axis else c for ax, c in enumerate(a.chunks)]
        if keepdims
        else None,
    )

    result = handle_out(out, result)
    return result
