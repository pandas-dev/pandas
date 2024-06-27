from __future__ import annotations

import math
import warnings
from collections import Counter
from functools import reduce
from itertools import product
from operator import mul

import numpy as np

from dask import config
from dask.array.core import Array, normalize_chunks
from dask.array.utils import meta_from_array
from dask.base import tokenize
from dask.core import flatten
from dask.highlevelgraph import HighLevelGraph
from dask.utils import M, parse_bytes

_not_implemented_message = """
Dask's reshape only supports operations that merge or split existing dimensions
evenly. For example:

>>> x = da.ones((6, 5, 4), chunks=(3, 2, 2))
>>> x.reshape((3, 2, 5, 4))  # supported, splits 6 into 3 & 2
>>> x.reshape((30, 4))       # supported, merges 6 & 5 into 30
>>> x.reshape((4, 5, 6))     # unsupported, existing dimensions split unevenly

To work around this you may call reshape in multiple passes, or (if your data
is small enough) call ``compute`` first and handle reshaping in ``numpy``
directly.
"""


def reshape_rechunk(inshape, outshape, inchunks):
    assert all(isinstance(c, tuple) for c in inchunks)
    ii = len(inshape) - 1
    oi = len(outshape) - 1
    result_inchunks = [None for i in range(len(inshape))]
    result_outchunks = [None for i in range(len(outshape))]

    while ii >= 0 or oi >= 0:
        if inshape[ii] == outshape[oi]:
            result_inchunks[ii] = inchunks[ii]
            result_outchunks[oi] = inchunks[ii]
            ii -= 1
            oi -= 1
            continue
        din = inshape[ii]
        dout = outshape[oi]
        if din == 1:
            result_inchunks[ii] = (1,)
            ii -= 1
        elif dout == 1:
            result_outchunks[oi] = (1,)
            oi -= 1
        elif din < dout:  # (4, 4, 4) -> (64,)
            ileft = ii - 1
            while (
                ileft >= 0 and reduce(mul, inshape[ileft : ii + 1]) < dout
            ):  # 4 < 64, 4*4 < 64, 4*4*4 == 64
                ileft -= 1
            if reduce(mul, inshape[ileft : ii + 1]) != dout:
                raise NotImplementedError(_not_implemented_message)
            # Special case to avoid intermediate rechunking:
            # When all the lower axis are completely chunked (chunksize=1) then
            # we're simply moving around blocks.
            if all(len(inchunks[i]) == inshape[i] for i in range(ii)):
                for i in range(ii + 1):
                    result_inchunks[i] = inchunks[i]
                result_outchunks[oi] = inchunks[ii] * math.prod(
                    map(len, inchunks[ileft:ii])
                )
            else:
                for i in range(ileft + 1, ii + 1):  # need single-shape dimensions
                    result_inchunks[i] = (inshape[i],)  # chunks[i] = (4,)

                chunk_reduction = reduce(mul, map(len, inchunks[ileft + 1 : ii + 1]))
                result_inchunks[ileft] = expand_tuple(inchunks[ileft], chunk_reduction)

                prod = reduce(mul, inshape[ileft + 1 : ii + 1])  # 16
                result_outchunks[oi] = tuple(
                    prod * c for c in result_inchunks[ileft]
                )  # (1, 1, 1, 1) .* 16

            oi -= 1
            ii = ileft - 1
        elif din > dout:  # (64,) -> (4, 4, 4)
            oleft = oi - 1
            while oleft >= 0 and reduce(mul, outshape[oleft : oi + 1]) < din:
                oleft -= 1
            if reduce(mul, outshape[oleft : oi + 1]) != din:
                raise NotImplementedError(_not_implemented_message)
            # TODO: don't coalesce shapes unnecessarily
            cs = reduce(mul, outshape[oleft + 1 : oi + 1])

            result_inchunks[ii] = contract_tuple(inchunks[ii], cs)  # (16, 16, 16, 16)

            for i in range(oleft + 1, oi + 1):
                result_outchunks[i] = (outshape[i],)

            result_outchunks[oleft] = tuple(c // cs for c in result_inchunks[ii])

            oi = oleft - 1
            ii -= 1

    return tuple(result_inchunks), tuple(result_outchunks)


def expand_tuple(chunks, factor):
    """

    >>> expand_tuple((2, 4), 2)
    (1, 1, 2, 2)

    >>> expand_tuple((2, 4), 3)
    (1, 1, 1, 1, 2)

    >>> expand_tuple((3, 4), 2)
    (1, 2, 2, 2)

    >>> expand_tuple((7, 4), 3)
    (2, 2, 3, 1, 1, 2)
    """
    if factor == 1:
        return chunks

    out = []
    for c in chunks:
        x = c
        part = max(x / factor, 1)
        while x >= 2 * part:
            out.append(int(part))
            x -= int(part)
        if x:
            out.append(x)
    assert sum(chunks) == sum(out)
    return tuple(out)


def contract_tuple(chunks, factor):
    """Return simple chunks tuple such that factor divides all elements

    Examples
    --------

    >>> contract_tuple((2, 2, 8, 4), 4)
    (4, 8, 4)
    """
    assert sum(chunks) % factor == 0

    out = []
    residual = 0
    for chunk in chunks:
        chunk += residual
        div = chunk // factor
        residual = chunk % factor
        good = factor * div
        if good:
            out.append(good)
    return tuple(out)


def reshape(x, shape, merge_chunks=True, limit=None):
    """Reshape array to new shape

    Parameters
    ----------
    shape : int or tuple of ints
        The new shape should be compatible with the original shape. If
        an integer, then the result will be a 1-D array of that length.
        One shape dimension can be -1. In this case, the value is
        inferred from the length of the array and remaining dimensions.
    merge_chunks : bool, default True
        Whether to merge chunks using the logic in :meth:`dask.array.rechunk`
        when communication is necessary given the input array chunking and
        the output shape. With ``merge_chunks==False``, the input array will
        be rechunked to a chunksize of 1, which can create very many tasks.
    limit: int (optional)
        The maximum block size to target in bytes. If no limit is provided,
        it defaults to using the ``array.chunk-size`` Dask config value.

    Notes
    -----
    This is a parallelized version of the ``np.reshape`` function with the
    following limitations:

    1.  It assumes that the array is stored in `row-major order`_
    2.  It only allows for reshapings that collapse or merge dimensions like
        ``(1, 2, 3, 4) -> (1, 6, 4)`` or ``(64,) -> (4, 4, 4)``

    .. _`row-major order`: https://en.wikipedia.org/wiki/Row-_and_column-major_order

    When communication is necessary this algorithm depends on the logic within
    rechunk.  It endeavors to keep chunk sizes roughly the same when possible.

    See :ref:`array-chunks.reshaping` for a discussion the tradeoffs of
    ``merge_chunks``.

    See Also
    --------
    dask.array.rechunk
    numpy.reshape
    """
    # Sanitize inputs, look for -1 in shape
    from dask.array.core import PerformanceWarning
    from dask.array.slicing import sanitize_index

    shape = tuple(map(sanitize_index, shape))
    known_sizes = [s for s in shape if s != -1]
    if len(known_sizes) < len(shape):
        if len(shape) - len(known_sizes) > 1:
            raise ValueError("can only specify one unknown dimension")
        # Fastpath for x.reshape(-1) on 1D arrays, allows unknown shape in x
        # for this case only.
        if len(shape) == 1 and x.ndim == 1:
            return x
        missing_size = sanitize_index(x.size / reduce(mul, known_sizes, 1))
        shape = tuple(missing_size if s == -1 else s for s in shape)

    if np.isnan(sum(x.shape)):
        raise ValueError(
            "Array chunk size or shape is unknown. shape: %s\n\n"
            "Possible solution with x.compute_chunk_sizes()" % str(x.shape)
        )

    if reduce(mul, shape, 1) != x.size:
        raise ValueError("total size of new array must be unchanged")

    if x.shape == shape:
        return x

    meta = meta_from_array(x, len(shape))

    name = "reshape-" + tokenize(x, shape)

    if x.npartitions == 1:
        key = next(flatten(x.__dask_keys__()))
        dsk = {(name,) + (0,) * len(shape): (M.reshape, key, shape)}
        chunks = tuple((d,) for d in shape)
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[x])
        return Array(graph, name, chunks, meta=meta)

    # Logic or how to rechunk
    din = len(x.shape)
    dout = len(shape)
    if not merge_chunks and din > dout:
        x = x.rechunk({i: 1 for i in range(din - dout)})

    inchunks, outchunks = reshape_rechunk(x.shape, shape, x.chunks)
    # Check output chunks are not too large
    max_chunksize_in_bytes = reduce(mul, [max(i) for i in outchunks]) * x.dtype.itemsize

    if limit is None:
        limit = parse_bytes(config.get("array.chunk-size"))
        split = config.get("array.slicing.split-large-chunks", None)
    else:
        limit = parse_bytes(limit)
        split = True

    if max_chunksize_in_bytes > limit:
        if split is None:
            msg = (
                "Reshaping is producing a large chunk. To accept the large\n"
                "chunk and silence this warning, set the option\n"
                "    >>> with dask.config.set(**{'array.slicing.split_large_chunks': False}):\n"
                "    ...     array.reshape(shape)\n\n"
                "To avoid creating the large chunks, set the option\n"
                "    >>> with dask.config.set(**{'array.slicing.split_large_chunks': True}):\n"
                "    ...     array.reshape(shape)"
                "Explicitly passing ``limit`` to ``reshape`` will also silence this warning\n"
                "    >>> array.reshape(shape, limit='128 MiB')"
            )
            warnings.warn(msg, PerformanceWarning, stacklevel=6)
        elif split:
            # Leave chunk sizes unaltered where possible
            matching_chunks = Counter(inchunks) & Counter(outchunks)
            chunk_plan = []
            for out in outchunks:
                if matching_chunks[out] > 0:
                    chunk_plan.append(out)
                    matching_chunks[out] -= 1
                else:
                    chunk_plan.append("auto")
            outchunks = normalize_chunks(
                chunk_plan,
                shape=shape,
                limit=limit,
                dtype=x.dtype,
                previous_chunks=inchunks,
            )

    x2 = x.rechunk(inchunks)

    # Construct graph
    in_keys = list(product([x2.name], *[range(len(c)) for c in inchunks]))
    out_keys = list(product([name], *[range(len(c)) for c in outchunks]))
    shapes = list(product(*outchunks))
    dsk = {a: (M.reshape, b, shape) for a, b, shape in zip(out_keys, in_keys, shapes)}

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[x2])
    return Array(graph, name, outchunks, meta=meta)
