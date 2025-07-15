from __future__ import annotations

import math
from functools import reduce
from itertools import count, product
from operator import mul
from typing import Literal

import numpy as np

from dask import config
from dask._task_spec import DataNode, List, Task, TaskRef
from dask.array.chunk import getitem
from dask.array.core import Array, unknown_chunk_message
from dask.array.dispatch import concatenate_lookup, take_lookup
from dask.base import tokenize
from dask.highlevelgraph import HighLevelGraph


def shuffle(x, indexer: list[list[int]], axis: int, chunks: Literal["auto"] = "auto"):
    """
    Reorders one dimensions of a Dask Array based on an indexer.

    The indexer defines a list of positional groups that will end up in the same chunk
    together. A single group is in at most one chunk on this dimension, but a chunk
    might contain multiple groups to avoid fragmentation of the array.

    The algorithm tries to balance the chunksizes as much as possible to ideally keep the
    number of chunks consistent or at least manageable.

    Parameters
    ----------
    x: dask array
        Array to be shuffled.
    indexer:  list[list[int]]
        The indexer that determines which elements along the dimension will end up in the
        same chunk. Multiple groups can be in the same chunk to avoid fragmentation, but
        each group will end up in exactly one chunk.
    axis: int
        The axis to shuffle along.
    chunks: "auto"
        Hint on how to rechunk if single groups are becoming too large. The default is
        to split chunks along the other dimensions evenly to keep the chunksize
        consistent. The rechunking is done in a way that ensures that non all-to-all
        network communication is necessary, chunks are only split and not combined with
        other chunks.

    Examples
    --------
    >>> import dask.array as da
    >>> import numpy as np
    >>> arr = np.array([[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11, 12, 13, 14, 15, 16]])
    >>> x = da.from_array(arr, chunks=(2, 4))

    Separate the elements in different groups.

    >>> y = x.shuffle([[6, 5, 2], [4, 1], [3, 0, 7]], axis=1)

    The shuffle algorihthm will combine the first 2 groups into a single chunk to keep
    the number of chunks small.

    The tolerance of increasing the chunk size is controlled by the configuration
    "array.chunk-size-tolerance". The default value is 1.25.

    >>> y.chunks
    ((2,), (5, 3))

    The array was reordered along axis 1 according to the positional indexer that was given.

    >>> y.compute()
    array([[ 7,  6,  3,  5,  2,  4,  1,  8],
           [15, 14, 11, 13, 10, 12,  9, 16]])
    """
    if np.isnan(x.shape).any():
        raise ValueError(
            f"Shuffling only allowed with known chunk sizes. {unknown_chunk_message}"
        )
    assert isinstance(axis, int), "axis must be an integer"
    _validate_indexer(x.chunks, indexer, axis)

    x = _rechunk_other_dimensions(x, max(map(len, indexer)), axis, chunks)

    token = tokenize(x, indexer, axis)
    out_name = f"shuffle-{token}"

    chunks, layer = _shuffle(x.chunks, indexer, axis, x.name, out_name, token)
    if len(layer) == 0:
        return Array(x.dask, x.name, x.chunks, meta=x)

    graph = HighLevelGraph.from_collections(out_name, layer, dependencies=[x])

    return Array(graph, out_name, chunks, meta=x)


def _calculate_new_chunksizes(
    input_chunks, new_chunks, changeable_dimensions: set, maximum_chunk: int
):

    chunksize_tolerance = config.get("array.chunk-size-tolerance")
    maximum_chunk = max(maximum_chunk, 1)

    # iterate until we distributed the increase in chunksize across all dimensions
    # or every non-shuffle dimension is all 1
    while changeable_dimensions:
        n_changeable_dimensions = len(changeable_dimensions)
        chunksize_inc_factor = reduce(mul, map(max, new_chunks)) / maximum_chunk
        if chunksize_inc_factor <= 1:
            break

        for i in list(changeable_dimensions):
            new_chunksizes = []
            # calculate what the max chunk size in this dimension is and split every
            # chunk that is larger than that. We split the increase factor evenly
            # between all dimensions that are not shuffled.
            up_chunksize_limit_for_dim = max(new_chunks[i]) / (
                chunksize_inc_factor ** (1 / n_changeable_dimensions)
            )
            for c in input_chunks[i]:
                if c > chunksize_tolerance * up_chunksize_limit_for_dim:
                    factor = math.ceil(c / up_chunksize_limit_for_dim)

                    # Ensure that we end up at least with chunksize 1
                    factor = min(factor, c)

                    chunksize, remainder = divmod(c, factor)
                    nc = [chunksize] * factor
                    for ii in range(remainder):
                        # Add remainder parts to the first few chunks
                        nc[ii] += 1
                    new_chunksizes.extend(nc)

                else:
                    new_chunksizes.append(c)

            if tuple(new_chunksizes) == new_chunks[i] or max(new_chunksizes) == 1:
                changeable_dimensions.remove(i)

            new_chunks[i] = tuple(new_chunksizes)
    return new_chunks


def _rechunk_other_dimensions(
    x: Array, longest_group: int, axis: int, chunks: Literal["auto"]
) -> Array:
    assert chunks == "auto", "Only auto is supported for now"
    chunksize_tolerance = config.get("array.chunk-size-tolerance")

    if longest_group <= max(x.chunks[axis]) * chunksize_tolerance:
        # We are staying below our threshold, so don't rechunk
        return x

    changeable_dimensions = set(range(len(x.chunks))) - {axis}
    new_chunks = list(x.chunks)
    new_chunks[axis] = (longest_group,)

    # How large is the largest chunk in the input
    maximum_chunk = reduce(mul, map(max, x.chunks))

    new_chunks = _calculate_new_chunksizes(
        x.chunks, new_chunks, changeable_dimensions, maximum_chunk
    )
    new_chunks[axis] = x.chunks[axis]
    return x.rechunk(tuple(new_chunks))


def _validate_indexer(chunks, indexer, axis):
    if not isinstance(indexer, list) or not all(isinstance(i, list) for i in indexer):
        raise ValueError("indexer must be a list of lists of positional indices")

    if not axis <= len(chunks):
        raise ValueError(
            f"Axis {axis} is out of bounds for array with {len(chunks)} axes"
        )

    if max(map(max, indexer)) >= sum(chunks[axis]):
        raise IndexError(
            f"Indexer contains out of bounds index. Dimension only has {sum(chunks[axis])} elements."
        )


def _shuffle(chunks, indexer, axis, in_name, out_name, token):
    _validate_indexer(chunks, indexer, axis)

    if len(indexer) == len(chunks[axis]):
        # check if the array is already shuffled the way we want
        ctr = 0
        for idx, c in zip(indexer, chunks[axis]):
            if idx != list(range(ctr, ctr + c)):
                break
            ctr += c
        else:
            return chunks, {}

    chunksize_tolerance = config.get("array.chunk-size-tolerance")
    chunk_size_limit = int(sum(chunks[axis]) / len(chunks[axis]) * chunksize_tolerance)

    # Figure out how many groups we can put into one chunk
    current_chunk, new_chunks = [], []
    for idx in indexer:
        if len(current_chunk) + len(idx) > chunk_size_limit and len(current_chunk) > 0:
            new_chunks.append(current_chunk)
            current_chunk = idx.copy()
        else:
            current_chunk.extend(idx)
            if len(current_chunk) > chunk_size_limit / chunksize_tolerance:
                new_chunks.append(current_chunk)
                current_chunk = []
    if len(current_chunk) > 0:
        new_chunks.append(current_chunk)

    chunk_boundaries = np.cumsum(chunks[axis])

    # Get existing chunk tuple locations
    chunk_tuples = list(
        product(*(range(len(c)) for i, c in enumerate(chunks) if i != axis))
    )

    intermediates = dict()
    merges = dict()
    dtype = np.min_scalar_type(max(*chunks[axis], chunk_size_limit))
    split_name = f"shuffle-split-{token}"
    slices = [slice(None)] * len(chunks)
    split_name_suffixes = count()
    sorter_name = "shuffle-sorter-"
    taker_name = "shuffle-taker-"

    old_blocks = {
        old_index: (in_name,) + old_index
        for old_index in np.ndindex(tuple([len(c) for c in chunks]))
    }
    for new_chunk_idx, new_chunk_taker in enumerate(new_chunks):
        new_chunk_taker = np.array(new_chunk_taker)
        sorter = np.argsort(new_chunk_taker).astype(dtype)
        sorter_key = None

        sorted_array = new_chunk_taker[sorter]
        source_chunk_nr, taker_boundary = np.unique(
            np.searchsorted(chunk_boundaries, sorted_array, side="right"),
            return_index=True,
        )
        taker_boundary = taker_boundary.tolist()
        taker_boundary.append(len(new_chunk_taker))

        taker_cache = {}
        for chunk_tuple in chunk_tuples:
            merge_keys = []

            for c, b_start, b_end in zip(
                source_chunk_nr, taker_boundary[:-1], taker_boundary[1:]
            ):
                # insert our axis chunk id into the chunk_tuple
                chunk_key = convert_key(chunk_tuple, c, axis)
                name = (split_name, next(split_name_suffixes))
                this_slice = slices.copy()

                # Cache the takers to allow de-duplication when serializing
                # Ugly!
                if c in taker_cache:
                    taker_key = taker_cache[c]
                else:
                    this_slice[axis] = (
                        sorted_array[b_start:b_end]
                        - (chunk_boundaries[c - 1] if c > 0 else 0)
                    ).astype(dtype)
                    if len(source_chunk_nr) == 1:
                        this_slice[axis] = this_slice[axis][np.argsort(sorter)]

                    taker_key = taker_name + tokenize(this_slice)
                    # low level fusion can't deal with arrays on first position
                    intermediates[taker_key] = DataNode(
                        taker_key, (1, tuple(this_slice))
                    )
                    taker_cache[c] = taker_key

                intermediates[name] = Task(
                    name, _getitem, TaskRef(old_blocks[chunk_key]), TaskRef(taker_key)
                )
                merge_keys.append(name)

            merge_suffix = convert_key(chunk_tuple, new_chunk_idx, axis)
            out_name_merge = (out_name,) + merge_suffix
            if len(merge_keys) > 1:
                if sorter_key is None:
                    sorter_key = sorter_name + tokenize(sorter)
                    # low level fusion can't deal with arrays on first position
                    merges[sorter_key] = DataNode(sorter_key, (1, sorter))
                merges[out_name_merge] = Task(
                    out_name_merge,
                    concatenate_arrays,
                    List(*(TaskRef(m) for m in merge_keys)),
                    TaskRef(sorter_key),
                    axis,
                )
            elif len(merge_keys) == 1:
                t = intermediates.pop(merge_keys[0])
                t.key = out_name_merge
                merges[out_name_merge] = t
            else:
                raise NotImplementedError

    output_chunks = []
    for i, c in enumerate(chunks):
        if i == axis:
            output_chunks.append(tuple(map(len, new_chunks)))
        else:
            output_chunks.append(c)

    layer = {**merges, **intermediates}
    return tuple(output_chunks), layer


def _getitem(obj, index):
    return getitem(obj, index[1])


def concatenate_arrays(arrs, sorter, axis):
    return take_lookup(
        concatenate_lookup.dispatch(type(arrs[0]))(arrs, axis=axis),
        np.argsort(sorter[1]),
        axis=axis,
    )


def convert_key(key, chunk, axis):
    key = list(key)
    key.insert(axis, chunk)
    return tuple(key)
