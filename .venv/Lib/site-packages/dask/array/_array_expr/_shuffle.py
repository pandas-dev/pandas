from __future__ import annotations

import copy
import functools
from itertools import count, product

import numpy as np

from dask import config
from dask._task_spec import DataNode, List, Task, TaskRef
from dask.array._array_expr._expr import ArrayExpr
from dask.array.chunk import getitem
from dask.array.dispatch import concatenate_lookup, take_lookup
from dask.base import tokenize


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


def _shuffle(x, indexer, axis, name):
    _validate_indexer(x.chunks, indexer, axis)

    if len(indexer) == len(x.chunks[axis]):
        # check if the array is already shuffled the way we want
        ctr = 0
        for idx, c in zip(indexer, x.chunks[axis]):
            if idx != list(range(ctr, ctr + c)):
                break
            ctr += c
        else:
            return x
    return Shuffle(x, indexer, axis, name)


class Shuffle(ArrayExpr):
    _parameters = ["array", "indexer", "axis", "name"]

    @functools.cached_property
    def _meta(self):
        return self.array._meta

    @functools.cached_property
    def _name(self):
        return f"{self.operand('name')}-{self._token}"

    @functools.cached_property
    def chunks(self):
        output_chunks = []
        for i, c in enumerate(self.array.chunks):
            if i == self.axis:
                output_chunks.append(tuple(map(len, self._new_chunks)))
            else:
                output_chunks.append(c)
        return tuple(output_chunks)

    @functools.cached_property
    def _chunksize_tolerance(self):
        return config.get("array.chunk-size-tolerance")

    @functools.cached_property
    def _chunk_size_limit(self):
        return int(
            sum(self.array.chunks[self.axis])
            / len(self.array.chunks[self.axis])
            * self._chunksize_tolerance
        )

    @functools.cached_property
    def _new_chunks(self):
        current_chunk, new_chunks = [], []
        for idx in copy.deepcopy(self.indexer):
            if (
                len(current_chunk) + len(idx) > self._chunk_size_limit
                and len(current_chunk) > 0
            ):
                new_chunks.append(current_chunk)
                current_chunk = idx.copy()
            else:
                current_chunk.extend(idx)
                if (
                    len(current_chunk)
                    > self._chunk_size_limit / self._chunksize_tolerance
                ):
                    new_chunks.append(current_chunk)
                    current_chunk = []
        if len(current_chunk) > 0:
            new_chunks.append(current_chunk)
        return new_chunks

    def _layer(self) -> dict:
        chunks = self.array.chunks
        axis = self.axis

        chunk_boundaries = np.cumsum(chunks[axis])

        # Get existing chunk tuple locations
        chunk_tuples = list(
            product(*(range(len(c)) for i, c in enumerate(chunks) if i != axis))
        )

        intermediates: dict = dict()
        merges: dict = dict()
        dtype = np.min_scalar_type(max(max(chunks[axis]), self._chunk_size_limit))
        split_name = f"shuffle-split-{self.deterministic_token}"
        slices = [slice(None)] * len(chunks)
        split_name_suffixes = count()
        sorter_name = "shuffle-sorter-"
        taker_name = "shuffle-taker-"

        old_blocks = {
            old_index: (self.array._name,) + old_index
            for old_index in np.ndindex(tuple([len(c) for c in chunks]))
        }

        for new_chunk_idx, new_chunk_taker in enumerate(self._new_chunks):
            new_chunk_taker = np.array(new_chunk_taker)
            sorter = np.argsort(new_chunk_taker).astype(dtype)
            sorter_key = sorter_name + tokenize(sorter)
            # low level fusion can't deal with arrays on first position
            merges[sorter_key] = DataNode(sorter_key, (1, sorter))

            sorted_array = new_chunk_taker[sorter]
            source_chunk_nr, taker_boundary = np.unique(
                np.searchsorted(chunk_boundaries, sorted_array, side="right"),
                return_index=True,
            )
            taker_boundary = taker_boundary.tolist()
            taker_boundary.append(len(new_chunk_taker))

            taker_cache: dict = {}
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
                        name,
                        _getitem,
                        TaskRef(old_blocks[chunk_key]),
                        TaskRef(taker_key),
                    )
                    merge_keys.append(name)

                merge_suffix = convert_key(chunk_tuple, new_chunk_idx, axis)
                out_name_merge = (self._name,) + merge_suffix
                if len(merge_keys) > 1:
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

        return {**merges, **intermediates}


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
