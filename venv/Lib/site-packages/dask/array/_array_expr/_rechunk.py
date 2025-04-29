from __future__ import annotations

import itertools
import operator

import numpy as np
import toolz

from dask.array._array_expr._expr import ArrayExpr
from dask.array.core import concatenate3
from dask.array.rechunk import (
    _balance_chunksizes,
    _choose_rechunk_method,
    _validate_rechunk,
    intersect_chunks,
    normalize_chunks,
    plan_rechunk,
    tokenize,
    validate_axis,
)
from dask.utils import cached_property


class Rechunk(ArrayExpr):
    _parameters = [
        "array",
        "_chunks",
        "threshold",
        "block_size_limit",
        "balance",
        "method",
    ]

    _defaults = {
        "_chunks": "auto",
        "threshold": None,
        "block_size_limit": None,
        "balance": None,
        "method": None,
    }

    @property
    def _meta(self):
        return self.array._meta

    @property
    def _name(self):
        return "rechunk-merge-" + tokenize(*self.operands)

    @cached_property
    def chunks(self):
        x = self.array
        chunks = self.operand("_chunks")

        # don't rechunk if array is empty
        if x.ndim > 0 and all(s == 0 for s in x.shape):
            return x.chunks

        if isinstance(chunks, dict):
            chunks = {validate_axis(c, x.ndim): v for c, v in chunks.items()}
            for i in range(x.ndim):
                if i not in chunks:
                    chunks[i] = x.chunks[i]
                elif chunks[i] is None:
                    chunks[i] = x.chunks[i]
        if isinstance(chunks, (tuple, list)):
            chunks = tuple(
                lc if lc is not None else rc for lc, rc in zip(chunks, x.chunks)
            )
        chunks = normalize_chunks(
            chunks,
            x.shape,
            limit=self.block_size_limit,
            dtype=x.dtype,
            previous_chunks=x.chunks,
        )

        if not len(chunks) == x.ndim:
            raise ValueError("Provided chunks are not consistent with shape")

        if self.balance:
            chunks = tuple(_balance_chunksizes(chunk) for chunk in chunks)

        _validate_rechunk(x.chunks, chunks)

        return chunks

    def _lower(self):

        if not self.balance and (self.chunks == self.array.chunks):
            return self.array

        method = self.method or _choose_rechunk_method(
            self.array.chunks, self.chunks, threshold=self.threshold
        )
        if method == "p2p":
            raise NotImplementedError
        else:
            return TasksRechunk(
                self.array, self.chunks, self.threshold, self.block_size_limit
            )


class TasksRechunk(Rechunk):
    _parameters = ["array", "_chunks", "threshold", "block_size_limit"]

    @cached_property
    def chunks(self):
        return self.operand("_chunks")

    def _lower(self):
        return

    def _layer(self):
        steps = plan_rechunk(
            self.array.chunks,
            self.chunks,
            self.array.dtype.itemsize,
            self.threshold,
            self.block_size_limit,
        )
        name = self.array.name
        old_chunks = self.array.chunks
        layers = []
        for i, c in enumerate(steps):
            level = len(steps) - i - 1
            name, old_chunks, layer = _compute_rechunk(
                name, old_chunks, c, level, self.name
            )
            layers.append(layer)

        return toolz.merge(*layers)


def _compute_rechunk(old_name, old_chunks, chunks, level, name):
    """Compute the rechunk of *x* to the given *chunks*."""
    # TODO: redo this logic
    # if x.size == 0:
    #     # Special case for empty array, as the algorithm below does not behave correctly
    #     return empty(x.shape, chunks=chunks, dtype=x.dtype)

    ndim = len(old_chunks)
    crossed = intersect_chunks(old_chunks, chunks)
    x2 = dict()
    intermediates = dict()
    # token = tokenize(old_name, chunks)
    if level != 0:
        merge_name = name.replace("rechunk-merge-", f"rechunk-merge-{level}-")
        split_name = name.replace("rechunk-merge-", f"rechunk-split-{level}-")
    else:
        merge_name = name.replace("rechunk-merge-", "rechunk-merge-")
        split_name = name.replace("rechunk-merge-", "rechunk-split-")
    split_name_suffixes = itertools.count()

    # Pre-allocate old block references, to allow re-use and reduce the
    # graph's memory footprint a bit.
    old_blocks = np.empty([len(c) for c in old_chunks], dtype="O")
    for index in np.ndindex(old_blocks.shape):
        old_blocks[index] = (old_name,) + index

    # Iterate over all new blocks
    new_index = itertools.product(*(range(len(c)) for c in chunks))

    for new_idx, cross1 in zip(new_index, crossed):
        key = (merge_name,) + new_idx
        old_block_indices = [[cr[i][0] for cr in cross1] for i in range(ndim)]
        subdims1 = [len(set(old_block_indices[i])) for i in range(ndim)]

        rec_cat_arg = np.empty(subdims1, dtype="O")
        rec_cat_arg_flat = rec_cat_arg.flat

        # Iterate over the old blocks required to build the new block
        for rec_cat_index, ind_slices in enumerate(cross1):
            old_block_index, slices = zip(*ind_slices)
            name = (split_name, next(split_name_suffixes))
            old_index = old_blocks[old_block_index][1:]
            if all(
                slc.start == 0 and slc.stop == old_chunks[i][ind]
                for i, (slc, ind) in enumerate(zip(slices, old_index))
            ):
                rec_cat_arg_flat[rec_cat_index] = old_blocks[old_block_index]
            else:
                intermediates[name] = (
                    operator.getitem,
                    old_blocks[old_block_index],
                    slices,
                )
                rec_cat_arg_flat[rec_cat_index] = name

        assert rec_cat_index == rec_cat_arg.size - 1

        # New block is formed by concatenation of sliced old blocks
        if all(d == 1 for d in rec_cat_arg.shape):
            x2[key] = rec_cat_arg.flat[0]
        else:
            x2[key] = (concatenate3, rec_cat_arg.tolist())

    del old_blocks, new_index

    return name, chunks, {**x2, **intermediates}
