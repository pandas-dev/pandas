"""
Functions for handling chunked arrays.
"""

from __future__ import annotations

import itertools
from collections.abc import Hashable, Mapping
from functools import lru_cache
from numbers import Number
from typing import TYPE_CHECKING, Any, Literal, TypeVar, Union, overload

from xarray.core import utils
from xarray.core.utils import emit_user_level_warning
from xarray.core.variable import IndexVariable, Variable
from xarray.namedarray.parallelcompat import (
    ChunkManagerEntrypoint,
    get_chunked_array_type,
    guess_chunkmanager,
)

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import T_ChunkDim

    MissingCoreDimOptions = Literal["raise", "copy", "drop"]


@lru_cache(maxsize=512)
def _get_breaks_cached(
    *,
    size: int,
    chunk_sizes: tuple[int, ...],
    preferred_chunk_sizes: int | tuple[int, ...],
) -> int | None:
    if isinstance(preferred_chunk_sizes, int) and preferred_chunk_sizes == 1:
        # short-circuit for the trivial case
        return None
    # Determine the stop indices of the preferred chunks, but omit the last stop
    # (equal to the dim size).  In particular, assume that when a sequence
    # expresses the preferred chunks, the sequence sums to the size.
    preferred_stops = (
        range(preferred_chunk_sizes, size, preferred_chunk_sizes)
        if isinstance(preferred_chunk_sizes, int)
        else set(itertools.accumulate(preferred_chunk_sizes[:-1]))
    )

    # Gather any stop indices of the specified chunks that are not a stop index
    # of a preferred chunk. Again, omit the last stop, assuming that it equals
    # the dim size.
    actual_stops = itertools.accumulate(chunk_sizes[:-1])
    # This copy is required for parallel iteration
    actual_stops_2 = itertools.accumulate(chunk_sizes[:-1])

    disagrees = itertools.compress(
        actual_stops_2, (a not in preferred_stops for a in actual_stops)
    )
    try:
        return next(disagrees)
    except StopIteration:
        return None


def _get_chunk(var: Variable, chunks, chunkmanager: ChunkManagerEntrypoint):
    """
    Return map from each dim to chunk sizes, accounting for backend's preferred chunks.
    """
    if isinstance(var, IndexVariable):
        return {}
    dims = var.dims
    shape = var.shape

    # Determine the explicit requested chunks.
    preferred_chunks = var.encoding.get("preferred_chunks", {})
    preferred_chunk_shape = tuple(
        itertools.starmap(preferred_chunks.get, zip(dims, shape, strict=True))
    )
    if isinstance(chunks, Number) or (chunks == "auto"):
        chunks = dict.fromkeys(dims, chunks)
    chunk_shape = tuple(
        chunks.get(dim, None) or preferred_chunk_sizes
        for dim, preferred_chunk_sizes in zip(dims, preferred_chunk_shape, strict=True)
    )

    chunk_shape = chunkmanager.normalize_chunks(
        chunk_shape, shape=shape, dtype=var.dtype, previous_chunks=preferred_chunk_shape
    )

    # Warn where requested chunks break preferred chunks, provided that the variable
    # contains data.
    if var.size:
        for dim, size, chunk_sizes in zip(dims, shape, chunk_shape, strict=True):
            try:
                preferred_chunk_sizes = preferred_chunks[dim]
            except KeyError:
                continue
            disagreement = _get_breaks_cached(
                size=size,
                chunk_sizes=chunk_sizes,
                preferred_chunk_sizes=preferred_chunk_sizes,
            )
            if disagreement:
                emit_user_level_warning(
                    "The specified chunks separate the stored chunks along "
                    f'dimension "{dim}" starting at index {disagreement}. This could '
                    "degrade performance. Instead, consider rechunking after loading.",
                )

    return dict(zip(dims, chunk_shape, strict=True))


def _maybe_chunk(
    name: Hashable,
    var: Variable,
    chunks: Mapping[Any, T_ChunkDim] | None,
    token=None,
    lock=None,
    name_prefix: str = "xarray-",
    overwrite_encoded_chunks: bool = False,
    inline_array: bool = False,
    chunked_array_type: str | ChunkManagerEntrypoint | None = None,
    from_array_kwargs=None,
) -> Variable:
    from xarray.namedarray.daskmanager import DaskManager

    if chunks is not None:
        chunks = {dim: chunks[dim] for dim in var.dims if dim in chunks}

    if var.ndim:
        chunked_array_type = guess_chunkmanager(
            chunked_array_type
        )  # coerce string to ChunkManagerEntrypoint type
        if isinstance(chunked_array_type, DaskManager):
            from dask.base import tokenize

            # when rechunking by different amounts, make sure dask names change
            # by providing chunks as an input to tokenize.
            # subtle bugs result otherwise. see GH3350
            # we use str() for speed, and use the name for the final array name on the next line
            token2 = tokenize(token or var._data, str(chunks))
            name2 = f"{name_prefix}{name}-{token2}"

            from_array_kwargs = utils.consolidate_dask_from_array_kwargs(
                from_array_kwargs,
                name=name2,
                lock=lock,
                inline_array=inline_array,
            )

        var = var.chunk(
            chunks,
            chunked_array_type=chunked_array_type,
            from_array_kwargs=from_array_kwargs,
        )

        if overwrite_encoded_chunks and var.chunks is not None:
            var.encoding["chunks"] = tuple(x[0] for x in var.chunks)
        return var
    else:
        return var


_T = TypeVar("_T", bound=Union["Dataset", "DataArray"])
_U = TypeVar("_U", bound=Union["Dataset", "DataArray"])
_V = TypeVar("_V", bound=Union["Dataset", "DataArray"])


@overload
def unify_chunks(obj: _T, /) -> tuple[_T]: ...


@overload
def unify_chunks(obj1: _T, obj2: _U, /) -> tuple[_T, _U]: ...


@overload
def unify_chunks(obj1: _T, obj2: _U, obj3: _V, /) -> tuple[_T, _U, _V]: ...


@overload
def unify_chunks(*objects: Dataset | DataArray) -> tuple[Dataset | DataArray, ...]: ...


def unify_chunks(*objects: Dataset | DataArray) -> tuple[Dataset | DataArray, ...]:
    """
    Given any number of Dataset and/or DataArray objects, returns
    new objects with unified chunk size along all chunked dimensions.

    Returns
    -------
    unified (DataArray or Dataset) – Tuple of objects with the same type as
    *objects with consistent chunk sizes for all dask-array variables

    See Also
    --------
    dask.array.core.unify_chunks
    """
    from xarray.core.dataarray import DataArray

    # Convert all objects to datasets
    datasets = [
        obj._to_temp_dataset() if isinstance(obj, DataArray) else obj.copy()
        for obj in objects
    ]

    # Get arguments to pass into dask.array.core.unify_chunks
    unify_chunks_args = []
    sizes: dict[Hashable, int] = {}
    for ds in datasets:
        for v in ds._variables.values():
            if v.chunks is not None:
                # Check that sizes match across different datasets
                for dim, size in v.sizes.items():
                    try:
                        if sizes[dim] != size:
                            raise ValueError(
                                f"Dimension {dim!r} size mismatch: {sizes[dim]} != {size}"
                            )
                    except KeyError:
                        sizes[dim] = size
                unify_chunks_args += [v._data, v._dims]

    # No dask arrays: Return inputs
    if not unify_chunks_args:
        return objects

    chunkmanager = get_chunked_array_type(*list(unify_chunks_args))
    _, chunked_data = chunkmanager.unify_chunks(*unify_chunks_args)
    chunked_data_iter = iter(chunked_data)
    out: list[Dataset | DataArray] = []
    for obj, ds in zip(objects, datasets, strict=True):
        for k, v in ds._variables.items():
            if v.chunks is not None:
                ds._variables[k] = v.copy(data=next(chunked_data_iter))
        out.append(obj._from_temp_dataset(ds) if isinstance(obj, DataArray) else ds)

    return tuple(out)
