from __future__ import annotations

import importlib
import itertools
import sys
import warnings
from collections.abc import Hashable, Iterable, Iterator, Mapping
from functools import lru_cache
from numbers import Number
from typing import TYPE_CHECKING, Any, TypeVar, cast

import numpy as np
from packaging.version import Version

from xarray.namedarray._typing import ErrorOptionsWithWarn, _DimsLike

if TYPE_CHECKING:
    from typing import TypeGuard

    from numpy.typing import NDArray

    try:
        from dask.array.core import Array as DaskArray
        from dask.typing import DaskCollection
    except ImportError:
        DaskArray = NDArray  # type: ignore[assignment, misc]
        DaskCollection: Any = NDArray  # type: ignore[no-redef]

    from xarray.core.types import T_ChunkDim
    from xarray.namedarray._typing import DuckArray, _Dim, duckarray
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint


K = TypeVar("K")
V = TypeVar("V")
T = TypeVar("T")


@lru_cache
def module_available(module: str, minversion: str | None = None) -> bool:
    """Checks whether a module is installed without importing it.

    Use this for a lightweight check and lazy imports.

    Parameters
    ----------
    module : str
        Name of the module.
    minversion : str, optional
        Minimum version of the module

    Returns
    -------
    available : bool
        Whether the module is installed.
    """
    if importlib.util.find_spec(module) is None:
        return False

    if minversion is not None:
        version = importlib.metadata.version(module)

        return Version(version) >= Version(minversion)

    return True


def is_dask_collection(x: object) -> TypeGuard[DaskCollection]:
    if module_available("dask"):
        from dask.base import is_dask_collection

        # use is_dask_collection function instead of dask.typing.DaskCollection
        # see https://github.com/pydata/xarray/pull/8241#discussion_r1476276023
        return is_dask_collection(x)
    return False


def is_duck_array(value: Any) -> TypeGuard[duckarray[Any, Any]]:
    # TODO: replace is_duck_array with runtime checks via _arrayfunction_or_api protocol on
    # python 3.12 and higher (see https://github.com/pydata/xarray/issues/8696#issuecomment-1924588981)
    if isinstance(value, np.ndarray):
        return True
    return (
        hasattr(value, "ndim")
        and hasattr(value, "shape")
        and hasattr(value, "dtype")
        and (
            (hasattr(value, "__array_function__") and hasattr(value, "__array_ufunc__"))
            or hasattr(value, "__array_namespace__")
        )
    )


def is_duck_dask_array(x: duckarray[Any, Any]) -> TypeGuard[DaskArray]:
    return is_duck_array(x) and is_dask_collection(x)


def to_0d_object_array(
    value: object,
) -> NDArray[np.object_]:
    """Given a value, wrap it in a 0-D numpy.ndarray with dtype=object."""
    result = np.empty((), dtype=object)
    result[()] = value
    return result


def is_dict_like(value: Any) -> TypeGuard[Mapping[Any, Any]]:
    return hasattr(value, "keys") and hasattr(value, "__getitem__")


def drop_missing_dims(
    supplied_dims: Iterable[_Dim],
    dims: Iterable[_Dim],
    missing_dims: ErrorOptionsWithWarn,
) -> _DimsLike:
    """Depending on the setting of missing_dims, drop any dimensions from supplied_dims that
    are not present in dims.

    Parameters
    ----------
    supplied_dims : Iterable of Hashable
    dims : Iterable of Hashable
    missing_dims : {"raise", "warn", "ignore"}
    """

    if missing_dims == "raise":
        supplied_dims_set = {val for val in supplied_dims if val is not ...}
        if invalid := supplied_dims_set - set(dims):
            raise ValueError(
                f"Dimensions {invalid} do not exist. Expected one or more of {dims}"
            )

        return supplied_dims

    elif missing_dims == "warn":
        if invalid := set(supplied_dims) - set(dims):
            warnings.warn(
                f"Dimensions {invalid} do not exist. Expected one or more of {dims}",
                stacklevel=2,
            )

        return [val for val in supplied_dims if val in dims or val is ...]

    elif missing_dims == "ignore":
        return [val for val in supplied_dims if val in dims or val is ...]

    else:
        raise ValueError(
            f"Unrecognised option {missing_dims} for missing_dims argument"
        )


def infix_dims(
    dims_supplied: Iterable[_Dim],
    dims_all: Iterable[_Dim],
    missing_dims: ErrorOptionsWithWarn = "raise",
) -> Iterator[_Dim]:
    """
    Resolves a supplied list containing an ellipsis representing other items, to
    a generator with the 'realized' list of all items
    """
    if ... in dims_supplied:
        dims_all_list = list(dims_all)
        if len(set(dims_all)) != len(dims_all_list):
            raise ValueError("Cannot use ellipsis with repeated dims")
        if list(dims_supplied).count(...) > 1:
            raise ValueError("More than one ellipsis supplied")
        other_dims = [d for d in dims_all if d not in dims_supplied]
        existing_dims = drop_missing_dims(dims_supplied, dims_all, missing_dims)
        for d in existing_dims:
            if d is ...:
                yield from other_dims
            else:
                yield d
    else:
        existing_dims = drop_missing_dims(dims_supplied, dims_all, missing_dims)
        if set(existing_dims) ^ set(dims_all):
            raise ValueError(
                f"{dims_supplied} must be a permuted list of {dims_all}, unless `...` is included"
            )
        yield from existing_dims


def either_dict_or_kwargs(
    pos_kwargs: Mapping[Any, T] | None,
    kw_kwargs: Mapping[str, T],
    func_name: str,
) -> Mapping[Hashable, T]:
    if pos_kwargs is None or pos_kwargs == {}:
        # Need an explicit cast to appease mypy due to invariance; see
        # https://github.com/python/mypy/issues/6228
        return cast(Mapping[Hashable, T], kw_kwargs)

    if not is_dict_like(pos_kwargs):
        raise ValueError(f"the first argument to .{func_name} must be a dictionary")
    if kw_kwargs:
        raise ValueError(
            f"cannot specify both keyword and positional arguments to .{func_name}"
        )
    return pos_kwargs


def _get_chunk(  # type: ignore[no-untyped-def]
    data: DuckArray[Any],
    chunks,
    chunkmanager: ChunkManagerEntrypoint[Any],
    *,
    preferred_chunks,
    dims=None,
) -> Mapping[Any, T_ChunkDim]:
    """
    Return map from each dim to chunk sizes, accounting for backend's preferred chunks.
    """
    from xarray.core.common import _contains_cftime_datetimes
    from xarray.core.utils import emit_user_level_warning
    from xarray.structure.chunks import _get_breaks_cached

    dims = chunks.keys() if dims is None else dims
    shape = data.shape

    # Determine the explicit requested chunks.
    preferred_chunk_shape = tuple(
        itertools.starmap(preferred_chunks.get, zip(dims, shape, strict=True))
    )
    if isinstance(chunks, Number) or (chunks == "auto"):
        chunks = dict.fromkeys(dims, chunks)
    chunk_shape = tuple(
        chunks.get(dim, None) or preferred_chunk_sizes
        for dim, preferred_chunk_sizes in zip(dims, preferred_chunk_shape, strict=True)
    )

    limit: int | None
    if _contains_cftime_datetimes(data):
        limit, dtype = fake_target_chunksize(data, chunkmanager.get_auto_chunk_size())
    else:
        limit = None
        dtype = data.dtype

    chunk_shape = chunkmanager.normalize_chunks(
        chunk_shape,
        shape=shape,
        dtype=dtype,
        limit=limit,
        previous_chunks=preferred_chunk_shape,
    )

    # Warn where requested chunks break preferred chunks, provided that the variable
    # contains data.
    if data.size:  # type: ignore[unused-ignore,attr-defined]  # DuckArray protocol doesn't include 'size' - should it?
        for dim, size, chunk_sizes in zip(dims, shape, chunk_shape, strict=True):
            if preferred_chunk_sizes := preferred_chunks.get(dim):
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


def fake_target_chunksize(
    data: DuckArray[Any],
    limit: int,
) -> tuple[int, np.dtype[Any]]:
    """
    The `normalize_chunks` algorithm takes a size `limit` in bytes, but will not
    work for object dtypes.  So we rescale the `limit` to an appropriate one based
    on `float64` dtype, and pass that to `normalize_chunks`.

    Arguments
    ---------
    data : Variable or ChunkedArray
        The data for which we want to determine chunk sizes.
    limit : int
        The target chunk size in bytes. Passed to the chunk manager's `normalize_chunks` method.
    """

    # Short circuit for non-object dtypes
    from xarray.core.common import _contains_cftime_datetimes

    if not _contains_cftime_datetimes(data):
        return limit, data.dtype

    from xarray.core.formatting import first_n_items

    output_dtype = np.dtype(np.float64)

    nbytes_approx: int = sys.getsizeof(first_n_items(data, 1))  # type: ignore[no-untyped-call]

    f64_nbytes = output_dtype.itemsize

    limit = int(limit * (f64_nbytes / nbytes_approx))

    return limit, output_dtype


class ReprObject:
    """Object that prints as the given value, for use with sentinel values."""

    __slots__ = ("_value",)

    _value: str

    def __init__(self, value: str):
        self._value = value

    def __repr__(self) -> str:
        return self._value

    def __eq__(self, other: ReprObject | Any) -> bool:
        # TODO: What type can other be? ArrayLike?
        return self._value == other._value if isinstance(other, ReprObject) else False

    def __hash__(self) -> int:
        return hash((type(self), self._value))

    def __dask_tokenize__(self) -> object:
        from dask.base import normalize_token

        return normalize_token((type(self), self._value))
