from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from xarray.core.indexing import ImplicitToExplicitIndexingAdapter
from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint, T_ChunkedArray
from xarray.namedarray.utils import is_duck_dask_array, module_available

if TYPE_CHECKING:
    from xarray.namedarray._typing import (
        T_Chunks,
        _DType_co,
        _NormalizedChunks,
        duckarray,
    )

    try:
        from dask.array import Array as DaskArray
    except ImportError:
        DaskArray = np.ndarray[Any, Any]


dask_available = module_available("dask")


class DaskManager(ChunkManagerEntrypoint["DaskArray"]):
    array_cls: type[DaskArray]
    available: bool = dask_available

    def __init__(self) -> None:
        # TODO can we replace this with a class attribute instead?

        from dask.array import Array

        self.array_cls = Array

    def is_chunked_array(self, data: duckarray[Any, Any]) -> bool:
        return is_duck_dask_array(data)

    def chunks(self, data: Any) -> _NormalizedChunks:
        return data.chunks  # type: ignore[no-any-return]

    def normalize_chunks(
        self,
        chunks: T_Chunks | _NormalizedChunks,
        shape: tuple[int, ...] | None = None,
        limit: int | None = None,
        dtype: _DType_co | None = None,
        previous_chunks: _NormalizedChunks | None = None,
    ) -> Any:
        """Called by open_dataset"""
        from dask.array.core import normalize_chunks

        return normalize_chunks(
            chunks,
            shape=shape,
            limit=limit,
            dtype=dtype,
            previous_chunks=previous_chunks,
        )  # type: ignore[no-untyped-call]

    def from_array(
        self, data: Any, chunks: T_Chunks | _NormalizedChunks, **kwargs: Any
    ) -> DaskArray | Any:
        import dask.array as da

        if isinstance(data, ImplicitToExplicitIndexingAdapter):
            # lazily loaded backend array classes should use NumPy array operations.
            kwargs["meta"] = np.ndarray

        return da.from_array(
            data,
            chunks,
            **kwargs,
        )  # type: ignore[no-untyped-call]

    def compute(
        self, *data: Any, **kwargs: Any
    ) -> tuple[np.ndarray[Any, _DType_co], ...]:
        from dask.array import compute

        return compute(*data, **kwargs)  # type: ignore[no-untyped-call, no-any-return]

    def persist(self, *data: Any, **kwargs: Any) -> tuple[DaskArray | Any, ...]:
        from dask import persist

        return persist(*data, **kwargs)  # type: ignore[no-untyped-call, no-any-return]

    @property
    def array_api(self) -> Any:
        from dask import array as da

        return da

    def reduction(
        self,
        arr: T_ChunkedArray,
        func: Callable[..., Any],
        combine_func: Callable[..., Any] | None = None,
        aggregate_func: Callable[..., Any] | None = None,
        axis: int | Sequence[int] | None = None,
        dtype: _DType_co | None = None,
        keepdims: bool = False,
    ) -> DaskArray | Any:
        from dask.array import reduction

        return reduction(
            arr,
            chunk=func,
            combine=combine_func,
            aggregate=aggregate_func,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
        )  # type: ignore[no-untyped-call]

    def scan(
        self,
        func: Callable[..., Any],
        binop: Callable[..., Any],
        ident: float,
        arr: T_ChunkedArray,
        axis: int | None = None,
        dtype: _DType_co | None = None,
        **kwargs: Any,
    ) -> DaskArray | Any:
        from dask.array.reductions import cumreduction

        return cumreduction(
            func,
            binop,
            ident,
            arr,
            axis=axis,
            dtype=dtype,
            **kwargs,
        )  # type: ignore[no-untyped-call]

    def apply_gufunc(
        self,
        func: Callable[..., Any],
        signature: str,
        *args: Any,
        axes: Sequence[tuple[int, ...]] | None = None,
        axis: int | None = None,
        keepdims: bool = False,
        output_dtypes: Sequence[_DType_co] | None = None,
        output_sizes: dict[str, int] | None = None,
        vectorize: bool | None = None,
        allow_rechunk: bool = False,
        meta: tuple[np.ndarray[Any, _DType_co], ...] | None = None,
        **kwargs: Any,
    ) -> Any:
        from dask.array.gufunc import apply_gufunc

        return apply_gufunc(
            func,
            signature,
            *args,
            axes=axes,
            axis=axis,
            keepdims=keepdims,
            output_dtypes=output_dtypes,
            output_sizes=output_sizes,
            vectorize=vectorize,
            allow_rechunk=allow_rechunk,
            meta=meta,
            **kwargs,
        )  # type: ignore[no-untyped-call]

    def map_blocks(
        self,
        func: Callable[..., Any],
        *args: Any,
        dtype: _DType_co | None = None,
        chunks: tuple[int, ...] | None = None,
        drop_axis: int | Sequence[int] | None = None,
        new_axis: int | Sequence[int] | None = None,
        **kwargs: Any,
    ) -> Any:
        from dask.array import map_blocks

        # pass through name, meta, token as kwargs
        return map_blocks(
            func,
            *args,
            dtype=dtype,
            chunks=chunks,
            drop_axis=drop_axis,
            new_axis=new_axis,
            **kwargs,
        )  # type: ignore[no-untyped-call]

    def blockwise(
        self,
        func: Callable[..., Any],
        out_ind: Iterable[Any],
        *args: Any,
        # can't type this as mypy assumes args are all same type, but dask blockwise args alternate types
        name: str | None = None,
        token: Any | None = None,
        dtype: _DType_co | None = None,
        adjust_chunks: dict[Any, Callable[..., Any]] | None = None,
        new_axes: dict[Any, int] | None = None,
        align_arrays: bool = True,
        concatenate: bool | None = None,
        meta: tuple[np.ndarray[Any, _DType_co], ...] | None = None,
        **kwargs: Any,
    ) -> DaskArray | Any:
        from dask.array import blockwise

        return blockwise(
            func,
            out_ind,
            *args,
            name=name,
            token=token,
            dtype=dtype,
            adjust_chunks=adjust_chunks,
            new_axes=new_axes,
            align_arrays=align_arrays,
            concatenate=concatenate,
            meta=meta,
            **kwargs,
        )  # type: ignore[no-untyped-call]

    def unify_chunks(
        self,
        *args: Any,  # can't type this as mypy assumes args are all same type, but dask unify_chunks args alternate types
        **kwargs: Any,
    ) -> tuple[dict[str, _NormalizedChunks], list[DaskArray]]:
        from dask.array.core import unify_chunks

        return unify_chunks(*args, **kwargs)  # type: ignore[no-any-return, no-untyped-call]

    def store(
        self,
        sources: Any | Sequence[Any],
        targets: Any,
        **kwargs: Any,
    ) -> Any:
        from dask.array import store

        return store(
            sources=sources,
            targets=targets,
            **kwargs,
        )

    def shuffle(
        self, x: DaskArray, indexer: list[list[int]], axis: int, chunks: T_Chunks
    ) -> DaskArray:
        import dask.array

        if not module_available("dask", minversion="2024.08.1"):
            raise ValueError(
                "This method is very inefficient on dask<2024.08.1. Please upgrade."
            )
        if chunks is None:
            chunks = "auto"
        if chunks != "auto":
            raise NotImplementedError("Only chunks='auto' is supported at present.")
        return dask.array.shuffle(x, indexer, axis, chunks="auto")

    def get_auto_chunk_size(self) -> int:
        from dask import config as dask_config
        from dask.utils import parse_bytes

        return parse_bytes(dask_config.get("array.chunk-size"))
