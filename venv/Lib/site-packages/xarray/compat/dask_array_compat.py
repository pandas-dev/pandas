from typing import Any

from xarray.namedarray.utils import module_available


def reshape_blockwise(
    x: Any,
    shape: int | tuple[int, ...],
    chunks: tuple[tuple[int, ...], ...] | None = None,
):
    if module_available("dask", "2024.08.2"):
        from dask.array import reshape_blockwise

        return reshape_blockwise(x, shape=shape, chunks=chunks)
    else:
        return x.reshape(shape)


def sliding_window_view(
    x, window_shape, axis=None, *, automatic_rechunk=True, **kwargs
):
    # Backcompat for handling `automatic_rechunk`, delete when dask>=2024.11.0
    # Note that subok, writeable are unsupported by dask, so we ignore those in kwargs
    from dask.array.lib.stride_tricks import sliding_window_view

    if module_available("dask", "2024.11.0"):
        return sliding_window_view(
            x, window_shape=window_shape, axis=axis, automatic_rechunk=automatic_rechunk
        )
    else:
        # automatic_rechunk is not supported
        return sliding_window_view(x, window_shape=window_shape, axis=axis)
