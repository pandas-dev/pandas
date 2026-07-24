from typing import Any

from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array
from xarray.namedarray.utils import module_available


def reshape_blockwise(
    x: Any,
    shape: int | tuple[int, ...],
    chunks: tuple[tuple[int, ...], ...] | None = None,
):
    try:
        array_api = get_chunked_array_type(x).array_api
    except TypeError:
        if is_chunked_array(x):
            raise
        array_api = None

    if array_api is not None and hasattr(array_api, "reshape_blockwise"):
        return array_api.reshape_blockwise(x, shape=shape, chunks=chunks)
    elif module_available("dask", "2024.08.2"):
        from dask.array import reshape_blockwise as dask_reshape_blockwise

        return dask_reshape_blockwise(x, shape=shape, chunks=chunks)
    else:
        return x.reshape(shape)


def sliding_window_view(
    x, window_shape, axis=None, *, automatic_rechunk=True, **kwargs
):
    # Backcompat for handling `automatic_rechunk`, delete when dask>=2024.11.0
    # Note that subok, writeable are unsupported by dask, so we ignore those in kwargs
    try:
        array_api = get_chunked_array_type(x).array_api
    except TypeError:
        if is_chunked_array(x):
            raise
        array_api = None

    if array_api is not None:
        array_sliding_window_view = getattr(array_api, "sliding_window_view", None)
        if array_sliding_window_view is not None:
            return array_sliding_window_view(
                x,
                window_shape=window_shape,
                axis=axis,
                automatic_rechunk=automatic_rechunk,
            )

    from dask.array.lib.stride_tricks import sliding_window_view

    if module_available("dask", "2024.11.0"):
        return sliding_window_view(
            x, window_shape=window_shape, axis=axis, automatic_rechunk=automatic_rechunk
        )
    else:
        # automatic_rechunk is not supported
        return sliding_window_view(x, window_shape=window_shape, axis=axis)
