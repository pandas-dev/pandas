from array_api_compat import (
    common as common,
    device as device,  # pyright: ignore[reportUnknownVariableType]
    get_namespace as get_namespace,
    is_array_api_obj as is_array_api_obj,  # pyright: ignore[reportUnknownVariableType]
    is_array_api_strict_namespace as is_array_api_strict_namespace,
    is_cupy_array as is_cupy_array,
    is_cupy_namespace as is_cupy_namespace,
    is_dask_array as is_dask_array,  # pyright: ignore[reportUnknownVariableType]
    is_dask_namespace as is_dask_namespace,
    is_jax_array as is_jax_array,  # pyright: ignore[reportUnknownVariableType]
    is_jax_namespace as is_jax_namespace,
    is_lazy_array as is_lazy_array,
    is_ndonnx_array as is_ndonnx_array,  # pyright: ignore[reportUnknownVariableType]
    is_ndonnx_namespace as is_ndonnx_namespace,
    is_numpy_array as is_numpy_array,
    is_numpy_namespace as is_numpy_namespace,
    is_pydata_sparse_array as is_pydata_sparse_array,  # pyright: ignore[reportUnknownVariableType]
    is_pydata_sparse_namespace as is_pydata_sparse_namespace,
    is_torch_array as is_torch_array,  # pyright: ignore[reportUnknownVariableType]
    is_torch_namespace as is_torch_namespace,
    is_writeable_array as is_writeable_array,
    size as size,
    to_device as to_device,
)

from ._array_api import array_namespace as scipy_array_namespace

array_namespace = scipy_array_namespace
