import numpy as np

from xarray.namedarray.pycompat import array_type


def is_weak_scalar_type(t):
    return isinstance(t, bool | int | float | complex | str | bytes)


def _future_array_api_result_type(*arrays_and_dtypes, xp):
    # fallback implementation for `xp.result_type` with python scalars. Can be removed once a
    # version of the Array API that includes https://github.com/data-apis/array-api/issues/805
    # can be required
    strongly_dtyped = [t for t in arrays_and_dtypes if not is_weak_scalar_type(t)]
    weakly_dtyped = [t for t in arrays_and_dtypes if is_weak_scalar_type(t)]

    if not strongly_dtyped:
        strongly_dtyped = [
            xp.asarray(x) if not isinstance(x, type) else x for x in weakly_dtyped
        ]
        weakly_dtyped = []

    dtype = xp.result_type(*strongly_dtyped)
    if not weakly_dtyped:
        return dtype

    possible_dtypes = {
        complex: "complex64",
        float: "float32",
        int: "int8",
        bool: "bool",
        str: "str",
        bytes: "bytes",
    }
    dtypes = [possible_dtypes.get(type(x), "object") for x in weakly_dtyped]

    return xp.result_type(dtype, *dtypes)


def result_type(*arrays_and_dtypes, xp) -> np.dtype:
    if xp is np or any(
        isinstance(getattr(t, "dtype", t), np.dtype) for t in arrays_and_dtypes
    ):
        return xp.result_type(*arrays_and_dtypes)
    else:
        return _future_array_api_result_type(*arrays_and_dtypes, xp=xp)


def get_array_namespace(*values):
    def _get_single_namespace(x):
        if hasattr(x, "__array_namespace__"):
            return x.__array_namespace__()
        elif isinstance(x, array_type("cupy")):
            # cupy is fully compliant from xarray's perspective, but will not expose
            # __array_namespace__ until at least v14. Special case it for now
            import cupy as cp

            return cp
        else:
            return np

    namespaces = {_get_single_namespace(t) for t in values}
    non_numpy = namespaces - {np}

    if len(non_numpy) > 1:
        names = [module.__name__ for module in non_numpy]
        raise TypeError(f"Mixed array types {names} are not supported.")
    elif non_numpy:
        [xp] = non_numpy
    else:
        xp = np

    return xp


def to_like_array(array, like):
    # Mostly for cupy compatibility, because cupy binary ops require all cupy arrays
    xp = get_array_namespace(like)
    if xp is not np:
        return xp.asarray(array)
    # avoid casting things like pint quantities to numpy arrays
    return array
