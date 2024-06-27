import numpy as np


def is_weak_scalar_type(t):
    return isinstance(t, (bool, int, float, complex, str, bytes))


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
