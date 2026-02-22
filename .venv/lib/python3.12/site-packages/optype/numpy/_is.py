import sys
from typing import Any

if sys.version_info >= (3, 13):
    from typing import TypeIs, TypeVar
else:
    from typing_extensions import TypeIs, TypeVar

import numpy as np

from ._array import Array0D, Array1D, Array2D, Array3D, ArrayND
from ._dtype import ToDType

__all__ = [
    "is_array_0d",
    "is_array_1d",
    "is_array_2d",
    "is_array_3d",
    "is_array_nd",
    "is_dtype",
    "is_sctype",
]


def __dir__() -> list[str]:
    return __all__


DTypeT = TypeVar("DTypeT", bound=np.dtype[Any])
ScalarT = TypeVar("ScalarT", bound=np.generic, default=Any)


def _issubdtype(
    x: ToDType[Any],
    /,
    dtype: ToDType[ScalarT] | None,
) -> TypeIs[np.dtype[ScalarT]]:
    """Checks if `x` is a subdtype of `dtype`, or if `dtype` is None."""

    # NOTE(jorenham): The `bool()` works around an the awkward return type annotation
    # of `numpy.issubdtype` on `numpy==2.2.*` (mea culpa).
    return dtype is None or bool(np.issubdtype(x, dtype))


def is_dtype(
    x: object,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[np.dtype[ScalarT]]:
    return isinstance(x, np.dtype) and _issubdtype(x, dtype)


def is_sctype(
    x: object,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[type[ScalarT]]:
    return isinstance(x, type) and issubclass(x, np.generic) and _issubdtype(x, dtype)


def is_array_nd(
    a: Any,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[ArrayND[ScalarT]]:
    """Checks if `a` is a `ndarray` of the given dtype (defaults to `generic`)."""
    return isinstance(a, np.ndarray) and _issubdtype(a.dtype, dtype)


def is_array_0d(
    a: Any,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[Array0D[ScalarT]]:
    """Checks if `a` is a 0-d `ndarray` of the given dtype (defaults to `generic`)."""
    return is_array_nd(a, dtype) and a.ndim == 0


def is_array_1d(
    a: Any,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[Array1D[ScalarT]]:
    """Checks if `a` is a 1-d `ndarray` of the given dtype (defaults to `generic`)."""
    return is_array_nd(a, dtype) and a.ndim == 1


def is_array_2d(
    a: Any,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[Array2D[ScalarT]]:
    """Checks if `a` is a 2-d `ndarray` of the given dtype (defaults to `generic`)."""
    return is_array_nd(a, dtype) and a.ndim == 2


def is_array_3d(
    a: Any,
    /,
    dtype: ToDType[ScalarT] | None = None,
) -> TypeIs[Array3D[ScalarT]]:
    """Checks if `a` is a 3-d `ndarray` of the given dtype (defaults to `generic`)."""
    return is_array_nd(a, dtype) and a.ndim == 3
