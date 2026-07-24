# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from .basic import get_lapack_funcs  # pyrefly:ignore[deprecated] # ty:ignore[deprecated]

__all__ = ["LinAlgError", "LinAlgWarning", "get_blas_funcs", "get_lapack_funcs", "norm"]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgWarning(RuntimeWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def norm(a: object, ord: object = None, axis: object = None, keepdims: object = False, check_finite: object = True) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def get_blas_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = "preferred") -> object: ...
