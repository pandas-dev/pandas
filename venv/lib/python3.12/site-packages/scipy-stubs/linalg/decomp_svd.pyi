# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = ["LinAlgError", "diagsvd", "get_lapack_funcs", "null_space", "orth", "subspace_angles", "svd", "svdvals"]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(np.linalg.LinAlgError): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = ..., dtype: object = ..., ilp64: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def diagsvd(s: object, M: object, N: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def null_space(
    A: object, rcond: object = ..., *, overwrite_a: bool = False, check_finite: bool = True, lapack_driver: str = "gesdd"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def orth(A: object, rcond: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def subspace_angles(A: object, B: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def svd(
    a: object,
    full_matrices: object = ...,
    compute_uv: object = ...,
    overwrite_a: object = ...,
    check_finite: object = ...,
    lapack_driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def svdvals(a: object, overwrite_a: object = ..., check_finite: object = ...) -> object: ...
