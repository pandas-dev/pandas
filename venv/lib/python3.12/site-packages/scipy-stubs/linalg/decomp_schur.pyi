# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = ["LinAlgError", "eigvals", "get_lapack_funcs", "norm", "rsf2csf", "schur"]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(np.linalg.LinAlgError): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = ..., dtype: object = ..., ilp64: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals(
    a: object, b: object = ..., overwrite_a: object = ..., check_finite: object = ..., homogeneous_eigvals: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def norm(x: object, ord: object = ..., axis: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def schur(
    a: object,
    output: object = ...,
    lwork: object = ...,
    overwrite_a: object = ...,
    sort: object = ...,
    check_finite: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rsf2csf(T: object, Z: object, check_finite: object = ...) -> object: ...
