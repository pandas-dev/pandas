# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = ["LinAlgError", "eigvals", "get_lapack_funcs", "norm", "rsf2csf", "schur"]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(np.linalg.LinAlgError): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals(
    a: object, b: object = None, overwrite_a: object = False, check_finite: object = True, homogeneous_eigvals: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def norm(x: object, ord: object = None, axis: object = None, keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def schur(
    a: object,
    output: object = "real",
    lwork: object = None,
    overwrite_a: object = False,
    sort: object = None,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rsf2csf(T: object, Z: object, check_finite: object = True) -> object: ...
