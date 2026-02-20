# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = [
    "LinAlgError",
    "cdf2rdf",
    "eig",
    "eig_banded",
    "eigh",
    "eigh_tridiagonal",
    "eigvals",
    "eigvals_banded",
    "eigvalsh",
    "eigvalsh_tridiagonal",
    "get_lapack_funcs",
    "hessenberg",
    "norm",
]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(np.linalg.LinAlgError): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def norm(a: object, ord: object = None, axis: object = None, keepdims: object = False, check_finite: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cdf2rdf(w: object, v: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eig(
    a: object,
    b: object = None,
    left: object = False,
    right: object = True,
    overwrite_a: object = False,
    overwrite_b: object = False,
    check_finite: object = True,
    homogeneous_eigvals: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigh(
    a: object,
    b: object = None,
    *,
    lower: object = True,
    eigvals_only: object = False,
    overwrite_a: object = False,
    overwrite_b: object = False,
    type: object = 1,
    check_finite: object = True,
    subset_by_index: object = None,
    subset_by_value: object = None,
    driver: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eig_banded(
    a_band: object,
    lower: object = False,
    eigvals_only: object = False,
    overwrite_a_band: object = False,
    select: object = "a",
    select_range: object = None,
    max_ev: int = 0,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigh_tridiagonal(
    d: object,
    e: object,
    eigvals_only: object = False,
    select: object = "a",
    select_range: object = None,
    check_finite: object = True,
    tol: object = 0.0,
    lapack_driver: object = "auto",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals(
    a: object, b: object = None, overwrite_a: object = False, check_finite: object = True, homogeneous_eigvals: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals_banded(
    a_band: object,
    lower: object = False,
    overwrite_a_band: object = False,
    select: object = "a",
    select_range: object = None,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvalsh_tridiagonal(
    d: object,
    e: object,
    select: object = "a",
    select_range: object = None,
    check_finite: object = True,
    tol: object = 0.0,
    lapack_driver: object = "auto",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvalsh(
    a: object,
    b: object = None,
    *,
    lower: object = True,
    overwrite_a: object = False,
    overwrite_b: object = False,
    type: object = 1,
    check_finite: object = True,
    subset_by_index: object = None,
    subset_by_value: object = None,
    driver: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hessenberg(a: object, calc_q: object = False, overwrite_a: object = False, check_finite: object = True) -> object: ...
