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
def get_lapack_funcs(names: object, arrays: object = ..., dtype: object = ..., ilp64: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def norm(a: object, ord: object = ..., axis: object = ..., keepdims: object = ..., check_finite: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cdf2rdf(w: object, v: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eig(
    a: object,
    b: object = ...,
    left: object = ...,
    right: object = ...,
    overwrite_a: object = ...,
    overwrite_b: object = ...,
    check_finite: object = ...,
    homogeneous_eigvals: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigh(
    a: object,
    b: object = ...,
    *,
    lower: object = ...,
    eigvals_only: object = ...,
    overwrite_a: object = ...,
    overwrite_b: object = ...,
    type: object = ...,
    check_finite: object = ...,
    subset_by_index: object = ...,
    subset_by_value: object = ...,
    driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eig_banded(
    a_band: object,
    lower: object = ...,
    eigvals_only: object = ...,
    overwrite_a_band: object = ...,
    select: object = ...,
    select_range: object = ...,
    max_ev: int = ...,
    check_finite: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigh_tridiagonal(
    d: object,
    e: object,
    eigvals_only: object = ...,
    select: object = ...,
    select_range: object = ...,
    check_finite: object = ...,
    tol: object = ...,
    lapack_driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals(
    a: object, b: object = ..., overwrite_a: object = ..., check_finite: object = ..., homogeneous_eigvals: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvals_banded(
    a_band: object,
    lower: object = ...,
    overwrite_a_band: object = ...,
    select: object = ...,
    select_range: object = ...,
    check_finite: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvalsh_tridiagonal(
    d: object,
    e: object,
    select: object = ...,
    select_range: object = ...,
    check_finite: object = ...,
    tol: object = ...,
    lapack_driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigvalsh(
    a: object,
    b: object = ...,
    *,
    lower: object = ...,
    overwrite_a: object = ...,
    overwrite_b: object = ...,
    type: object = ...,
    check_finite: object = ...,
    subset_by_index: object = ...,
    subset_by_value: object = ...,
    driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hessenberg(a: object, calc_q: object = ..., overwrite_a: object = ..., check_finite: object = ...) -> object: ...
