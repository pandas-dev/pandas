# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = [
    "coshm",
    "cosm",
    "expm",
    "expm_cond",
    "expm_frechet",
    "fractional_matrix_power",
    "funm",
    "inv",
    "khatri_rao",
    "logm",
    "norm",
    "rsf2csf",
    "schur",
    "signm",
    "sinhm",
    "sinm",
    "solve",
    "sqrtm",
    "svd",
    "tanhm",
    "tanm",
]

@deprecated("will be removed in SciPy v2.0.0")
def norm(a: object, ord: object = ..., axis: object = ..., keepdims: object = ..., check_finite: object = ...) -> object: ...
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
def rsf2csf(T: object, Z: object, check_finite: object = ...) -> object: ...
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
def inv(a: object, overwrite_a: object = ..., check_finite: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve(
    a: object,
    b: object,
    lower: object = ...,
    overwrite_a: object = ...,
    overwrite_b: object = ...,
    check_finite: object = ...,
    assume_a: object = ...,
    transposed: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fractional_matrix_power(A: object, t: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def logm(A: object, disp: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def expm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cosm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sinm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tanm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def coshm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sinhm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tanhm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def funm(A: object, func: object, disp: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def signm(A: object, disp: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def khatri_rao(a: object, b: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sqrtm(A: object, disp: object = ..., blocksize: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def expm_frechet(
    A: object, E: object, method: object = ..., compute_expm: object = ..., check_finite: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def expm_cond(A: object, check_finite: object = ...) -> object: ...
