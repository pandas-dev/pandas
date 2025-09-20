# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["get_lapack_funcs", "qr", "qr_multiply", "rq"]

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = ..., dtype: object = ..., ilp64: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qr(
    a: object,
    overwrite_a: object = ...,
    lwork: object = ...,
    mode: object = ...,
    pivoting: object = ...,
    check_finite: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qr_multiply(
    a: object,
    c: object,
    mode: object = ...,
    pivoting: object = ...,
    conjugate: object = ...,
    overwrite_a: object = ...,
    overwrite_c: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rq(a: object, overwrite_a: object = ..., lwork: object = ..., mode: object = ..., check_finite: object = ...) -> object: ...
