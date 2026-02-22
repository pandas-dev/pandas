# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["get_lapack_funcs", "qr", "qr_multiply", "rq"]

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qr(
    a: object,
    overwrite_a: object = False,
    lwork: object = None,
    mode: object = "full",
    pivoting: object = False,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qr_multiply(
    a: object,
    c: object,
    mode: object = "right",
    pivoting: object = False,
    conjugate: object = False,
    overwrite_a: object = False,
    overwrite_c: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rq(
    a: object, overwrite_a: object = False, lwork: object = None, mode: object = "full", check_finite: object = True
) -> object: ...
