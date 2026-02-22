# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = [
    "LinAlgError",
    "LinAlgWarning",
    "det",
    "get_lapack_funcs",
    "inv",
    "lstsq",
    "matmul_toeplitz",
    "matrix_balance",
    "pinv",
    "pinvh",
    "solve",
    "solve_banded",
    "solve_circulant",
    "solve_toeplitz",
    "solve_triangular",
    "solveh_banded",
]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgWarning(RuntimeWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve(
    a: object,
    b: object,
    lower: object = False,
    overwrite_a: object = False,
    overwrite_b: object = False,
    check_finite: object = True,
    assume_a: object = None,
    transposed: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_triangular(
    a: object,
    b: object,
    trans: object = 0,
    lower: object = False,
    unit_diagonal: object = False,
    overwrite_b: object = False,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_banded(
    l_and_u: object, ab: object, b: object, overwrite_ab: object = False, overwrite_b: object = False, check_finite: object = True
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solveh_banded(
    ab: object,
    b: object,
    overwrite_ab: object = False,
    overwrite_b: object = False,
    lower: object = False,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_toeplitz(c_or_cr: object | tuple[object, object], b: object, check_finite: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_circulant(
    c: object,
    b: object,
    singular: object = "raise",
    tol: object = None,
    caxis: object = -1,
    baxis: object = 0,
    outaxis: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def inv(
    a: object, overwrite_a: object = False, check_finite: object = True, *, assume_a: str | None = None, lower: bool = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def det(a: object, overwrite_a: object = False, check_finite: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lstsq(
    a: object,
    b: object,
    cond: object = None,
    overwrite_a: object = False,
    overwrite_b: object = False,
    check_finite: object = True,
    lapack_driver: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pinv(
    a: object, *, atol: object = None, rtol: object = None, return_rank: object = False, check_finite: object = True
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pinvh(
    a: object,
    atol: object = None,
    rtol: object = None,
    lower: object = True,
    return_rank: object = False,
    check_finite: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def matrix_balance(
    A: object, permute: object = True, scale: object = True, separate: object = False, overwrite_a: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def matmul_toeplitz(
    c_or_cr: object | tuple[object, object], x: object, check_finite: object = False, workers: object = None
) -> object: ...
