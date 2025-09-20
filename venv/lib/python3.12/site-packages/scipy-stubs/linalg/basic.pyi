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
def get_lapack_funcs(names: object, arrays: object = ..., dtype: object = ..., ilp64: object = ...) -> object: ...
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
def solve_triangular(
    a: object,
    b: object,
    trans: object = ...,
    lower: object = ...,
    unit_diagonal: object = ...,
    overwrite_b: object = ...,
    check_finite: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_banded(
    l_and_u: object, ab: object, b: object, overwrite_ab: object = ..., overwrite_b: object = ..., check_finite: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solveh_banded(
    ab: object, b: object, overwrite_ab: object = ..., overwrite_b: object = ..., lower: object = ..., check_finite: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_toeplitz(c_or_cr: object | tuple[object, object], b: object, check_finite: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def solve_circulant(
    c: object,
    b: object,
    singular: object = ...,
    tol: object = ...,
    caxis: object = ...,
    baxis: object = ...,
    outaxis: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def inv(a: object, overwrite_a: object = ..., check_finite: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def det(a: object, overwrite_a: object = ..., check_finite: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lstsq(
    a: object,
    b: object,
    cond: object = ...,
    overwrite_a: object = ...,
    overwrite_b: object = ...,
    check_finite: object = ...,
    lapack_driver: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pinv(
    a: object, *, atol: object = ..., rtol: object = ..., return_rank: object = ..., check_finite: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pinvh(
    a: object, atol: object = ..., rtol: object = ..., lower: object = ..., return_rank: object = ..., check_finite: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def matrix_balance(
    A: object, permute: object = ..., scale: object = ..., separate: object = ..., overwrite_a: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def matmul_toeplitz(
    c_or_cr: object | tuple[object, object], x: object, check_finite: object = ..., workers: object = ...
) -> object: ...
