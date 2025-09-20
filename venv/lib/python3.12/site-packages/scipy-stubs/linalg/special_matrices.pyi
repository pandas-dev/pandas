# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = [
    "block_diag",
    "circulant",
    "companion",
    "convolution_matrix",
    "dft",
    "fiedler",
    "fiedler_companion",
    "hadamard",
    "hankel",
    "helmert",
    "hilbert",
    "invhilbert",
    "invpascal",
    "kron",
    "leslie",
    "pascal",
    "toeplitz",
]

@deprecated("will be removed in SciPy v2.0.0")
def toeplitz(c: object, r: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circulant(c: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hankel(c: object, r: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hadamard(n: object, dtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def leslie(f: object, s: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kron(a: object, b: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def block_diag(*arrs: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def companion(a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def helmert(n: object, full: bool = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hilbert(n: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invhilbert(n: object, exact: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pascal(n: object, kind: object = ..., exact: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invpascal(n: object, kind: object = ..., exact: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def dft(n: object, scale: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fiedler(a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fiedler_companion(a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolution_matrix(a: object, n: object, mode: object = ...) -> object: ...
