# This module is not meant for public use and will be removed in SciPy v2.0.0.

from types import ModuleType
from typing_extensions import deprecated

from . import _eigen

__all__ = ["ArpackError", "ArpackNoConvergence", "eigs", "eigsh", "lobpcg", "svds", "test"]

test: ModuleType

@deprecated("will be removed in SciPy v2.0.0")
class ArpackError(_eigen.ArpackError): ...

@deprecated("will be removed in SciPy v2.0.0")
class ArpackNoConvergence(_eigen.ArpackNoConvergence): ...

@deprecated("will be removed in SciPy v2.0.0")
def eigs(
    A: object,
    k: object = ...,
    M: object = ...,
    sigma: object = ...,
    which: object = ...,
    v0: object = ...,
    ncv: object = ...,
    maxiter: object = ...,
    tol: object = ...,
    return_eigenvectors: object = ...,
    Minv: object = ...,
    OPinv: object = ...,
    OPpart: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigsh(
    A: object,
    k: object = ...,
    M: object = ...,
    sigma: object = ...,
    which: object = ...,
    v0: object = ...,
    ncv: object = ...,
    maxiter: object = ...,
    tol: object = ...,
    return_eigenvectors: object = ...,
    Minv: object = ...,
    OPinv: object = ...,
    mode: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lobpcg(
    A: object,
    X: object,
    B: object = ...,
    M: object = ...,
    Y: object = ...,
    tol: object = ...,
    maxiter: object = ...,
    largest: object = ...,
    verbosityLevel: object = ...,
    retLambdaHistory: object = ...,
    retResidualNormsHistory: object = ...,
    restartControl: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def svds(
    A: object,
    k: object = ...,
    ncv: object = ...,
    tol: object = ...,
    which: object = ...,
    v0: object = ...,
    maxiter: object = ...,
    return_singular_vectors: object = ...,
    solver: object = ...,
    rng: object = None,
    options: object = None,
    *,
    random_state: object = None,
) -> object: ...
