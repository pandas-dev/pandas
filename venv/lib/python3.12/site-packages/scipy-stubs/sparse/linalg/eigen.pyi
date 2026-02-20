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
    k: object = 6,
    M: object = None,
    sigma: object = None,
    which: object = "LM",
    v0: object = None,
    ncv: object = None,
    maxiter: object = None,
    tol: object = 0,
    return_eigenvectors: object = True,
    Minv: object = None,
    OPinv: object = None,
    OPpart: object = None,
    rng: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def eigsh(
    A: object,
    k: object = 6,
    M: object = None,
    sigma: object = None,
    which: object = "LM",
    v0: object = None,
    ncv: object = None,
    maxiter: object = None,
    tol: object = 0,
    return_eigenvectors: object = True,
    Minv: object = None,
    OPinv: object = None,
    mode: object = "normal",
    rng: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lobpcg(
    A: object,
    X: object,
    B: object = None,
    M: object = None,
    Y: object = None,
    tol: object = None,
    maxiter: object = None,
    largest: object = True,
    verbosityLevel: object = 0,
    retLambdaHistory: object = False,
    retResidualNormsHistory: object = False,
    restartControl: object = 20,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def svds(
    A: object,
    k: object = 6,
    ncv: object = None,
    tol: object = 0,
    which: object = "LM",
    v0: object = None,
    maxiter: object = None,
    return_singular_vectors: object = True,
    solver: object = "arpack",
    rng: object = None,
    options: object = None,
    *,
    random_state: object = None,
) -> object: ...
