# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = [
    "BroydenFirst",
    "InverseJacobian",
    "KrylovJacobian",
    "anderson",
    "broyden1",
    "broyden2",
    "diagbroyden",
    "excitingmixing",
    "linearmixing",
    "newton_krylov",
]

@deprecated("will be removed in SciPy v2.0.0")
class BroydenFirst:
    def __init__(self, /, alpha: object = ..., reduction_method: object = ..., max_rank: object = ...) -> None: ...
    def todense(self, /) -> object: ...
    def solve(self, /, f: object, tol: object = ...) -> object: ...
    def matvec(self, /, f: object) -> object: ...
    def rsolve(self, /, f: object, tol: object = ...) -> object: ...
    def rmatvec(self, /, f: object) -> object: ...
    def setup(self, /, x: object, F: object, func: object) -> None: ...

@deprecated("will be removed in SciPy v2.0.0")
class InverseJacobian:
    def __init__(self, /, jacobian: object) -> None: ...
    @property
    def shape(self, /) -> object: ...
    @property
    def dtype(self, /) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class KrylovJacobian:
    def __init__(
        self,
        /,
        rdiff: object = ...,
        method: object = ...,
        inner_maxiter: object = ...,
        inner_M: object = ...,
        outer_k: object = ...,
        **kw: object,
    ) -> None: ...
    def matvec(self, /, v: object) -> object: ...
    def solve(self, /, rhs: object, tol: object = ...) -> object: ...
    def update(self, /, x: object, f: object) -> object: ...
    def setup(self, /, x: object, f: object, func: object) -> object: ...

# Deprecated
broyden1: object
broyden2: object
anderson: object
linearmixing: object
diagbroyden: object
excitingmixing: object
newton_krylov: object
