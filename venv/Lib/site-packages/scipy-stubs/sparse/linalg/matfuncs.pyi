# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import _interface

__all__ = ["LinearOperator", "expm", "inv", "spsolve"]

@deprecated("will be removed in SciPy v2.0.0")
class LinearOperator(_interface.LinearOperator): ...

@deprecated("will be removed in SciPy v2.0.0")
def spsolve(A: object, b: object, permc_spec: object = None, use_umfpack: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def expm(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def inv(A: object) -> object: ...
