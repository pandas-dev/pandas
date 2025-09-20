# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import _interface

__all__ = ["LinearOperator", "aslinearoperator"]

@deprecated("will be removed in SciPy v2.0.0")
class LinearOperator(_interface.LinearOperator): ...

@deprecated("will be removed in SciPy v2.0.0")
def aslinearoperator(A: object) -> object: ...
