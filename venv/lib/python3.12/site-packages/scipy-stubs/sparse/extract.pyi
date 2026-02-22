# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _coo

__all__ = ["coo_matrix", "find", "tril", "triu"]

@deprecated("will be removed in SciPy v2.0.0")
class coo_matrix(_coo.coo_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def find(A: object) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def tril(A: object, k: object = 0, format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def triu(A: object, k: object = 0, format: object = None) -> Any: ...
