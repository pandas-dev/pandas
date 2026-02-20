# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _base, _csc

__all__ = ["csc_matrix", "isspmatrix_csc", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_base.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class csc_matrix(_csc.csc_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_csc(x: object) -> Any: ...
