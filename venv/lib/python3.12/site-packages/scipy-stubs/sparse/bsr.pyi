# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _bsr, _matrix

__all__ = ["bsr_matrix", "isspmatrix_bsr", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class bsr_matrix(_bsr.bsr_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_bsr(x: object) -> Any: ...
