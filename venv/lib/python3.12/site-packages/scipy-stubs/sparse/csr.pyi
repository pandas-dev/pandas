# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _csr, _matrix

__all__ = ["csr_matrix", "isspmatrix_csr", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class csr_matrix(_csr.csr_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_csr(x: object) -> Any: ...
