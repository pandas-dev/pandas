# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _base, _coo, _matrix

__all__ = ["SparseEfficiencyWarning", "coo_matrix", "isspmatrix_coo", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class SparseEfficiencyWarning(_base.SparseEfficiencyWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class coo_matrix(_coo.coo_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_coo(x: object) -> Any: ...
