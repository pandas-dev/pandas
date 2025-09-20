# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _dok, _matrix

__all__ = ["dok_matrix", "isspmatrix_dok", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class dok_matrix(_dok.dok_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_dok(x: object) -> Any: ...
