# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _dia, _matrix

__all__ = ["dia_matrix", "isspmatrix_dia", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class dia_matrix(_dia.dia_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_dia(x: object) -> Any: ...
