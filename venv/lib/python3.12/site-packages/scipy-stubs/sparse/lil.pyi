# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _lil

__all__ = ["isspmatrix_lil", "lil_array", "lil_matrix"]

@deprecated("will be removed in SciPy v2.0.0")
class lil_array(_lil.lil_array): ...

@deprecated("will be removed in SciPy v2.0.0")
class lil_matrix(_lil.lil_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix_lil(x: object) -> Any: ...
