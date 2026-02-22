# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from . import _base, _matrix

__all__ = ["SparseEfficiencyWarning", "SparseWarning", "issparse", "isspmatrix", "spmatrix"]

@deprecated("will be removed in SciPy v2.0.0")
class SparseWarning(_base.SparseWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
class SparseEfficiencyWarning(_base.SparseEfficiencyWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
class spmatrix(_matrix.spmatrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def issparse(x: object) -> bool: ...
@deprecated("will be removed in SciPy v2.0.0")
def isspmatrix(x: object) -> bool: ...
