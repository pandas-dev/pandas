# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from . import _base

__all__ = ["SparseEfficiencyWarning"]

@deprecated("will be removed in SciPy v2.0.0")
class SparseEfficiencyWarning(_base.SparseEfficiencyWarning): ...
