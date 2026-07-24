# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import _rbf

__all__ = ["Rbf"]

@deprecated("will be removed in SciPy v2.0.0")
class Rbf(_rbf.Rbf): ...
