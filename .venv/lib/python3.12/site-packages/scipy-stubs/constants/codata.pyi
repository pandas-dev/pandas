# This file is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

# Constants can't be tagged with `deprecated` so no point to stub these out again.
from ._codata import c, k, physical_constants

__all__ = ["ConstantWarning", "c", "find", "k", "physical_constants", "precision", "unit", "value"]

@deprecated("will be removed in SciPy v2.0.0")
class ConstantWarning(DeprecationWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def find(sub: str | None = None, disp: bool = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def precision(key: str) -> float: ...
@deprecated("will be removed in SciPy v2.0.0")
def unit(key: str) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def value(key: str) -> float: ...
