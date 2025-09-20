# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["sepfir2d"]

@deprecated("will be removed in SciPy v2.0.0")
def sepfir2d(input: object, hrow: object, hcol: object) -> object: ...
