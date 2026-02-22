# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _polyint

__all__ = [
    "BarycentricInterpolator",
    "KroghInterpolator",
    "approximate_taylor_polynomial",
    "barycentric_interpolate",
    "krogh_interpolate",
]

@deprecated("will be removed in SciPy v2.0.0")
class KroghInterpolator(_polyint.KroghInterpolator): ...

@deprecated("will be removed in SciPy v2.0.0")
class BarycentricInterpolator(_polyint.BarycentricInterpolator): ...

@deprecated("will be removed in SciPy v2.0.0")
def krogh_interpolate(xi: object, yi: object, x: object, der: object = 0, axis: object = 0) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def approximate_taylor_polynomial(f: object, x: object, degree: object, scale: object, order: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def barycentric_interpolate(xi: object, yi: object, x: object, axis: int = 0, *, der: object = 0, rng: object = None) -> Any: ...
