# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _bsplines, _fitpack2, _interpolate, _rgi

__all__ = [
    "BPoly",
    "BSpline",
    "NdPPoly",
    "PPoly",
    "RectBivariateSpline",
    "RegularGridInterpolator",
    "interp1d",
    "interp2d",
    "interpn",
    "lagrange",
    "make_interp_spline",
]

# _bsplines

@deprecated("will be removed in SciPy v2.0.0")
class BSpline(_bsplines.BSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
def make_interp_spline(
    x: object, y: object, k: object = 3, t: object = None, bc_type: object = None, axis: object = 0, check_finite: object = True
) -> Any: ...

# _fitpack2

@deprecated("will be removed in SciPy v2.0.0")
class RectBivariateSpline(_fitpack2.RectBivariateSpline): ...

# _interpolate

@deprecated("will be removed in SciPy v2.0.0")
class BPoly(_interpolate.BPoly): ...

@deprecated("will be removed in SciPy v2.0.0")
class NdPPoly(_interpolate.NdPPoly): ...

@deprecated("will be removed in SciPy v2.0.0")
class PPoly(_interpolate.PPoly): ...

@deprecated("will be removed in SciPy v2.0.0")
class interp1d(_interpolate.interp1d): ...

@deprecated("will be removed in SciPy v2.0.0")
class interp2d(_interpolate.interp2d): ...  # pyright: ignore[reportDeprecated]  # ty: ignore[deprecated]

@deprecated("will be removed in SciPy v2.0.0")
def lagrange(x: object, w: object) -> Any: ...

# _rgi

@deprecated("will be removed in SciPy v2.0.0")
class RegularGridInterpolator(_rgi.RegularGridInterpolator): ...

@deprecated("will be removed in SciPy v2.0.0")
def interpn(
    points: object, values: object, xi: object, method: object = "linear", bounds_error: object = True, fill_value: object = ...
) -> Any: ...
