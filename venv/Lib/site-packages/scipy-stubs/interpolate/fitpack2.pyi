# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import _fitpack2

__all__ = [
    "BivariateSpline",
    "InterpolatedUnivariateSpline",
    "LSQBivariateSpline",
    "LSQSphereBivariateSpline",
    "LSQUnivariateSpline",
    "RectBivariateSpline",
    "RectSphereBivariateSpline",
    "SmoothBivariateSpline",
    "SmoothSphereBivariateSpline",
    "UnivariateSpline",
]

@deprecated("will be removed in SciPy v2.0.0")
class UnivariateSpline(_fitpack2.UnivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class InterpolatedUnivariateSpline(_fitpack2.InterpolatedUnivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class LSQUnivariateSpline(_fitpack2.LSQUnivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class BivariateSpline(_fitpack2.BivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class SmoothBivariateSpline(_fitpack2.SmoothBivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class LSQBivariateSpline(_fitpack2.LSQBivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class RectBivariateSpline(_fitpack2.RectBivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class SmoothSphereBivariateSpline(_fitpack2.SmoothSphereBivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class LSQSphereBivariateSpline(_fitpack2.LSQSphereBivariateSpline): ...

@deprecated("will be removed in SciPy v2.0.0")
class RectSphereBivariateSpline(_fitpack2.RectSphereBivariateSpline): ...
