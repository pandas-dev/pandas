from typing_extensions import deprecated

from . import _interpnd

__all__ = ["CloughTocher2DInterpolator", "LinearNDInterpolator"]

@deprecated(
    "Please import `CloughTocher2DInterpolator` from the `scipy.interpolate` namespace; "
    "the `scipy.interpolate.interpnd` namespace is deprecated and will be removed in SciPy 2.0.0."
)
class CloughTocher2DInterpolator(_interpnd.CloughTocher2DInterpolator): ...

@deprecated(
    "Please import `LinearNDInterpolator` from the `scipy.interpolate` namespace; "
    "the `scipy.interpolate.interpnd` namespace is deprecated and will be removed in SciPy 2.0.0."
)
class LinearNDInterpolator(_interpnd.LinearNDInterpolator): ...
