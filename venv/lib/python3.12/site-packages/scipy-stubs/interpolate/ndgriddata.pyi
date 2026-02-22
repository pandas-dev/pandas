# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from . import _interpnd, _ndgriddata

__all__ = ["CloughTocher2DInterpolator", "LinearNDInterpolator", "NearestNDInterpolator", "griddata"]

# interpnd

@deprecated("will be removed in SciPy v2.0.0")
class LinearNDInterpolator(_interpnd.LinearNDInterpolator): ...

# _ndgriddata

@deprecated("will be removed in SciPy v2.0.0")
class CloughTocher2DInterpolator(_ndgriddata.CloughTocher2DInterpolator): ...

@deprecated("will be removed in SciPy v2.0.0")
class NearestNDInterpolator(_ndgriddata.NearestNDInterpolator): ...

@deprecated("will be removed in SciPy v2.0.0")
def griddata(
    points: object, values: object, xi: object, method: object = "linear", fill_value: float = ..., rescale: object = False
) -> object: ...
