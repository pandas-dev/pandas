# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.
from typing_extensions import deprecated

from ._qhull import (
    ConvexHull as _ConvexHull,
    Delaunay as _Delaunay,
    HalfspaceIntersection as _HalfspaceIntersection,
    QhullError as _QhullError,
    Voronoi as _Voronoi,
)

__all__ = ["ConvexHull", "Delaunay", "HalfspaceIntersection", "QhullError", "Voronoi", "tsearch"]

@deprecated("will be removed in SciPy v2.0.0")
class QhullError(_QhullError): ...

@deprecated("will be removed in SciPy v2.0.0")
class ConvexHull(_ConvexHull): ...

@deprecated("will be removed in SciPy v2.0.0")
class Delaunay(_Delaunay): ...

@deprecated("will be removed in SciPy v2.0.0")
class HalfspaceIntersection(_HalfspaceIntersection): ...

@deprecated("will be removed in SciPy v2.0.0")
class Voronoi(_Voronoi): ...

@deprecated("will be removed in SciPy v2.0.0")
def tsearch(tri: _Delaunay, xi: object) -> object: ...
