from typing import Generic, Protocol, TypeAlias, final, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

__all__ = ["ConvexHull", "Delaunay", "HalfspaceIntersection", "QhullError", "Voronoi", "tsearch"]

_Bool1D: TypeAlias = onp.Array1D[np.bool_]

_IntC1D: TypeAlias = onp.Array1D[np.int32]
_IntC2D: TypeAlias = onp.Array2D[np.int32]
_IntCND: TypeAlias = onp.ArrayND[np.int32]

_Int: TypeAlias = int | np.int32 | np.int64 | np.intp | np.int_
_Int1D: TypeAlias = onp.Array1D[np.int32 | np.int64 | np.intp]
_Int2D: TypeAlias = onp.Array2D[np.int32 | np.int64 | np.intp]

_Float: TypeAlias = float | np.float64
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_Float3D: TypeAlias = onp.Array3D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]

@type_check_only
class DelaunayInfo_t(Protocol):
    # interface (for internal use) of `ctypedef struct DelaunayInfo_t` in `scipy/spatial/_qhull.pxd`
    @property
    def ndim(self, /) -> _Int: ...
    @property
    def npoints(self, /) -> _Int: ...
    @property
    def nsimplex(self, /) -> _Int: ...
    @property
    def points(self, /) -> _Float2D: ...
    @property
    def simplices(self, /) -> _Int2D: ...
    @property
    def neighbors(self, /) -> _Int2D: ...
    @property
    def equations(self, /) -> _Float2D: ...
    @property
    def transform(self, /) -> _Float3D: ...
    @property
    def vertex_to_simplex(self, /) -> _Int1D: ...
    @property
    def paraboloid_scale(self, /) -> _Float: ...
    @property
    def paraboloid_shift(self, /) -> _Float: ...
    @property
    def max_bound(self, /) -> _Float1D: ...
    @property
    def min_bound(self, /) -> _Float1D: ...
    @property
    def vertex_neighbors_indices(self, /) -> _Int1D: ...
    @property
    def vertex_neighbors_indptr(self, /) -> _Int1D: ...

###

class QhullError(RuntimeError): ...

@final
class _Qhull:
    mode_option: bytes
    options: bytes
    furthest_site: bool

    @property
    def ndim(self, /) -> int: ...
    def __init__(
        self,
        /,
        mode_option: bytes,
        points: onp.ToFloat2D,
        options: bytes | None = None,
        required_options: bytes | None = None,
        furthest_site: bool = False,
        incremental: bool = False,
        interior_point: onp.ToFloat1D | None = None,
    ) -> None: ...
    def check_active(self, /) -> None: ...
    def close(self, /) -> None: ...
    def get_points(self, /) -> _Float2D: ...
    def add_points(self, /, points: onp.ToFloat2D, interior_point: onp.ToFloat1D | None = None) -> None: ...
    def get_paraboloid_shift_scale(self, /) -> tuple[float, float]: ...
    def volume_area(self, /) -> tuple[float, float]: ...
    def triangulate(self, /) -> None: ...
    def get_simplex_facet_array(self, /) -> tuple[_IntC2D, _IntC2D, _Float2D, _IntC2D, _IntC1D]: ...
    def get_hull_points(self, /) -> _Float2D: ...
    def get_hull_facets(self, /) -> tuple[list[list[int]], _Float2D]: ...
    def get_voronoi_diagram(self, /) -> tuple[_Float2D, _IntC2D, list[list[int]], list[list[int]], _Int1D]: ...
    def get_extremes_2d(self, /) -> _IntC1D: ...

_QT = TypeVar("_QT", bound=_Qhull, default=_Qhull)

class _QhullUser(Generic[_QT]):
    _qhull: _QT | None = None
    ndim: int
    npoints: int
    min_bound: _Float1D
    max_bound: _Float1D

    def __init__(self, /, qhull: _QT, incremental: bool = False) -> None: ...
    def __del__(self, /) -> None: ...
    def _update(self, /, qhull: _QT) -> None: ...
    def _add_points(
        self, /, points: onp.ToFloat2D, restart: bool = False, interior_point: onp.ToFloat1D | None = None
    ) -> None: ...
    def close(self, /) -> None: ...

class Delaunay(_QhullUser[_Qhull]):
    furthest_site: bool
    paraboloid_scale: float
    paraboloid_shift: float
    simplices: _IntC2D
    neighbors: _IntC2D
    equations: _Float2D
    coplanar: _IntC2D
    good: _IntC1D
    nsimplex: int
    vertices: _Float2D

    @property
    def points(self, /) -> _Float2D: ...
    @property
    def transform(self, /) -> _Float3D: ...
    @property
    def vertex_to_simplex(self, /) -> _IntC1D: ...
    @property
    def vertex_neighbor_vertices(self, /) -> tuple[_IntC1D, _IntC1D]: ...
    @property
    def convex_hull(self, /) -> _IntC2D: ...

    #
    def __init__(
        self, /, points: onp.ToFloat2D, furthest_site: bool = False, incremental: bool = False, qhull_options: str | None = None
    ) -> None: ...
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...
    def find_simplex(self, /, xi: onp.ToFloatND, bruteforce: bool = False, tol: float | None = None) -> _IntCND: ...
    def plane_distance(self, /, xi: onp.ToFloatND) -> _Float2D: ...
    def lift_points(self, /, x: onp.ToFloatND) -> _FloatND: ...

class ConvexHull(_QhullUser[_Qhull]):
    simplices: _IntC2D
    neighbors: _IntC2D
    equations: _Float2D
    coplanar: _IntC2D
    good: _Bool1D | None
    volume: float
    area: float
    nsimplex: int

    @property
    def points(self, /) -> _Float2D: ...
    @property
    def vertices(self, /) -> _IntC2D: ...

    #
    def __init__(self, /, points: onp.ToFloat2D, incremental: bool = False, qhull_options: str | None = None) -> None: ...
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...

class Voronoi(_QhullUser):
    vertices: _Float2D
    ridge_points: _IntC2D
    ridge_vertices: list[list[int]]
    regions: list[list[int]]
    point_region: _Int1D
    furthest_site: bool

    @property
    def points(self, /) -> _Float2D: ...
    @property
    def ridge_dict(self, /) -> dict[tuple[int, int], list[int]]: ...

    #
    def __init__(
        self, /, points: onp.ToFloat2D, furthest_site: bool = False, incremental: bool = False, qhull_options: str | None = None
    ) -> None: ...
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...

class HalfspaceIntersection(_QhullUser):
    interior_point: _Float1D
    intersections: _Float2D
    dual_facets: list[list[int]]
    dual_equations: _Float2D
    dual_points: _Float2D
    dual_volume: float
    dual_area: float
    ndim: int
    nineq: int

    @property
    def halfspaces(self, /) -> _Float2D: ...
    @property
    def dual_vertices(self, /) -> _Int1D: ...

    #
    def __init__(
        self,
        /,
        halfspaces: onp.ToFloat2D,
        interior_point: onp.ToFloat1D,
        incremental: bool = False,
        qhull_options: str | None = None,
    ) -> None: ...
    def add_halfspaces(self, /, halfspaces: onp.ToFloat2D, restart: bool = False) -> None: ...

def tsearch(tri: Delaunay, xi: onp.ToFloatND) -> _IntCND: ...
