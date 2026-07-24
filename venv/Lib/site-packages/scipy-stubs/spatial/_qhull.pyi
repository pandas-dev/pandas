from typing import Never, Protocol, final, overload, type_check_only
from typing_extensions import TypeVar, deprecated

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ConvexHull", "Delaunay", "HalfspaceIntersection", "QhullError", "Voronoi", "tsearch"]

###

# on numpy<2.1 pyright reports false positive incompatible overload errors for `plane_distance` and `lift_points` in `Delaunay`
# pyright: reportOverlappingOverload=false

# workaround for type-checkers that don't comply to the overload spec
type _JustAnyShape = tuple[Never, Never, Never, Never]
type _ToArrayStrictND = onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

@type_check_only
class DelaunayInfo_t(Protocol):
    # interface (for internal use) of `ctypedef struct DelaunayInfo_t` in `scipy/spatial/_qhull.pxd`
    @property
    def ndim(self, /) -> int: ...
    @property
    def npoints(self, /) -> int: ...
    @property
    def nsimplex(self, /) -> int: ...
    @property
    def points(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def simplices(self, /) -> onp.Array2D[np.int32]: ...
    @property
    def neighbors(self, /) -> onp.Array2D[np.int32]: ...
    @property
    def equations(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def transform(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def vertex_to_simplex(self, /) -> onp.Array1D[np.int32]: ...
    @property
    def paraboloid_scale(self, /) -> float: ...
    @property
    def paraboloid_shift(self, /) -> float: ...
    @property
    def max_bound(self, /) -> onp.Array1D[np.float64]: ...
    @property
    def min_bound(self, /) -> onp.Array1D[np.float64]: ...

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
    def get_points(self, /) -> onp.Array2D[np.float64]: ...
    def add_points(self, /, points: onp.ToFloat2D, interior_point: onp.ToFloat1D | None = None) -> None: ...
    def get_paraboloid_shift_scale(self, /) -> tuple[float, float]: ...
    def volume_area(self, /) -> tuple[float, float]: ...
    def triangulate(self, /) -> None: ...
    def get_simplex_facet_array(
        self, /
    ) -> tuple[
        onp.Array2D[np.int32],  # facets
        onp.Array2D[np.int32],  # neighbors
        onp.Array2D[np.float64],  # equations
        onp.Array2D[np.int32],  # coplanar
        onp.Array1D[np.int32],  # good
    ]: ...
    def get_hull_points(self, /) -> onp.Array2D[np.float64]: ...
    def get_hull_facets(
        self, /
    ) -> tuple[
        list[list[int]],  # facets
        onp.Array2D[np.float64],  # equations
    ]: ...
    def get_voronoi_diagram(
        self, /
    ) -> tuple[
        onp.Array2D[np.float64],  # vertices
        onp.Array2D[np.int32],  # ridge_points
        list[list[int]],  # ridge_vertices
        list[list[int]],  # regions
        onp.Array1D[np.intp],  # point_region
    ]: ...
    def get_extremes_2d(self, /) -> onp.Array1D[np.int32]: ...

class _QhullUser:
    _qhull: _Qhull | None = None
    ndim: int
    npoints: int
    min_bound: onp.Array1D[np.float64]
    max_bound: onp.Array1D[np.float64]

    def __init__(self, /, qhull: _Qhull, incremental: bool = False) -> None: ...
    def __del__(self, /) -> None: ...
    def _update(self, /, qhull: _Qhull) -> None: ...
    def _add_points(
        self, /, points: onp.ToFloat2D, restart: bool = False, interior_point: onp.ToFloat1D | None = None
    ) -> None: ...
    def close(self, /) -> None: ...

class Delaunay(_QhullUser):
    furthest_site: bool
    paraboloid_scale: float
    paraboloid_shift: float
    simplices: onp.Array2D[np.int32]
    neighbors: onp.Array2D[np.int32]
    equations: onp.Array2D[np.float64]
    coplanar: onp.Array2D[np.int32]
    good: onp.Array1D[np.int32]
    nsimplex: int

    @property
    def points(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def transform(self, /) -> onp.Array3D[np.float64]: ...
    @property
    def vertex_to_simplex(self, /) -> onp.Array1D[np.int32]: ...
    @property
    def vertex_neighbor_vertices(self, /) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
    @property
    def convex_hull(self, /) -> onp.Array2D[np.int32]: ...

    #
    def __init__(
        self, /, points: onp.ToFloat2D, furthest_site: bool = False, incremental: bool = False, qhull_options: str | None = None
    ) -> None: ...

    #
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...

    #
    @overload
    def find_simplex(
        self, /, xi: _ToArrayStrictND, bruteforce: bool = False, tol: float | None = None
    ) -> onp.ArrayND[np.int32]: ...
    @overload
    def find_simplex(
        self, /, xi: onp.ToFloatStrict1D, bruteforce: bool = False, tol: float | None = None
    ) -> onp.Array0D[np.int32]: ...
    @overload
    def find_simplex(
        self, /, xi: onp.ToFloatStrict2D, bruteforce: bool = False, tol: float | None = None
    ) -> onp.Array1D[np.int32]: ...
    @overload
    def find_simplex(
        self, /, xi: onp.ToFloatStrict3D, bruteforce: bool = False, tol: float | None = None
    ) -> onp.Array2D[np.int32]: ...
    @overload
    def find_simplex(self, /, xi: onp.ToFloatND, bruteforce: bool = False, tol: float | None = None) -> onp.ArrayND[np.int32]: ...

    #
    @overload
    def plane_distance(self, /, xi: onp.ArrayND[npc.floating | npc.integer, _ShapeT]) -> onp.ArrayND[np.float64, _ShapeT]: ...
    @overload
    def plane_distance(self, /, xi: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def lift_points(self, /, x: onp.ArrayND[npc.floating | npc.integer, _ShapeT]) -> onp.ArrayND[np.float64, _ShapeT]: ...
    @overload
    def lift_points(self, /, x: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

class ConvexHull(_QhullUser):
    simplices: onp.Array2D[np.int32]
    neighbors: onp.Array2D[np.int32]
    equations: onp.Array2D[np.float64]
    coplanar: onp.Array2D[np.int32]
    good: onp.Array1D[np.bool] | None
    volume: float
    area: float
    nsimplex: int

    @property
    def points(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def vertices(self, /) -> onp.Array1D[np.int32]: ...

    #
    def __init__(self, /, points: onp.ToFloat2D, incremental: bool = False, qhull_options: str | None = None) -> None: ...
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...

class Voronoi(_QhullUser):
    vertices: onp.Array2D[np.float64]
    ridge_points: onp.Array2D[np.int32]
    ridge_vertices: list[list[int]]
    regions: list[list[int]]
    point_region: onp.Array1D[np.intp]
    furthest_site: bool

    @property
    def points(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def ridge_dict(self, /) -> dict[tuple[int, int], list[int]]: ...

    #
    def __init__(
        self, /, points: onp.ToFloat2D, furthest_site: bool = False, incremental: bool = False, qhull_options: str | None = None
    ) -> None: ...
    def add_points(self, /, points: onp.ToFloat2D, restart: bool = False) -> None: ...

class HalfspaceIntersection(_QhullUser):
    interior_point: onp.Array1D[np.float64]
    intersections: onp.Array2D[np.float64]
    dual_facets: list[list[int]]
    dual_equations: onp.Array2D[np.float64]
    dual_points: onp.Array2D[np.float64]
    dual_volume: float
    dual_area: float
    ndim: int
    nineq: int

    @property
    def halfspaces(self, /) -> onp.Array2D[np.float64]: ...
    @property
    def dual_vertices(self, /) -> onp.Array1D[np.int_]: ...

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

@overload
@deprecated("`tsearch` is deprecated in favor of `Delaunay.find_simplex` and will be removed in SciPy 1.22.0.")
def tsearch(tri: Delaunay, xi: _ToArrayStrictND) -> onp.ArrayND[np.int32]: ...
@overload
@deprecated("`tsearch` is deprecated in favor of `Delaunay.find_simplex` and will be removed in SciPy 1.22.0.")
def tsearch(tri: Delaunay, xi: onp.ToFloatStrict1D) -> onp.Array0D[np.int32]: ...
@overload
@deprecated("`tsearch` is deprecated in favor of `Delaunay.find_simplex` and will be removed in SciPy 1.22.0.")
def tsearch(tri: Delaunay, xi: onp.ToFloatStrict2D) -> onp.Array1D[np.int32]: ...
@overload
@deprecated("`tsearch` is deprecated in favor of `Delaunay.find_simplex` and will be removed in SciPy 1.22.0.")
def tsearch(tri: Delaunay, xi: onp.ToFloatStrict3D) -> onp.Array2D[np.int32]: ...
@overload
@deprecated("`tsearch` is deprecated in favor of `Delaunay.find_simplex` and will be removed in SciPy 1.22.0.")
def tsearch(tri: Delaunay, xi: onp.ToFloatND) -> onp.ArrayND[np.int32]: ...
