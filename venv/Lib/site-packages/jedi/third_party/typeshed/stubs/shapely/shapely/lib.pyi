from typing import Literal, SupportsIndex, final, overload
from typing_extensions import Never, Self, disjoint_base

import numpy as np
from numpy.typing import NDArray

area: np.ufunc
boundary: np.ufunc
bounds: np.ufunc
box: np.ufunc
buffer: np.ufunc
build_area: np.ufunc
centroid: np.ufunc
clip_by_rect: np.ufunc
concave_hull: np.ufunc
constrained_delaunay_triangles: np.ufunc
contains: np.ufunc
contains_properly: np.ufunc
contains_xy: np.ufunc
convex_hull: np.ufunc
coverage_invalid_edges: np.ufunc
coverage_is_valid: np.ufunc
coverage_simplify: np.ufunc
coverage_union: np.ufunc
covered_by: np.ufunc
covers: np.ufunc
create_collection: np.ufunc
crosses: np.ufunc
delaunay_triangles: np.ufunc
destroy_prepared: np.ufunc
difference: np.ufunc
difference_prec: np.ufunc
disjoint: np.ufunc
disjoint_subset_union: np.ufunc
distance: np.ufunc
dwithin: np.ufunc
envelope: np.ufunc
equals: np.ufunc
equals_exact: np.ufunc
equals_identical: np.ufunc
extract_unique_points: np.ufunc
force_2d: np.ufunc
force_3d: np.ufunc
frechet_distance: np.ufunc
frechet_distance_densify: np.ufunc
from_geojson: np.ufunc
from_wkb: np.ufunc
from_wkt: np.ufunc
get_coordinate_dimension: np.ufunc
get_dimensions: np.ufunc
get_exterior_ring: np.ufunc
get_geometry: np.ufunc
get_interior_ring: np.ufunc
get_m: np.ufunc
get_num_coordinates: np.ufunc
get_num_geometries: np.ufunc
get_num_interior_rings: np.ufunc
get_num_points: np.ufunc
get_point: np.ufunc
get_precision: np.ufunc
get_srid: np.ufunc
get_type_id: np.ufunc
has_m: np.ufunc
get_x: np.ufunc
get_y: np.ufunc
get_z: np.ufunc
has_z: np.ufunc
hausdorff_distance: np.ufunc
hausdorff_distance_densify: np.ufunc
intersection: np.ufunc
intersection_all: np.ufunc
intersection_prec: np.ufunc
intersects: np.ufunc
intersects_xy: np.ufunc
is_ccw: np.ufunc
is_closed: np.ufunc
is_empty: np.ufunc
is_geometry: np.ufunc
is_missing: np.ufunc
is_prepared: np.ufunc
is_ring: np.ufunc
is_simple: np.ufunc
is_valid: np.ufunc
is_valid_input: np.ufunc
is_valid_reason: np.ufunc
length: np.ufunc
line_interpolate_point: np.ufunc
line_interpolate_point_normalized: np.ufunc
line_locate_point: np.ufunc
line_locate_point_normalized: np.ufunc
line_merge: np.ufunc
line_merge_directed: np.ufunc
linearrings: np.ufunc
linestrings: np.ufunc
make_valid: np.ufunc
make_valid_with_params: np.ufunc
maximum_inscribed_circle: np.ufunc
minimum_bounding_circle: np.ufunc
minimum_bounding_radius: np.ufunc
minimum_clearance: np.ufunc
minimum_clearance_line: np.ufunc
node: np.ufunc
normalize: np.ufunc
offset_curve: np.ufunc
orient_polygons: np.ufunc
oriented_envelope: np.ufunc
overlaps: np.ufunc
point_on_surface: np.ufunc
points: np.ufunc
polygonize: np.ufunc
polygonize_full: np.ufunc
polygons: np.ufunc
prepare: np.ufunc
relate: np.ufunc
relate_pattern: np.ufunc
remove_repeated_points: np.ufunc
reverse: np.ufunc
segmentize: np.ufunc
set_precision: np.ufunc
set_srid: np.ufunc
shared_paths: np.ufunc
shortest_line: np.ufunc
simplify: np.ufunc
simplify_preserve_topology: np.ufunc
snap: np.ufunc
symmetric_difference: np.ufunc
symmetric_difference_all: np.ufunc
symmetric_difference_prec: np.ufunc
to_geojson: np.ufunc
to_wkb: np.ufunc
to_wkt: np.ufunc
touches: np.ufunc
unary_union: np.ufunc
unary_union_prec: np.ufunc
union: np.ufunc
union_prec: np.ufunc
voronoi_polygons: np.ufunc
within: np.ufunc
geos_capi_version: tuple[int, int, int]
geos_capi_version_string: str
geos_version: tuple[int, int, int]
geos_version_string: str
registry: list[type[Geometry]]

@disjoint_base
class Geometry:
    def __hash__(self) -> int: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __ne__(self, other: object, /) -> bool: ...
    def __ge__(self, other: Never, /) -> bool: ...
    def __gt__(self, other: Never, /) -> bool: ...
    def __le__(self, other: Never, /) -> bool: ...
    def __lt__(self, other: Never, /) -> bool: ...

@final
class STRtree:
    count: int
    def __new__(cls, geoms: NDArray[np.object_], node_capacity: SupportsIndex, /, **kwargs: object) -> Self: ...
    def dwithin(self, geoms: NDArray[np.object_], distances: NDArray[np.float64], /) -> NDArray[np.int64]: ...
    def nearest(self, geoms: NDArray[np.object_], /) -> NDArray[np.int64]: ...
    def query(self, geoms: NDArray[np.object_], predicate: SupportsIndex, /) -> NDArray[np.int64]: ...
    def query_nearest(
        self, geoms: NDArray[np.object_], max_distance: float, exclusive: SupportsIndex, all_matches: SupportsIndex, /
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]: ...

class ShapelyError(Exception): ...
class GEOSException(ShapelyError): ...

def count_coordinates(geoms: NDArray[np.object_], /) -> int: ...
@overload
def get_coordinates(arr: NDArray[np.object_], include_z: bool, return_index: Literal[False], /) -> NDArray[np.float64]: ...
@overload
def get_coordinates(
    arr: NDArray[np.object_], include_z: bool, return_index: Literal[True], /
) -> tuple[NDArray[np.float64], NDArray[np.int64]]: ...
@overload
def get_coordinates(
    arr: NDArray[np.object_], include_z: bool, return_index: bool, /
) -> NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.int64]]: ...
def set_coordinates(geoms: NDArray[np.object_], coords: NDArray[np.float64], /) -> NDArray[np.object_]: ...
