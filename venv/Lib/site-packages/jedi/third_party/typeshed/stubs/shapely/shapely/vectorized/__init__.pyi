from typing import overload
from typing_extensions import deprecated

import numpy as np
from numpy.typing import NDArray

from .._typing import ArrayLike, ArrayLikeSeq
from ..lib import Geometry
from ..prepared import PreparedGeometry

@overload
@deprecated("Use 'shapely.contains_xy' instead (available since shapely 2.0.0).")
def contains(geometry: Geometry | PreparedGeometry[Geometry], x: float, y: float) -> np.bool_: ...
@overload
@deprecated("Use 'shapely.contains_xy' instead (available since shapely 2.0.0).")
def contains(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLikeSeq[float], y: ArrayLike[float]
) -> NDArray[np.bool_]: ...
@overload
@deprecated("Use 'shapely.contains_xy' instead (available since shapely 2.0.0).")
def contains(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLike[float], y: ArrayLikeSeq[float]
) -> NDArray[np.bool_]: ...
@overload
@deprecated("Use 'shapely.contains_xy' instead (available since shapely 2.0.0).")
def contains(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLike[float], y: ArrayLike[float]
) -> np.bool_ | NDArray[np.bool_]: ...
@overload
@deprecated("Use 'shapely.intersects_xy' instead (available since shapely 2.0.0).")
def touches(geometry: Geometry | PreparedGeometry[Geometry], x: float, y: float) -> np.bool_: ...
@overload
@deprecated("Use 'shapely.intersects_xy' instead (available since shapely 2.0.0).")
def touches(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLikeSeq[float], y: ArrayLike[float]
) -> NDArray[np.bool_]: ...
@overload
@deprecated("Use 'shapely.intersects_xy' instead (available since shapely 2.0.0).")
def touches(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLike[float], y: ArrayLikeSeq[float]
) -> NDArray[np.bool_]: ...
@overload
@deprecated("Use 'shapely.intersects_xy' instead (available since shapely 2.0.0).")
def touches(
    geometry: Geometry | PreparedGeometry[Geometry], x: ArrayLike[float], y: ArrayLike[float]
) -> np.bool_ | NDArray[np.bool_]: ...
