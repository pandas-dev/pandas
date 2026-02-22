import numpy as np
import optype.numpy as onp

__all__ = ["SphericalVoronoi"]

def calculate_solid_angles(R: onp.ToComplexND) -> onp.ArrayND[np.float64]: ...

class SphericalVoronoi:
    points: onp.Array2D[np.float64]
    center: onp.Array1D[np.float64]
    radius: float

    def __init__(
        self, /, points: onp.ToFloat2D, radius: float = 1, center: onp.ToFloat1D | None = None, threshold: float = 1e-06
    ) -> None: ...
    def sort_vertices_of_regions(self, /) -> None: ...
    def calculate_areas(self, /) -> onp.Array1D[np.float64]: ...
