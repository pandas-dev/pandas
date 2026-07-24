# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.spatial` namespace for importing the functions
# included below.
from typing_extensions import deprecated

from ._ckdtree import cKDTree as _cKDTree
from ._kdtree import KDTree as _KDTree, Rectangle as _Rectangle

__all__ = ["KDTree", "Rectangle", "cKDTree", "distance_matrix", "minkowski_distance", "minkowski_distance_p"]

@deprecated("will be removed in SciPy v2.0.0")
class cKDTree(_cKDTree): ...

@deprecated("will be removed in SciPy v2.0.0")
class KDTree(_KDTree): ...

@deprecated("will be removed in SciPy v2.0.0")
class Rectangle(_Rectangle): ...

def minkowski_distance_p(x: object, y: object, p: object = 2.0) -> object: ...
def minkowski_distance(x: object, y: object, p: object = 2.0) -> object: ...
def distance_matrix(x: object, y: object, p: object = 2.0, threshold: object = 1_000_000) -> object: ...
