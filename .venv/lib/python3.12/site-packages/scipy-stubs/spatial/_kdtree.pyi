from typing import Final, Generic, Self, TypeAlias, overload
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp

from ._ckdtree import cKDTree, cKDTreeNode

__all__ = ["KDTree", "Rectangle", "distance_matrix", "minkowski_distance", "minkowski_distance_p"]

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_BoxSizeT_co = TypeVar("_BoxSizeT_co", bound=_Float2D | None, default=_Float2D | None, covariant=True)
_BoxSizeDataT_co = TypeVar("_BoxSizeDataT_co", bound=_Float1D | None, default=_Float1D | None, covariant=True)

###

class Rectangle:
    maxes: Final[onp.Array1D[np.float64]]
    mins: Final[onp.Array1D[np.float64]]
    def __init__(self, /, maxes: onp.ToFloat1D, mins: onp.ToFloat1D) -> None: ...
    def volume(self, /) -> np.float64: ...
    def split(self, /, d: op.CanIndex, split: onp.ToFloat) -> tuple[Self, Self]: ...
    def min_distance_point(self, /, x: onp.ToFloat | onp.ToFloatND, p: onp.ToFloat = 2.0) -> onp.ArrayND[np.float64]: ...
    def max_distance_point(self, /, x: onp.ToFloat | onp.ToFloatND, p: onp.ToFloat = 2.0) -> onp.ArrayND[np.float64]: ...
    def min_distance_rectangle(self, /, other: Rectangle, p: onp.ToFloat = 2.0) -> onp.ArrayND[np.float64]: ...
    def max_distance_rectangle(self, /, other: Rectangle, p: onp.ToFloat = 2.0) -> onp.ArrayND[np.float64]: ...

class KDTree(cKDTree[_BoxSizeT_co, _BoxSizeDataT_co], Generic[_BoxSizeT_co, _BoxSizeDataT_co]):
    class node:
        @staticmethod
        def _create(ckdtree_node: cKDTreeNode | None = None) -> KDTree.leafnode | KDTree.innernode: ...
        def __init__(self, /, ckdtree_node: cKDTreeNode | None = None) -> None: ...
        def __lt__(self, other: object, /) -> bool: ...
        def __gt__(self, other: object, /) -> bool: ...
        def __le__(self, other: object, /) -> bool: ...
        def __ge__(self, other: object, /) -> bool: ...

    class leafnode(node):
        @property
        def idx(self, /) -> onp.ArrayND[np.intp]: ...
        @property
        def children(self, /) -> int: ...

    class innernode(node):
        less: Final[KDTree.innernode | KDTree.leafnode]
        greater: Final[KDTree.innernode | KDTree.leafnode]

        def __init__(self, /, ckdtreenode: cKDTreeNode) -> None: ...
        @property
        def split_dim(self, /) -> int: ...
        @property
        def split(self, /) -> float: ...
        @property
        def children(self, /) -> int: ...

    @overload
    def __init__(
        self: KDTree[None, None],
        /,
        data: onp.ToComplexND,
        leafsize: onp.ToInt = 10,
        compact_nodes: bool = True,
        copy_data: bool = False,
        balanced_tree: bool = True,
        boxsize: None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: KDTree[_Float2D, _Float1D],
        /,
        data: onp.ToComplexND,
        leafsize: onp.ToInt,
        compact_nodes: bool,
        copy_data: bool,
        balanced_tree: bool,
        boxsize: onp.ToFloat2D,
    ) -> None: ...
    @overload
    def __init__(
        self: KDTree[_Float2D, _Float1D],
        /,
        data: onp.ToComplexND,
        leafsize: onp.ToInt = 10,
        compact_nodes: bool = True,
        copy_data: bool = False,
        balanced_tree: bool = True,
        *,
        boxsize: onp.ToFloat2D,
    ) -> None: ...

    #
    @override
    def query(
        self,
        /,
        x: onp.ToFloat1D,
        k: onp.ToInt | onp.ToInt1D = 1,
        eps: onp.ToFloat = 0.0,
        p: onp.ToFloat = 2.0,
        distance_upper_bound: float = float("inf"),  # noqa: PYI011
        workers: int | None = 1,
    ) -> tuple[float, np.intp] | tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...

@overload
def minkowski_distance_p(x: onp.ToFloatND, y: onp.ToFloatND, p: float = 2.0) -> onp.ArrayND[np.float64]: ...
@overload
def minkowski_distance_p(x: onp.ToComplexND, y: onp.ToComplexND, p: float = 2.0) -> onp.ArrayND[np.float64 | np.complex128]: ...

#
@overload
def minkowski_distance(x: onp.ToFloatND, y: onp.ToFloatND, p: float = 2.0) -> onp.ArrayND[np.float64]: ...
@overload
def minkowski_distance(x: onp.ToComplexND, y: onp.ToComplexND, p: float = 2.0) -> onp.ArrayND[np.float64 | np.complex128]: ...

#
@overload
def distance_matrix(
    x: onp.ToFloatND, y: onp.ToFloatND, p: float = 2.0, threshold: int = 1_000_000
) -> onp.Array2D[np.float64]: ...
@overload
def distance_matrix(
    x: onp.ToComplexND, y: onp.ToComplexND, p: float = 2.0, threshold: int = 1_000_000
) -> onp.Array2D[np.float64 | np.complex128]: ...
