from collections.abc import Callable
from types import ModuleType
from typing import Final, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._disjoint_set import DisjointSet
from scipy.spatial.distance import _Metric

__all__ = [
    "ClusterNode",
    "DisjointSet",
    "average",
    "centroid",
    "complete",
    "cophenet",
    "correspond",
    "cut_tree",
    "dendrogram",
    "fcluster",
    "fclusterdata",
    "from_mlab_linkage",
    "inconsistent",
    "is_isomorphic",
    "is_monotonic",
    "is_valid_im",
    "is_valid_linkage",
    "leaders",
    "leaves_list",
    "linkage",
    "maxRstat",
    "maxdists",
    "maxinconsts",
    "median",
    "num_obs_linkage",
    "optimal_leaf_ordering",
    "set_link_color_palette",
    "single",
    "to_mlab_linkage",
    "to_tree",
    "ward",
    "weighted",
]

_T = TypeVar("_T")
_SCT = TypeVar("_SCT", bound=npc.number, default=np.float64)

_LinkageArray: TypeAlias = onp.Array2D[_SCT]
_LinkageMethod: TypeAlias = Literal["single", "complete", "average", "weighted", "centroid", "median", "ward"]
_ClusterCriterion: TypeAlias = Literal["inconsistent", "distance", "maxclust", "monocrit", "maxclust_monocrit"]
_SortOrder: TypeAlias = Literal["ascending", "descending"]
_TruncateMode: TypeAlias = Literal["lastp", "level"]
_Orientation: TypeAlias = Literal["top", "bottom", "left", "right"]

# for the lack of better types
_MatplotlibAxes: TypeAlias = object
_ArrayAPINamespace: TypeAlias = ModuleType

@type_check_only
class _DendrogramResult(TypedDict):
    color_list: list[str]
    icoord: list[list[int]]
    dcoord: list[list[int]]
    ivl: list[str]
    leaves: list[int] | None
    leaves_color_list: list[str]

###

class ClusterWarning(UserWarning): ...

# NOTE: this can't be made generic, because mypy doesn't support cyclical generic types (classic mypy...)
class ClusterNode:
    id: Final[int]
    left: Final[ClusterNode | None]
    right: Final[ClusterNode | None]
    dist: Final[float]
    count: Final[int]

    # NOTE: either both `left` and `right` are None, or both are `ClusterNode`
    @overload
    def __init__(
        self: ClusterNode, /, id: int, left: None = None, right: None = None, dist: float = 0.0, count: int = 1
    ) -> None: ...
    @overload
    def __init__(self, /, id: int, left: ClusterNode, right: ClusterNode, dist: float = 0, count: int = 1) -> None: ...

    # NOTE: These raise a `ValueError` if passed anything other than `ClusterNode`
    @override
    def __eq__(self, node: ClusterNode, /) -> bool: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    def __lt__(self, node: ClusterNode, /) -> bool: ...
    def __gt__(self, node: ClusterNode, /) -> bool: ...

    # NOTE: These getters are basically redundant, as the attributes they (directly) return are public anyway
    def get_id(self, /) -> int: ...
    def get_count(self, /) -> int: ...
    def get_left(self, /) -> ClusterNode | None: ...
    def get_right(self, /) -> ClusterNode | None: ...

    # NOTE: True iff `left` (and therefore `right`) is `None`
    def is_leaf(self, /) -> bool: ...

    # NOTE: `func` defaults to `(x) -> x.id`
    @overload
    def pre_order(self, /, func: Callable[[ClusterNode], int] = ...) -> list[int]: ...
    @overload
    def pre_order(self, /, func: Callable[[ClusterNode], _T]) -> list[_T]: ...

#
def int_floor(arr: onp.ToArrayND, xp: _ArrayAPINamespace) -> int: ...

#
def single(y: onp.ToArrayND) -> _LinkageArray: ...
def complete(y: onp.ToArrayND) -> _LinkageArray: ...
def average(y: onp.ToArrayND) -> _LinkageArray: ...
def weighted(y: onp.ToArrayND) -> _LinkageArray: ...
def centroid(y: onp.ToArrayND) -> _LinkageArray: ...
def median(y: onp.ToArrayND) -> _LinkageArray: ...
def ward(y: onp.ToArrayND) -> _LinkageArray: ...
def linkage(
    y: onp.ToArrayND, method: _LinkageMethod = "single", metric: _Metric = "euclidean", optimal_ordering: bool = False
) -> _LinkageArray[np.int_ | np.float64 | np.complex128]: ...

#
def cut_tree(
    Z: onp.ToArray2D, n_clusters: onp.ToInt1D | None = None, height: onp.ToFloat1D | None = None
) -> onp.Array2D[np.int64]: ...

#
@overload
def to_tree(Z: onp.ToArray2D, rd: onp.ToFalse = False) -> ClusterNode: ...
@overload
def to_tree(Z: onp.ToArray2D, rd: onp.ToTrue) -> tuple[ClusterNode, list[ClusterNode]]: ...

#
def optimal_leaf_ordering(Z: onp.ToArray2D, y: onp.ToArrayND, metric: _Metric = "euclidean") -> _LinkageArray: ...

#
@overload
def cophenet(Z: onp.ToArray2D, Y: None = None) -> onp.Array1D[np.float64]: ...
@overload
def cophenet(Z: onp.ToArray2D, Y: onp.ToArrayND) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64]]: ...

#
def inconsistent(Z: onp.ToArray2D, d: int = 2) -> _LinkageArray: ...

#
def from_mlab_linkage(Z: onp.ToArray2D) -> _LinkageArray: ...
def to_mlab_linkage(Z: onp.ToArray2D) -> _LinkageArray: ...

#
def is_monotonic(Z: onp.ToArray2D) -> bool: ...
def is_valid_im(R: onp.ToArrayND, warning: bool = False, throw: bool = False, name: str | None = None) -> bool: ...
def is_valid_linkage(Z: onp.ToArray2D, warning: bool = False, throw: bool = False, name: str | None = None) -> bool: ...
def is_isomorphic(T1: onp.ToArrayND, T2: onp.ToArrayND) -> bool: ...

#
def num_obs_linkage(Z: onp.ToArray2D) -> int: ...

#
def correspond(Z: onp.ToArray2D, Y: onp.ToArrayND) -> bool: ...

#
def fcluster(
    Z: onp.ToArray2D,
    t: onp.ToFloat,
    criterion: _ClusterCriterion = "inconsistent",
    depth: op.JustInt = 2,
    R: onp.ToArrayND | None = None,
    monocrit: onp.ToArrayND | None = None,
) -> onp.Array1D[np.int32]: ...

#
def fclusterdata(
    X: onp.ToArrayND,
    t: onp.ToFloat,
    criterion: _ClusterCriterion = "inconsistent",
    metric: _Metric = "euclidean",
    depth: op.JustInt = 2,
    method: _LinkageMethod = "single",
    R: onp.ToArrayND | None = None,
) -> onp.Array1D[np.int32]: ...

#
def leaves_list(Z: onp.ToArray2D) -> onp.Array1D[np.int32]: ...

#
def set_link_color_palette(palette: list[str] | tuple[str, ...] | None) -> None: ...

#
def dendrogram(
    Z: onp.ToArray2D,
    p: int = 30,
    truncate_mode: _TruncateMode | None = None,
    color_threshold: float | npc.floating | None = None,
    get_leaves: bool = True,
    orientation: _Orientation = "top",
    labels: onp.ToArrayND | None = None,
    count_sort: _SortOrder | bool = False,
    distance_sort: _SortOrder | bool = False,
    show_leaf_counts: bool = True,
    no_plot: bool = False,
    no_labels: bool = False,
    leaf_font_size: float | npc.floating | None = None,
    leaf_rotation: float | npc.floating | None = None,
    leaf_label_func: Callable[[int], str] | None = None,
    show_contracted: bool = False,
    link_color_func: Callable[[int], str] | None = None,
    ax: _MatplotlibAxes | None = None,
    above_threshold_color: str = "C0",
) -> _DendrogramResult: ...

#
def maxdists(Z: onp.ToArray2D) -> onp.Array1D[np.float64]: ...
def maxinconsts(Z: onp.ToArray2D, R: onp.ToArrayND) -> onp.Array1D[np.float64]: ...
def maxRstat(Z: onp.ToArray2D, R: onp.ToArrayND, i: int) -> onp.Array1D[np.float64]: ...
def leaders(Z: onp.ToArray2D, T: onp.ToArrayND) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
