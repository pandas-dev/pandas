from collections.abc import Callable
from types import ModuleType
from typing import Final, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype.numpy as onp

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
    def __init__(self, /, id: int, left: ClusterNode, right: ClusterNode, dist: float = 0.0, count: int = 1) -> None: ...

    # NOTE: These raise a `ValueError` if passed anything other than `ClusterNode`
    @override
    def __eq__(self, node: ClusterNode, /) -> bool: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]
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
def int_floor(arr: onp.ToArrayND, xp: _ArrayAPINamespace) -> int: ...  # undocumented

#
def linkage(
    y: onp.ToFloat1D | onp.ToFloat2D,
    method: _LinkageMethod = "single",
    metric: _Metric = "euclidean",
    optimal_ordering: bool = False,
) -> onp.Array2D[np.float64]: ...
def single(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def complete(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def average(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def weighted(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def centroid(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def median(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def ward(y: onp.ToFloat1D | onp.ToFloat2D) -> onp.Array2D[np.float64]: ...

#
@overload  # Y: None  (default)
def cophenet(Z: onp.ToFloat2D, Y: None = None) -> onp.Array1D[np.float64]: ...
@overload  # Y: 1d +f64
def cophenet(Z: onp.ToFloat2D, Y: onp.ToFloat64_1D) -> tuple[np.float64, onp.Array1D[np.float64]]: ...
@overload  # Y: 1d ~f80
def cophenet(Z: onp.ToFloat2D, Y: onp.ToJustLongDouble1D) -> tuple[np.longdouble, onp.Array1D[np.float64]]: ...
@overload  # Y: 1d ~c128
def cophenet(Z: onp.ToFloat2D, Y: onp.ToJustComplex128_1D) -> tuple[np.complex128, onp.Array1D[np.float64]]: ...
@overload  # Y: 1d ~c160
def cophenet(Z: onp.ToFloat2D, Y: onp.ToJustCLongDouble1D) -> tuple[np.clongdouble, onp.Array1D[np.float64]]: ...

#
def from_mlab_linkage(Z: onp.ToFloat2D) -> onp.Array2D[np.float64]: ...
def to_mlab_linkage(Z: onp.ToFloat2D) -> onp.Array2D[np.float64]: ...

#
def inconsistent(Z: onp.ToFloat2D, d: int = 2) -> onp.Array2D[np.float64]: ...
def maxinconsts(Z: onp.ToFloat2D, R: onp.ToFloat2D) -> onp.Array1D[np.float64]: ...
def maxdists(Z: onp.ToFloat2D) -> onp.Array1D[np.float64]: ...
def maxRstat(Z: onp.ToFloat2D, R: onp.ToFloat2D, i: int) -> onp.Array1D[np.float64]: ...

#
def dendrogram(
    Z: onp.ToFloat2D,
    p: int = 30,
    truncate_mode: _TruncateMode | None = None,
    color_threshold: float | None = None,
    get_leaves: bool = True,
    orientation: _Orientation = "top",
    labels: onp.ToArray1D[str, np.str_] | None = None,
    count_sort: _SortOrder | bool = False,
    distance_sort: _SortOrder | bool = False,
    show_leaf_counts: bool = True,
    no_plot: bool = False,
    no_labels: bool = False,
    leaf_font_size: float | None = None,
    leaf_rotation: float | None = None,
    leaf_label_func: Callable[[int], str] | None = None,
    show_contracted: bool = False,
    link_color_func: Callable[[int], str] | None = None,
    ax: _MatplotlibAxes | None = None,
    above_threshold_color: str = "C0",
) -> _DendrogramResult: ...

#
def set_link_color_palette(palette: list[str] | tuple[str, ...] | None) -> None: ...

#
def leaves_list(Z: onp.ToFloat2D) -> onp.Array1D[np.int32]: ...

#
@overload
def to_tree(Z: onp.ToFloat2D, rd: onp.ToFalse = False) -> ClusterNode: ...
@overload
def to_tree(Z: onp.ToFloat2D, rd: onp.ToTrue) -> tuple[ClusterNode, list[ClusterNode]]: ...

#
def cut_tree(
    Z: onp.ToFloat2D, n_clusters: onp.ToInt1D | None = None, height: onp.ToFloat1D | None = None
) -> onp.Array2D[np.int64]: ...

#
def optimal_leaf_ordering(Z: onp.ToFloat2D, y: onp.ToFloat1D, metric: _Metric = "euclidean") -> onp.Array2D[np.float64]: ...

# keep in sync with `is_valid_linkage`
@overload
def is_valid_im(R: onp.ToArrayND, warning: bool = False, throw: Literal[False] = False, name: str | None = None) -> bool: ...
@overload
def is_valid_im(R: onp.ToJustFloat64_2D, warning: bool = False, *, throw: Literal[True], name: str | None = None) -> bool: ...

# keep in sync with `is_valid_im`
@overload
def is_valid_linkage(Z: onp.ToArrayND, warning: bool = False, throw: Literal[False] = False, name: str | None = None) -> bool: ...
@overload
def is_valid_linkage(
    Z: onp.ToJustFloat64_2D, warning: bool = False, *, throw: Literal[True], name: str | None = None
) -> bool: ...

#
def is_isomorphic(T1: onp.ToArray1D, T2: onp.ToArray1D) -> bool: ...
def is_monotonic(Z: onp.ToJustFloat64_2D) -> bool: ...
def correspond(Z: onp.ToJustFloat64_2D, Y: onp.ToFloat1D) -> bool: ...
def num_obs_linkage(Z: onp.ToJustFloat64_2D) -> int: ...

#
@overload
def fcluster(
    Z: onp.ToFloat2D,
    t: onp.ToFloat,
    criterion: _ClusterCriterion = "inconsistent",
    depth: int = 2,
    R: None = None,
    monocrit: onp.ToFloat1D | None = None,
) -> onp.Array1D[np.int32]: ...
@overload  # criterion="inconsistent"  (default)
def fcluster(
    Z: onp.ToFloat2D,
    t: onp.ToFloat,
    criterion: Literal["inconsistent"] = "inconsistent",
    depth: int = 2,
    R: onp.ToFloat2D | None = None,
    monocrit: onp.ToFloat1D | None = None,
) -> onp.Array1D[np.int32]: ...

# keep in sync with `fcluster`
@overload
def fclusterdata(
    X: onp.ToFloat2D,
    t: onp.ToFloat,
    criterion: _ClusterCriterion = "inconsistent",
    metric: _Metric = "euclidean",
    depth: int = 2,
    method: _LinkageMethod = "single",
    R: None = None,
) -> onp.Array1D[np.int32]: ...
@overload
def fclusterdata(
    X: onp.ToFloat2D,
    t: onp.ToFloat,
    criterion: Literal["inconsistent"] = "inconsistent",
    metric: _Metric = "euclidean",
    depth: int = 2,
    method: _LinkageMethod = "single",
    R: onp.ToFloat1D | None = None,
) -> onp.Array1D[np.int32]: ...

#
def leaders(Z: onp.ToFloat2D, T: onp.Array1D[np.int32]) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
