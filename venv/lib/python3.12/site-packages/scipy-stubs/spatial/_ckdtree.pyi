from typing import Generic, Literal as L, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.sparse import coo_matrix, dok_matrix

__all__ = ["cKDTree"]

_Weights: TypeAlias = onp.ToFloatND | tuple[onp.ToFloatND, onp.ToFloatND]
_Indices: TypeAlias = onp.Array1D[np.intp]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_NodeT_co = TypeVar("_NodeT_co", bound=_KDTreeNode | None, default=_KDTreeNode | None, covariant=True)
_BoxSizeT_co = TypeVar("_BoxSizeT_co", bound=_Float2D | None, default=_Float2D | None, covariant=True)
_BoxSizeDataT_co = TypeVar("_BoxSizeDataT_co", bound=_Float1D | None, default=_Float1D | None, covariant=True)

@type_check_only
class _CythonMixin:
    def __setstate_cython__(self, pyx_state: object, /) -> None: ...
    def __reduce_cython__(self, /) -> None: ...

# workaround for mypy's lack of cyclical TypeVar support
@type_check_only
class _KDTreeNode(Protocol):
    @property
    def level(self, /) -> int: ...
    @property
    def split_dim(self, /) -> int: ...
    @property
    def split(self, /) -> float: ...
    @property
    def children(self, /) -> int: ...
    @property
    def data_points(self, /) -> _Float2D: ...
    @property
    def indices(self, /) -> _Indices: ...
    @property
    def start_idx(self, /) -> int: ...
    @property
    def end_idx(self, /) -> int: ...
    @property
    def lesser(self, /) -> _KDTreeNode | None: ...
    @property
    def greater(self, /) -> _KDTreeNode | None: ...

###

class cKDTreeNode(_CythonMixin, _KDTreeNode, Generic[_NodeT_co]):
    @property
    @override
    def lesser(self, /) -> _NodeT_co: ...
    @property
    @override
    def greater(self, /) -> _NodeT_co: ...

class cKDTree(_CythonMixin, Generic[_BoxSizeT_co, _BoxSizeDataT_co]):
    @property
    def data(self, /) -> _Float2D: ...
    @property
    def leafsize(self, /) -> int: ...
    @property
    def m(self, /) -> int: ...
    @property
    def n(self, /) -> int: ...
    @property
    def maxes(self, /) -> _Float1D: ...
    @property
    def mins(self, /) -> _Float1D: ...
    @property
    def tree(self, /) -> cKDTreeNode: ...
    @property
    def size(self, /) -> int: ...
    @property
    def indices(self, /) -> _Indices: ...
    @property
    def boxsize(self, /) -> _BoxSizeT_co: ...
    boxsize_data: _BoxSizeDataT_co

    #
    @overload
    def __init__(
        self: cKDTree[None, None],
        /,
        data: onp.ToFloat2D,
        leafsize: int = 16,
        compact_nodes: bool = True,
        copy_data: bool = False,
        balanced_tree: bool = True,
        boxsize: None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: cKDTree[_Float2D, _Float1D],
        /,
        data: onp.ToFloat2D,
        leafsize: int,
        compact_nodes: bool,
        copy_data: bool,
        balanced_tree: bool,
        boxsize: onp.ToFloat2D,
    ) -> None: ...
    @overload
    def __init__(
        self: cKDTree[_Float2D, _Float1D],
        /,
        data: onp.ToFloat2D,
        leafsize: int = 16,
        compact_nodes: bool = True,
        copy_data: bool = False,
        balanced_tree: bool = True,
        *,
        boxsize: onp.ToFloat2D,
    ) -> None: ...

    #
    def query(
        self,
        /,
        x: onp.ToFloat1D,
        k: onp.ToInt | onp.ToInt1D = 1,
        eps: onp.ToFloat = 0.0,
        p: onp.ToFloat = 2.0,
        distance_upper_bound: float = float("inf"),  # noqa: PYI011
        workers: int | None = None,
    ) -> tuple[float, np.intp] | tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...

    # NOTE: The parameters `eps` and `p` default to `0.0` and `2.0` in `cKDTree`, but are overridden in KDTree to default to
    # `0` and `2` (or `2.0`) respectively. Filling in these defaults would therefore require us to override these methods in
    # `KDTree`, which otherwise wouldn't be necessary. Hence, we leave these parameters without defaults here, so that we avoid a
    # lot of duplicated code.
    # In scipy 1.17.0 this will no longer be necessary (scipy/scipy#23727).

    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatStrict1D,
        r: onp.ToFloat,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        return_length: onp.ToFalse = False,
    ) -> list[int]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatStrict1D,
        r: onp.ToFloat,
        p: onp.ToFloat,
        eps: onp.ToFloat,
        workers: op.CanIndex | None,
        return_sorted: onp.ToBool | None,
        return_length: onp.ToTrue,
    ) -> np.intp: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatStrict1D,
        r: onp.ToFloat,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        *,
        return_length: onp.ToTrue,
    ) -> np.intp: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloatND,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        return_length: onp.ToFalse = False,
    ) -> onp.ArrayND[np.object_]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloatND,
        p: onp.ToFloat,
        eps: onp.ToFloat,
        workers: op.CanIndex | None,
        return_sorted: onp.ToBool | None,
        return_length: onp.ToTrue,
    ) -> onp.ArrayND[np.intp]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloatND,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        *,
        return_length: onp.ToTrue,
    ) -> onp.ArrayND[np.intp]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloat | onp.ToFloatND,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        return_length: onp.ToFalse = False,
    ) -> list[int] | onp.ArrayND[np.object_]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloat | onp.ToFloatND,
        p: onp.ToFloat,
        eps: onp.ToFloat,
        workers: op.CanIndex | None,
        return_sorted: onp.ToBool | None,
        return_length: onp.ToTrue,
    ) -> np.intp | onp.ArrayND[np.intp]: ...
    @overload
    def query_ball_point(
        self,
        /,
        x: onp.ToFloatND,
        r: onp.ToFloat | onp.ToFloatND,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # # stubdefaulter: ignore[missing-default]
        workers: op.CanIndex | None = None,
        return_sorted: onp.ToBool | None = None,
        *,
        return_length: onp.ToTrue,
    ) -> np.intp | onp.ArrayND[np.intp]: ...

    #
    def query_ball_tree(
        self,
        /,
        other: cKDTree,
        r: onp.ToFloat,
        p: onp.ToFloat = 2.0,
        eps: onp.ToFloat = ...,  # stubdefaulter: ignore[missing-default]
    ) -> list[list[int]]: ...

    #
    @overload
    def query_pairs(
        self, /, r: onp.ToFloat, p: onp.ToFloat = 2.0, eps: onp.ToFloat = 0.0, output_type: L["set"] = "set"
    ) -> set[tuple[int, int]]: ...
    @overload
    def query_pairs(
        self, /, r: onp.ToFloat, p: onp.ToFloat, eps: onp.ToFloat, output_type: L["ndarray"]
    ) -> onp.ArrayND[np.intp]: ...
    @overload
    def query_pairs(
        self, /, r: onp.ToFloat, p: onp.ToFloat = 2.0, eps: onp.ToFloat = 0.0, *, output_type: L["ndarray"]
    ) -> onp.ArrayND[np.intp]: ...

    #
    @overload
    def count_neighbors(
        self,
        /,
        other: cKDTree,
        r: onp.ToFloat,
        p: onp.ToFloat = 2.0,
        weights: tuple[None, None] | None = None,
        cumulative: bool = True,
    ) -> np.intp: ...
    @overload
    def count_neighbors(
        self, /, other: cKDTree, r: onp.ToFloat, p: onp.ToFloat, weights: _Weights, cumulative: bool = True
    ) -> np.float64: ...
    @overload
    def count_neighbors(
        self, /, other: cKDTree, r: onp.ToFloat, p: onp.ToFloat = 2.0, *, weights: _Weights, cumulative: bool = True
    ) -> np.float64: ...
    @overload
    def count_neighbors(
        self,
        /,
        other: cKDTree,
        r: onp.ToFloat | onp.ToFloat1D,
        p: onp.ToFloat = 2.0,
        weights: tuple[None, None] | None = None,
        cumulative: bool = True,
    ) -> np.intp | onp.Array1D[np.intp]: ...
    @overload
    def count_neighbors(
        self, /, other: cKDTree, r: onp.ToFloat | onp.ToFloat1D, p: onp.ToFloat, weights: _Weights, cumulative: bool = True
    ) -> np.float64 | onp.Array1D[np.float64]: ...
    @overload
    def count_neighbors(
        self,
        /,
        other: cKDTree,
        r: onp.ToFloat | onp.ToFloat1D,
        p: onp.ToFloat = 2.0,
        *,
        weights: _Weights,
        cumulative: bool = True,
    ) -> np.float64 | onp.Array1D[np.float64]: ...

    #
    @overload
    def sparse_distance_matrix(
        self, /, other: cKDTree, max_distance: onp.ToFloat, p: onp.ToFloat = 2.0, output_type: L["dok_matrix"] = "dok_matrix"
    ) -> dok_matrix[np.float64]: ...
    @overload
    def sparse_distance_matrix(
        self, /, other: cKDTree, max_distance: onp.ToFloat, p: onp.ToFloat = 2.0, *, output_type: L["coo_matrix"]
    ) -> coo_matrix[np.float64]: ...
    @overload
    def sparse_distance_matrix(
        self, /, other: cKDTree, max_distance: onp.ToFloat, p: onp.ToFloat = 2.0, *, output_type: L["dict"]
    ) -> dict[tuple[int, int], float]: ...
    @overload
    def sparse_distance_matrix(
        self, /, other: cKDTree, max_distance: onp.ToFloat, p: onp.ToFloat = 2.0, *, output_type: L["ndarray"]
    ) -> onp.ArrayND[np.void]: ...
