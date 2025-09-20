from typing import Any, Generic, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar, Unpack, override

import numpy as np
import optype.numpy as onp

from ._interpnd import CloughTocher2DInterpolator, LinearNDInterpolator, NDInterpolatorBase
from scipy.spatial._ckdtree import cKDTree

__all__ = ["CloughTocher2DInterpolator", "LinearNDInterpolator", "NearestNDInterpolator", "griddata"]

###

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

_Method: TypeAlias = Literal["nearest", "linear", "cubic"]
_ToXi: TypeAlias = onp.ToFloat2D | tuple[onp.ToFloat1D | onp.ToFloat2D, ...]

@type_check_only
class _TreeOptions(TypedDict, total=False):
    leafsize: onp.ToJustInt
    compact_nodes: onp.ToBool
    copy_data: onp.ToBool
    balanced_tree: onp.ToBool
    boxsize: onp.ToFloatND | None

@type_check_only
class _QueryOptions(TypedDict, total=False):
    eps: onp.ToFloat
    p: onp.ToFloat
    distance_upper_bound: onp.ToFloat
    workers: int

###

class NearestNDInterpolator(NDInterpolatorBase[_CT_co], Generic[_CT_co]):
    tree: cKDTree

    @overload
    def __init__(
        self: NearestNDInterpolator[np.float64],
        /,
        x: onp.ToFloat2D,
        y: onp.ToFloat1D,
        rescale: bool = False,
        tree_options: _TreeOptions | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: NearestNDInterpolator[np.complex128],
        /,
        x: onp.ToFloat2D,
        y: onp.ToJustComplex1D,
        rescale: bool = False,
        tree_options: _TreeOptions | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: NearestNDInterpolator[Any],
        /,
        x: onp.ToFloat2D,
        y: onp.ToComplex1D,
        rescale: bool = False,
        tree_options: _TreeOptions | None = None,
    ) -> None: ...

    #
    @override
    def __call__(self, /, *args: onp.ToFloatND, **query_options: Unpack[_QueryOptions]) -> onp.Array[onp.AtLeast1D, _CT_co]: ...

#
@overload
def griddata(
    points: onp.ToFloat1D | onp.ToFloat2D,
    values: onp.ToFloat1D,
    xi: _ToXi,
    method: _Method = "linear",
    fill_value: onp.ToFloat = ...,  # np.nan
    rescale: onp.ToBool = False,
) -> onp.Array[onp.AtLeast1D, np.float64]: ...
@overload
def griddata(
    points: onp.ToFloat1D | onp.ToFloat2D,
    values: onp.ToJustComplex1D,
    xi: _ToXi,
    method: _Method = "linear",
    fill_value: onp.ToComplex = ...,  # np.nan
    rescale: onp.ToBool = False,
) -> onp.Array[onp.AtLeast1D, np.complex128]: ...
@overload
def griddata(
    points: onp.ToFloat1D | onp.ToFloat2D,
    values: onp.ToComplex1D,
    xi: _ToXi,
    method: _Method = "linear",
    fill_value: onp.ToComplex = ...,  # np.nan
    rescale: onp.ToBool = False,
) -> onp.Array[onp.AtLeast1D, Any]: ...
