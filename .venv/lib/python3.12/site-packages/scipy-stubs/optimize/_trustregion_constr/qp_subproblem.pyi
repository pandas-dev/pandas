from collections.abc import Sequence
from typing import Literal, NotRequired, TypeAlias, TypedDict, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

__all__ = [
    "box_intersections",
    "box_sphere_intersections",
    "eqp_kktfact",
    "inside_box_boundaries",
    "modified_dogleg",
    "projected_cg",
    "sphere_intersections",
]

_ScalarB1: TypeAlias = bool | np.bool_
_ScalarF8: TypeAlias = float | np.float64
_VectorF8: TypeAlias = onp.Array1D[np.float64]

_ScalarInt_co: TypeAlias = npc.integer
_ScalarFloat_co: TypeAlias = npc.floating | _ScalarInt_co

_ScalarLikeInt_co: TypeAlias = int | _ScalarInt_co
_ScalarLikeFloat_co: TypeAlias = float | _ScalarFloat_co
_VectorLikeFloat_co: TypeAlias = Sequence[_ScalarLikeFloat_co] | onp.CanArray1D[_ScalarFloat_co]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_SCT_float = TypeVar("_SCT_float", bound=_ScalarFloat_co)

_SparseArray: TypeAlias = sparray | spmatrix

class _SphereInfoDict(TypedDict):
    ta: _ScalarF8
    tb: _ScalarF8
    intersect: bool

class _ProjectedCGDict(TypedDict):
    niter: int
    stop_cond: Literal[1, 2, 3, 4]
    hits_boundary: bool
    allvecs: NotRequired[Sequence[_VectorF8]]

def eqp_kktfact(
    H: _SparseArray, c: _VectorLikeFloat_co, A: _SparseArray, b: _VectorLikeFloat_co
) -> tuple[_VectorF8, _VectorF8]: ...
def sphere_intersections(
    z: _VectorLikeFloat_co, d: _VectorLikeFloat_co, trust_radius: _ScalarLikeFloat_co, entire_line: _ScalarB1 = False
) -> tuple[_ScalarF8, _ScalarF8, _ScalarB1]: ...
def box_intersections(
    z: _VectorLikeFloat_co,
    d: _VectorLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
    entire_line: _ScalarB1 = False,
) -> tuple[_ScalarF8, _ScalarF8, _ScalarB1]: ...
@overload
def box_sphere_intersections(
    z: _VectorLikeFloat_co,
    d: _VectorLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co,
    entire_line: _ScalarB1 = False,
    extra_info: onp.ToFalse | None = False,
) -> tuple[_ScalarF8, _ScalarF8, _ScalarB1]: ...
@overload
def box_sphere_intersections(
    z: _VectorLikeFloat_co,
    d: _VectorLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co,
    entire_line: _ScalarB1,
    extra_info: onp.ToTrue,
) -> tuple[_ScalarF8, _ScalarF8, _ScalarB1, _SphereInfoDict, _SphereInfoDict]: ...
@overload
def box_sphere_intersections(
    z: _VectorLikeFloat_co,
    d: _VectorLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co,
    entire_line: _ScalarB1 = False,
    *,
    extra_info: onp.ToTrue,
) -> tuple[_ScalarF8, _ScalarF8, _ScalarB1, _SphereInfoDict, _SphereInfoDict]: ...
def inside_box_boundaries(
    x: onp.Array[_ShapeT, _ScalarFloat_co], lb: onp.Array[_ShapeT, _ScalarFloat_co], ub: onp.Array[_ShapeT, _ScalarFloat_co]
) -> np.bool_: ...
def reinforce_box_boundaries(
    x: onp.Array[_ShapeT, _SCT_float], lb: onp.Array[_ShapeT, _SCT_float], ub: onp.Array[_ShapeT, _SCT_float]
) -> onp.Array[_ShapeT, _SCT_float]: ...
def modified_dogleg(
    A: LinearOperator | _SparseArray | onp.ArrayND[_ScalarFloat_co],
    Y: LinearOperator | _SparseArray | onp.ArrayND[_ScalarFloat_co],
    b: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
) -> _VectorF8: ...
def projected_cg(
    H: LinearOperator | _SparseArray | onp.ArrayND[_ScalarFloat_co],
    c: _VectorLikeFloat_co,
    Z: LinearOperator | _SparseArray | onp.ArrayND[_ScalarFloat_co],
    Y: LinearOperator | _SparseArray | onp.ArrayND[_ScalarFloat_co],
    b: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co = ...,
    lb: _ScalarLikeFloat_co | None = None,
    ub: _ScalarLikeFloat_co | None = None,
    tol: _ScalarLikeFloat_co | None = None,
    max_iter: _ScalarLikeInt_co | None = None,
    max_infeasible_iter: _ScalarLikeInt_co | npc.integer | None = None,
    return_all: _ScalarB1 = False,
) -> tuple[_VectorF8, _ProjectedCGDict]: ...
