from collections.abc import Sequence
from typing import Literal, NotRequired, TypedDict, overload, type_check_only

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

###

type _ScalarB1 = bool | np.bool
type _ScalarF8 = float | np.float64
type _VectorF8 = onp.Array1D[np.float64]

type _CoFloating = npc.floating | npc.integer

type _ScalarLikeInt_co = int | npc.integer
type _ScalarLikeFloat_co = float | _CoFloating
type _VectorLikeFloat_co = Sequence[_ScalarLikeFloat_co] | onp.CanArray1D[_CoFloating]

type _SparseArray = sparray | spmatrix

@type_check_only
class _SphereInfoDict(TypedDict):
    ta: _ScalarF8
    tb: _ScalarF8
    intersect: bool

@type_check_only
class _ProjectedCGDict(TypedDict):
    niter: int
    stop_cond: Literal[1, 2, 3, 4]
    hits_boundary: bool
    allvecs: NotRequired[Sequence[_VectorF8]]

###

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
def inside_box_boundaries[ShapeT: tuple[int, ...]](
    x: onp.Array[ShapeT, _CoFloating], lb: onp.Array[ShapeT, _CoFloating], ub: onp.Array[ShapeT, _CoFloating]
) -> np.bool: ...
def reinforce_box_boundaries[ShapeT: tuple[int, ...], ScalarT: _CoFloating](
    x: onp.Array[ShapeT, ScalarT], lb: onp.Array[ShapeT, ScalarT], ub: onp.Array[ShapeT, ScalarT]
) -> onp.Array[ShapeT, ScalarT]: ...
def modified_dogleg(
    A: LinearOperator | _SparseArray | onp.ArrayND[_CoFloating],
    Y: LinearOperator | _SparseArray | onp.ArrayND[_CoFloating],
    b: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co,
    lb: _VectorLikeFloat_co,
    ub: _VectorLikeFloat_co,
) -> _VectorF8: ...
def projected_cg(
    H: LinearOperator | _SparseArray | onp.ArrayND[_CoFloating],
    c: _VectorLikeFloat_co,
    Z: LinearOperator | _SparseArray | onp.ArrayND[_CoFloating],
    Y: LinearOperator | _SparseArray | onp.ArrayND[_CoFloating],
    b: _VectorLikeFloat_co,
    trust_radius: _ScalarLikeFloat_co = ...,
    lb: _ScalarLikeFloat_co | None = None,
    ub: _ScalarLikeFloat_co | None = None,
    tol: _ScalarLikeFloat_co | None = None,
    max_iter: _ScalarLikeInt_co | None = None,
    max_infeasible_iter: _ScalarLikeInt_co | npc.integer | None = None,
    return_all: _ScalarB1 = False,
) -> tuple[_VectorF8, _ProjectedCGDict]: ...
