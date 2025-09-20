from collections.abc import Sequence
from typing import Literal, LiteralString, TypeAlias, TypedDict, type_check_only
from typing_extensions import TypeAliasType

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._constraints import Bounds, LinearConstraint
from ._optimize import OptimizeResult as _OptimizeResult
from scipy.sparse._base import _spbase

###

_Max3: TypeAlias = Literal[0, 1, 2, 3]
_Max4: TypeAlias = Literal[_Max3, 4]

_Float: TypeAlias = float | np.float64
_ToFloat2D: TypeAlias = onp.ToFloat2D | _spbase[np.bool_ | npc.integer | np.float32 | np.float64 | npc.floating80]
# The `TypeAliasType` helps make error messages less unreadable
_ToLinearConstraint = TypeAliasType(
    "_ToLinearConstraint",
    LinearConstraint
    | tuple[_ToFloat2D]
    | tuple[_ToFloat2D, onp.ToFloat | onp.ToFloat1D]
    | tuple[_ToFloat2D, onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D],
)

@type_check_only
class _OptionsMILP(TypedDict, total=False):
    disp: onp.ToBool  # default: False
    presolve: onp.ToBool  # default: True
    node_limit: int  # default: no limit
    time_limit: onp.ToFloat  # default: no limit
    mip_rel_gap: onp.ToFloat  # default: ?

@type_check_only
class _OptimizeResultSucess(_OptimizeResult):
    x: onp.Array1D[np.float64]
    fun: _Float
    status: Literal[0]
    message: LiteralString
    success: Literal[True]
    mip_node_count: int
    mip_dual_bound: _Float
    mip_gap: _Float

@type_check_only
class _OptimizeResultNoSolution(_OptimizeResult):
    x: None
    fun: None
    status: Literal[1, 2, 3, 4]
    message: LiteralString
    success: Literal[False]
    mip_node_count: int | None
    mip_dual_bound: _Float | None
    mip_gap: _Float | None

###

# only used as implicit export
class OptimizeResult(_OptimizeResult):
    x: onp.Array1D[np.float64] | None
    fun: _Float | None
    status: _Max4
    message: LiteralString
    success: bool  # status == 0
    mip_node_count: int | None
    mip_dual_bound: _Float | None
    mip_gap: _Float | None

def milp(
    c: onp.ToFloat1D,
    *,
    integrality: Sequence[_Max3 | npc.integer] | onp.CanArrayND[npc.integer] | None = None,
    bounds: Bounds | None = None,
    constraints: _ToLinearConstraint | Sequence[_ToLinearConstraint] | None = None,
    options: _OptionsMILP | None = None,
) -> _OptimizeResultSucess | _OptimizeResultNoSolution: ...
