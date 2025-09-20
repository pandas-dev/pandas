from collections.abc import Callable, Sequence
from typing import Final, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import deprecated

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._optimize import OptimizeResult as _OptimizeResult
from ._typing import Bound, MethodLinprog, MethodLinprogLegacy

__all__ = ["linprog", "linprog_terse_callback", "linprog_verbose_callback"]

###

_Ignored: TypeAlias = object
_Max3: TypeAlias = Literal[0, 1, 2, 3]

_Int: TypeAlias = int | np.int32 | np.int64
_Float: TypeAlias = float | np.float64

@type_check_only
class _OptionsCommon(TypedDict, total=False):
    maxiter: _Int  # default: method-specific
    disp: onp.ToBool  # default: False
    presolve: onp.ToBool  # default: True

# highs-ds
@type_check_only
class _OptionsHighsDS(_OptionsCommon, TypedDict, total=False):
    time_limit: _Float  # default: np.finfo(float).max
    dual_feasibility_tolerance: _Float  # default: 1e-7
    primal_feasibility_tolerance: _Float  # default: 1e-7
    simplex_dual_edge_weight_strategy: Literal["dantzig", "devex", "steepest", "steepest-devex"] | None  # default: None

# highs-ips
@type_check_only
class _OptionsHighsIPM(_OptionsHighsDS, TypedDict, total=False):
    ipm_optimality_tolerance: _Float  # default: 1e-8

# highs
@type_check_only
class _OptionsHighs(_OptionsHighsIPM, TypedDict, total=False):
    min_rel_gap: _Float | None  # default: None

@type_check_only
class _OptionsCommonLegacy(_OptionsCommon, TypedDict, total=False):
    tol: _Float
    autoscale: onp.ToBool  # default: False
    rr: onp.ToBool  # default: True
    rr_method: Literal["SVD", "pivot", "ID", "None"] | None  # default: None

# interior-point (legacy, see https://github.com/scipy/scipy/issues/15707)
@type_check_only
class _OptionsInteriorPoint(_OptionsCommonLegacy, TypedDict, total=False):
    alpha0: _Float  # default: 0.99995
    beta: _Float  # default: 0.1
    sparse: onp.ToBool  # default: False
    lstq: onp.ToBool  # default: False
    sym_pos: onp.ToBool  # default: True
    cholsky: onp.ToBool  # default: True
    pc: onp.ToBool  # default: True
    ip: onp.ToBool  # default: False
    perm_spec: Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"] | None  # default: "MMD_AT_PLUS_A"

# revised simplex (legacy, see https://github.com/scipy/scipy/issues/15707)
@type_check_only
class _OptionsRevisedSimplex(_OptionsCommonLegacy, TypedDict, total=False):
    maxupdate: _Int  # default: 10
    mast: onp.ToBool  # default: False
    pivot: Literal["mrc", "bland"]

# simplex (legacy, see https://github.com/scipy/scipy/issues/15707)
@type_check_only
class _OptionsSimplex(_OptionsCommonLegacy, TypedDict, total=False):
    bland: onp.ToBool  # default: False

###

class OptimizeResult(_OptimizeResult):
    x: onp.Array1D[np.float64] | None
    fun: float | None
    slack: onp.Array1D[np.float64] | None
    con: onp.Array1D[np.float64] | None
    success: bool
    status: Literal[0, 1, 2, 3, 4]
    message: str
    nit: int

@type_check_only
class _OptimizeResultSensitivity(_OptimizeResult):
    residual: onp.Array1D[np.float64]
    marginals: onp.Array1D[np.float64]

@type_check_only
class _OptimizeResultHighs(OptimizeResult):
    crossover_nit: int
    lower: _OptimizeResultSensitivity
    upper: _OptimizeResultSensitivity
    eqlin: _OptimizeResultSensitivity
    ineqlin: _OptimizeResultSensitivity
    mip_node_count: int  # only exists if `success=True`
    mip_dual_bound: float  # only exists if `success=True`
    mip_gap: float  # only exists if `success=True`

###

__docformat__: Final = "restructuredtext en"  # undocumented
LINPROG_METHODS: Final[Sequence[MethodLinprog | MethodLinprogLegacy]] = ...  # undocumented

def linprog_verbose_callback(res: OptimizeResult) -> None: ...
def linprog_terse_callback(res: OptimizeResult) -> None: ...

#
@overload  # highs (default)
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    method: Literal["highs"] = "highs",
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsHighs | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> _OptimizeResultHighs: ...
@overload  # highs-ds
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    *,
    method: Literal["highs-ds"],
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsHighsDS | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> _OptimizeResultHighs: ...
@overload  # highs-ipm
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    *,
    method: Literal["highs-ipm"],
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsHighsIPM | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> _OptimizeResultHighs: ...
@overload  # interior-point (legacy, see https://github.com/scipy/scipy/issues/15707)
@deprecated("`method='interior-point'` is deprecated and will be removed in SciPy 1.17. Please use one of the HIGHS solvers.")
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    *,
    method: Literal["interior-point"],
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsInteriorPoint | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> OptimizeResult: ...
@overload  # revised simplex (legacy, see https://github.com/scipy/scipy/issues/15707)
@deprecated("`method='revised simplex'` is deprecated and will be removed in SciPy 1.17. Please use one of the HIGHS solvers.")
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    *,
    method: Literal["revised simplex"],
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsRevisedSimplex | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> OptimizeResult: ...
@overload  # simplex (legacy, see https://github.com/scipy/scipy/issues/15707)
@deprecated("`method='simplex'` is deprecated and will be removed in SciPy 1.17. Please use one of the HIGHS solvers.")
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    *,
    method: Literal["simplex"],
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsSimplex | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> OptimizeResult: ...
@overload  # any "highs"
def linprog(
    c: onp.ToFloat1D,
    A_ub: onp.ToFloat2D | None = None,
    b_ub: onp.ToFloat1D | None = None,
    A_eq: onp.ToFloat2D | None = None,
    b_eq: onp.ToFloat1D | None = None,
    bounds: Bound | Sequence[Bound] = (0, None),
    method: MethodLinprog = "highs",
    callback: Callable[[OptimizeResult], _Ignored] | None = None,
    options: _OptionsHighs | None = None,
    x0: onp.ToFloat1D | None = None,
    integrality: _Max3 | Sequence[_Max3] | onp.CanArrayND[npc.integer] | None = None,
) -> _OptimizeResultHighs: ...
