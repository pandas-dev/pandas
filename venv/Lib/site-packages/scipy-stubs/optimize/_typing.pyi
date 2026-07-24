# type-check-only helper types for internal use

from collections.abc import Callable, Sequence
from typing import Concatenate, Literal, NotRequired, Required, type_check_only
from typing_extensions import TypedDict

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds as _Bounds, LinearConstraint, NonlinearConstraint
from ._hessian_update_strategy import HessianUpdateStrategy
from ._minimize import _MinimizeOptions

__all__ = [
    "Bound",
    "Bounds",
    "Brack",
    "Constraint",
    "Constraints",
    "MethodAll",
    "MethodLinprog",
    "MethodLinprogLegacy",
    "MethodMimimize",
    "MethodMinimizeScalar",
    "MethodRootScalar",
    "MinimizerKwargs",
    "MinimizerKwargsJac",
    "Solver",
    "TRSolver",
]

type _Float1D = onp.Array1D[np.float64]
type _Args = tuple[object, ...]

# bounds
type Bound = tuple[onp.ToFloat | None, onp.ToFloat | None]
type Bounds = Sequence[Bound] | _Bounds

# constaints
@type_check_only
class _ConstraintDict(TypedDict):
    type: Literal["eq", "ineq"]
    fun: Callable[Concatenate[_Float1D, ...], onp.ToFloat]
    jac: NotRequired[Callable[Concatenate[_Float1D, ...], onp.ToFloat1D]]
    args: NotRequired[_Args]

type Constraint = LinearConstraint | NonlinearConstraint | _ConstraintDict
type Constraints = Constraint | Sequence[Constraint]

type Brack = tuple[onp.ToFloat, onp.ToFloat] | tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat]

type Solver = Literal["minimize", "minimize_scalar", "root", "root_salar", "linprog", "quadratic_assignment"]
type TRSolver = Literal["exact", "lsmr"]

type MethodMimimize = Literal[
    "Nelder-Mead", "nelder-mead",
    "Powell", "powell",
    "CG", "cg",
    "BFGS", "bfgs",
    "Newton-CG", "newton-cg",
    "L-BFGS-B", "l-bfgs-b",
    "TNC", "tnc",
    "COBYLA", "cobyla",
    "COBYQA", "cobyqa",
    "SLSQP", "slsqp",
    "Trust-Constr", "trust-constr",
    "Dogleg", "dogleg",
    "Trust-NCG", "trust-ncg",
    "Trust-Exact", "trust-exact",
    "Trust-Krylov", "trust-krylov",
]  # fmt: skip
type MethodMinimizeScalar = Literal["brent", "golden", "bounded"]
type MethodLinprog = Literal["highs", "highs-ds", "highs-ipm"]
type MethodLinprogLegacy = Literal["interior-point", "revised simplex", "simplex"]
type _MethodRoot = Literal[
    "hybr",
    "lm",
    "broyden1",
    "broyden2",
    "anderson",
    "linearmixing",
    "diagbroyden",
    "excitingmixing",
    "krylov",
    "df-sane",
]  # fmt: skip
type MethodRootScalar = Literal["bisect", "brentq", "brenth", "ridder", "toms748", "newton", "secant", "halley"]
type _MethodQuadraticAssignment = Literal["faq", "2opt"]
type MethodAll = Literal[
    MethodMimimize,
    _MethodRoot,
    MethodMinimizeScalar,
    MethodRootScalar,
    MethodLinprog,
    _MethodQuadraticAssignment,
]  # fmt: skip

type _FDMethod = Literal["2-point", "3-point", "cs"]

@type_check_only
class _MinimizerKwargsBase(TypedDict, total=False):
    args: _Args
    method: MethodMimimize
    hess: Callable[Concatenate[_Float1D, ...], onp.ToFloat2D] | _FDMethod | HessianUpdateStrategy
    hessp: Callable[Concatenate[_Float1D, _Float1D, ...], onp.ToFloat1D]
    bounds: Bounds
    constraints: Constraints
    tol: onp.ToFloat
    options: _MinimizeOptions

@type_check_only
class MinimizerKwargs(_MinimizerKwargsBase):
    jac: NotRequired[Callable[Concatenate[_Float1D, ...], onp.ToFloat1D] | _FDMethod | Literal[False]]

@type_check_only
class MinimizerKwargsJac(_MinimizerKwargsBase):
    jac: Required[Literal[True]]
