from collections.abc import Callable
from typing import Final, Generic, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.interpolate import PPoly
from scipy.sparse import csc_matrix

###

_SCT_fc = TypeVar("_SCT_fc", bound=npc.inexact, default=np.float64 | np.complex128)

_FunRHS: TypeAlias = Callable[[onp.Array1D, onp.Array2D[_SCT_fc]], onp.ArrayND[_SCT_fc]]
_FunRHS_p: TypeAlias = Callable[[onp.Array1D, onp.Array2D[_SCT_fc], onp.Array1D], onp.ArrayND[_SCT_fc]]
_FunRHS_x: TypeAlias = Callable[[onp.Array1D, onp.Array2D[_SCT_fc], onp.Array1D], onp.Array2D[_SCT_fc]]
_FunBCR: TypeAlias = Callable[[onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc]], onp.ArrayND[_SCT_fc]]
_FunBCR_p: TypeAlias = Callable[[onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc], onp.Array1D], onp.ArrayND[_SCT_fc]]
_FunBCR_x: TypeAlias = Callable[[onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc], onp.Array1D], onp.Array1D[_SCT_fc]]
_Funs_x: TypeAlias = tuple[_FunRHS_x[_SCT_fc], _FunBCR_x[_SCT_fc], _FunRHS_jac_x[_SCT_fc], _FunBCR_jac_x[_SCT_fc]]

_FunRHS_jac: TypeAlias = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_SCT_fc]],
    onp.ArrayND[_SCT_fc],
]  # fmt: skip
_FunRHS_jac_p: TypeAlias = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_SCT_fc], onp.Array1D[np.float64]],
    tuple[onp.ArrayND[_SCT_fc], onp.ArrayND[_SCT_fc]],
]  # fmt: skip
_FunRHS_jac_x: TypeAlias = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_SCT_fc], onp.Array1D[np.float64]],
    tuple[onp.Array3D[_SCT_fc], onp.Array3D[_SCT_fc] | None],
]  # fmt: skip

_FunBCR_jac: TypeAlias = Callable[
    [onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc]],
    tuple[onp.ArrayND[_SCT_fc], onp.ArrayND[_SCT_fc]]
]  # fmt: skip
_FunBCR_jac_p: TypeAlias = Callable[
    [onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc], onp.Array1D[np.float64]],
    tuple[onp.ArrayND[_SCT_fc], onp.ArrayND[_SCT_fc], onp.ArrayND[_SCT_fc]],
]
_FunBCR_jac_x: TypeAlias = Callable[
    [onp.Array1D[_SCT_fc], onp.Array1D[_SCT_fc], onp.Array1D[np.float64]],
    tuple[onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc] | None],
]

_FunCol: TypeAlias = Callable[
    [onp.Array2D[_SCT_fc], onp.Array1D[np.float64]],
    tuple[onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc]],
]
_FunCol_jac: TypeAlias = Callable[
    [
        onp.Array1D[_SCT_fc],
        onp.Array1D[_SCT_fc],
        onp.Array2D[_SCT_fc],
        onp.Array2D[_SCT_fc],
        onp.Array2D[_SCT_fc],
        onp.Array1D[_SCT_fc],
    ],
    csc_matrix,
]

###

EPS: Final[float] = ...
TERMINATION_MESSAGES: Final[dict[Literal[0, 1, 2, 3], str]] = ...

# NOTE: this inherits from `scipy.optimize.OptimizeResult` at runtime.
# But because `BVPResult` doesn't share all members (and optional attributes
# still aren't a thing), it was omitted as a base class here.
class BVPResult(Generic[_SCT_fc]):
    sol: Final[PPoly]
    p: Final[onp.Array1D[np.float64] | None]
    x: Final[onp.Array1D[np.float64]]
    rms_residuals: Final[onp.Array1D[np.float64]]
    niter: Final[int]
    status: Final[Literal[0, 1, 2]]
    message: Final[str]
    success: Final[bool]

    y: onp.Array2D[_SCT_fc]
    yp: onp.Array2D[_SCT_fc]

def estimate_fun_jac(
    fun: _FunRHS_x[_SCT_fc],
    x: onp.Array1D[np.float64],
    y: onp.Array2D[_SCT_fc],
    p: onp.Array1D[np.float64],
    f0: onp.Array2D[_SCT_fc] | None = None,
) -> tuple[onp.Array3D[_SCT_fc], onp.Array3D[_SCT_fc] | None]: ...  # undocumented
def estimate_bc_jac(
    bc: _FunBCR_x[_SCT_fc],
    ya: onp.Array1D[_SCT_fc],
    yb: onp.Array1D[_SCT_fc],
    p: onp.Array1D[np.float64],
    bc0: onp.Array1D[_SCT_fc] | None = None,
) -> tuple[onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc] | None]: ...  # undocumented
def compute_jac_indices(n: int, m: int, k: int) -> tuple[onp.Array1D[np.intp], onp.Array1D[np.intp]]: ...  # undocumented
def stacked_matmul(a: onp.ArrayND[_SCT_fc], b: onp.ArrayND[_SCT_fc]) -> onp.ArrayND[_SCT_fc]: ...  # undocumented
def construct_global_jac(
    n: int,
    m: int,
    k: int,
    i_jac: onp.Array1D[np.intp],
    j_jac: onp.Array1D[np.intp],
    h: float,
    df_dy: onp.Array3D[_SCT_fc],
    df_dy_middle: onp.Array3D[_SCT_fc],
    df_dp: onp.Array3D[_SCT_fc] | None,
    df_dp_middle: onp.Array3D[_SCT_fc] | None,
    dbc_dya: onp.Array2D[_SCT_fc],
    dbc_dyb: onp.Array2D[_SCT_fc],
    dbc_dp: onp.Array2D[_SCT_fc] | None,
) -> csc_matrix: ...  # undocumented
def collocation_fun(
    fun: _FunRHS_x[_SCT_fc], y: onp.Array2D[_SCT_fc], p: onp.Array1D[np.float64], x: onp.Array1D[np.float64], h: float
) -> tuple[onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc], onp.Array2D[_SCT_fc]]: ...  # undocumented
def prepare_sys(
    n: int,
    m: int,
    k: int,
    fun: _FunRHS_x[_SCT_fc],
    bc: _FunBCR_x[_SCT_fc],
    fun_jac: _FunRHS_jac_x[_SCT_fc] | None,
    bc_jac: _FunBCR_jac_x[_SCT_fc] | None,
    x: onp.Array1D[np.float64],
    h: float,
) -> tuple[_FunCol[_SCT_fc], _FunCol_jac[_SCT_fc]]: ...  # undocumented
def solve_newton(
    n: int,
    m: int,
    h: int,
    col_fun: _FunCol[_SCT_fc],
    bc: _FunBCR_x[_SCT_fc],
    jac: _FunCol_jac[_SCT_fc],
    y: onp.Array2D[_SCT_fc],
    p: onp.Array1D[np.float64],
    B: onp.Array2D[np.float64] | None,
    bvp_tol: float,
    bc_tol: float,
) -> tuple[onp.Array2D[_SCT_fc], onp.Array1D[np.float64], bool]: ...  # undocumented
def print_iteration_header() -> None: ...  # undocumented
def print_iteration_progress(
    iteration: int, residual: complex, bc_residual: complex, total_nodes: int, nodes_added: int
) -> None: ...  # undocumented
def estimate_rms_residuals(
    fun: _FunRHS_x[_SCT_fc],
    sol: PPoly,
    x: onp.Array1D,
    h: float,
    p: onp.Array1D,
    r_middle: onp.Array2D[_SCT_fc],
    f_middle: onp.Array2D[_SCT_fc],
) -> onp.Array1D: ...  # undocumented
def create_spline(
    y: onp.Array2D[_SCT_fc], yp: onp.Array2D[_SCT_fc], x: onp.Array1D[np.float64], h: float
) -> PPoly: ...  # undocumented
def modify_mesh(
    x: onp.Array1D, insert_1: onp.Array1D[np.intp], insert_2: onp.Array1D[np.intp]
) -> onp.Array1D: ...  # undocumented

#
@overload
def wrap_functions(
    fun: _FunRHS[_SCT_fc],
    bc: _FunBCR[_SCT_fc],
    fun_jac: _FunRHS_jac[_SCT_fc] | None,
    bc_jac: _FunBCR_jac[_SCT_fc] | None,
    k: onp.ToFalse,
    a: onp.ToFloat,
    S: onp.Array2D[np.float64] | None,
    D: onp.Array2D[np.float64] | None,
    dtype: type[float | complex],
) -> _Funs_x[_SCT_fc]: ...  # undocumented
@overload
def wrap_functions(
    fun: _FunRHS_p[_SCT_fc],
    bc: _FunBCR_p[_SCT_fc],
    fun_jac: _FunRHS_jac_p[_SCT_fc] | None,
    bc_jac: _FunBCR_jac_p[_SCT_fc] | None,
    k: onp.ToTrue,
    a: onp.ToFloat,
    S: onp.Array2D[np.float64] | None,
    D: onp.Array2D[np.float64] | None,
    dtype: type[float | complex],
) -> _Funs_x[_SCT_fc]: ...  # undocumented

#
@overload
def solve_bvp(
    fun: _FunRHS[_SCT_fc],
    bc: _FunBCR[_SCT_fc],
    x: onp.ToFloat1D,
    y: onp.ToComplex2D,
    p: None = None,
    S: onp.ToFloat2D | None = None,
    fun_jac: _FunRHS_jac[_SCT_fc] | None = None,
    bc_jac: _FunBCR_jac[_SCT_fc] | None = None,
    tol: float = 0.001,
    max_nodes: int = 1_000,
    verbose: Literal[0, 1, 2] = 0,
    bc_tol: float | None = None,
) -> BVPResult[_SCT_fc]: ...
@overload
def solve_bvp(
    fun: _FunRHS_p[_SCT_fc],
    bc: _FunBCR_p[_SCT_fc],
    x: onp.ToFloat1D,
    y: onp.ToComplex2D,
    p: onp.ToFloat1D,
    S: onp.ToFloat2D | None = None,
    fun_jac: _FunRHS_jac_p[_SCT_fc] | None = None,
    bc_jac: _FunBCR_jac_p[_SCT_fc] | None = None,
    tol: float = 0.001,
    max_nodes: int = 1_000,
    verbose: Literal[0, 1, 2] = 0,
    bc_tol: float | None = None,
) -> BVPResult[_SCT_fc]: ...
