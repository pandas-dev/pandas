from collections.abc import Callable
from typing import Final, Generic, Literal, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.interpolate import PPoly
from scipy.sparse import csc_matrix

###

_InexactT = TypeVar("_InexactT", bound=npc.inexact, default=np.float64 | np.complex128)

type _FunRHS[_InexactT: npc.inexact] = Callable[[onp.Array1D, onp.Array2D[_InexactT]], onp.ArrayND[_InexactT]]
type _FunRHS_p[_InexactT: npc.inexact] = Callable[[onp.Array1D, onp.Array2D[_InexactT], onp.Array1D], onp.ArrayND[_InexactT]]
type _FunRHS_x[_InexactT: npc.inexact] = Callable[[onp.Array1D, onp.Array2D[_InexactT], onp.Array1D], onp.Array2D[_InexactT]]
type _FunBCR[_InexactT: npc.inexact] = Callable[[onp.Array1D[_InexactT], onp.Array1D[_InexactT]], onp.ArrayND[_InexactT]]
type _FunBCR_p[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[_InexactT], onp.Array1D[_InexactT], onp.Array1D],
    onp.ArrayND[_InexactT],
]  # fmt: skip
type _FunBCR_x[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[_InexactT], onp.Array1D[_InexactT], onp.Array1D],
    onp.Array1D[_InexactT],
]  # fmt: skip
type _Funs_x[_InexactT: npc.inexact] = tuple[
    _FunRHS_x[_InexactT],
    _FunBCR_x[_InexactT],
    _FunRHS_jac_x[_InexactT],
    _FunBCR_jac_x[_InexactT],
]  # fmt: skip

type _FunRHS_jac[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_InexactT]],
    onp.ArrayND[_InexactT],
]  # fmt: skip
type _FunRHS_jac_p[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_InexactT], onp.Array1D[np.float64]],
    tuple[onp.ArrayND[_InexactT], onp.ArrayND[_InexactT]],
]  # fmt: skip
type _FunRHS_jac_x[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[np.float64], onp.Array2D[_InexactT], onp.Array1D[np.float64]],
    tuple[onp.Array3D[_InexactT], onp.Array3D[_InexactT] | None],
]  # fmt: skip

type _FunBCR_jac[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[_InexactT], onp.Array1D[_InexactT]],
    tuple[onp.ArrayND[_InexactT], onp.ArrayND[_InexactT]]
]  # fmt: skip
type _FunBCR_jac_p[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[_InexactT], onp.Array1D[_InexactT], onp.Array1D[np.float64]],
    tuple[onp.ArrayND[_InexactT], onp.ArrayND[_InexactT], onp.ArrayND[_InexactT]],
]
type _FunBCR_jac_x[_InexactT: npc.inexact] = Callable[
    [onp.Array1D[_InexactT], onp.Array1D[_InexactT], onp.Array1D[np.float64]],
    tuple[onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT] | None],
]

type _FunCol[_InexactT: npc.inexact] = Callable[
    [onp.Array2D[_InexactT], onp.Array1D[np.float64]],
    tuple[onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT]],
]
type _FunCol_jac[_InexactT: npc.inexact] = Callable[
    [
        onp.Array1D[_InexactT],
        onp.Array1D[_InexactT],
        onp.Array2D[_InexactT],
        onp.Array2D[_InexactT],
        onp.Array2D[_InexactT],
        onp.Array1D[_InexactT],
    ],
    csc_matrix,
]

###

EPS: Final[float] = ...
TERMINATION_MESSAGES: Final[dict[Literal[0, 1, 2, 3], str]] = ...

# NOTE: this inherits from `scipy.optimize.OptimizeResult` at runtime.
# But because `BVPResult` doesn't share all members (and optional attributes
# still aren't a thing), it was omitted as a base class here.
class BVPResult(Generic[_InexactT]):
    sol: Final[PPoly]
    p: Final[onp.Array1D[np.float64] | None]
    x: Final[onp.Array1D[np.float64]]
    rms_residuals: Final[onp.Array1D[np.float64]]
    niter: Final[int]
    status: Final[Literal[0, 1, 2]]
    message: Final[str]
    success: Final[bool]

    y: onp.Array2D[_InexactT]
    yp: onp.Array2D[_InexactT]

def estimate_fun_jac(
    fun: _FunRHS_x[_InexactT],
    x: onp.Array1D[np.float64],
    y: onp.Array2D[_InexactT],
    p: onp.Array1D[np.float64],
    f0: onp.Array2D[_InexactT] | None = None,
) -> tuple[onp.Array3D[_InexactT], onp.Array3D[_InexactT] | None]: ...  # undocumented
def estimate_bc_jac(
    bc: _FunBCR_x[_InexactT],
    ya: onp.Array1D[_InexactT],
    yb: onp.Array1D[_InexactT],
    p: onp.Array1D[np.float64],
    bc0: onp.Array1D[_InexactT] | None = None,
) -> tuple[onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT] | None]: ...  # undocumented
def compute_jac_indices(n: int, m: int, k: int) -> tuple[onp.Array1D[np.intp], onp.Array1D[np.intp]]: ...  # undocumented
def stacked_matmul(a: onp.ArrayND[_InexactT], b: onp.ArrayND[_InexactT]) -> onp.ArrayND[_InexactT]: ...  # undocumented
def construct_global_jac(
    n: int,
    m: int,
    k: int,
    i_jac: onp.Array1D[np.intp],
    j_jac: onp.Array1D[np.intp],
    h: float,
    df_dy: onp.Array3D[_InexactT],
    df_dy_middle: onp.Array3D[_InexactT],
    df_dp: onp.Array3D[_InexactT] | None,
    df_dp_middle: onp.Array3D[_InexactT] | None,
    dbc_dya: onp.Array2D[_InexactT],
    dbc_dyb: onp.Array2D[_InexactT],
    dbc_dp: onp.Array2D[_InexactT] | None,
) -> csc_matrix: ...  # undocumented
def collocation_fun(
    fun: _FunRHS_x[_InexactT], y: onp.Array2D[_InexactT], p: onp.Array1D[np.float64], x: onp.Array1D[np.float64], h: float
) -> tuple[onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT]]: ...  # undocumented
def prepare_sys(
    n: int,
    m: int,
    k: int,
    fun: _FunRHS_x[_InexactT],
    bc: _FunBCR_x[_InexactT],
    fun_jac: _FunRHS_jac_x[_InexactT] | None,
    bc_jac: _FunBCR_jac_x[_InexactT] | None,
    x: onp.Array1D[np.float64],
    h: float,
) -> tuple[_FunCol[_InexactT], _FunCol_jac[_InexactT]]: ...  # undocumented
def solve_newton(
    n: int,
    m: int,
    h: int,
    col_fun: _FunCol[_InexactT],
    bc: _FunBCR_x[_InexactT],
    jac: _FunCol_jac[_InexactT],
    y: onp.Array2D[_InexactT],
    p: onp.Array1D[np.float64],
    B: onp.Array2D[np.float64] | None,
    bvp_tol: float,
    bc_tol: float,
) -> tuple[onp.Array2D[_InexactT], onp.Array1D[np.float64], bool]: ...  # undocumented
def print_iteration_header() -> None: ...  # undocumented
def print_iteration_progress(
    iteration: int, residual: complex, bc_residual: complex, total_nodes: int, nodes_added: int
) -> None: ...  # undocumented
def estimate_rms_residuals(
    fun: _FunRHS_x[_InexactT],
    sol: PPoly,
    x: onp.Array1D,
    h: float,
    p: onp.Array1D,
    r_middle: onp.Array2D[_InexactT],
    f_middle: onp.Array2D[_InexactT],
) -> onp.Array1D: ...  # undocumented
def create_spline(
    y: onp.Array2D[_InexactT], yp: onp.Array2D[_InexactT], x: onp.Array1D[np.float64], h: float
) -> PPoly: ...  # undocumented
def modify_mesh(
    x: onp.Array1D, insert_1: onp.Array1D[np.intp], insert_2: onp.Array1D[np.intp]
) -> onp.Array1D: ...  # undocumented

#
@overload
def wrap_functions(
    fun: _FunRHS[_InexactT],
    bc: _FunBCR[_InexactT],
    fun_jac: _FunRHS_jac[_InexactT] | None,
    bc_jac: _FunBCR_jac[_InexactT] | None,
    k: onp.ToFalse,
    a: onp.ToFloat,
    S: onp.Array2D[np.float64] | None,
    D: onp.Array2D[np.float64] | None,
    dtype: type[float | complex],
) -> _Funs_x[_InexactT]: ...  # undocumented
@overload
def wrap_functions(
    fun: _FunRHS_p[_InexactT],
    bc: _FunBCR_p[_InexactT],
    fun_jac: _FunRHS_jac_p[_InexactT] | None,
    bc_jac: _FunBCR_jac_p[_InexactT] | None,
    k: onp.ToTrue,
    a: onp.ToFloat,
    S: onp.Array2D[np.float64] | None,
    D: onp.Array2D[np.float64] | None,
    dtype: type[float | complex],
) -> _Funs_x[_InexactT]: ...  # undocumented

#
@overload
def solve_bvp(
    fun: _FunRHS[_InexactT],
    bc: _FunBCR[_InexactT],
    x: onp.ToFloat1D,
    y: onp.ToComplex2D,
    p: None = None,
    S: onp.ToFloat2D | None = None,
    fun_jac: _FunRHS_jac[_InexactT] | None = None,
    bc_jac: _FunBCR_jac[_InexactT] | None = None,
    tol: float = 0.001,
    max_nodes: int = 1_000,
    verbose: Literal[0, 1, 2] = 0,
    bc_tol: float | None = None,
) -> BVPResult[_InexactT]: ...
@overload
def solve_bvp(
    fun: _FunRHS_p[_InexactT],
    bc: _FunBCR_p[_InexactT],
    x: onp.ToFloat1D,
    y: onp.ToComplex2D,
    p: onp.ToFloat1D,
    S: onp.ToFloat2D | None = None,
    fun_jac: _FunRHS_jac_p[_InexactT] | None = None,
    bc_jac: _FunBCR_jac_p[_InexactT] | None = None,
    tol: float = 0.001,
    max_nodes: int = 1_000,
    verbose: Literal[0, 1, 2] = 0,
    bc_tol: float | None = None,
) -> BVPResult[_InexactT]: ...
