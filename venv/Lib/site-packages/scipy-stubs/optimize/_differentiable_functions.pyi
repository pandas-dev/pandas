from collections.abc import Callable, Iterable, Sequence
from types import ModuleType
from typing import Concatenate, Final, Generic, Literal, Protocol, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._hessian_update_strategy import HessianUpdateStrategy
from scipy.sparse import csr_array, sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

###

type _ToFloat64Vec = Sequence[float | np.float64 | npc.integer | np.bool] | onp.CanArrayND[np.float64]
type _ToJac = onp.ToFloat2D | spmatrix | sparray
type _ToHess = _ToJac | LinearOperator

type _Vec[FloatT: npc.floating] = onp.Array1D[FloatT]
type _Jac[FloatT: npc.floating] = onp.Array2D[FloatT] | csr_array[FloatT]
type _Hess[FloatT: npc.floating] = _Jac[FloatT] | LinearOperator

type _ScalarFun[FloatT: npc.floating] = Callable[Concatenate[onp.Array1D[FloatT], ...], onp.ToFloat]
type _VectorFun[FloatT: npc.floating] = Callable[Concatenate[onp.Array1D[FloatT], ...], onp.ToFloat1D]
type _JacFun[FloatT: npc.floating] = Callable[Concatenate[onp.Array1D[FloatT], ...], _ToJac]
type _HessFun[FloatT: npc.floating] = Callable[Concatenate[onp.Array1D[FloatT], ...], _ToHess]

type _FDMethod = Literal["2-point", "3-point", "cs"]
type _FDBounds = onp.ToFloat1D | onp.ToFloat2D  # len-2 array-like of scalar- or vector-likes

type _ToGradFun[FloatT: npc.floating] = _VectorFun[FloatT] | _FDMethod
type _ToJacFun[FloatT: npc.floating] = _JacFun[FloatT] | _FDMethod
type _ToHessFun[FloatT: npc.floating] = _HessFun[FloatT] | _FDMethod | HessianUpdateStrategy

@type_check_only
class _DoesMap(Protocol):
    def __call__[VT, RT](self, func: Callable[[VT], RT], iterable: Iterable[VT], /) -> Iterable[RT]: ...

type _Workers = int | _DoesMap

_XT_contra = TypeVar("_XT_contra", bound=npc.floating, default=npc.floating, contravariant=True)

###

FD_METHODS: Final = "2-point", "3-point", "cs"

class ScalarFunction(Generic[_XT_contra]):
    xp: Final[ModuleType]

    _args: Final[tuple[object, ...]]

    n: Final[int]
    x_dtype: _XT_contra  # readonly
    x: _Vec[_XT_contra]
    x_prev: _Vec[_XT_contra] | None
    _lowest_x: _Vec[_XT_contra] | None

    _orig_fun: _ScalarFun[_XT_contra]  # readonly
    _wrapped_fun: Callable[[onp.Array1D[_XT_contra]], onp.ToFloat]  # readonly
    _nfev: Final[list[int]]  # size 1
    _lowest_f: onp.ToFloat
    f_updated: bool

    _orig_grad: _ToGradFun[_XT_contra]  # readonly
    _wrapped_grad: Callable[[onp.Array1D[_XT_contra]], _Vec[npc.floating]]  # readonly
    _ngev: Final[list[int]]  # size 1
    g_prev: _Vec[npc.floating] | None
    g_updated: bool

    _orig_hess: _ToHessFun[_XT_contra]  # readonly
    _wrapped_hess: Callable[[onp.Array1D[_XT_contra]], _Hess[npc.floating]]  # readonly
    _nhev: Final[list[int]]  # size 1
    H: _Hess[npc.floating]
    H_updated: bool

    #
    @property
    def nfev(self, /) -> int: ...
    @property
    def ngev(self, /) -> int: ...
    @property
    def nhev(self, /) -> int: ...

    #
    @overload
    def __init__(
        self: ScalarFunction[np.float64],
        /,
        fun: _ScalarFun[np.float64],
        x0: _ToFloat64Vec,
        args: tuple[object, ...],
        grad: _ToGradFun[np.float64],
        hess: _ToHessFun[np.float64],
        finite_diff_rel_step: onp.ToFloat | onp.ToFloat1D | None = None,
        finite_diff_bounds: _FDBounds = ...,
        epsilon: onp.ToFloat | onp.ToFloat1D | None = None,
        workers: int | _DoesMap | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        fun: _ScalarFun[_XT_contra],
        x0: _Vec[_XT_contra] | onp.ToFloat1D,
        args: tuple[object, ...],
        grad: _ToGradFun[_XT_contra],
        hess: _ToHessFun[_XT_contra],
        finite_diff_rel_step: onp.ToFloat | onp.ToFloat1D | None = None,
        finite_diff_bounds: _FDBounds = ...,
        epsilon: onp.ToFloat | onp.ToFloat1D | None = None,
        workers: int | _DoesMap | None = None,
    ) -> None: ...

    #
    def _update_x(self, /, x: onp.ToFloat1D) -> None: ...
    def _update_fun(self, /) -> None: ...
    def _update_grad(self, /) -> None: ...
    def _update_hess(self, /) -> None: ...

    #
    def fun(self, /, x: onp.ToFloat1D) -> float | npc.floating: ...
    def grad(self, /, x: onp.ToFloat1D) -> _Vec[npc.floating]: ...
    def hess(self, /, x: onp.ToFloat1D) -> _Hess[npc.floating]: ...
    def fun_and_grad(self, /, x: onp.ToFloat1D) -> tuple[float | npc.floating, _Vec[npc.floating]]: ...

class VectorFunction(Generic[_XT_contra]):
    xp: Final[ModuleType]

    n: Final[int]
    x_dtype: _XT_contra  # readonly
    x: _Vec[_XT_contra]
    x_diff: _Vec[_XT_contra]
    x_prev: _Vec[_XT_contra] | None

    m: Final[int]
    f: _Vec[npc.floating]
    v: _Vec[npc.floating]
    f_updated: bool

    sparse_jacobian: Final[bool]
    J: _Jac[npc.floating]
    J_prev: _Jac[npc.floating] | None
    J_updated: bool

    H: _Hess[npc.floating]
    H_updated: bool

    @overload
    def __init__(
        self: VectorFunction[np.float64],
        /,
        fun: _VectorFun[np.float64],
        x0: _ToFloat64Vec,
        jac: _ToJacFun[np.float64],
        hess: _ToHessFun[np.float64],
        finite_diff_rel_step: onp.ToFloat | onp.ToFloat1D | None = None,
        finite_diff_jac_sparsity: _ToJac | None = None,
        finite_diff_bounds: _FDBounds = ...,
        sparse_jacobian: bool | None = None,
        workers: _Workers | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        fun: _VectorFun[_XT_contra],
        x0: _Vec[_XT_contra] | onp.ToFloat1D,
        jac: _ToJacFun[_XT_contra],
        hess: _ToHessFun[_XT_contra],
        finite_diff_rel_step: onp.ToFloat | onp.ToFloat1D | None = None,
        finite_diff_jac_sparsity: _ToJac | None = None,
        finite_diff_bounds: _FDBounds = ...,
        sparse_jacobian: bool | None = None,
        workers: int | _DoesMap | None = None,
    ) -> None: ...

    #
    @property
    def nfev(self, /) -> int: ...
    @property
    def njev(self, /) -> int: ...
    @property
    def nhev(self, /) -> int: ...

    #
    def _update_v(self, /, v: onp.ToFloat1D) -> None: ...
    def _update_x(self, /, x: onp.ToFloat1D) -> None: ...
    def _update_fun(self, /) -> None: ...
    def _update_jac(self, /) -> None: ...
    def _update_hess(self, /) -> None: ...

    #
    def fun(self, /, x: onp.ToFloat1D) -> _Vec[npc.floating]: ...
    def jac(self, /, x: onp.ToFloat1D) -> _Jac[npc.floating]: ...
    def hess(self, /, x: onp.ToFloat1D, v: onp.ToFloat1D) -> _Hess[npc.floating]: ...

class LinearVectorFunction(Generic[_XT_contra]):
    xp: Final[ModuleType]

    n: Final[int]
    x_dtype: _XT_contra  # readonly
    x: _Vec[_XT_contra]

    m: Final[int]
    f: _Vec[npc.floating]
    f_updated: bool

    sparse_jacobian: Final[bool]
    J: Final[_Jac[npc.floating]]

    H: Final[csr_array]
    v: _Vec[np.float64]

    @overload
    def __init__(
        self: LinearVectorFunction[np.float64],
        /,
        A: onp.ToFloat2D | spmatrix | sparray,
        x0: _ToFloat64Vec,
        sparse_jacobian: bool | None,
    ) -> None: ...
    @overload
    def __init__(
        self, /, A: onp.ToFloat2D | spmatrix | sparray, x0: _Vec[_XT_contra] | onp.ToFloat1D, sparse_jacobian: bool | None
    ) -> None: ...

    #
    def _update_x(self, /, x: onp.ToFloat1D) -> None: ...

    #
    def fun(self, /, x: onp.ToFloat1D) -> _Vec[npc.floating]: ...
    def jac(self, /, x: onp.ToFloat1D) -> _Jac[npc.floating]: ...
    def hess(self, /, x: onp.ToFloat1D, v: _Vec[np.float64]) -> csr_array: ...

class IdentityVectorFunction(LinearVectorFunction[_XT_contra]):
    @overload
    def __init__(self: IdentityVectorFunction[np.float64], /, x0: _ToFloat64Vec, sparse_jacobian: bool | None) -> None: ...
    @overload
    def __init__(self, /, x0: onp.CanArrayND[_XT_contra] | onp.ToFloat1D, sparse_jacobian: bool | None) -> None: ...
