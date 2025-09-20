from collections.abc import Callable
from typing import Any, ClassVar, Final, Generic, Literal, Self, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar, TypeVarTuple, Unpack, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["complex_ode", "ode"]

_Ts = TypeVarTuple("_Ts", default=Unpack[tuple[Any, ...]])
_Inexact64T_co = TypeVar("_Inexact64T_co", bound=npc.inexact, default=np.float64 | np.complex128, covariant=True)

_IntegratorReal: TypeAlias = Literal["vode", "dopri5", "dop853", "lsoda"]
_IntegratorComplex: TypeAlias = Literal["vode", "zvode"]

@type_check_only
class _IntegratorParams(TypedDict, total=False):
    with_jacobian: bool
    rtol: float
    atol: float
    lband: float | None
    uband: float | None
    order: int
    nsteps: int
    max_step: float
    min_step: float
    first_step: float
    ixpr: int
    max_hnil: int
    max_order_ns: int
    max_order_s: int
    method: Literal["adams", "bds"] | None
    safety: float
    ifactor: float
    dfactor: float
    beta: float
    verbosity: int

###

class IntegratorConcurrencyError(RuntimeError):
    def __init__(self, /, name: str) -> None: ...

class ode(Generic[_Inexact64T_co, *_Ts]):
    f: Callable[[float, onp.Array1D[_Inexact64T_co], *_Ts], complex | onp.ToComplex1D]
    jac: Callable[[float, onp.Array1D[_Inexact64T_co], *_Ts], complex | onp.ToComplex2D] | None
    f_params: tuple[*_Ts]
    jac_params: tuple[*_Ts]
    stiff: Literal[0, 1]
    t: float

    def __init__(
        self,
        /,
        f: Callable[[float, onp.Array1D[_Inexact64T_co], *_Ts], complex | onp.ToComplex1D],
        jac: Callable[[float, onp.Array1D[_Inexact64T_co], *_Ts], complex | onp.ToComplex2D] | None = None,
    ) -> None: ...

    #
    @property
    def y(self, /) -> onp.Array1D[_Inexact64T_co]: ...

    #
    @overload
    def set_initial_value(
        self: ode[np.float64, *_Ts], /, y: float | onp.ToFloat1D, t: float = 0.0
    ) -> ode[_Inexact64T_co, *_Ts]: ...
    @overload
    def set_initial_value(
        self: ode[np.complex128, *_Ts], /, y: complex | onp.ToComplex1D, t: float = 0.0
    ) -> ode[_Inexact64T_co, *_Ts]: ...

    #
    @overload
    def set_integrator(
        self: ode[np.float64, *_Ts], /, name: _IntegratorReal, **integrator_params: Unpack[_IntegratorParams]
    ) -> ode[_Inexact64T_co, *_Ts]: ...
    @overload
    def set_integrator(
        self: ode[np.complex128, *_Ts], /, name: _IntegratorComplex, **integrator_params: Unpack[_IntegratorParams]
    ) -> ode[_Inexact64T_co, *_Ts]: ...

    #
    def integrate(self, /, t: float, step: bool = False, relax: bool = False) -> onp.Array1D[_Inexact64T_co]: ...
    def successful(self, /) -> bool: ...
    def get_return_code(self, /) -> Literal[-7, -6, -5, -4, -3, -2, -1, 1, 2]: ...
    def set_f_params(self, /, *args: *_Ts) -> Self: ...
    def set_jac_params(self, /, *args: *_Ts) -> Self: ...
    def set_solout(self, /, solout: Callable[[float, onp.Array1D[_Inexact64T_co]], Literal[-1, 0] | None]) -> None: ...

class complex_ode(ode[np.complex128, *_Ts], Generic[*_Ts]):
    cf: Callable[[float, onp.Array1D[np.complex128], *_Ts], complex | onp.ToComplex1D]
    cjac: Callable[[float, onp.Array1D[np.complex128], *_Ts], complex | onp.ToComplex2D] | None
    tmp: onp.Array1D[np.float64]

    @override
    def set_integrator(self, /, name: _IntegratorReal, **integrator_params: Unpack[_IntegratorParams]) -> Self: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @override
    def set_initial_value(self, /, y: complex | onp.ToComplex1D, t: float = 0.0) -> Self: ...

class IntegratorBase(Generic[_Inexact64T_co]):
    runner: ClassVar[Callable[..., tuple[Any, ...]] | None]  # fortran function or unavailable
    supports_run_relax: ClassVar[Literal[0, 1] | None] = None
    supports_step: ClassVar[Literal[0, 1] | None] = None
    supports_solout: ClassVar[bool] = ...
    scalar: ClassVar[type] = ...

    handle: int
    success: Literal[0, 1] | bool | None = None
    integrator_classes: list[type[IntegratorBase]]
    istate: int | None = None

    def acquire_new_handle(self, /) -> None: ...
    def check_handle(self, /) -> None: ...
    def reset(self, /, n: int, has_jac: bool) -> None: ...
    def run(
        self,
        /,
        f: Callable[..., _Inexact64T_co],
        jac: Callable[..., onp.ArrayND[_Inexact64T_co]] | None,
        y0: complex,
        t0: float,
        t1: float,
        f_params: tuple[object, ...],
        jac_params: tuple[object, ...],
    ) -> tuple[_Inexact64T_co, float]: ...
    def step(
        self,
        /,
        f: Callable[..., _Inexact64T_co],
        jac: Callable[..., onp.ArrayND[_Inexact64T_co]],
        y0: complex,
        t0: float,
        t1: float,
        f_params: tuple[object, ...],
        jac_params: tuple[object, ...],
    ) -> tuple[_Inexact64T_co, float]: ...
    def run_relax(
        self,
        /,
        f: Callable[..., _Inexact64T_co],
        jac: Callable[..., onp.ArrayND[_Inexact64T_co]],
        y0: complex,
        t0: float,
        t1: float,
        f_params: tuple[object, ...],
        jac_params: tuple[object, ...],
    ) -> tuple[_Inexact64T_co, float]: ...

class vode(IntegratorBase[_Inexact64T_co], Generic[_Inexact64T_co]):
    messages: ClassVar[dict[int, str]] = ...

    active_global_handle: int
    meth: int
    with_jacobian: bool
    rtol: float
    atol: float
    mu: float
    ml: float
    order: int
    nsteps: int
    max_step: float
    min_step: float
    first_step: float
    initialized: bool
    rwork: onp.Array1D[np.float64]
    iwork: onp.Array1D[np.int32]
    call_args: list[float | onp.ArrayND[np.float64] | onp.ArrayND[np.int32]]

    def __init__(
        self,
        /,
        method: Literal["adams", "bdf"] = "adams",
        with_jacobian: bool = False,
        rtol: float = 1e-06,
        atol: float = 1e-12,
        lband: float | None = None,
        uband: float | None = None,
        order: int = 12,
        nsteps: int = 500,
        max_step: float = 0.0,
        min_step: float = 0.0,
        first_step: float = 0.0,
    ) -> None: ...

class zvode(vode[np.complex128]):
    active_global_handle: int
    zwork: onp.Array1D[np.complex128]
    call_args: list[float | onp.ArrayND[np.complex128] | onp.ArrayND[np.float64] | onp.ArrayND[np.int32]]  # type: ignore[assignment] # pyright: ignore[reportIncompatibleVariableOverride]
    initialized: bool

class dopri5(IntegratorBase[np.float64]):
    name: ClassVar[str] = "dopri5"
    messages: ClassVar[dict[int, str]] = ...

    rtol: Final[float]
    atol: Final[float]
    nsteps: Final[int]
    max_step: Final[float]
    first_step: Final[float]
    safety: Final[float]
    ifactor: Final[float]
    dfactor: Final[float]
    beta: Final[float]
    verbosity: Final[int]
    solout: Callable[[float, onp.Array1D[npc.inexact]], Literal[0, -1]] | None
    solout_cmplx: bool
    iout: int
    work: onp.Array1D[np.float64]
    iwork: onp.Array1D[np.int32]
    call_args: list[float | Callable[..., Literal[0, -1, 1]] | onp.ArrayND[np.float64] | onp.ArrayND[np.int32]]

    def __init__(
        self,
        /,
        rtol: float = 1e-06,
        atol: float = 1e-12,
        nsteps: int = 500,
        max_step: float = 0.0,
        first_step: float = 0.0,
        safety: float = 0.9,
        ifactor: float = 10.0,
        dfactor: float = 0.2,
        beta: float = 0.0,
        method: None = None,  # unused
        verbosity: int = -1,
    ) -> None: ...
    def set_solout(
        self, /, solout: Callable[[float, onp.Array1D[np.float64]], Literal[0, -1]] | None, complex: bool = False
    ) -> None: ...
    def _solout(
        self,
        /,
        nr: int,  # unused
        xold: object,  # unused
        x: float,
        y: onp.Array1D[np.float64],
        nd: int,  # unused
        icomp: int,  # unused
        con: object,  # unused
    ) -> Literal[0, -1, 1]: ...

class dop853(dopri5):
    name: ClassVar[str] = "dop853"
    def __init__(
        self,
        /,
        rtol: float = 1e-06,
        atol: float = 1e-12,
        nsteps: int = 500,
        max_step: float = 0.0,
        first_step: float = 0.0,
        safety: float = 0.9,
        ifactor: float = 6.0,
        dfactor: float = 0.3,
        beta: float = 0.0,
        method: None = None,  # ignored
        verbosity: int = -1,
    ) -> None: ...

class lsoda(IntegratorBase[np.float64]):
    active_global_handle: ClassVar[int] = 0
    messages: ClassVar[dict[int, str]] = ...

    with_jacobian: Final[bool]
    rtol: Final[float]
    atol: Final[float]
    mu: Final[float | None]
    ml: Final[float | None]
    max_order_ns: Final[int]
    max_order_s: Final[int]
    nsteps: Final[int]
    max_step: Final[float]
    min_step: Final[float]
    first_step: Final[float]
    ixpr: Final[int]
    max_hnil: Final[int]
    initialized: Final[bool]
    rwork: onp.Array1D[np.float64]
    iwork: onp.Array1D[np.int32]
    call_args: list[float | onp.Array1D[np.float64] | onp.Array1D[np.int32]]
    def __init__(
        self,
        /,
        with_jacobian: bool = False,
        rtol: float = 1e-06,
        atol: float = 1e-12,
        lband: float | None = None,
        uband: float | None = None,
        nsteps: int = 500,
        max_step: float = 0.0,
        min_step: float = 0.0,
        first_step: float = 0.0,
        ixpr: int = 0,
        max_hnil: int = 0,
        max_order_ns: int = 12,
        max_order_s: int = 5,
        method: None = None,  # ignored
    ) -> None: ...

@overload
def find_integrator(name: Literal["vode"]) -> type[vode]: ...
@overload
def find_integrator(name: Literal["zvode"]) -> type[zvode]: ...
@overload
def find_integrator(name: Literal["dopri5"]) -> type[dopri5]: ...
@overload
def find_integrator(name: Literal["dop853"]) -> type[dop853]: ...
@overload
def find_integrator(name: Literal["lsoda"]) -> type[lsoda]: ...
@overload
def find_integrator(name: str) -> type[IntegratorBase] | None: ...
