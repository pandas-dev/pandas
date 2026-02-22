from collections.abc import Callable
from typing import Any, Final, Literal as L, LiteralString, Protocol, TypeAlias, final, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

_Float1D: TypeAlias = onp.Array1D[np.float64]

_NameT = TypeVar("_NameT", bound=str)
_FuncT_co = TypeVar("_FuncT_co", bound=Callable[..., object], default=Callable[..., Any], covariant=True)

@final
@type_check_only
class _TypesType:
    intvar: Final[onp.Array0D[np.int32]]

@final
@type_check_only
class _FortranFunction(Protocol[_NameT, _FuncT_co]):
    __name__: _NameT
    @property
    def __call__(self, /) -> _FuncT_co: ...

@type_check_only
class _Func_fpchec(Protocol):
    def __call__(self, /, x: onp.ToFloat1D, t: onp.ToFloat1D, k: onp.ToJustInt) -> int: ...

@type_check_only
class _Func_splev(Protocol):
    def __call__(
        self,
        /,
        t: onp.ToFloat1D,
        c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
        k: onp.ToJustInt,  # 0 <= k < len(c)
        x: onp.ToFloat1D,
        e: L[0, 1, 2, 3] = 0,
    ) -> tuple[_Float1D, int]: ...

@type_check_only
class _Func_splder(Protocol):
    def __call__(
        self,
        /,
        t: onp.ToFloat1D,
        c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
        k: onp.ToJustInt,  # 0 <= k < len(c)
        x: onp.ToFloat1D,
        nu: onp.ToJustInt = 1,  # 0 <= nu <= k < len(c)
        e: L[0, 1, 2, 3] = 0,
    ) -> tuple[_Float1D, int]: ...

@type_check_only
class _Func_splint(Protocol):
    def __call__(
        self,
        /,
        t: onp.ToFloat1D,
        c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
        k: onp.ToJustInt,  # 0 <= k < len(c)
        a: onp.ToFloat | onp.ToFloatND,
        b: onp.ToFloat | onp.ToFloatND,
    ) -> tuple[float, _Float1D]: ...

@type_check_only
class _Func_sproot(Protocol):
    def __call__(
        self,
        /,
        t: onp.ToFloat1D,  # len(t) >= 8
        c: onp.ToFloat1D,  # len(c) >= len(t) - 4
        mest: onp.ToJustInt = ...,  # mest = 3 * (len(t) - 7)
    ) -> tuple[_Float1D, int, int]: ...

@type_check_only
class _Func_spalde(Protocol):
    def __call__(
        self,
        /,
        t: onp.ToFloat1D,  # len(t) >= 8
        c: onp.ToFloat1D,  # len(c) >= len(t) - 4
        k1: onp.ToJustInt,
        x: onp.ToFloat | onp.ToFloatND,
    ) -> tuple[_Float1D, int]: ...

@type_check_only
class _Func_bispe_(Protocol):
    def __call__(
        self,
        /,
        tx: onp.ToFloat1D,
        ty: onp.ToFloat1D,
        c: onp.ToFloat1D,  # len(c) == (len(x) - kx - 1) * (len(y) - ky - 1)
        kx: onp.ToJustInt,  # 0 <= kx < len(x)
        ky: onp.ToJustInt,  # 0 <= ky < len(y)
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
    ) -> tuple[_Float1D, int]: ...

###

__f2py_numpy_version__: Final[LiteralString] = ...
types: Final[_TypesType] = ...

# (x, t, k) -> ier
fpchec: _FortranFunction[L["function fpchec"], _Func_fpchec] = ...

# univariate spline

# (t, c, k, x, e=0) -> (y, ier)
splev: _FortranFunction[L["function splev"], _Func_splev] = ...
# (t, c, k, x, nu=1, e=0) -> (dy, ier)
splder: _FortranFunction[L["function splder"], _Func_splder] = ...
# (t, c, k, a, b) -> (iy, wrk)
splint: _FortranFunction[L["splint"], _Func_splint] = ...
# (t, c, mest=...) -> (zero, m, ier)
sproot: _FortranFunction[L["function sproot"], _Func_sproot] = ...
# (t, c, k1, x) -> (d, ier)
spalde: _FortranFunction[L["function spalde"], _Func_spalde] = ...

_Func_cur: TypeAlias = Callable[..., tuple[int, _Float1D, float, int]]

# (iopt, x, y, w, t, wrk, iwrk, xb?, xe?, k?, s?) -> (n, c, fp, ier)
curfit: _FortranFunction[L["function curfit"], _Func_cur] = ...
# (iopt, x, y, w, t, wrk, iwrk, k?, s?) -> (n, c, fp, ier)
percur: _FortranFunction[L["function percur"], _Func_cur] = ...
# (iopt, ipar, idim, u, x, w, ub, ue, t, wrk, iwrk, k?, s?) -> (n, c, fp, ier)
parcur: _FortranFunction[L["function parcur"], _Func_cur] = ...

# these have ridculously large signatures...
fpcurf0: _FortranFunction[L["function fpcurf0"]] = ...
fpcurf1: _FortranFunction[L["function fpcurf1"]] = ...
fpcurfm1: _FortranFunction[L["function fpcurfm1"]] = ...

# biariate spline

# (tx, ty, c, kx, ky, x, y) -> (z, ier)
bispev: _FortranFunction[L["function bispev"], _Func_bispe_] = ...
# (tx, ty, c, kx, ky, x, y) -> (z, ier)
bispeu: _FortranFunction[L["function bispeu"], _Func_bispe_] = ...

# TODO(jorenham)
parder: _FortranFunction[L["function parder"]] = ...
pardtc: _FortranFunction[L["function pardtc"]] = ...
pardeu: _FortranFunction[L["function pardeu"]] = ...
surfit_smth: _FortranFunction[L["function surfit_smth"]] = ...
surfit_lsq: _FortranFunction[L["function surfit_lsq"]] = ...
spherfit_smth: _FortranFunction[L["function spherfit_smth"]] = ...
spherfit_lsq: _FortranFunction[L["function spherfit_lsq"]] = ...
regrid_smth: _FortranFunction[L["function regrid_smth"]] = ...
regrid_smth_spher: _FortranFunction[L["function regrid_smth_spher"]] = ...
dblint: _FortranFunction[L["dblint"]] = ...
