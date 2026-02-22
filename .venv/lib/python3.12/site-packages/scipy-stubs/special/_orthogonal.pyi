from collections.abc import Callable, Sequence
from typing import Final, TypeAlias, overload
from typing_extensions import TypeVar, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_PointsWeights: TypeAlias = tuple[onp.Array1D[np.float64], onp.Array1D[np.float64]]
_PointsWeightsMu: TypeAlias = tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], np.float64]

__all__ = [
    "c_roots",
    "cg_roots",
    "chebyc",
    "chebys",
    "chebyt",
    "chebyu",
    "gegenbauer",
    "genlaguerre",
    "h_roots",
    "he_roots",
    "hermite",
    "hermitenorm",
    "j_roots",
    "jacobi",
    "js_roots",
    "l_roots",
    "la_roots",
    "laguerre",
    "legendre",
    "p_roots",
    "ps_roots",
    "roots_chebyc",
    "roots_chebys",
    "roots_chebyt",
    "roots_chebyu",
    "roots_gegenbauer",
    "roots_genlaguerre",
    "roots_hermite",
    "roots_hermitenorm",
    "roots_jacobi",
    "roots_laguerre",
    "roots_legendre",
    "roots_sh_chebyt",
    "roots_sh_chebyu",
    "roots_sh_jacobi",
    "roots_sh_legendre",
    "s_roots",
    "sh_chebyt",
    "sh_chebyu",
    "sh_jacobi",
    "sh_legendre",
    "t_roots",
    "ts_roots",
    "u_roots",
    "us_roots",
]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

###

# we need this to avoid false positives on numpy<2.1
# pyright: reportOverlappingOverload = false

class orthopoly1d(np.poly1d):  # undocumented
    limits: Final[tuple[float, float]]
    weights: Final[onp.Array2D[np.float64]]
    weight_func: Final[Callable[[float], float]]
    normcoef: Final[float]

    def __init__(
        self,
        /,
        roots: onp.ToComplex1D,
        weights: onp.ToFloat1D | None = None,
        hn: float = 1.0,
        kn: float = 1.0,
        wfunc: Callable[[float], float] | None = None,
        limits: tuple[float, float] | None = None,
        monic: bool = False,
        eval_func: np.ufunc | None = None,
    ) -> None: ...

    #
    @override  # type: ignore[override]
    @overload
    def __call__(self, /, v: np.poly1d) -> np.poly1d: ...
    @overload
    def __call__(self, /, v: onp.ToFloat) -> np.float64: ...
    @overload
    def __call__(self, /, v: onp.ToJustComplex) -> np.complex128: ...
    @overload
    def __call__(self, /, v: onp.ArrayND[npc.floating | npc.integer | np.bool_, _ShapeT]) -> onp.ArrayND[np.float64, _ShapeT]: ...
    @overload
    def __call__(self, /, v: onp.ArrayND[npc.complexfloating, _ShapeT]) -> onp.ArrayND[np.complex128, _ShapeT]: ...
    @overload
    def __call__(self, /, v: Sequence[float]) -> onp.Array1D[np.float64]: ...
    @overload
    def __call__(self, /, v: Sequence[Sequence[float]]) -> onp.Array2D[np.float64]: ...
    @overload
    def __call__(self, /, v: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...
    @overload
    def __call__(self, /, v: list[complex]) -> onp.Array1D[np.complex128]: ...
    @overload
    def __call__(self, /, v: list[list[complex]]) -> onp.Array2D[np.complex128]: ...
    @overload
    def __call__(self, /, v: onp.ToJustComplexND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def __call__(self, /, v: onp.ToComplexND) -> onp.ArrayND[np.float64 | np.complex128]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]

#
@overload
def roots_jacobi(n: onp.ToInt, alpha: onp.ToFloat, beta: onp.ToFloat, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_jacobi(n: onp.ToInt, alpha: onp.ToFloat, beta: onp.ToFloat, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_sh_jacobi(n: onp.ToInt, p1: onp.ToFloat, q1: onp.ToFloat, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_sh_jacobi(n: onp.ToInt, p1: onp.ToFloat, q1: onp.ToFloat, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_genlaguerre(n: onp.ToInt, alpha: onp.ToFloat, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_genlaguerre(n: onp.ToInt, alpha: onp.ToFloat, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_laguerre(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_laguerre(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_hermite(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_hermite(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_hermitenorm(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_hermitenorm(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_gegenbauer(n: onp.ToInt, alpha: onp.ToFloat, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_gegenbauer(n: onp.ToInt, alpha: onp.ToFloat, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_chebyt(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_chebyt(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_chebyu(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_chebyu(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_chebyc(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_chebyc(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_chebys(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_chebys(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_sh_chebyt(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_sh_chebyt(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_sh_chebyu(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_sh_chebyu(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_legendre(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_legendre(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
@overload
def roots_sh_legendre(n: onp.ToInt, mu: onp.ToFalse = False) -> _PointsWeights: ...
@overload
def roots_sh_legendre(n: onp.ToInt, mu: onp.ToTrue) -> _PointsWeightsMu: ...

#
def legendre(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def chebyt(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def chebyu(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def chebyc(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def chebys(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def jacobi(n: onp.ToInt, alpha: onp.ToFloat, beta: onp.ToFloat, monic: bool = False) -> orthopoly1d: ...
def laguerre(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def genlaguerre(n: onp.ToInt, alpha: onp.ToFloat, monic: bool = False) -> orthopoly1d: ...
def hermite(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def hermitenorm(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def gegenbauer(n: onp.ToInt, alpha: onp.ToFloat, monic: bool = False) -> orthopoly1d: ...
def sh_legendre(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def sh_chebyt(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def sh_chebyu(n: onp.ToInt, monic: bool = False) -> orthopoly1d: ...
def sh_jacobi(n: onp.ToInt, p: onp.ToFloat, q: onp.ToFloat, monic: bool = False) -> orthopoly1d: ...

#
def _roots_hermite_asy(n: onp.ToInt) -> _PointsWeights: ...  # undocumented

p_roots = roots_legendre
t_roots = roots_chebyt
u_roots = roots_chebyu
c_roots = roots_chebyc
s_roots = roots_chebys
j_roots = roots_jacobi
l_roots = roots_laguerre
la_roots = roots_genlaguerre
h_roots = roots_hermite
he_roots = roots_hermitenorm
cg_roots = roots_gegenbauer
ps_roots = roots_sh_legendre
ts_roots = roots_sh_chebyt
us_roots = roots_sh_chebyu
js_roots = roots_sh_jacobi
