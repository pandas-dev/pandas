# pyright: reportIncompatibleMethodOverride=false

from collections.abc import Callable, Iterable
from typing import Any, Final, Literal as L, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "assoc_legendre_p",
    "assoc_legendre_p_all",
    "legendre_p",
    "legendre_p_all",
    "sph_harm_y",
    "sph_harm_y_all",
    "sph_legendre_p",
    "sph_legendre_p_all",
]

###

type _Complex = np.complex64 | np.complex128  # `clongdouble` isn't supported
type _ToJustComplex = op.JustComplex | _Complex
type _ToJustComplexND = onp.CanArrayND[_Complex] | onp.SequenceND[onp.CanArrayND[_Complex]] | onp.SequenceND[_ToJustComplex]
type _ToJustComplex_D = _ToJustComplex | _ToJustComplexND

type _ToInt_D = onp.ToInt | onp.ToIntND
type _ToFloat_D = onp.ToFloat | onp.ToFloatND
type _ToComplex_D = onp.ToComplex | onp.ToComplexND

type _Float1D = onp.Array1D[np.float64]
type _Float3D = onp.Array3D[np.float64]
type _Float1_D = onp.Array[onp.AtLeast1D[Any], np.float64]
type _Float2_D = onp.Array[onp.AtLeast2D[Any], np.float64]
type _Float3_D = onp.Array[onp.AtLeast3D[Any], np.float64]

type _Complex0D = onp.Array0D[np.complex128]
type _Complex1D = onp.Array1D[np.complex128]
type _Complex2D = onp.Array2D[np.complex128]
type _Complex3D = onp.Array3D[np.complex128]
type _Complex4D = onp.ArrayND[np.complex128, tuple[int, int, int, int]]
type _Complex1_D = onp.Array[onp.AtLeast1D[Any], np.complex128]
type _Complex2_D = onp.Array[onp.AtLeast2D[Any], np.complex128]
type _Complex3_D = onp.Array[onp.AtLeast3D[Any], np.complex128]

type _Complex01D = tuple[_Complex0D, _Complex1D]
type _Complex012D = tuple[_Complex0D, _Complex1D, _Complex2D]
type _Complex23D = tuple[_Complex2D, _Complex3D]
type _Complex234D = tuple[_Complex2D, _Complex3D, _Complex4D]
type _Complex12_D = tuple[_Complex1_D, _Complex2_D]
type _Complex123_D = tuple[_Complex1_D, _Complex2_D, _Complex3_D]
type _Complex33_D = tuple[_Complex3_D, _Complex3_D]
type _Complex333_D = tuple[_Complex3_D, _Complex3_D, _Complex3_D]

type _Branch = L[2, 3]
type _Branch_D = _Branch | onp.SequenceND[_Branch] | onp.CanArrayND[npc.integer]

type _D0 = L[False, 0]
type _D1 = L[True, 1]
type _D2 = L[2]
type _Dn = L[0, 1, 2] | bool | np.bool

_UFuncT_co = TypeVar("_UFuncT_co", bound=Callable[..., object], default=Callable[..., Any], covariant=True)

###

class MultiUFunc:  # undocumented
    @property
    @override
    # pyrefly: ignore [bad-override]
    def __doc__(self, /) -> str | None: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleVariableOverride]

    #
    def __init__(
        self,
        /,
        ufunc_or_ufuncs: _UFuncT_co | Iterable[_UFuncT_co],
        name: str | None = None,
        doc: str | None = None,
        *,
        force_complex_output: bool = False,
        **default_kwargs: object,
    ) -> None: ...
    def __call__(self, /, *args: Any, **kwargs: Any) -> Any: ...

@type_check_only
class _LegendreP(MultiUFunc):
    @overload  # 0-d, 0-d
    def __call__(self, /, n: int, z: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float1D: ...
    @overload  # 0-d, >0-d
    def __call__(self, /, n: int, z: onp.ToFloatND, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >0-d, >=0-d
    def __call__(self, /, n: onp.ToIntND, z: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...

@type_check_only
class _LegendrePAll(MultiUFunc):
    @overload  # 0d float
    def __call__(self, /, n: int, z: onp.ToFloat, *, diff_n: _Dn = 0) -> onp.Array2D[np.float64]: ...
    @overload  # 0d complex
    def __call__(self, /, n: int, z: onp.ToJustComplex, *, diff_n: _Dn = 0) -> onp.Array2D[np.complex128]: ...
    @overload  # 0d float/complex
    def __call__(
        self, /, n: int, z: onp.ToComplex, *, diff_n: _Dn = 0
    ) -> onp.Array2D[np.float64] | onp.Array2D[np.complex128]: ...
    @overload  # Nd float
    def __call__(self, /, n: int, z: onp.ToFloatND, *, diff_n: _Dn = 0) -> _Float3_D: ...
    @overload  # Nd complex
    def __call__(self, /, n: int, z: onp.ToJustComplexND, *, diff_n: _Dn = 0) -> _Complex3_D: ...
    @overload  # Nd float/complex
    def __call__(self, /, n: int, z: onp.ToComplexND, *, diff_n: _Dn = 0) -> _Float3_D | _Complex3_D: ...

@type_check_only
class _AssocLegendreP(MultiUFunc):
    @overload  # float
    def __call__(
        self, /, n: _ToInt_D, m: _ToInt_D, z: _ToFloat_D, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Float1_D: ...
    @overload  # complex
    def __call__(
        self, /, n: _ToInt_D, m: _ToInt_D, z: _ToJustComplex_D, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Complex1_D: ...
    @overload  # float or complex
    def __call__(
        self, /, n: _ToInt_D, m: _ToInt_D, z: _ToComplex_D, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Float1_D | _Complex1_D: ...

@type_check_only
class _AssocLegendrePAll(MultiUFunc):
    @overload  # z: 0-d float
    def __call__(
        self, /, n: int, m: int, z: onp.ToFloat, *, branch_cut: _Branch = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Float3D: ...
    @overload  # z: >=0-d float
    def __call__(
        self, /, n: int, m: int, z: _ToFloat_D, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Float3_D: ...
    @overload  # z: 0-d complex
    def __call__(
        self, /, n: int, m: int, z: _ToJustComplex, *, branch_cut: _Branch = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Complex3D: ...
    @overload  # z: >=0-d complex
    def __call__(
        self, /, n: int, m: int, z: _ToJustComplexND, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Complex3_D: ...
    @overload  # z: >=0-d float or complex
    def __call__(
        self, /, n: int, m: int, z: onp.ToComplexND, *, branch_cut: _Branch_D = 2, norm: bool = False, diff_n: _Dn = 0
    ) -> _Float3_D | _Complex3_D: ...

@type_check_only
class _SphLegendreP(MultiUFunc):
    @overload  # 0-d, 0-d, 0-d
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float1D: ...
    @overload  # >=0-d, >=0-d, >0-d
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >=0-d, >0-d, >=0-d
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >0-d, >=0-d, >=0-d
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...

@type_check_only
class _SphLegendrePAll(MultiUFunc):
    @overload  # 0-d, 0-d, 0-d
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float3D: ...
    @overload  # 0-d, 0-d, >=0-d
    def __call__(self, /, n: int, m: int, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float3_D: ...

@type_check_only
class _SphHarmY(MultiUFunc):
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 0
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D0 = 0) -> _Complex0D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 1
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D1) -> _Complex01D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 2
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D2) -> _Complex012D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...

@type_check_only
class _SphHarmYAll(MultiUFunc):
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 0
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D0 = 0) -> _Complex2D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 0
    def __call__(self, /, n: int, m: int, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D0 = 0) -> _Complex3_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 0
    def __call__(self, /, n: int, m: int, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex3_D: ...
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 1
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D1) -> _Complex23D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 1
    def __call__(self, /, n: int, m: int, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D1) -> _Complex33_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 1
    def __call__(self, /, n: int, m: int, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex33_D: ...
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 2
    def __call__(self, /, n: int, m: int, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D2) -> _Complex234D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 2
    def __call__(self, /, n: int, m: int, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D2) -> _Complex333_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 2
    def __call__(self, /, n: int, m: int, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex333_D: ...

###

legendre_p: Final[_LegendreP] = ...
legendre_p_all: Final[_LegendrePAll] = ...

assoc_legendre_p: Final[_AssocLegendreP] = ...
assoc_legendre_p_all: Final[_AssocLegendrePAll] = ...

sph_legendre_p: Final[_SphLegendreP] = ...
sph_legendre_p_all: Final[_SphLegendrePAll] = ...

sph_harm_y: Final[_SphHarmY] = ...
sph_harm_y_all: Final[_SphHarmYAll] = ...
