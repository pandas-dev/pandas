from collections.abc import Callable, Iterable
from typing import Any, Final, Generic, Literal as L, ParamSpec, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, override

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

_Tss = ParamSpec("_Tss")
_RT = TypeVar("_RT")
_UFuncT_co = TypeVar("_UFuncT_co", bound=Callable[..., object], default=Callable[..., Any], covariant=True)

_Complex: TypeAlias = np.complex64 | np.complex128  # `clongdouble` isn't supported
_ToJustComplex: TypeAlias = op.JustComplex | _Complex
_ToJustComplexND: TypeAlias = onp.CanArrayND[_Complex] | onp.SequenceND[onp.CanArrayND[_Complex]] | onp.SequenceND[_ToJustComplex]
_ToJustComplex_D: TypeAlias = _ToJustComplex | _ToJustComplexND

_ToInt_D: TypeAlias = onp.ToInt | onp.ToIntND
_ToFloat_D: TypeAlias = onp.ToFloat | onp.ToFloatND
_ToComplex_D: TypeAlias = onp.ToComplex | onp.ToComplexND

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float3D: TypeAlias = onp.Array3D[np.float64]
_Float1_D: TypeAlias = onp.Array[onp.AtLeast1D, np.float64]
_Float2_D: TypeAlias = onp.Array[onp.AtLeast2D, np.float64]
_Float3_D: TypeAlias = onp.Array[onp.AtLeast3D, np.float64]

_Complex0D: TypeAlias = onp.Array0D[np.complex128]
_Complex1D: TypeAlias = onp.Array1D[np.complex128]
_Complex2D: TypeAlias = onp.Array2D[np.complex128]
_Complex3D: TypeAlias = onp.Array3D[np.complex128]
_Complex4D: TypeAlias = onp.ArrayND[np.complex128, tuple[int, int, int, int]]
_Complex1_D: TypeAlias = onp.Array[onp.AtLeast1D, np.complex128]
_Complex2_D: TypeAlias = onp.Array[onp.AtLeast2D, np.complex128]
_Complex3_D: TypeAlias = onp.Array[onp.AtLeast3D, np.complex128]

_Complex01D: TypeAlias = tuple[_Complex0D, _Complex1D]
_Complex012D: TypeAlias = tuple[_Complex0D, _Complex1D, _Complex2D]
_Complex23D: TypeAlias = tuple[_Complex2D, _Complex3D]
_Complex234D: TypeAlias = tuple[_Complex2D, _Complex3D, _Complex4D]
_Complex12_D: TypeAlias = tuple[_Complex1_D, _Complex2_D]
_Complex123_D: TypeAlias = tuple[_Complex1_D, _Complex2_D, _Complex3_D]
_Complex33_D: TypeAlias = tuple[_Complex3_D, _Complex3_D]
_Complex333_D: TypeAlias = tuple[_Complex3_D, _Complex3_D, _Complex3_D]

_Branch: TypeAlias = L[2, 3]
_Branch_D: TypeAlias = _Branch | onp.SequenceND[_Branch] | onp.CanArrayND[npc.integer]

_D0: TypeAlias = L[False, 0]
_D1: TypeAlias = L[True, 1]
_D2: TypeAlias = L[2]
_Dn: TypeAlias = L[0, 1, 2] | bool | np.bool_

###

@type_check_only
class _LegendreP(Protocol):
    @overload  # 0-d, 0-d
    def __call__(self, /, n: onp.ToInt, z: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float1D: ...
    @overload  # 0-d, >0-d
    def __call__(self, /, n: onp.ToInt, z: onp.ToFloatND, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >0-d, >=0-d
    def __call__(self, /, n: onp.ToIntND, z: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...

@type_check_only
class _LegendrePAll(Protocol):
    @overload  # float
    def __call__(self, /, n: onp.ToInt, z: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float3_D: ...
    @overload  # complex
    def __call__(self, /, n: onp.ToInt, z: _ToJustComplex_D, *, diff_n: _Dn = 0) -> _Complex3_D: ...
    @overload  # float or complex
    def __call__(self, /, n: onp.ToInt, z: _ToComplex_D, *, diff_n: _Dn = 0) -> _Float3_D | _Complex3_D: ...

@type_check_only
class _AssocLegendreP(Protocol):
    @overload  # float
    def __call__(
        self, /, n: _ToInt_D, m: _ToInt_D, z: _ToFloat_D, *, branch_cut: _Branch_D = 2, norm: onp.ToBool = False, diff_n: _Dn = 0
    ) -> _Float1_D: ...
    @overload  # complex
    def __call__(
        self,
        /,
        n: _ToInt_D,
        m: _ToInt_D,
        z: _ToJustComplex_D,
        *,
        branch_cut: _Branch_D = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Complex1_D: ...
    @overload  # float or complex
    def __call__(
        self,
        /,
        n: _ToInt_D,
        m: _ToInt_D,
        z: _ToComplex_D,
        *,
        branch_cut: _Branch_D = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Float1_D | _Complex1_D: ...

@type_check_only
class _AssocLegendrePAll(Protocol):
    @overload  # z: 0-d float
    def __call__(
        self, /, n: onp.ToInt, m: onp.ToInt, z: onp.ToFloat, *, branch_cut: _Branch = 2, norm: onp.ToBool = False, diff_n: _Dn = 0
    ) -> _Float3D: ...
    @overload  # z: >=0-d float
    def __call__(
        self,
        /,
        n: onp.ToInt,
        m: onp.ToInt,
        z: _ToFloat_D,
        *,
        branch_cut: _Branch_D = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Float3_D: ...
    @overload  # z: 0-d complex
    def __call__(
        self,
        /,
        n: onp.ToInt,
        m: onp.ToInt,
        z: _ToJustComplex,
        *,
        branch_cut: _Branch = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Complex3D: ...
    @overload  # z: >=0-d complex
    def __call__(
        self,
        /,
        n: onp.ToInt,
        m: onp.ToInt,
        z: _ToJustComplexND,
        *,
        branch_cut: _Branch_D = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Complex3_D: ...
    @overload  # z: >=0-d float or complex
    def __call__(
        self,
        /,
        n: onp.ToInt,
        m: onp.ToInt,
        z: onp.ToComplexND,
        *,
        branch_cut: _Branch_D = 2,
        norm: onp.ToBool = False,
        diff_n: _Dn = 0,
    ) -> _Float3_D | _Complex3_D: ...

@type_check_only
class _SphLegendreP(Protocol):
    @overload  # 0-d, 0-d, 0-d
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float1D: ...
    @overload  # >=0-d, >=0-d, >0-d
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >=0-d, >0-d, >=0-d
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...
    @overload  # >0-d, >=0-d, >=0-d
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float2_D: ...

@type_check_only
class _SphLegendrePAll(Protocol):
    @overload  # 0-d, 0-d, 0-d
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, *, diff_n: _Dn = 0) -> _Float3D: ...
    @overload  # 0-d, 0-d, >=0-d
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: _ToFloat_D, *, diff_n: _Dn = 0) -> _Float3_D: ...

@type_check_only
class _SphHarmY(Protocol):
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 0
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D0 = 0) -> _Complex0D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 0
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D0 = 0) -> _Complex1_D: ...
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 1
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D1) -> _Complex01D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 1
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex12_D: ...
    @overload  # 0-d,     0-d,   0-d,   0-d, diff_n == 2
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D2) -> _Complex012D: ...
    @overload  # >=0-d, >=0-d, >=0-d, > 0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # >=0-d, >=0-d, > 0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: _ToInt_D, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # >=0-d, > 0-d, >=0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: _ToInt_D, m: onp.ToIntND, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...
    @overload  # > 0-d, >=0-d, >=0-d, >=0-d, diff_n == 2
    def __call__(self, /, n: onp.ToIntND, m: _ToInt_D, theta: _ToFloat_D, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex123_D: ...

@type_check_only
class _SphHarmYAll(Protocol):
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 0
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D0 = 0) -> _Complex2D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 0
    def __call__(
        self, /, n: onp.ToInt, m: onp.ToInt, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D0 = 0
    ) -> _Complex3_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 0
    def __call__(
        self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D0 = 0
    ) -> _Complex3_D: ...
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 1
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D1) -> _Complex23D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 1
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D1) -> _Complex33_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 1
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D1) -> _Complex33_D: ...
    @overload  # theta: 0-d, phi: 0-d,    diff_n == 2
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloat, phi: onp.ToFloat, *, diff_n: _D2) -> _Complex234D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 2
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: _ToFloat_D, phi: onp.ToFloatND, *, diff_n: _D2) -> _Complex333_D: ...
    @overload  # theta: >=0-d, phi: >0-d, diff_n == 2
    def __call__(self, /, n: onp.ToInt, m: onp.ToInt, theta: onp.ToFloatND, phi: _ToFloat_D, *, diff_n: _D2) -> _Complex333_D: ...

###

class MultiUFunc(Generic[_UFuncT_co]):
    @property
    @override
    def __doc__(self, /) -> str | None: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleVariableOverride]

    #
    def __init__(
        self,
        /,
        ufunc_or_ufuncs: _UFuncT_co | Iterable[_UFuncT_co],
        doc: str | None = None,
        *,
        force_complex_output: bool = False,
        **default_kwargs: object,
    ) -> None: ...
    def __call__(self: MultiUFunc[Callable[_Tss, _RT]], /, *args: _Tss.args, **kwargs: _Tss.kwargs) -> _RT: ...

###

legendre_p: Final[MultiUFunc[_LegendreP]] = ...
legendre_p_all: Final[MultiUFunc[_LegendrePAll]] = ...

assoc_legendre_p: Final[MultiUFunc[_AssocLegendreP]] = ...
assoc_legendre_p_all: Final[MultiUFunc[_AssocLegendrePAll]] = ...

sph_legendre_p: Final[MultiUFunc[_SphLegendreP]] = ...
sph_legendre_p_all: Final[MultiUFunc[_SphLegendrePAll]] = ...

sph_harm_y: Final[MultiUFunc[_SphHarmY]] = ...
sph_harm_y_all: Final[MultiUFunc[_SphHarmYAll]] = ...
