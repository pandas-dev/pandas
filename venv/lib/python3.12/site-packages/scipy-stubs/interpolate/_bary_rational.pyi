from collections.abc import Sequence
from typing import Generic, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["AAA", "FloaterHormannInterpolator"]

_SCT = TypeVar("_SCT", bound=npc.inexact)
_SCT_co = TypeVar("_SCT_co", bound=npc.inexact, default=npc.inexact, covariant=True)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[int, ...], covariant=True)

_ToFloat16: TypeAlias = np.bool_ | npc.integer | np.float16
_ToFloat64: TypeAlias = _ToFloat16 | np.float32 | np.float64
_ToComplex128: TypeAlias = _ToFloat64 | np.complex64 | np.complex128

###

class _BarycentricRational(Generic[_SCT_co, _ShapeT_co]):
    def __init__(self, /, x: onp.ToComplex1D, y: onp.ToComplexND) -> None: ...

    #
    @overload
    def __call__(
        self: _BarycentricRational[_SCT, tuple[int]], /, z: onp.ArrayND[_SCT, _ShapeT]
    ) -> onp.ArrayND[_SCT, _ShapeT]: ...
    @overload
    def __call__(self, /, z: onp.ToFloat) -> onp.ArrayND[_SCT_co, _ShapeT_co]: ...
    @overload
    def __call__(self, /, z: onp.ToComplex) -> onp.ArrayND[npc.inexact, _ShapeT_co]: ...
    @overload
    def __call__(self, /, z: onp.ToComplexND) -> onp.ArrayND[npc.inexact]: ...

    #
    def poles(self, /) -> onp.Array1D[npc.complexfloating]: ...
    def roots(self, /) -> onp.Array1D[npc.complexfloating]: ...
    def residues(self, /) -> onp.ArrayND[npc.inexact]: ...

class AAA(_BarycentricRational[_SCT_co, tuple[int]], Generic[_SCT_co]):
    weights: onp.Array1D[_SCT_co]

    @property
    def support_points(self, /) -> onp.Array1D[_SCT_co]: ...
    @property
    def support_values(self, /) -> onp.Array1D[_SCT_co]: ...

    #
    @overload
    def __init__(
        self,
        /,
        x: onp.CanArrayND[_SCT_co] | Sequence[_SCT_co],
        y: onp.CanArrayND[_SCT_co | _ToFloat16] | Sequence[_SCT_co | _ToFloat16],
        *,
        rtol: float | None = None,
        max_terms: int = 100,
        clean_up: bool = True,
        clean_up_tol: float = 1e-13,
    ) -> None: ...
    @overload
    def __init__(
        self: AAA[np.float64],
        /,
        x: Sequence[float],
        y: onp.CanArrayND[_ToFloat64] | Sequence[float | _ToFloat64],
        *,
        rtol: float | None = None,
        max_terms: int = 100,
        clean_up: bool = True,
        clean_up_tol: float = 1e-13,
    ) -> None: ...
    @overload
    def __init__(
        self: AAA[np.float64],
        /,
        x: onp.CanArrayND[_ToFloat64] | Sequence[float | _ToFloat64],
        y: Sequence[float],
        *,
        rtol: float | None = None,
        max_terms: int = 100,
        clean_up: bool = True,
        clean_up_tol: float = 1e-13,
    ) -> None: ...
    @overload
    def __init__(
        self: AAA[np.complex128],
        /,
        x: Sequence[complex],
        y: onp.CanArrayND[_ToFloat64] | Sequence[complex | _ToComplex128],
        *,
        rtol: float | None = None,
        max_terms: int = 100,
        clean_up: bool = True,
        clean_up_tol: float = 1e-13,
    ) -> None: ...
    @overload
    def __init__(
        self: AAA[np.complex128],
        /,
        x: onp.CanArrayND[_ToComplex128] | Sequence[complex | _ToComplex128],
        y: Sequence[complex],
        *,
        rtol: float | None = None,
        max_terms: int = 100,
        clean_up: bool = True,
        clean_up_tol: float = 1e-13,
    ) -> None: ...

    #
    def clean_up(self, /, cleanup_tol: float = 1e-13) -> int: ...

class FloaterHormannInterpolator(_BarycentricRational[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    @overload
    def __init__(
        self,
        /,
        points: onp.CanArrayND[_SCT_co | _ToFloat16] | Sequence[_SCT_co | _ToFloat16 | int],
        values: onp.CanArrayND[_SCT_co, _ShapeT_co],
        *,
        d: int = 3,
    ) -> None: ...
    @overload
    def __init__(
        self: FloaterHormannInterpolator[npc.floating, tuple[int]],
        /,
        points: onp.ToFloat1D,
        values: onp.ToFloatStrict1D,
        *,
        d: int = 3,
    ) -> None: ...
    @overload
    def __init__(
        self: FloaterHormannInterpolator[npc.floating, tuple[int, int]],
        /,
        points: onp.ToFloat1D,
        values: onp.ToFloatStrict2D,
        *,
        d: int = 3,
    ) -> None: ...
    @overload
    def __init__(
        self: FloaterHormannInterpolator[npc.inexact, tuple[int]],
        /,
        points: onp.ToComplex1D,
        values: onp.ToComplexStrict1D,
        *,
        d: int = 3,
    ) -> None: ...
    @overload
    def __init__(
        self: FloaterHormannInterpolator[npc.inexact, tuple[int, int]],
        /,
        points: onp.ToComplex1D,
        values: onp.ToComplexStrict2D,
        *,
        d: int = 3,
    ) -> None: ...
    @overload
    def __init__(self, /, points: onp.ToComplex1D, values: onp.ToComplexND, *, d: int = 3) -> None: ...
