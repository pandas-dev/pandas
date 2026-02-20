from typing import Any, Generic, Self, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

###

_ScalarT = TypeVar("_ScalarT", bound=np.float64 | np.complex128)
_UncensoredT_co = TypeVar("_UncensoredT_co", bound=np.float64 | np.complex128, default=np.float64 | Any, covariant=True)
_LeftT = TypeVar("_LeftT", bound=np.float64 | np.complex128, default=np.float64)
_LeftT_co = TypeVar("_LeftT_co", bound=np.float64 | np.complex128, default=np.float64 | Any, covariant=True)
_RightT = TypeVar("_RightT", bound=np.float64 | np.complex128, default=np.float64)
_RightT_co = TypeVar("_RightT_co", bound=np.float64 | np.complex128, default=np.float64 | Any, covariant=True)
_IntervalT = TypeVar("_IntervalT", bound=np.float64 | np.complex128, default=np.float64)
_IntervalT_co = TypeVar("_IntervalT_co", bound=np.float64 | np.complex128, default=np.float64 | Any, covariant=True)

_LeftCensored: TypeAlias = CensoredData[_ScalarT, _ScalarT, np.float64, np.float64]
_RightCensored: TypeAlias = CensoredData[_ScalarT, np.float64, _ScalarT, np.float64]
_IntervalCensored: TypeAlias = CensoredData[_ScalarT, np.float64, np.float64, _ScalarT]

###

# NOTE: We ignore the `[c]longdouble` dtypes for now, as they are not commonly used and would complicate the overloads a lot.

class CensoredData(Generic[_UncensoredT_co, _LeftT_co, _RightT_co, _IntervalT_co]):
    _uncensored: onp.Array1D[_UncensoredT_co]
    _left: onp.Array1D[_LeftT_co]
    _right: onp.Array1D[_RightT_co]
    _interval: onp.Array2D[_IntervalT_co]

    #
    @overload
    @classmethod
    def left_censored(cls, /, x: onp.ToFloat1D, censored: onp.ToBool1D) -> _LeftCensored[np.float64]: ...
    @overload
    @classmethod
    def left_censored(cls, /, x: onp.ToJustComplex1D, censored: onp.ToBool1D) -> _LeftCensored[np.complex128]: ...

    #
    @overload
    @classmethod
    def right_censored(cls, /, x: onp.ToFloat1D, censored: onp.ToBool1D) -> _RightCensored[np.float64]: ...
    @overload
    @classmethod
    def right_censored(cls, /, x: onp.ToJustComplex1D, censored: onp.ToBool1D) -> _RightCensored[np.complex128]: ...

    #
    @overload
    @classmethod
    def interval_censored(cls, /, low: onp.ToFloat1D, high: onp.ToFloat1D) -> _IntervalCensored[np.float64]: ...
    @overload
    @classmethod
    def interval_censored(cls, /, low: onp.ToJustComplex1D, high: onp.ToComplex1D) -> _IntervalCensored[np.complex128]: ...
    @overload
    @classmethod
    def interval_censored(cls, /, low: onp.ToComplex1D, high: onp.ToJustComplex1D) -> _IntervalCensored[np.complex128]: ...

    #
    @overload
    def __init__(
        self: CensoredData[_ScalarT, _LeftT, _RightT, _IntervalT],
        /,
        uncensored: onp.ToArray1D[float, _ScalarT] | None = None,
        *,
        left: onp.ToArray1D[float, _LeftT] | None = None,
        right: onp.ToArray1D[float, _RightT] | None = None,
        interval: onp.ToArray2D[float, _IntervalT] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CensoredData[np.complex128, _LeftT, _RightT, _IntervalT],
        /,
        uncensored: onp.ToJustComplex1D,
        *,
        left: onp.ToArray1D[float, _LeftT] | None = None,
        right: onp.ToArray1D[float, _RightT] | None = None,
        interval: onp.ToArray2D[float, _IntervalT] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CensoredData[np.complex128, np.complex128, _RightT, _IntervalT],
        /,
        uncensored: onp.ToJustComplex1D,
        *,
        left: onp.ToJustComplex1D,
        right: onp.ToArray1D[float, _RightT] | None = None,
        interval: onp.ToArray2D[float, _IntervalT] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CensoredData[np.complex128, _LeftT, np.complex128, _IntervalT],
        /,
        uncensored: onp.ToJustComplex1D,
        *,
        left: onp.ToArray1D[float, _LeftT] | None = None,
        right: onp.ToJustComplex1D,
        interval: onp.ToArray2D[float, _IntervalT] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CensoredData[np.complex128, _LeftT, _RightT, np.complex128],
        /,
        uncensored: onp.ToJustComplex1D,
        *,
        left: onp.ToArray1D[float, _LeftT] | None = None,
        right: onp.ToArray1D[float, _RightT] | None = None,
        interval: onp.ToJustComplex128_2D,
    ) -> None: ...

    #
    @overload
    def __sub__(self, other: onp.ToFloat, /) -> Self: ...
    @overload
    def __sub__(
        self, other: onp.ToJustComplex, /
    ) -> CensoredData[np.complex128, np.complex128, np.complex128, np.complex128]: ...

    #
    @overload
    def __truediv__(self, other: onp.ToFloat, /) -> Self: ...
    @overload
    def __truediv__(
        self, other: onp.ToJustComplex, /
    ) -> CensoredData[np.complex128, np.complex128, np.complex128, np.complex128]: ...

    #
    def __len__(self, /) -> int: ...
    def num_censored(self, /) -> int: ...

    #
    @overload
    def _uncensor(self: CensoredData[np.float64, np.float64, np.float64, np.float64], /) -> onp.Array1D[np.float64]: ...
    @overload
    def _uncensor(self: CensoredData[np.complex128, Any, Any, Any], /) -> onp.Array1D[np.complex128]: ...
    @overload
    def _uncensor(self: CensoredData[Any, np.complex128, Any, Any], /) -> onp.Array1D[np.complex128]: ...
    @overload
    def _uncensor(self: CensoredData[Any, Any, np.complex128, Any], /) -> onp.Array1D[np.complex128]: ...
    @overload
    def _uncensor(self: CensoredData[Any, Any, Any, np.complex128], /) -> onp.Array1D[np.complex128]: ...

    #
    def _supported(self, /, a: onp.ToComplex, b: onp.ToComplex) -> Self: ...
