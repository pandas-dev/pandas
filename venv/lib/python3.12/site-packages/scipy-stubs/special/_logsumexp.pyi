from typing import Any, Never, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

__all__ = ["log_softmax", "logsumexp", "softmax"]

_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_FloatingT = TypeVar("_FloatingT", bound=npc.floating)
_CFloatingT = TypeVar("_CFloatingT", bound=npc.complexfloating)
_InexactOrArrayT = TypeVar("_InexactOrArrayT", bound=npc.inexact | onp.ArrayND[npc.inexact])

###

# Mypy reports four false positive `overload-overlap` only with `numpy<2.1`
# mypy: disable-error-code="overload-overlap"

@overload  # 0d/nd T, axis=None (default), keepdims=False (default)
def logsumexp(
    a: _InexactT | onp.ToArrayND[Never, _InexactT],
    axis: None = None,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> _InexactT: ...
@overload  # 0d/nd +float , axis=None (default), keepdims=False (default)
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: None = None,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> np.float64: ...
@overload  # 0d/nd ~complex, axis=None (default), keepdims=False (default)
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: None = None,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> np.complex128: ...
@overload  # 0d/nd T, keepdims=True
def logsumexp(
    a: _InexactT | onp.ToArrayND[Never, _InexactT],
    axis: AnyShape | None = None,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[_InexactT]: ...
@overload  # 0d/nd +float, keepdims=True
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: AnyShape | None = None,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # 0d/nd ~complex, keepdims=True
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: AnyShape | None = None,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 0d/nd T, axis=<given>
def logsumexp(
    a: _InexactT | onp.ToArrayND[Never, _InexactT],
    axis: AnyShape,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    *,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[_InexactT] | Any: ...
@overload  # 0d/nd +float, axis=<given>
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: AnyShape,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # 0d/nd ~complex, axis=<given>
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: AnyShape,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    keepdims: onp.ToFalse = False,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # floating fallback, return_sign=False
def logsumexp(
    a: onp.ToFloat | onp.ToFloatND,
    axis: AnyShape | None = None,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    keepdims: bool = False,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.float64 | Any] | Any: ...
@overload  # complex fallback, return_sign=False
def logsumexp(
    a: onp.ToComplex | onp.ToComplexND,
    axis: AnyShape | None = None,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    keepdims: bool = False,
    return_sign: onp.ToFalse = False,
) -> onp.ArrayND[np.complex128 | Any] | Any: ...
@overload  # 0d/nd T@floating, axis=None (default), keepdims=False (default), return_sign=True
def logsumexp(
    a: _FloatingT | onp.ToArrayND[Never, _FloatingT],
    axis: None = None,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[_FloatingT, _FloatingT]: ...
@overload  # 0d/nd +float , axis=None (default), keepdims=False (default), return_sign=True
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: None = None,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[np.float64, np.float64]: ...
@overload  # 0d/nd ~complex, axis=None (default), keepdims=False (default), return_sign=True
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: None = None,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[np.float64, np.complex128]: ...
@overload  # 0d/nd T@complexfloating, axis=None (default), keepdims=False (default), return_sign=True
def logsumexp(
    a: _CFloatingT | onp.ToArrayND[Never, _CFloatingT],
    axis: None = None,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[npc.floating, _CFloatingT]: ...
@overload  # 0d/nd T@floatinv, keepdims=True, return_sign=True
def logsumexp(
    a: _FloatingT | onp.ToArrayND[Never, _FloatingT],
    axis: AnyShape | None = None,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[_FloatingT], onp.ArrayND[_FloatingT]]: ...
@overload  # 0d/nd +float, keepdims=True, return_sign=True
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: AnyShape | None = None,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # 0d/nd ~complex, keepdims=True, return_sign=True
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: AnyShape | None = None,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.complex128]]: ...
@overload  # 0d/nd T@complexfloating, keepdims=True, return_sign=True
def logsumexp(
    a: _CFloatingT | onp.ToArrayND[Never, _CFloatingT],
    axis: AnyShape | None = None,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    *,
    keepdims: onp.ToTrue,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[npc.floating], onp.ArrayND[_CFloatingT]]: ...
@overload  # 0d/nd T@floatinv, axis=<given>, return_sign=True
def logsumexp(
    a: _FloatingT | onp.ToArrayND[Never, _FloatingT],
    axis: AnyShape,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[_FloatingT] | Any, onp.ArrayND[_FloatingT] | Any]: ...
@overload  # 0d/nd +float, axis=<given>, return_sign=True
def logsumexp(
    a: onp.ToInt | onp.ToIntND | onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: AnyShape,
    b: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64] | Any, onp.ArrayND[np.float64] | Any]: ...
@overload  # 0d/nd ~complex, axis=<given>, return_sign=True
def logsumexp(
    a: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: AnyShape,
    b: onp.ToComplex128 | onp.ToComplex128_ND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64] | Any, onp.ArrayND[np.complex128] | Any]: ...
@overload  # 0d/nd T@complexfloating, axis=<given>, return_sign=True
def logsumexp(
    a: _CFloatingT | onp.ToArrayND[Never, _CFloatingT],
    axis: AnyShape,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    keepdims: onp.ToFalse = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[npc.floating] | Any, onp.ArrayND[_CFloatingT] | Any]: ...
@overload  # floating fallback, return_sign=True
def logsumexp(
    a: onp.ToFloat | onp.ToFloatND,
    axis: AnyShape | None = None,
    b: onp.ToFloat | onp.ToFloatND | None = None,
    keepdims: bool = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64 | Any] | Any, onp.ArrayND[np.float64 | Any] | Any]: ...
@overload  # complex fallback, return_sign=True
def logsumexp(
    a: onp.ToComplex | onp.ToComplexND,
    axis: AnyShape | None = None,
    b: onp.ToComplex | onp.ToComplexND | None = None,
    keepdims: bool = False,
    *,
    return_sign: onp.ToTrue,
) -> tuple[onp.ArrayND[np.float64 | Any] | Any, onp.ArrayND[np.complex128 | Any] | Any]: ...

# NOTE: keep in sync with `log_softmax`
@overload  # T
def softmax(x: _InexactOrArrayT, axis: AnyShape | None = None) -> _InexactOrArrayT: ...
@overload  # 0d +float64
def softmax(x: onp.ToInt | onp.ToJustFloat64, axis: AnyShape | None = None) -> np.float64: ...
@overload  # 0d ~complex128
def softmax(x: onp.ToJustComplex128, axis: AnyShape | None = None) -> np.complex128: ...
@overload  # nd T@inexact
def softmax(x: onp.ToArrayND[Never, _InexactT], axis: AnyShape | None = None) -> onp.ArrayND[_InexactT]: ...
@overload  # nd +float64
def softmax(x: onp.ToIntND | onp.ToJustFloat64_ND, axis: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~complex128
def softmax(x: onp.ToJustComplex128_ND, axis: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload  # 0d float fallback
def softmax(x: onp.ToFloat, axis: AnyShape | None = None) -> np.float64 | Any: ...
@overload  # 0d complex fallback
def softmax(x: onp.ToComplex, axis: AnyShape | None = None) -> np.complex128 | Any: ...
@overload  # nd float fallback
def softmax(x: onp.ToFloatND, axis: AnyShape | None = None) -> onp.ArrayND[np.float64 | Any]: ...
@overload  # nd complex fallback
def softmax(x: onp.ToComplexND, axis: AnyShape | None = None) -> onp.ArrayND[np.complex128 | Any]: ...

# NOTE: keep in sync with `softmax`
@overload  # T
def log_softmax(x: _InexactOrArrayT, axis: AnyShape | None = None) -> _InexactOrArrayT: ...
@overload  # 0d +float64
def log_softmax(x: onp.ToInt | onp.ToJustFloat64, axis: AnyShape | None = None) -> np.float64: ...
@overload  # 0d ~complex128
def log_softmax(x: onp.ToJustComplex128, axis: AnyShape | None = None) -> np.complex128: ...
@overload  # nd T@inexact
def log_softmax(x: onp.ToArrayND[Never, _InexactT], axis: AnyShape | None = None) -> onp.ArrayND[_InexactT]: ...
@overload  # nd +float64
def log_softmax(x: onp.ToIntND | onp.ToJustFloat64_ND, axis: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~complex128
def log_softmax(x: onp.ToJustComplex128_ND, axis: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload  # 0d float fallback
def log_softmax(x: onp.ToFloat, axis: AnyShape | None = None) -> np.float64 | Any: ...
@overload  # 0d complex fallback
def log_softmax(x: onp.ToComplex, axis: AnyShape | None = None) -> np.complex128 | Any: ...
@overload  # nd float fallback
def log_softmax(x: onp.ToFloatND, axis: AnyShape | None = None) -> onp.ArrayND[np.float64 | Any]: ...
@overload  # nd complex fallback
def log_softmax(x: onp.ToComplexND, axis: AnyShape | None = None) -> onp.ArrayND[np.complex128 | Any]: ...
