from typing import Literal, TypeAlias, overload
from typing_extensions import TypeAliasType, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._typing import NanPolicy

__all__ = ["differential_entropy", "entropy"]

_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)
_PyScalarT = TypeVar("_PyScalarT")

_AsF32: TypeAlias = np.float16 | np.float32
_AsF64: TypeAlias = np.float64 | npc.integer | np.bool_

_DifferentialMethod: TypeAlias = Literal["vasicek", "van es", "ebrahimi", "correa", "auto"]

_ToArrayMaxND = TypeAliasType(
    "_ToArrayMaxND", onp.ToArrayND[_PyScalarT, _ScalarT] | _PyScalarT | _ScalarT, type_params=(_ScalarT, _PyScalarT)
)
_ToArrayMax1D = TypeAliasType(
    "_ToArrayMax1D", onp.ToArrayStrict1D[_PyScalarT, _ScalarT] | _PyScalarT | _ScalarT, type_params=(_ScalarT, _PyScalarT)
)
_ToArrayMax2D = TypeAliasType(
    "_ToArrayMax2D",
    onp.ToArrayStrict2D[_PyScalarT, _ScalarT] | _ToArrayMax1D[_ScalarT, _PyScalarT],
    type_params=(_ScalarT, _PyScalarT),
)
_ToArrayMax3D = TypeAliasType(
    "_ToArrayMax3D",
    onp.ToArrayStrict3D[_PyScalarT, _ScalarT] | _ToArrayMax2D[_ScalarT, _PyScalarT],
    type_params=(_ScalarT, _PyScalarT),
)

###

# NOTE: The (many) [overload-overlap] mypy errors in `entropy` are false positives, so we instead rely on pyright for this.
# mypy: disable-error-code=overload-overlap

@overload  # nd float64 | int, axis=None (positional) -> 0d float64
def entropy(
    pk: _ToArrayMaxND[_AsF64, float],
    qk: onp.ToFloat64_ND | onp.ToFloat64 | None,
    base: float | None,
    axis: None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def entropy(
    pk: onp.ToFloat64_ND | onp.ToFloat64,
    qk: _ToArrayMaxND[_AsF64, float],
    base: float | None,
    axis: None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # nd float64 | int, axis=None (keyword) -> 0d float64
def entropy(
    pk: _ToArrayMaxND[_AsF64, float],
    qk: onp.ToFloat64_ND | onp.ToFloat64 | None = None,
    base: float | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def entropy(
    pk: onp.ToFloat64_ND | onp.ToFloat64,
    qk: _ToArrayMaxND[_AsF64, float],
    base: float | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # 0d or 1d float64 | int -> 0d float64
def entropy(
    pk: _ToArrayMax1D[_AsF64, float],
    qk: onp.ToFloat64Strict1D | onp.ToFloat64 | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def entropy(
    pk: onp.ToFloat64 | onp.ToFloat64Strict1D,
    qk: _ToArrayMax1D[_AsF64, float],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # 2d float64 | int -> 1d float64
def entropy(
    pk: onp.ToArrayStrict2D[float, _AsF64],
    qk: onp.ToFloat64Strict2D | onp.ToFloat64Strict1D | onp.ToFloat64 | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64Strict2D | onp.ToFloat64Strict1D | onp.ToFloat64,
    qk: onp.ToArrayStrict2D[float, _AsF64],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload
def entropy(
    pk: _ToArrayMax2D[_AsF64, float],
    qk: onp.ToFloat64Strict2D,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64Strict2D,
    qk: _ToArrayMax2D[_AsF64, float],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 3d float64 | int -> 2d float64
def entropy(
    pk: onp.ToArrayStrict3D[float, _AsF64],
    qk: onp.ToFloat64Strict3D | onp.ToFloat64Strict2D | onp.ToFloat64Strict1D | onp.ToFloat64 | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64Strict3D | onp.ToFloat64Strict2D | onp.ToFloat64Strict1D | onp.ToFloat64,
    qk: onp.ToArrayStrict3D[float, _AsF64],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload
def entropy(
    pk: _ToArrayMax3D[_AsF64, float],
    qk: onp.ToFloat64Strict3D,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64Strict3D,
    qk: _ToArrayMax3D[_AsF64, float],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd float64 | int, keepdims=True -> Nd float64
def entropy(
    pk: _ToArrayMaxND[_AsF64, float],
    qk: onp.ToFloat64_ND | onp.ToFloat64 | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64_ND | onp.ToFloat64,
    qk: _ToArrayMaxND[_AsF64, float],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float64 | int -> 0d float64 | Nd float64
def entropy(
    pk: _ToArrayMaxND[_AsF64, float],
    qk: onp.ToFloat64_ND | onp.ToFloat64 | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload
def entropy(
    pk: onp.ToFloat64_ND | onp.ToFloat64,
    qk: _ToArrayMaxND[_AsF64, float],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # Nd float32 | float16, axis=None (positional) -> 0d float32
def entropy(
    pk: _ToArrayMaxND[_AsF32, np.float32],
    qk: _ToArrayMaxND[_AsF32, np.float32] | None,
    base: float | None,
    axis: None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float32: ...
@overload  # Nd float32 | float16, axis=None (keyword) -> 0d float32
def entropy(
    pk: _ToArrayMaxND[_AsF32, np.float32],
    qk: _ToArrayMaxND[_AsF32, np.float32] | None = None,
    base: float | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float32: ...
@overload  # 1d float32 | float16 -> 0d float32
def entropy(
    pk: _ToArrayMax1D[_AsF32, np.float32],
    qk: _ToArrayMax1D[_AsF32, np.float32] | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float32: ...
@overload  # 2d float32 | float16 -> 1d float32
def entropy(
    pk: onp.ToArrayStrict2D[np.float32, _AsF32],
    qk: _ToArrayMax2D[_AsF32, np.float32] | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float32]: ...
@overload
def entropy(
    pk: _ToArrayMax2D[_AsF32, np.float32],
    qk: onp.ToArrayStrict2D[np.float32, _AsF32],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float32]: ...
@overload  # 3d float32 | float16 -> 2d float32
def entropy(
    pk: onp.ToArrayStrict3D[np.float32, _AsF32],
    qk: _ToArrayMax3D[_AsF32, np.float32] | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float32]: ...
@overload
def entropy(
    pk: _ToArrayMax3D[_AsF32, np.float32],
    qk: onp.ToArrayStrict3D[np.float32, _AsF32],
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float32]: ...
@overload  # Nd float32 | float16, keepdims=True -> Nd float32
def entropy(
    pk: _ToArrayMaxND[_AsF32, np.float32],
    qk: _ToArrayMaxND[_AsF32, np.float32] | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float32]: ...
@overload  # Nd float32 | float16 -> 0d float32 | Nd float32
def entropy(
    pk: _ToArrayMaxND[_AsF32, np.float32],
    qk: _ToArrayMaxND[_AsF32, np.float32] | None = None,
    base: float | None = None,
    axis: int = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float32 | onp.ArrayND[np.float32]: ...
@overload  # fallback
def entropy(
    pk: onp.ToFloat64 | onp.ToFloat64_ND,
    qk: onp.ToFloat64 | onp.ToFloat64_ND | None = None,
    base: float | None = None,
    axis: int | None = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> npc.floating | onp.ArrayND[npc.floating]: ...

#
@overload  # Nd known inexact dtype, axis=None
def differential_entropy(
    values: _ToArrayMaxND[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: None,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> _InexactT: ...
@overload  # 0d or 1d known inexact dtype
def differential_entropy(
    values: _ToArrayMax1D[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> _InexactT: ...
@overload  # 2d known inexact dtype
def differential_entropy(
    values: onp.ToArrayStrict2D[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[_InexactT]: ...
@overload  # 2d known inexact dtype
def differential_entropy(
    values: onp.ToArrayStrict3D[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[_InexactT]: ...
@overload  # Nd known inexact dtype, keepdims=True
def differential_entropy(
    values: onp.ToArrayND[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[_InexactT]: ...
@overload  # Nd known inexact dtype
def differential_entropy(
    values: onp.ToArrayND[_InexactT, _InexactT],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _InexactT | onp.ArrayND[_InexactT]: ...
@overload  # Nd float64 | +integer, axis=None
def differential_entropy(
    values: _ToArrayMaxND[_AsF64, float],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: None,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # 0d or 1d float64 | +integer
def differential_entropy(
    values: _ToArrayMax1D[_AsF64, float],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # 2d float64 | +integer
def differential_entropy(
    values: onp.ToArrayStrict2D[float, _AsF64],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 3d float64 | +integer
def differential_entropy(
    values: onp.ToArrayStrict3D[float, _AsF64],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd float64 | +integer, keepdims=True
def differential_entropy(
    values: _ToArrayMaxND[_AsF64, float],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float64 | +integer
def differential_entropy(
    values: _ToArrayMaxND[_AsF64, float],
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # Nd ~complex, axis=None
def differential_entropy(
    values: onp.ToJustComplex128_ND | onp.ToJustComplex128,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: None,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.complex128: ...
@overload  # 0d or 1d ~complex
def differential_entropy(
    values: onp.ToJustComplex128Strict1D | onp.ToJustComplex128,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.complex128: ...
@overload  # 2d ~complex
def differential_entropy(
    values: onp.ToJustComplex128Strict2D,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # 3d ~complex
def differential_entropy(
    values: onp.ToJustComplex128Strict3D,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex, keepdims=True
def differential_entropy(
    values: onp.ToJustComplex128_ND | onp.ToJustComplex128,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex
def differential_entropy(
    values: onp.ToJustComplex128_ND | onp.ToJustComplex128,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.complex128 | onp.ArrayND[np.complex128]: ...
@overload  # floating fallback, keepdims=True
def differential_entropy(
    values: onp.ToFloatND | onp.ToFloat,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int | None = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[npc.floating]: ...
@overload  # floating fallback
def differential_entropy(
    values: onp.ToFloatND | onp.ToFloat,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int | None = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> npc.floating | onp.ArrayND[npc.floating]: ...
@overload  # inexact fallback, keepdims=True
def differential_entropy(
    values: onp.ToComplexND | onp.ToComplex,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int | None = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[npc.inexact]: ...
@overload  # inexact fallback
def differential_entropy(
    values: onp.ToComplexND | onp.ToComplex,
    *,
    window_length: int | None = None,
    base: float | None = None,
    axis: int | None = 0,
    method: _DifferentialMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> npc.inexact | onp.ArrayND[npc.inexact]: ...
