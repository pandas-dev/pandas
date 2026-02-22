from typing import Any, Generic, Literal as L, NamedTuple, Never, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._resampling import PermutationMethod
from ._typing import Alternative, NanPolicy

###

# the `Any` shapes in the bounds are workarounds for a pyrefly bug
_StatisticT_co = TypeVar("_StatisticT_co", bound=npc.floating | onp.ArrayND[npc.floating, Any], default=Any, covariant=True)
_PValueT_co = TypeVar("_PValueT_co", bound=np.float64 | onp.ArrayND[np.float64, Any], default=Any, covariant=True)
_FloatT = TypeVar("_FloatT", bound=npc.floating)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[Any, ...])

_MannwhitneyuResult0D: TypeAlias = MannwhitneyuResult[_FloatT, np.float64]
_MannwhitneyuResultND: TypeAlias = MannwhitneyuResult[onp.ArrayND[_FloatT, _ShapeT], onp.ArrayND[np.float64, _ShapeT]]

_MannWhitneyUMethod: TypeAlias = L["auto", "asymptotic", "exact"] | PermutationMethod

_JustAnyShape: TypeAlias = tuple[Never, Never, Never, Never]  # workaround for https://github.com/microsoft/pyright/issues/10232

###

class MannwhitneyuResult(NamedTuple, Generic[_StatisticT_co, _PValueT_co]):
    statistic: _StatisticT_co
    pvalue: _PValueT_co

#
@overload  # ?d ~f64
def mannwhitneyu(
    x: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    y: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[np.float64] | _MannwhitneyuResultND[np.float64]: ...
@overload  # ?d ~T
def mannwhitneyu(
    x: onp.ArrayND[_FloatT, _JustAnyShape],
    y: onp.ArrayND[_FloatT, _JustAnyShape],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[_FloatT] | _MannwhitneyuResultND[_FloatT]: ...
@overload  # 1d ~f64
def mannwhitneyu(
    x: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[np.float64]: ...
@overload  # 1d ~T
def mannwhitneyu(
    x: onp.ToArrayStrict1D[_FloatT, _FloatT],
    y: onp.ToArrayStrict1D[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[_FloatT]: ...
@overload  # 2d ~f64
def mannwhitneyu(
    x: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResultND[np.float64, tuple[int]]: ...
@overload  # 2d ~T
def mannwhitneyu(
    x: onp.ToArrayStrict2D[_FloatT, _FloatT],
    y: onp.ToArrayStrict2D[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResultND[_FloatT, tuple[int]]: ...
@overload  # 3d ~f64
def mannwhitneyu(
    x: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResultND[np.float64, tuple[int, int]]: ...
@overload  # 3d ~T
def mannwhitneyu(
    x: onp.ToArrayStrict3D[_FloatT, _FloatT],
    y: onp.ToArrayStrict3D[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResultND[_FloatT, tuple[int, int]]: ...
@overload  # nd ~f64
def mannwhitneyu(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[np.float64] | _MannwhitneyuResultND[np.float64]: ...
@overload  # nd ~f64, axis=None  (keyword)
def mannwhitneyu(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    method: _MannWhitneyUMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[np.float64]: ...
@overload  # nd ~f64, keepdims=True
def mannwhitneyu(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    y: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> _MannwhitneyuResultND[np.float64]: ...
@overload  # nd ~T
def mannwhitneyu(
    x: onp.ToArrayND[_FloatT, _FloatT],
    y: onp.ToArrayND[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[_FloatT] | _MannwhitneyuResultND[_FloatT]: ...
@overload  # nd ~T, axis=None  (keyword)
def mannwhitneyu(
    x: onp.ToArrayND[_FloatT, _FloatT],
    y: onp.ToArrayND[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    method: _MannWhitneyUMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[_FloatT]: ...
@overload  # nd ~T, keepdims=True
def mannwhitneyu(
    x: onp.ToArrayND[_FloatT, _FloatT],
    y: onp.ToArrayND[_FloatT, _FloatT],
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> _MannwhitneyuResultND[_FloatT]: ...
@overload  # nd +floating
def mannwhitneyu(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[npc.floating] | _MannwhitneyuResultND[npc.floating]: ...
@overload  # nd +floating, axis=None  (keyword)
def mannwhitneyu(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    method: _MannWhitneyUMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _MannwhitneyuResult0D[npc.floating]: ...
@overload  # nd +floating, keepdims=True
def mannwhitneyu(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    use_continuity: bool = True,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    method: _MannWhitneyUMethod = "auto",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> _MannwhitneyuResultND[npc.floating]: ...
