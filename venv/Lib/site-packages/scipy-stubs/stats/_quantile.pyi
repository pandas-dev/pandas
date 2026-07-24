from typing import Any, Literal, Never, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._typing import NanPolicy

###

type _QuantileMethod = Literal[
    "inverted_cdf",
    "averaged_inverted_cdf",
    "closest_observation",
    "interpolated_inverted_cdf",
    "hazen",
    "weibull",
    "linear",  # default
    "median_unbiased",
    "normal_unbiased",
]

###

# NOTE: There is a false positive `overload-overlap` mypy error for `quantile` that only occurs with `numpy<2.2`
# mypy: disable-error-code=overload-overlap

# NOTE: And of course, Pyright reports a different overload error (obvously a false positive) for `estimated_cdf` on `numpy<2.1`
# pyright: reportOverlappingOverload=false

# TODO(@jorenham): propagate floating dtype
@overload  # 1d, 0d
def quantile(
    x: onp.ToFloatStrict1D,
    p: onp.ToJustFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
    weights: onp.ToJustFloatStrict1D | None = None,
) -> np.float64: ...
@overload  # 2d, 0d
def quantile(
    x: onp.ToFloatStrict2D,
    p: onp.ToJustFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
    weights: onp.ToJustFloatStrict1D | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # 3d, 0d
def quantile(
    x: onp.ToFloatStrict3D,
    p: onp.ToJustFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
    weights: onp.ToJustFloatStrict1D | None = None,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd, >0d
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
    weights: onp.ToJustFloatND | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # axis=None
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloat | onp.ToJustFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
    weights: onp.ToJustFloatND | None = None,
) -> np.float64: ...
@overload  # keepdims=True
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloat | onp.ToJustFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
    weights: onp.ToJustFloatND | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd, Nd  (fallback)
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloat | onp.ToJustFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
    weights: onp.ToJustFloatND | None = None,
) -> onp.ArrayND[np.float64] | np.float64: ...

#
@overload  # ?d T, 0d float (workaround)
def estimated_cdf[FloatT: npc.floating](
    x: onp.ArrayND[FloatT, tuple[Never, Never, Never, Never]],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> FloatT | onp.ArrayND[FloatT]: ...
@overload  # ?d +f64, 0d float (workaround)
def estimated_cdf(
    x: onp.ArrayND[npc.integer, tuple[Never, Never, Never, Never]],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T, 0d float
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayStrict1D[FloatT, FloatT],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> FloatT: ...
@overload  # 1d +f64, 0d float
def estimated_cdf(
    x: onp.ToArrayStrict1D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> np.float64: ...
@overload  # 1d +f64, 0d float, keepdims=True
def estimated_cdf(
    x: onp.ToArrayStrict1D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.Array1D[np.float64]: ...
@overload  # 2d T, 0d float
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayStrict2D[FloatT, FloatT],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.Array1D[FloatT]: ...
@overload  # 2d +f64, 0d float
def estimated_cdf(
    x: onp.ToArrayStrict2D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +f64, 0d float, keepdims=True
def estimated_cdf(
    x: onp.ToArrayStrict2D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.Array2D[np.float64]: ...
@overload  # 3d T, 0d float
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayStrict3D[FloatT, FloatT],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.Array2D[FloatT]: ...
@overload  # 3d +f64, 0d float
def estimated_cdf(
    x: onp.ToArrayStrict3D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d +f64, 0d float, keepdims=True
def estimated_cdf(
    x: onp.ToArrayStrict3D[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.Array3D[np.float64]: ...
@overload  # Nd T, 0d float
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayND[FloatT, FloatT],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.ArrayND[FloatT] | Any: ...
@overload  # Nd +f64, 0d float
def estimated_cdf(
    x: onp.ToArrayND[float, npc.integer],
    y: float,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] | None = None,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd T, 0d float, keepdims=True
def estimated_cdf[FloatT: npc.floating, ShapeT: tuple[int, ...]](
    x: onp.ArrayND[FloatT, ShapeT],
    y: float | onp.SequenceND[float],
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[FloatT, ShapeT]: ...
@overload  # Nd ^f64, Nd ^f64
def estimated_cdf(
    x: onp.ToArrayND[float, np.float64 | npc.integer],
    y: onp.ToArrayND[float, np.float64 | npc.integer],
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd +f64, Nd ~f64
def estimated_cdf(
    x: onp.ToFloat64_ND,
    y: onp.ArrayND[npc.floating64],
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd ~f64, Nd +f64
def estimated_cdf(
    x: onp.ArrayND[npc.floating64],
    y: onp.ToFloat64_ND,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float, Nd T
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayND[float, FloatT],
    y: onp.ToArrayND[FloatT, FloatT],
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[FloatT]: ...
@overload  # Nd ~f32, Nd +f32
def estimated_cdf[FloatT: npc.floating](
    x: onp.ToArrayND[FloatT, FloatT],
    y: onp.ToArrayND[float, FloatT],
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[FloatT]: ...
@overload  # Nd +floating, Nd +floating  (fallback)
def estimated_cdf(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
) -> onp.ArrayND[npc.floating]: ...
@overload  # Nd +floating, 0d +floating, keepdims=True  (fallback)
def estimated_cdf(
    x: onp.ToFloatND,
    y: onp.ToFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[npc.floating]: ...
