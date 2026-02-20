from typing import Literal, TypeAlias, overload

import numpy as np
import optype.numpy as onp

from ._typing import NanPolicy

_QuantileMethod: TypeAlias = Literal[
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

# NOTE: There is a false positive `overload-overlap` mypy error that only occurs with `numpy<2.2`
# mypy: disable-error-code=overload-overlap

@overload
def quantile(
    x: onp.ToFloatStrict1D,
    p: onp.ToJustFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: Literal[-1, 0] | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool | None = None,
    weights: onp.ToJustFloatStrict1D | None = None,
) -> np.float64: ...
@overload
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloat,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse | None = None,
    weights: onp.ToJustFloatND | None = None,
) -> np.float64: ...
@overload
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
@overload
def quantile(
    x: onp.ToFloatND,
    p: onp.ToJustFloat | onp.ToJustFloatND,
    *,
    method: _QuantileMethod = "linear",
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToTrue,
    weights: onp.ToJustFloatND | None = None,
) -> onp.ArrayND[np.float64]: ...
