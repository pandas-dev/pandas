__all__ = ["chatterjeexi", "spearmanrho"]

from typing import Literal as L, TypeAlias, overload

import numpy as np
import optype as op
import optype.numpy as onp

from . import _resampling
from ._stats_py import SignificanceResult
from ._typing import NanPolicy

_PermutationMethod: TypeAlias = L["asymptotic"] | _resampling.PermutationMethod
_Alternative: TypeAlias = L["two-sided", "less", "greater"]

###

@overload
def chatterjeexi(
    x: onp.ToComplexStrict1D,
    y: onp.ToComplexStrict1D,
    *,
    axis: op.CanIndex = 0,
    y_continuous: bool = False,
    method: _PermutationMethod = "asymptotic",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload
def chatterjeexi(
    x: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    y: onp.ToComplexStrict2D,
    *,
    axis: op.CanIndex = 0,
    y_continuous: bool = False,
    method: _PermutationMethod = "asymptotic",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def chatterjeexi(
    x: onp.ToComplexStrict2D,
    y: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    *,
    axis: op.CanIndex = 0,
    y_continuous: bool = False,
    method: _PermutationMethod = "asymptotic",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def chatterjeexi(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    *,
    axis: op.CanIndex = 0,
    y_continuous: bool = False,
    method: _PermutationMethod = "asymptotic",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64 | onp.ArrayND[np.float64]]: ...
@overload
def chatterjeexi(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    *,
    axis: op.CanIndex = 0,
    y_continuous: bool = False,
    method: _PermutationMethod = "asymptotic",
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...

# keep in sync with `chatterjeexi` above
@overload
def spearmanrho(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    /,
    *,
    alternative: _Alternative = "two-sided",
    method: _resampling.ResamplingMethod | None = None,
    axis: L[0] = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload
def spearmanrho(
    x: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    y: onp.ToComplexStrict2D,
    /,
    *,
    alternative: _Alternative = "two-sided",
    method: _resampling.ResamplingMethod | None = None,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def spearmanrho(
    x: onp.ToComplexStrict2D,
    y: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    /,
    *,
    alternative: _Alternative = "two-sided",
    method: _resampling.ResamplingMethod | None = None,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def spearmanrho(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    /,
    *,
    alternative: _Alternative = "two-sided",
    method: _resampling.ResamplingMethod | None = None,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64 | onp.ArrayND[np.float64]]: ...
@overload
def spearmanrho(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    /,
    *,
    alternative: _Alternative = "two-sided",
    method: _resampling.ResamplingMethod | None = None,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
