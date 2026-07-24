from typing import Any, Literal as L, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from . import _resampling
from ._stats_mstats_common import SiegelslopesResult, TheilslopesResult
from ._stats_py import SignificanceResult
from ._typing import NanPolicy

__all__ = ["chatterjeexi", "siegelslopes", "spearmanrho", "theilslopes"]

###

type _PermutationMethod = L["asymptotic"] | _resampling.PermutationMethod
type _Alternative = L["two-sided", "less", "greater"]

type _SlopesMethod = L["hierarchical", "separate"]

###

@overload
def chatterjeexi(
    x: onp.ToComplexStrict1D,
    y: onp.ToComplexStrict1D,
    *,
    axis: SupportsIndex = 0,
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
    axis: SupportsIndex = 0,
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
    axis: SupportsIndex = 0,
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
    axis: SupportsIndex = 0,
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
    axis: SupportsIndex = 0,
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

#
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: None = None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict1D,
    x: onp.ToFloatStrict1D | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict2D,
    x: onp.ToFloatStrict2D | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: int,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.Array1D[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict3D,
    x: onp.ToFloatStrict3D | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: int,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.Array2D[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.ArrayND[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _SlopesMethod = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64 | Any]: ...

#
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: None = None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict1D,
    x: onp.ToFloatStrict1D | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: int | None = None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict2D,
    x: onp.ToFloatStrict2D | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: int,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.Array1D[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict3D,
    x: onp.ToFloatStrict3D | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: int,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.Array2D[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: int | None = None,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.ArrayND[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: float | npc.floating = 0.95,
    method: _SlopesMethod = "separate",
    *,
    axis: int | None = None,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64 | Any]: ...
