# defined in scipy/stats/_stats.pyx

from collections.abc import Callable
from typing import Final, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import CapsuleType, ReadOnly

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

# matches the `ctypedef fused ordered`
_Ordered: TypeAlias = np.int32 | np.int64 | np.float32 | np.float64

# matches the `ctypedef fused real`
_Real: TypeAlias = np.float32 | np.float64 | npc.floating80

# castable to `_Real`
_AsReal: TypeAlias = npc.floating | npc.integer | np.bool_

# castable to a real distance matrix
_Dist2D: TypeAlias = onp.Array2D[_AsReal]

# the (presumed) type of the `global_corr` parameters
_GlobalCorr: TypeAlias = Literal["mgc", "mantel", "biased", "rank"]

###

@type_check_only
class _CApiDict(TypedDict):
    _geninvgauss_pdf: ReadOnly[CapsuleType]
    _studentized_range_cdf: ReadOnly[CapsuleType]
    _studentized_range_cdf_asymptotic: ReadOnly[CapsuleType]
    _studentized_range_pdf: ReadOnly[CapsuleType]
    _studentized_range_pdf_asymptotic: ReadOnly[CapsuleType]
    _studentized_range_moment: ReadOnly[CapsuleType]
    _genhyperbolic_pdf: ReadOnly[CapsuleType]
    _genhyperbolic_logpdf: ReadOnly[CapsuleType]

###

__pyx_capi__: Final[_CApiDict] = ...  # undocumented

def von_mises_cdf(k_obj: onp.ToFloatND, x_obj: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...  # undocumented
def _kendall_dis(x: onp.Array1D[np.intp], y: onp.Array1D[np.intp]) -> int: ...  # undocumented
def _toint64(x: onp.ToIntND) -> onp.Array1D[np.int64]: ...  # undocumented
def _weightedrankedtau(
    x: onp.Array1D[_Ordered],
    y: onp.Array1D[_Ordered],
    rank: onp.Array1D[np.intp],
    weigher: Callable[[float], onp.ToFloat],
    additive: bool,
) -> np.float64: ...  # undocumented
def _rank_distance_matrix(distx: onp.Array2D[npc.floating | npc.integer]) -> onp.Array2D[np.intp]: ...  # undocumented

#
@overload
def _center_distance_matrix(
    distx: _Dist2D, global_corr: _GlobalCorr = "mgc", is_ranked: onp.ToTrue = True
) -> tuple[onp.Array2D[np.float64], onp.Array2D[np.intp]]: ...
@overload
def _center_distance_matrix(
    distx: _Dist2D, global_corr: _GlobalCorr, is_ranked: onp.ToFalse
) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.float64]]: ...
@overload
def _center_distance_matrix(
    distx: _Dist2D, global_corr: _GlobalCorr = "mgc", *, is_ranked: onp.ToFalse
) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.float64]]: ...  # undocumented

#
@overload
def _transform_distance_matrix(
    distx: _Dist2D, disty: _Dist2D, global_corr: _GlobalCorr = "mgc", is_ranked: onp.ToTrue = True
) -> dict[str, onp.Array2D[np.float64] | onp.Array2D[np.intp]]: ...
@overload
def _transform_distance_matrix(
    distx: _Dist2D, disty: _Dist2D, global_corr: _GlobalCorr, is_ranked: onp.ToFalse
) -> dict[str, onp.Array2D[np.float64] | onp.Array1D[np.float64]]: ...
@overload
def _transform_distance_matrix(
    distx: _Dist2D, disty: _Dist2D, global_corr: _GlobalCorr = "mgc", *, is_ranked: onp.ToFalse
) -> dict[str, onp.Array2D[np.float64] | onp.Array1D[np.float64]]: ...  # undocumented

#
def _local_covariance(
    distx: _Dist2D, disty: _Dist2D, rank_distx: onp.ArrayND[_AsReal], rank_disty: onp.ArrayND[_AsReal]
) -> onp.Array2D[np.float64]: ...  # undocumented
def _local_correlations(
    distx: _Dist2D, disty: _Dist2D, global_corr: _GlobalCorr = "mgc"
) -> onp.Array2D[np.float64]: ...  # undocumented
def geninvgauss_logpdf(x: float, p: float, b: float) -> float: ...  # undocumented
def _studentized_range_cdf_logconst(k: float, df: float) -> float: ...  # undocumented
def _studentized_range_pdf_logconst(k: float, df: float) -> float: ...  # undocumented
def genhyperbolic_pdf(x: float, p: float, a: float, b: float) -> float: ...  # undocumented
def genhyperbolic_logpdf(x: float, p: float, a: float, b: float) -> float: ...  # undocumented

# NOTE: There are two false positive `overload-overlap` mypy errors that only occur with `numpy>=2.2`.
# mypy: disable-error-code=overload-overlap

# keep in sync with `gaussian_kernel_estimate_log`
@overload
def gaussian_kernel_estimate(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyFloat32DType,
    _: _Real | float = 0,
) -> onp.Array2D[np.float32]: ...
@overload
def gaussian_kernel_estimate(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyFloat64DType | None,
    _: _Real | float = 0,
) -> onp.Array2D[np.float64]: ...
@overload
def gaussian_kernel_estimate(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyLongDoubleDType,
    _: _Real | float = 0,
) -> onp.Array2D[np.longdouble]: ...  # undocumented

# keep in sync with `gaussian_kernel_estimate`
@overload
def gaussian_kernel_estimate_log(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyFloat32DType,
    _: _Real | float = 0,
) -> onp.Array2D[np.float32]: ...
@overload
def gaussian_kernel_estimate_log(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyFloat64DType | None,
    _: _Real | float = 0,
) -> onp.Array2D[np.float64]: ...
@overload
def gaussian_kernel_estimate_log(
    points: onp.Array2D[_AsReal],
    values: onp.Array2D[_Real],
    xi: onp.Array2D[_AsReal],
    cho_cov: onp.Array2D[_AsReal],
    dtype: onp.AnyLongDoubleDType,
    _: _Real | float = 0,
) -> onp.Array2D[np.longdouble]: ...  # undocumented
