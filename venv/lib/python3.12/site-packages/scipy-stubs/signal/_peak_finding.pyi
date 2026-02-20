from collections.abc import Callable, Sequence
from typing import Concatenate, Literal, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["argrelextrema", "argrelmax", "argrelmin", "find_peaks", "find_peaks_cwt", "peak_prominences", "peak_widths"]

###

_ScalarT = TypeVar("_ScalarT", bound=np.generic)

_Mode: TypeAlias = Literal["clip", "wrap"]
_PeakProminences: TypeAlias = tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp], onp.ArrayND[np.intp]]
_PeakCondition: TypeAlias = float | onp.ToFloatND | Sequence[float | None]
_FnWavelet: TypeAlias = (
    Callable[Concatenate[int, float, ...], onp.ToComplex1D]
    | Callable[Concatenate[np.intp, np.float64, ...], onp.ToComplex1D]
)  # fmt: skip

###

def argrelmin(data: onp.Array, axis: int = 0, order: int = 1, mode: _Mode = "clip") -> tuple[onp.ArrayND[np.intp], ...]: ...
def argrelmax(data: onp.Array, axis: int = 0, order: int = 1, mode: _Mode = "clip") -> tuple[onp.ArrayND[np.intp], ...]: ...
def argrelextrema(
    data: onp.ArrayND[_ScalarT],
    comparator: Callable[[onp.ArrayND[_ScalarT], onp.ArrayND[_ScalarT]], onp.ToBoolND],
    axis: int = 0,
    order: int = 1,
    mode: _Mode = "clip",
) -> tuple[onp.ArrayND[np.intp], ...]: ...

#
def peak_prominences(x: onp.ToArray1D, peaks: onp.ToIntND, wlen: float | None = None) -> _PeakProminences: ...
def peak_widths(
    x: onp.ToArray1D,
    peaks: onp.ToIntND,
    rel_height: float = 0.5,
    prominence_data: _PeakProminences | None = None,
    wlen: float | None = None,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int_], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...

#
def find_peaks_cwt(
    vector: onp.Array,
    widths: float | onp.ToFloatND,
    wavelet: _FnWavelet | None = None,
    max_distances: onp.ArrayND[npc.floating | npc.integer] | None = None,
    gap_thresh: float | None = None,
    min_length: int | None = None,
    min_snr: float = 1,
    noise_perc: float = 10,
    window_size: int | None = None,
) -> onp.Array1D[np.intp]: ...

# We need these 2^5=32 (combinations of) TypedDicts for each combination of 5 optional find_peaks parameters
# https://github.com/scipy/scipy-stubs/issues/944#issuecomment-3413406314

# 0 (5 choose 0 = 1)

@type_check_only  # {}
class _PeakProperties_0(TypedDict): ...

# 1 (5 choose 1 = 5)

@type_check_only  # {height}
class _PeakProperties_h(TypedDict):
    peak_heights: onp.Array1D[np.float64]

@type_check_only  # {threshold}
class _PeakProperties_t(TypedDict):
    left_thresholds: onp.Array1D[np.float64]
    right_thresholds: onp.Array1D[np.float64]

@type_check_only  # {prominence}
class _PeakProperties_p(TypedDict):
    prominences: onp.Array1D[np.float64]
    left_bases: onp.Array1D[np.intp]
    right_bases: onp.Array1D[np.intp]

@type_check_only  # {width}
class _PeakProperties_w(TypedDict):
    widths: onp.Array1D[np.float64]
    width_heights: onp.Array1D[np.float64]
    left_ips: onp.Array1D[np.float64]
    right_ips: onp.Array1D[np.float64]

@type_check_only  # {plateau_size}
class _PeakProperties_s(TypedDict):
    plateau_sizes: onp.Array1D[np.intp]
    left_edges: onp.Array1D[np.intp]
    right_edges: onp.Array1D[np.intp]

# 2 (5 choose 2 = 10)

@type_check_only  # {height, threshold}
class _PeakProperties_ht(_PeakProperties_h, _PeakProperties_t): ...

@type_check_only  # {height, prominence}
class _PeakProperties_hp(_PeakProperties_h, _PeakProperties_p): ...

@type_check_only  # {height, width}
class _PeakProperties_hw(_PeakProperties_h, _PeakProperties_w): ...

@type_check_only  # {height, plateau_size}
class _PeakProperties_hs(_PeakProperties_h, _PeakProperties_s): ...

@type_check_only  # {threshold, prominence}
class _PeakProperties_tp(_PeakProperties_t, _PeakProperties_p): ...

@type_check_only  # {threshold, width}
class _PeakProperties_tw(_PeakProperties_t, _PeakProperties_w): ...

@type_check_only  # {threshold, plateau_size}
class _PeakProperties_ts(_PeakProperties_t, _PeakProperties_s): ...

@type_check_only  # {prominence, width}
class _PeakProperties_pw(_PeakProperties_p, _PeakProperties_w): ...

@type_check_only  # {prominence, plateau_size}
class _PeakProperties_ps(_PeakProperties_p, _PeakProperties_s): ...

@type_check_only  # {width, plateau_size}
class _PeakProperties_ws(_PeakProperties_w, _PeakProperties_s): ...

# 3 (5 choose 3 = 10)

@type_check_only  # {height, threshold, prominence}
class _PeakProperties_htp(_PeakProperties_ht, _PeakProperties_p): ...

@type_check_only  # {height, threshold, width}
class _PeakProperties_htw(_PeakProperties_ht, _PeakProperties_w): ...

@type_check_only  # {height, threshold, plateau_size}
class _PeakProperties_hts(_PeakProperties_ht, _PeakProperties_s): ...

@type_check_only  # {height, prominence, width}
class _PeakProperties_hpw(_PeakProperties_hp, _PeakProperties_w): ...

@type_check_only  # {height, prominence, plateau_size}
class _PeakProperties_hps(_PeakProperties_hp, _PeakProperties_s): ...

@type_check_only  # {height, width, plateau_size}
class _PeakProperties_hws(_PeakProperties_hw, _PeakProperties_s): ...

@type_check_only  # {threshold, prominence, width}
class _PeakProperties_tpw(_PeakProperties_t, _PeakProperties_pw): ...

@type_check_only  # {threshold, prominence, plateau_size}
class _PeakProperties_tps(_PeakProperties_t, _PeakProperties_ps): ...

@type_check_only  # {threshold, width, plateau_size}
class _PeakProperties_tws(_PeakProperties_t, _PeakProperties_ws): ...

@type_check_only  # {prominence, width, plateau_size}
class _PeakProperties_pws(_PeakProperties_p, _PeakProperties_ws): ...

# 4 (5 choose 4 = 5)

@type_check_only  # {height, threshold, prominence, width}
class _PeakProperties_htpw(_PeakProperties_ht, _PeakProperties_pw): ...

@type_check_only  # {height, threshold, prominence, plateau_size}
class _PeakProperties_htps(_PeakProperties_ht, _PeakProperties_ps): ...

@type_check_only  # {height, threshold, width, plateau_size}
class _PeakProperties_htws(_PeakProperties_ht, _PeakProperties_ws): ...

@type_check_only  # {height, prominence, width, plateau_size}
class _PeakProperties_hpws(_PeakProperties_hp, _PeakProperties_ws): ...

@type_check_only  # {threshold, prominence, width, plateau_size}
class _PeakProperties_tpws(_PeakProperties_t, _PeakProperties_pws): ...

# 5 (5 choose 5 = 1)

@type_check_only  # {height, threshold, prominence, width, plateau_size}
class _PeakProperties_htpws(_PeakProperties_htp, _PeakProperties_ws): ...

# 0
@overload  # {}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_0]: ...

# 1
@overload  # {height}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_h]: ...
@overload  # {threshold}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_t]: ...
@overload  # {prominence}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_p]: ...
@overload  # {width}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    *,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_w]: ...
@overload  # {plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    *,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_s]: ...

# 2
@overload  # {height, threshold}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_ht]: ...
@overload  # {height, prominence}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hp]: ...
@overload  # {height, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    *,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hw]: ...
@overload  # {height, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    *,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hs]: ...
@overload  # {threshold, prominence}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tp]: ...
@overload  # {threshold, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tw]: ...
@overload  # {threshold, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_ts]: ...
@overload  # {prominence, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_pw]: ...
@overload  # {prominence, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_ps]: ...
@overload  # {width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    *,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_ws]: ...

# 3
@overload  # {height, threshold, prominence}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htp]: ...
@overload  # {height, threshold, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htw]: ...
@overload  # {height, threshold, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hts]: ...
@overload  # {height, prominence, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hpw]: ...
@overload  # {height, prominence, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hps]: ...
@overload  # {height, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    prominence: None = None,
    *,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hws]: ...
@overload  # {threshold, prominence, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tpw]: ...
@overload  # {threshold, prominence, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tps]: ...
@overload  # {threshold, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tws]: ...
@overload  # {prominence, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_pws]: ...

# 4
@overload  # {height, threshold, prominence, width}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: None = None,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htpw]: ...
@overload  # {height, threshold, prominence, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: None = None,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htps]: ...
@overload  # {height, threshold, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: None = None,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htws]: ...
@overload  # {height, prominence, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: None = None,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_hpws]: ...
@overload  # {threshold, prominence, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: None = None,
    *,
    threshold: _PeakCondition,
    distance: float | None = None,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_tpws]: ...

# 5
@overload  # {height, threshold, prominence, width, plateau_size}
def find_peaks(
    x: onp.ToFloat1D,
    height: _PeakCondition,
    threshold: _PeakCondition,
    distance: float | None = None,
    *,
    prominence: _PeakCondition,
    width: _PeakCondition,
    wlen: int | None = None,
    rel_height: float = 0.5,
    plateau_size: _PeakCondition,
) -> tuple[onp.Array1D[np.intp], _PeakProperties_htpws]: ...
