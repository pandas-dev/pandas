from collections.abc import Callable, Sequence
from typing import Concatenate, Literal, TypeAlias, TypeVar, TypedDict, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["argrelextrema", "argrelmax", "argrelmin", "find_peaks", "find_peaks_cwt", "peak_prominences", "peak_widths"]

_SCT = TypeVar("_SCT", bound=np.generic)

_Int1D: TypeAlias = onp.Array1D[np.intp]
_IntND: TypeAlias = onp.ArrayND[np.intp]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]

_Mode: TypeAlias = Literal["clip", "wrap"]

_ArgRel: TypeAlias = tuple[_IntND, ...]
_PeakProminences: TypeAlias = tuple[_FloatND, _IntND, _IntND]
_PeakWidths: TypeAlias = tuple[_FloatND, _FloatND, _FloatND, _FloatND]

_WaveletFunc: TypeAlias = (
    Callable[Concatenate[int, float, ...], onp.ToComplex1D] | Callable[Concatenate[np.intp, np.float64, ...], onp.ToComplex1D]
)

@type_check_only
class _FindPeaksResultsDict(TypedDict, total=False):
    plateau_sizes: _Int1D
    left_edges: _Int1D
    right_edges: _Int1D

    peak_heights: _Float1D

    left_thresholds: _Float1D
    right_thresholds: _Float1D

    prominences: _Float1D
    left_bases: _Int1D
    right_bases: _Int1D

    widths: _Float1D
    width_heights: _Float1D
    left_ips: _Float1D
    right_ips: _Float1D

###

def argrelmin(data: onp.Array, axis: op.CanIndex = 0, order: onp.ToInt = 1, mode: _Mode = "clip") -> _ArgRel: ...
def argrelmax(data: onp.Array, axis: op.CanIndex = 0, order: onp.ToInt = 1, mode: _Mode = "clip") -> _ArgRel: ...
def argrelextrema(
    data: onp.ArrayND[_SCT],
    comparator: Callable[[onp.ArrayND[_SCT], onp.ArrayND[_SCT]], onp.ToBoolND],
    axis: op.CanIndex = 0,
    order: onp.ToInt = 1,
    mode: _Mode = "clip",
) -> _ArgRel: ...

#
def peak_prominences(x: onp.ToArray1D, peaks: onp.ToIntND, wlen: onp.ToFloat | None = None) -> _PeakProminences: ...
def peak_widths(
    x: onp.ToArray1D,
    peaks: onp.ToIntND,
    rel_height: onp.ToFloat = 0.5,
    prominence_data: _PeakProminences | None = None,
    wlen: onp.ToFloat | None = None,
) -> _PeakWidths: ...

#
def find_peaks(
    x: onp.ToArray1D,
    height: onp.ToFloat | onp.ToFloatND | Sequence[onp.ToFloat | None] | None = None,
    threshold: onp.ToInt | onp.ToFloatND | Sequence[onp.ToFloat | None] | None = None,
    distance: onp.ToFloat | None = None,
    prominence: onp.ToFloat | onp.ToFloatND | Sequence[onp.ToFloat | None] | None = None,
    width: onp.ToFloat | onp.ToFloatND | Sequence[onp.ToFloat | None] | None = None,
    wlen: onp.ToFloat | None = None,
    rel_height: onp.ToFloat = 0.5,
    plateau_size: onp.ToInt | onp.ToIntND | Sequence[onp.ToInt | None] | None = None,
) -> tuple[_Int1D, _FindPeaksResultsDict]: ...

#
def find_peaks_cwt(
    vector: onp.Array,
    widths: onp.ToFloat | onp.ToFloatND,
    wavelet: _WaveletFunc | None = None,
    max_distances: onp.ArrayND[npc.floating | npc.integer] | None = None,
    gap_thresh: onp.ToFloat | None = None,
    min_length: onp.ToInt | None = None,
    min_snr: onp.ToFloat = 1,
    noise_perc: onp.ToFloat = 10,
    window_size: onp.ToInt | None = None,
) -> _Int1D: ...
