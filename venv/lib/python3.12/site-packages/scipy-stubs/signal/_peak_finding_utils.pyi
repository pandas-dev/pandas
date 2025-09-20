# scipy/signal/_peak_finding_utils.pyx
from typing import TypeAlias

import numpy as np
import optype.numpy as onp

__all__ = ["_local_maxima_1d", "_peak_prominences", "_peak_widths", "_select_by_peak_distance"]

_Array_b_1d: TypeAlias = onp.Array1D[np.bool_]
_Array_n_1d: TypeAlias = onp.Array1D[np.intp]
_Array_f8_1d: TypeAlias = onp.Array1D[np.float64]

class PeakPropertyWarning(RuntimeWarning): ...  # undocumented

def _local_maxima_1d(x: _Array_f8_1d) -> tuple[_Array_n_1d, _Array_n_1d, _Array_n_1d]: ...
def _select_by_peak_distance(peaks: _Array_n_1d, priority: _Array_f8_1d, distance: np.float64) -> _Array_b_1d: ...
def _peak_prominences(x: _Array_f8_1d, peaks: _Array_n_1d, wlen: np.intp) -> tuple[_Array_f8_1d, _Array_n_1d, _Array_n_1d]: ...
def _peak_widths(
    x: _Array_f8_1d,
    peaks: _Array_n_1d,
    rel_height: np.float64,
    prominences: _Array_f8_1d,
    left_bases: _Array_n_1d,
    right_bases: _Array_n_1d,
) -> tuple[_Array_f8_1d, _Array_f8_1d, _Array_f8_1d, _Array_f8_1d]: ...
