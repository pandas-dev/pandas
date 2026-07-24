from collections.abc import Callable
from typing import Any
from typing_extensions import TypeAlias, TypeVar

import numpy as np

__all__ = ["resample", "resample_nu"]

# np.floating[Any] because precision is not important
_FloatArray = TypeVar("_FloatArray", bound=np.ndarray[tuple[int, ...], np.dtype[np.floating[Any]]])
_FilterType: TypeAlias = str | Callable[[int], np.ndarray[tuple[int], np.dtype[np.float64]]]

def resample(
    x: _FloatArray,
    sr_orig: float,
    sr_new: float,
    axis: int = -1,
    filter: _FilterType = "kaiser_best",
    parallel: bool = False,
    *,
    num_zeros: int = 64,
    precision: int = 9,
    rolloff: float = 0.945,
) -> _FloatArray: ...
def resample_nu(
    x: _FloatArray,
    sr_orig: float,
    t_out: _FloatArray,
    axis: int = -1,
    filter: _FilterType = "kaiser_best",
    parallel: bool = False,
    *,
    num_zeros: int = 64,
    precision: int = 9,
    rolloff: float = 0.945,
) -> _FloatArray: ...
