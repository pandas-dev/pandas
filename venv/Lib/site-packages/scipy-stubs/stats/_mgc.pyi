from collections.abc import Callable, Iterable, Sequence
from typing import Self, TypedDict, override, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._typing import BaseBunch

__all__ = ["multiscale_graphcorr"]

###

@type_check_only
class _MGCDict(TypedDict):
    mgc_map: onp.Array2D[np.float64]
    opt_scale: Sequence[int | np.intp]  # list of size 2
    null_dist: onp.Array1D[np.float64]

###

class MGCResult(BaseBunch[np.float64, np.float64, _MGCDict]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    def pvalue(self, /) -> np.float64: ...
    @property
    def mgc_dict(self, /) -> _MGCDict: ...

    #
    @override
    def __new__(_cls, statistic: np.float64, pvalue: np.float64, mgc_dict: _MGCDict) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, statistic: np.float64, pvalue: np.float64, mgc_dict: _MGCDict) -> None: ...  # pyrefly:ignore[bad-override]

def multiscale_graphcorr[T, R](
    x: onp.ArrayND[npc.floating | npc.integer | np.bool],
    y: onp.ArrayND[npc.floating | npc.integer | np.bool],
    compute_distance: Callable[[onp.ArrayND[np.float64]], onp.ArrayND[npc.floating]] = ...,
    reps: int = 1000,
    workers: int | Callable[[Callable[[T], R], Iterable[T]], Sequence[R]] = 1,
    is_twosamp: bool = False,
    random_state: onp.random.ToRNG | None = None,
) -> MGCResult: ...
