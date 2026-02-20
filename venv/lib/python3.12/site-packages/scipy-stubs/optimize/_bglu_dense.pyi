from collections.abc import Callable
from typing import Any, Concatenate, Generic, Never
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_MethodT = TypeVar("_MethodT", bound=Callable[Concatenate[Any, ...], object])
_NumberT_co = TypeVar("_NumberT_co", bound=npc.number, default=Any, covariant=True)

###

__all__ = ["BGLU", "LU"]

def _consider_refactor(method: _MethodT) -> _MethodT: ...  # undocumented

class LU(Generic[_NumberT_co]):  # undocumented
    A: onp.Array2D[_NumberT_co]
    B: onp.Array2D[_NumberT_co]
    b: onp.Array1D[np.intp]
    m: int
    n: int

    def __init__(self, /, A: onp.Array2D[_NumberT_co], b: onp.Array1D[np.intp]) -> None: ...
    def __setstate_cython__(self, /, state: tuple[object, ...]) -> None: ...
    def __reduce_cython__(self) -> tuple[Any, ...]: ...
    def update(self, /, i: int, j: int) -> None: ...
    def solve(self, /, q: onp.ArrayND[npc.number], transposed: bool = False) -> onp.ArrayND[_NumberT_co]: ...

class BGLU(LU[_NumberT_co], Generic[_NumberT_co]):  # undocumented
    plu: tuple[onp.ArrayND[Any], ...]
    L: onp.Array2D[_NumberT_co]
    U: onp.Array2D[_NumberT_co]
    pi: onp.Array2D[_NumberT_co]
    pit: onp.Array1D[np.int_]
    ops_list: list[onp.Array2D[np.float64]]
    bglu_time: float
    solves: int
    updates: int
    max_updates: int
    average_solve_times: list[float]
    mast: bool

    def __init__(
        self, /, A: onp.Array2D[_NumberT_co], b: onp.Array1D[np.intp], max_updates: int = 10, mast: bool = False
    ) -> None: ...
    def refactor(self, /, *args: Never) -> None: ...  # *args are needed for stubtest
    def update_basis(self, /, i: int, j: int) -> None: ...
    def perform_perm(self, /, p: onp.Array1D[np.intp]) -> onp.Array1D[np.intp]: ...
