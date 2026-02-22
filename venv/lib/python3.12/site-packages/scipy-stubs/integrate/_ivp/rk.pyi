from collections.abc import Callable
from typing import ClassVar, Final, Generic, Never, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

from .base import DenseOutput, OdeSolver

_SCT_fc = TypeVar("_SCT_fc", bound=np.float64 | np.complex128, default=np.float64 | np.complex128)

###

SAFETY: Final = 0.9
MIN_FACTOR: Final = 0.2
MAX_FACTOR: Final = 10

class RungeKutta(OdeSolver, Generic[_SCT_fc]):
    C: ClassVar[onp.ArrayND[np.float64]]
    A: ClassVar[onp.ArrayND[np.float64]]
    B: ClassVar[onp.ArrayND[np.float64]]
    E: ClassVar[onp.ArrayND[np.float64]]
    P: ClassVar[onp.ArrayND[np.float64]]
    order: ClassVar[int]
    error_estimator_order: ClassVar[int]
    n_stages: ClassVar[int]

    y_old: onp.ArrayND[_SCT_fc] | None
    f: onp.ArrayND[_SCT_fc]
    K: onp.ArrayND[_SCT_fc]
    max_step: float
    h_abs: float
    error_exponent: float
    h_previous: float | None

    def __init__(
        self,
        /,
        fun: Callable[[float, onp.ArrayND[_SCT_fc]], onp.ArrayND[_SCT_fc]],
        t0: float,
        y0: onp.ArrayND[_SCT_fc],
        t_bound: float,
        max_step: float = ...,
        rtol: float = 0.001,
        atol: float = 1e-06,
        vectorized: bool = False,
        first_step: float | None = None,
        **extraneous: Never,
    ) -> None: ...

class RK23(RungeKutta[_SCT_fc], Generic[_SCT_fc]): ...
class RK45(RungeKutta[_SCT_fc], Generic[_SCT_fc]): ...

class DOP853(RungeKutta[_SCT_fc], Generic[_SCT_fc]):
    E3: ClassVar[onp.ArrayND[np.float64]]
    E5: ClassVar[onp.ArrayND[np.float64]]
    D: ClassVar[onp.ArrayND[np.float64]]
    A_EXTRA: ClassVar[onp.ArrayND[np.float64]]
    C_EXTRA: ClassVar[onp.ArrayND[np.float64]]

    K_extended: onp.ArrayND[_SCT_fc]

class RkDenseOutput(DenseOutput[_SCT_fc], Generic[_SCT_fc]):
    h: float
    order: int
    Q: onp.ArrayND[_SCT_fc]
    y_old: onp.ArrayND[_SCT_fc]

    def __init__(self, /, t_old: float, t: float, y_old: onp.ArrayND[_SCT_fc], Q: onp.ArrayND[_SCT_fc]) -> None: ...

class Dop853DenseOutput(DenseOutput[_SCT_fc], Generic[_SCT_fc]):
    h: float
    F: onp.ArrayND[_SCT_fc]
    y_old: onp.ArrayND[_SCT_fc]

    def __init__(self, /, t_old: float, t: float, y_old: onp.ArrayND[_SCT_fc], F: onp.ArrayND[_SCT_fc]) -> None: ...

@overload
def rk_step(
    fun: Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]],
    t: float,
    y: onp.ArrayND[np.float64],
    f: onp.ArrayND[np.float64],
    h: float,
    A: onp.ArrayND[np.float64],
    B: onp.ArrayND[np.float64],
    C: onp.ArrayND[np.float64],
    K: onp.ArrayND[np.float64],
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def rk_step(
    fun: Callable[[float, onp.ArrayND[np.complex128]], onp.ArrayND[np.complex128]],
    t: float,
    y: onp.ArrayND[np.complex128],
    f: onp.ArrayND[np.complex128],
    h: float,
    A: onp.ArrayND[np.float64],
    B: onp.ArrayND[np.float64],
    C: onp.ArrayND[np.float64],
    K: onp.ArrayND[np.float64],
) -> tuple[onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...
