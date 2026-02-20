from collections.abc import Callable
from typing import Final, TypeAlias

import optype.numpy as onp

_Fun: TypeAlias = Callable[[float], onp.ToFloat]

###

class DCSRCH:
    phi: Final[_Fun]
    derphi: Final[_Fun]
    ftol: Final[float]
    gtol: Final[float]
    xtol: Final[float]
    stpmin: Final[float]
    stpmax: Final[float]

    stage: int | None
    ginit: float | None
    gtest: float | None
    gx: float | None
    gy: float | None
    finit: float | None
    fx: float | None
    fy: float | None
    stx: float | None
    sty: float | None
    stmin: float | None
    stmax: float | None
    width: float | None
    width1: float | None

    def __init__(
        self, /, phi: _Fun, derphi: _Fun, ftol: float, gtol: float, xtol: float, stpmin: float, stpmax: float
    ) -> None: ...
    def __call__(
        self, /, alpha1: float, phi0: float | None = None, derphi0: float | None = None, maxiter: int = 100
    ) -> tuple[float | None, float, float, bytes]: ...  # alpha, phi(alpha), phi(0), task

def dcstep(
    stx: float,
    fx: float,
    dx: float,
    sty: float,
    fy: float,
    dy: float,
    stp: float,
    fp: float,
    dp: float,
    brackt: bool,
    stpmin: float,
    stpmax: float,
) -> tuple[float, float, float, float, float, float, float, bool]: ...  # stx, fx, dx, sty, fy, dy, stp, brackt
