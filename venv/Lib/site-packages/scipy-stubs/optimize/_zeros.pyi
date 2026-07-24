# defined in scipy/optimize/zeros.c

from collections.abc import Callable
from typing import Concatenate

def _bisect(
    f: Callable[Concatenate[float, ...], float],
    a: float,
    b: float,
    xtol: float,
    rtol: float,
    iter: int,
    xargs: tuple[object, ...],
    fulloutput: int,
    disp: int = 1,
) -> tuple[float, int, int, int]: ...
def _brenth(
    f: Callable[Concatenate[float, ...], float],
    a: float,
    b: float,
    xtol: float,
    rtol: float,
    iter: int,
    xargs: tuple[object, ...],
    fulloutput: int,
    disp: int = 1,
) -> tuple[float, int, int, int]: ...
def _brentq(
    f: Callable[Concatenate[float, ...], float],
    a: float,
    b: float,
    xtol: float,
    rtol: float,
    iter: int,
    xargs: tuple[object, ...],
    fulloutput: int,
    disp: int = 1,
) -> tuple[float, int, int, int]: ...
def _ridder(
    f: Callable[Concatenate[float, ...], float],
    a: float,
    b: float,
    xtol: float,
    rtol: float,
    iter: int,
    xargs: tuple[object, ...],
    fulloutput: int,
    disp: int = 1,
) -> tuple[float, int, int, int]: ...
