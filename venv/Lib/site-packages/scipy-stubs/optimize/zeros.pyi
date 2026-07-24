# This file is not meant for public use and will be removed in SciPy v2.0.0.

from collections.abc import Callable
from typing_extensions import deprecated

from ._zeros_py import RootResults as _RootResults

__all__ = ["RootResults", "bisect", "brenth", "brentq", "newton", "ridder", "toms748"]

@deprecated("will be removed in SciPy v2.0.0")
class RootResults(_RootResults): ...

@deprecated("will be removed in SciPy v2.0.0")
def newton(
    func: Callable[..., object],
    x0: object,
    fprime: Callable[..., object] | None = None,
    args: tuple[object, ...] = (),
    tol: float = 1.48e-8,
    maxiter: int = 50,
    fprime2: Callable[..., object] | None = None,
    x1: object = None,
    rtol: float = 0.0,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bisect(
    f: Callable[..., object],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    xtol: float = 2e-12,
    rtol: float = ...,
    maxiter: int = 100,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ridder(
    f: Callable[..., object],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    xtol: float = 2e-12,
    rtol: float = ...,
    maxiter: int = 100,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brentq(
    f: Callable[..., object],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    xtol: float = 2e-12,
    rtol: float = ...,
    maxiter: int = 100,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brenth(
    f: Callable[..., object],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    xtol: float = 2e-12,
    rtol: float = ...,
    maxiter: int = 100,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def toms748(
    f: Callable[..., object],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    k: object = 1,
    xtol: float = 2e-12,
    rtol: float = ...,
    maxiter: int = 100,
    full_output: bool = False,
    disp: bool = True,
) -> object: ...
