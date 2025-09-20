__all__ = ["cubature"]

import dataclasses
from collections.abc import Callable, Iterable, Sequence
from types import ModuleType
from typing import Any, Concatenate, Generic, Literal, TypeAlias
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

_VT = TypeVar("_VT")
_RT = TypeVar("_RT")
_ShapeT_co = TypeVar("_ShapeT_co", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any], covariant=True)

_ArrayAPINamespace: TypeAlias = ModuleType
_Rule: TypeAlias = Literal["gk21", "gk15", "gauss-kronrod", "genz-malik"]
_Status: TypeAlias = Literal["converged", "not_converged"]

###

@dataclasses.dataclass
class CubatureRegion(Generic[_ShapeT_co]):
    estimate: onp.ArrayND[np.float64, _ShapeT_co]
    error: onp.ArrayND[np.float64, _ShapeT_co]
    a: onp.Array1D[np.float64]
    b: onp.Array1D[np.float64]
    _xp: _ArrayAPINamespace

    def __lt__(self, other: CubatureRegion, /) -> bool: ...

@dataclasses.dataclass
class CubatureResult(Generic[_ShapeT_co]):
    estimate: onp.ArrayND[np.float64, _ShapeT_co]
    error: onp.ArrayND[np.float64, _ShapeT_co]
    status: _Status
    regions: list[CubatureRegion[_ShapeT_co]]
    subdivisions: int
    atol: float
    rtol: float

def cubature(
    f: Callable[Concatenate[onp.Array2D[np.float64], ...], onp.ToFloatND],
    a: onp.ToFloat1D,
    b: onp.ToFloat1D,
    *,
    rule: _Rule = "gk21",
    rtol: float = 1e-8,
    atol: float = 0,
    max_subdivisions: onp.ToJustInt = 10_000,
    args: tuple[object, ...] = (),
    workers: int | Callable[[Callable[[_VT], _RT], Iterable[_VT]], Sequence[_RT]] = 1,
    points: Sequence[onp.ToFloat1D] | None = None,
) -> CubatureResult: ...
