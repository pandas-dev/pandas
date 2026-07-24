from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Generic, Literal, Never, Protocol, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._resampling import BootstrapResult

__all__ = ["sobol_indices"]

###

_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

type _SobolFunc[T] = Callable[[onp.Array2D[np.float64]], T]

# workaround for mypy & pyright's failure to conform to the overload typing specification
type _JustAnyShape = tuple[Never, Never, Never]

type _SobolKey = Literal["f_A", "f_B", "f_AB"]
type _SobolMethod = Callable[
    [onp.Array2D[np.float64], onp.Array2D[np.float64], onp.Array3D[np.float64]],
    tuple[onp.ToFloat | onp.ToFloatND, onp.ToFloat | onp.ToFloatND],
]
type _ToSobolMethod = _SobolMethod | Literal["saltelli_2010"]

@type_check_only
class _HasPPF(Protocol):
    @property
    def ppf(self, /) -> Callable[..., np.float64]: ...

###

@dataclass
class BootstrapSobolResult:
    first_order: BootstrapResult
    total_order: BootstrapResult

@dataclass
class SobolResult(Generic[_ShapeT_co]):
    first_order: onp.ArrayND[np.float64, _ShapeT_co]
    total_order: onp.ArrayND[np.float64, _ShapeT_co]

    _indices_method: _SobolMethod
    _f_A: onp.Array2D[np.float64]
    _f_B: onp.Array2D[np.float64]
    _f_AB: onp.Array3D[np.float64]
    _A: onp.Array2D[np.float64] | None = None
    _B: onp.Array2D[np.float64] | None = None
    _AB: onp.Array3D[np.float64] | None = None
    _bootstrap_result: BootstrapResult | None = None

    def bootstrap(self, /, confidence_level: onp.ToFloat = 0.95, n_resamples: onp.ToInt = 999) -> BootstrapSobolResult: ...

#
def f_ishigami(x: onp.ToFloat2D) -> onp.Array1D[npc.floating]: ...  # undocumented

#
def sample_A_B(
    n: onp.ToInt, dists: Sequence[_HasPPF], rng: onp.random.ToRNG | None = None
) -> onp.Array2D[np.float64]: ...  # undocumented
def sample_AB(A: onp.Array2D[np.float64], B: onp.Array2D[np.float64]) -> onp.Array3D[np.float64]: ...  # undocumented

#
def saltelli_2010(
    f_A: onp.Array2D[np.float64], f_B: onp.Array2D[np.float64], f_AB: onp.Array3D[np.float64]
) -> tuple[onp.Array2D[np.float64], onp.Array2D[np.float64]]: ...  # undocumented

#
@overload  # (2d f64) -> ?d +complex  (mypy & pyright workaround)
def sobol_indices(
    *,
    func: _SobolFunc[onp.ArrayND[npc.number | np.bool, _JustAnyShape]],
    n: onp.ToJustInt,
    dists: Sequence[_HasPPF] | None = None,
    method: _ToSobolMethod = "saltelli_2010",
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> SobolResult: ...
@overload  # (2d f64) -> 1d +complex  (single output)
def sobol_indices(
    *,
    func: _SobolFunc[onp.ToComplexStrict1D],
    n: onp.ToJustInt,
    dists: Sequence[_HasPPF] | None = None,
    method: _ToSobolMethod = "saltelli_2010",
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> SobolResult[tuple[int]]: ...
@overload  # (2d f64) -> 2d +complex  (multiple outputs)
def sobol_indices(
    *,
    func: _SobolFunc[onp.ToComplexStrict2D],
    n: onp.ToJustInt,
    dists: Sequence[_HasPPF] | None = None,
    method: _ToSobolMethod = "saltelli_2010",
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> SobolResult[tuple[int, int]]: ...
@overload  # fallback
def sobol_indices(
    *,
    func: _SobolFunc[onp.ToComplex2D] | Mapping[_SobolKey, onp.ArrayND[npc.number]],
    n: onp.ToJustInt,
    dists: Sequence[_HasPPF] | None = None,
    method: _ToSobolMethod = "saltelli_2010",
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> SobolResult: ...
