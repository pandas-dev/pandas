from collections.abc import Callable, Sequence
from types import ModuleType
from typing import Final, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["kmeans", "kmeans2", "vq", "whiten"]

_InexactT = TypeVar("_InexactT", bound=npc.inexact)

_InitMethod: TypeAlias = Literal["random", "points", "++", "matrix"]
_MissingMethod: TypeAlias = Literal["warn", "raise"]

_ToFloat32_2D: TypeAlias = onp.ToArray2D[int, np.float32 | np.float16 | npc.integer16 | npc.integer8]
_AsFloat64_2D: TypeAlias = onp.ToArray2D[float, npc.floating64 | npc.integer]
_PyFloatMax2D: TypeAlias = Sequence[float] | Sequence[Sequence[float]]

###
# NOTE: DO NOT RE-ORDER THE OVERLOADS, otherwise it'll trigger a pernicious bug in pyright (1.1.403).

class ClusterError(Exception): ...

#
@overload
def whiten(obs: onp.ArrayND[np.bool_ | npc.integer], check_finite: bool | None = None) -> onp.Array2D[np.float64]: ...
@overload
def whiten(obs: onp.ArrayND[_InexactT], check_finite: bool | None = None) -> onp.Array2D[_InexactT]: ...

#
@overload  # float32
def vq(
    obs: onp.CanArrayND[np.float32], code_book: _ToFloat32_2D, check_finite: bool = True
) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.float32]]: ...
@overload  # float64
def vq(
    obs: onp.ToJustFloat64_2D, code_book: _AsFloat64_2D, check_finite: bool = True
) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.float64]]: ...
@overload  # floating
def vq(
    obs: onp.ToJustFloat2D, code_book: onp.ToFloat2D, check_finite: bool = True
) -> tuple[onp.Array1D[np.int32], onp.Array1D[npc.floating]]: ...

#
@overload  # float64
def py_vq(
    obs: onp.ToFloat64_2D, code_book: onp.ToFloat64_2D, check_finite: bool = True
) -> tuple[onp.Array1D[np.intp], onp.Array1D[np.float64]]: ...
@overload  # floating
def py_vq(
    obs: onp.ToFloat2D, code_book: onp.ToFloat2D, check_finite: bool = True
) -> tuple[onp.Array1D[np.intp], onp.Array1D[npc.floating]]: ...

#
@overload  # float32
def kmeans(
    obs: onp.CanArrayND[np.float32],
    k_or_guess: int | _ToFloat32_2D,
    iter: int = 20,
    thresh: float = 1e-5,
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[np.float32], np.float32]: ...
@overload  # float64
def kmeans(
    obs: onp.ToJustFloat64_2D,
    k_or_guess: int | _AsFloat64_2D,
    iter: int = 20,
    thresh: float = 1e-5,
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[np.float64], np.float64]: ...
@overload  # floating
def kmeans(
    obs: onp.ToJustFloat2D,
    k_or_guess: int | onp.ToFloat2D,
    iter: int = 20,
    thresh: float = 1e-5,
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[npc.floating], npc.floating]: ...

#
def _kpoints(
    data: onp.ArrayND[_InexactT], k: int, rng: onp.random.ToRNG, xp: ModuleType
) -> onp.Array2D[_InexactT]: ...  # undocumented
def _krandinit(
    data: onp.ArrayND[npc.inexact], k: int, rng: onp.random.ToRNG, xp: ModuleType
) -> onp.Array2D[np.float64]: ...  # undocumented
def _kpp(
    data: onp.ArrayND[npc.inexact], k: int, rng: onp.random.ToRNG, xp: ModuleType
) -> onp.Array2D[np.float64]: ...  # undocumented

_valid_init_meth: Final[
    dict[str, Callable[[onp.ArrayND[npc.inexact], int, onp.random.ToRNG, ModuleType], onp.Array2D[npc.inexact]]]
] = ...  # undocumented

def _missing_warn() -> None: ...  # undocumented
def _missing_raise() -> None: ...  # undocumented

_valid_miss_meth: Final[dict[str, Callable[[], None]]] = ...  # undocumented

# NOTE: There is a false positive `overload-overlap` mypy error that only occurs with `numpy<2.2`
# mypy: disable-error-code=overload-overlap

#
@overload  # float32
def kmeans2(
    data: onp.CanArrayND[np.float32],
    k: int | _ToFloat32_2D,
    iter: int = 10,
    thresh: float = 1e-5,
    minit: _InitMethod = "random",
    missing: _MissingMethod = "warn",
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[np.float32], onp.Array1D[np.int32]]: ...
@overload  # float64
def kmeans2(
    data: onp.CanArrayND[np.float64] | _PyFloatMax2D,
    k: int | _AsFloat64_2D,
    iter: int = 10,
    thresh: float = 1e-5,
    minit: _InitMethod = "random",
    missing: _MissingMethod = "warn",
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.int32]]: ...
@overload  # floating
def kmeans2(
    data: onp.CanArrayND[npc.floating] | _PyFloatMax2D,
    k: int | onp.ToFloat2D,
    iter: int = 10,
    thresh: float = 1e-5,
    minit: _InitMethod = "random",
    missing: _MissingMethod = "warn",
    check_finite: bool = True,
    *,
    seed: onp.random.ToRNG | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[npc.floating], onp.Array1D[np.int32]]: ...
