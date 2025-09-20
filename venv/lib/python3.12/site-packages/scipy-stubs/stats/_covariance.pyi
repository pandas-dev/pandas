from typing import Final, Generic, Protocol, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["Covariance"]

_ScalarT = TypeVar("_ScalarT", bound=npc.floating | npc.integer)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.floating | npc.integer, default=np.float64, covariant=True)

class Covariance(Generic[_ScalarT_co]):
    @staticmethod
    @overload
    def from_diagonal(diagonal: onp.ToJustFloat64_1D) -> CovViaDiagonal[np.float64]: ...
    @staticmethod
    @overload
    def from_diagonal(diagonal: onp.ToJustInt64_1D) -> CovViaDiagonal[np.int_]: ...
    @staticmethod
    @overload
    def from_diagonal(diagonal: onp.ToArray1D[_ScalarT, _ScalarT]) -> CovViaDiagonal[_ScalarT]: ...

    #
    @staticmethod
    def from_precision(precision: onp.ToFloat2D, covariance: onp.ToFloat2D | None = None) -> CovViaPrecision: ...
    @staticmethod
    def from_cholesky(cholesky: onp.ToFloat2D) -> CovViaCholesky: ...
    @staticmethod
    def from_eigendecomposition(eigendecomposition: tuple[onp.ToFloat1D, onp.ToFloat2D]) -> CovViaEigendecomposition: ...

    #
    @property
    def log_pdet(self, /) -> np.float64: ...
    @property
    def rank(self, /) -> np.int_: ...
    @property
    def covariance(self, /) -> onp.Array2D[_ScalarT_co]: ...
    @property
    def shape(self, /) -> tuple[int, int]: ...

    #
    def whiten(self, /, x: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...
    def colorize(self, /, x: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...

class CovViaDiagonal(Covariance[_ScalarT_co], Generic[_ScalarT_co]):
    @overload
    def __init__(self: CovViaDiagonal[np.float64], /, diagonal: onp.ToJustFloat64_1D) -> None: ...
    @overload
    def __init__(self: CovViaDiagonal[np.int_], /, diagonal: onp.ToJustInt64_1D) -> None: ...
    @overload
    def __init__(self, /, diagonal: onp.ToArray1D[_ScalarT_co, _ScalarT_co]) -> None: ...

class CovViaPrecision(Covariance[np.float64]):
    def __init__(self, /, precision: onp.ToFloat2D, covariance: onp.ToFloat2D | None = None) -> None: ...

class CovViaCholesky(Covariance[np.float64]):
    def __init__(self, /, cholesky: onp.ToFloat2D) -> None: ...

class CovViaEigendecomposition(Covariance[np.float64]):
    def __init__(self, /, eigendecomposition: tuple[onp.ToFloat1D, onp.ToFloat2D]) -> None: ...

@type_check_only
class _PSD(Protocol):
    _M: onp.ArrayND[np.float64]
    V: onp.ArrayND[np.float64]
    U: onp.ArrayND[np.float64]
    eps: float
    log_pdet: float
    cond: float
    rank: int

    @property
    def pinv(self, /) -> onp.ArrayND[npc.floating]: ...

class CovViaPSD(Covariance[np.float64]):
    _LP: Final[onp.ArrayND[np.float64]]
    _log_pdet: Final[float]
    _rank: Final[int]
    _covariance: Final[onp.ArrayND[np.float64]]
    _shape: tuple[int, int]
    _psd: Final[_PSD]
    _allow_singular: Final = False

    def __init__(self, /, psd: _PSD) -> None: ...
