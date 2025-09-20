from typing import Any, Final, Generic, Literal, Self, TypeAlias
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase

__all__ = ["expm", "inv", "matrix_power"]

_SCT_co = TypeVar("_SCT_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_SparseT = TypeVar("_SparseT", bound=_spbase)

_Structure: TypeAlias = Literal["upper_triangular"]

###

UPPER_TRIANGULAR: Final[_Structure] = "upper_triangular"

class MatrixPowerOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    @property
    @override
    def T(self, /) -> Self: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    def __init__(self, /, A: onp.Array2D[_SCT_co] | _spbase, p: int, structure: _Structure | None = None) -> None: ...

class ProductOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    @property
    @override
    def T(self, /) -> Self: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    def __init__(self, /, *args: onp.Array2D[_SCT_co] | _spbase, structure: _Structure | None = None) -> None: ...

def inv(A: _SparseT) -> _SparseT: ...
def expm(A: _SparseT) -> _SparseT: ...
def matrix_power(A: _SparseT, power: op.CanIndex) -> _SparseT: ...
