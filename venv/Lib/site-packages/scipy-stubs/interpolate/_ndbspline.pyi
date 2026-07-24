import types
from collections.abc import Callable
from typing import Any, Concatenate, Generic, Self, SupportsIndex, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array

__all__ = ["NdBSpline"]

###

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

type _ToKnots = tuple[onp.ToFloat1D, ...]
type _ToDegrees = SupportsIndex | tuple[SupportsIndex, ...]

type _DesignMatrix = csr_array[np.float64, tuple[int, int]]
type _SolverFunc[ScalarT: np.generic] = Callable[Concatenate[_DesignMatrix, onp.Array2D[np.float64], ...], onp.ArrayND[ScalarT]]

###

class NdBSpline(Generic[_CT_co]):
    extrapolate: bool

    @property
    def k(self, /) -> tuple[np.int64, ...]: ...
    @property
    def t(self, /) -> tuple[onp.Array1D[np.float64], ...]: ...
    @property
    def c(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @classmethod
    def __class_getitem__(cls, arg: type | object, /) -> types.GenericAlias: ...
    @classmethod
    def design_matrix(cls, xvals: onp.ToFloat2D, t: _ToKnots, k: _ToDegrees, extrapolate: bool = True) -> _DesignMatrix: ...

    #
    @overload
    def __init__(
        self: NdBSpline[np.float64], /, t: _ToKnots, c: onp.ToFloatND, k: _ToDegrees, *, extrapolate: bool | None = None
    ) -> None: ...
    @overload
    def __init__(
        self: NdBSpline[np.complex128], /, t: _ToKnots, c: onp.ToJustComplexND, k: _ToDegrees, *, extrapolate: bool | None = None
    ) -> None: ...
    @overload
    def __init__(
        self: NdBSpline[Any], /, t: _ToKnots, c: onp.ToComplexND, k: _ToDegrees, *, extrapolate: bool | None = None
    ) -> None: ...

    #
    def __call__(
        self, /, xi: onp.ToFloatND, *, nu: onp.ToFloat1D | None = None, extrapolate: bool | None = None
    ) -> onp.ArrayND[_CT_co]: ...

    #
    def derivative(self, /, nu: onp.ToInt1D) -> Self: ...

#
@overload
def make_ndbspl(
    points: _ToKnots,
    values: onp.ToFloatND,
    k: _ToDegrees = 3,
    *,
    solver: _SolverFunc[npc.floating | npc.integer] = ...,  # scipy.sparse.linalg.gcrotmk
    **solver_args: object,
) -> NdBSpline[np.float64]: ...
@overload
def make_ndbspl(
    points: _ToKnots, values: onp.ToFloatND, k: _ToDegrees = 3, *, solver: _SolverFunc[npc.complexfloating], **solver_args: object
) -> NdBSpline[np.complex128]: ...
