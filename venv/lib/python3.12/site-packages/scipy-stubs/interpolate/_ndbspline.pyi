from collections.abc import Callable
from typing import Any, Concatenate, Generic, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array

__all__ = ["NdBSpline"]

###

_ScalarT = TypeVar("_ScalarT", bound=np.generic)
_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

_ToKnots: TypeAlias = tuple[onp.ToFloat1D, ...]
_ToDegrees: TypeAlias = op.CanIndex | tuple[op.CanIndex, ...]

_DesignMatrix: TypeAlias = csr_array[np.float64, tuple[int, int]]
_SolverFunc: TypeAlias = Callable[Concatenate[_DesignMatrix, onp.Array2D[np.float64], ...], onp.ArrayND[_ScalarT]]

###

class NdBSpline(Generic[_CT_co]):
    c: onp.ArrayND[np.float64]
    extrapolate: bool

    @property
    def k(self, /) -> tuple[np.int32, ...]: ...
    @property
    def t(self, /) -> tuple[onp.Array1D[np.float64], ...]: ...

    #
    @overload
    def __init__(
        self: NdBSpline[np.float64], /, t: _ToKnots, c: onp.ToFloatND, k: _ToDegrees, *, extrapolate: onp.ToBool | None = None
    ) -> None: ...
    @overload
    def __init__(
        self: NdBSpline[np.complex128],
        /,
        t: _ToKnots,
        c: onp.ToJustComplexND,
        k: _ToDegrees,
        *,
        extrapolate: onp.ToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: NdBSpline[Any], /, t: _ToKnots, c: onp.ToComplexND, k: _ToDegrees, *, extrapolate: onp.ToBool | None = None
    ) -> None: ...

    #
    def __call__(
        self, /, xi: onp.ToFloatND, *, nu: onp.ToFloat1D | None = None, extrapolate: onp.ToBool | None = None
    ) -> onp.ArrayND[_CT_co]: ...

    #
    @classmethod
    def design_matrix(cls, xvals: onp.ToFloat2D, t: _ToKnots, k: _ToDegrees, extrapolate: onp.ToBool = True) -> _DesignMatrix: ...

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
