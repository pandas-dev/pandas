from _typeshed import Incomplete
from collections.abc import Callable
from types import GenericAlias, ModuleType
from typing import Any, ClassVar, Final, Generic, Protocol, Self, SupportsIndex, final, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
import numpy_typing_compat as nptc
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

__all__ = ["LinearOperator", "aslinearoperator"]

###

type _Real = npc.floating | npc.integer | np.bool
type _Scalar = npc.number | np.bool

type _Shape = tuple[int, int, *tuple[int, ...]]
type _AnyShape = tuple[int, int, *tuple[Any, ...]]

type _FunMatVec = Callable[[onp.ArrayND[Any]], onp.ToComplexND]
type _FunMatMat = Callable[[onp.Array2D[Any]], onp.ToComplexND]

_SCT_co = TypeVar("_SCT_co", bound=npc.number | np.bool, default=Any, covariant=True)
_SCT1_co = TypeVar("_SCT1_co", bound=npc.number | np.bool, default=Any, covariant=True)
_SCT2_co = TypeVar("_SCT2_co", bound=npc.number | np.bool, default=_SCT1_co, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=_Shape, default=_AnyShape, covariant=True)
_LinearOperatorT_co = TypeVar("_LinearOperatorT_co", bound=LinearOperator[Any, Any], covariant=True)

@type_check_only
class _CanAdjoint(Protocol[_LinearOperatorT_co]):
    def _adjoint(self, /) -> _LinearOperatorT_co: ...

@type_check_only
class _HasShapeAndMatVec(Protocol[_SCT_co, _ShapeT_co]):
    @property
    def shape(self, /) -> _ShapeT_co: ...

    #
    @overload
    def matvec(self, /, x: onp.CanArrayND[np.float64]) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def matvec(self, /, x: onp.CanArrayND[np.complex128]) -> onp.ToComplexND: ...

@type_check_only
class _HasShapeAndDTypeAndMatVec(Protocol[_SCT_co, _ShapeT_co]):
    @property
    def dtype(self, /) -> np.dtype[_SCT_co]: ...
    @property
    def shape(self, /) -> _ShapeT_co: ...

    #
    def matvec(self, /, x: onp.CanArrayND[np.float64] | onp.CanArrayND[np.complex128]) -> onp.ToComplexND: ...

###

class LinearOperator(Generic[_SCT_co, _ShapeT_co]):
    __array_ufunc__: ClassVar[None] = None

    @classmethod
    def __class_getitem__(cls, arg: object, /) -> GenericAlias: ...

    dtype: np.dtype[_SCT_co]
    shape: _ShapeT_co
    ndim: Final[int]
    _xp: Final[ModuleType]

    # keep in sync with `_CustomLinearOperator.__init__`
    @overload  # no dtype
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: None = None,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.int8 | Any, ShapeT]: ...
    @overload  # dtype known (positional)
    def __new__[SCT: _Scalar, ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.ToDType[SCT],
        rmatmat: _FunMatMat | None = None,
        xp: ModuleType | None = None,
    ) -> _CustomLinearOperator[SCT, ShapeT]: ...
    @overload  # dtype known (keyword)
    def __new__[SCT: _Scalar, ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.ToDType[SCT],
        rmatmat: _FunMatMat | None = None,
        xp: ModuleType | None = None,
    ) -> _CustomLinearOperator[SCT, ShapeT]: ...
    @overload  # dtype-like int_ (positional)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.int_, ShapeT]: ...
    @overload  # dtype-like int_ (keyword)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.int_, ShapeT]: ...
    @overload  # dtype-like float64 (positional)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.float64, ShapeT]: ...
    @overload  # dtype-like float64 (positional)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.float64, ShapeT]: ...
    @overload  # dtype-like complex128 (positional)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.complex128, ShapeT]: ...
    @overload  # dtype-like complex128 (keyword)
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[np.complex128, ShapeT]: ...
    @overload  # unknown dtype
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: type | str | None = None,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> _CustomLinearOperator[Any, ShapeT]: ...
    @overload  # xp given
    def __new__[ShapeT: _Shape](
        cls,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: Incomplete | None = None,
        rmatmat: _FunMatMat | None = None,
        *,
        xp: ModuleType,
    ) -> _CustomLinearOperator[Any, ShapeT]: ...

    # NOTE: the `__init__` method cannot be annotated, because it will cause mypy to ignore `__new__`:
    # https://github.com/python/mypy/issues/17251

    # ruff: noqa: ERA001

    # @overload
    # def __init__(self, /, dtype: onp.ToDType[_SCT_co], shape: _ShapeT_co) -> None: ...
    # @overload
    # def __init__[ShapeT: _Shape](self: LinearOperator[np.int_, ShapeT], /, dtype: onp.AnyIntDType, shape: ShapeT) -> None: ...
    # @overload
    # def __init__[ShapeT: _Shape](
    #     self: LinearOperator[np.float64, ShapeT], /, dtype: onp.AnyFloat64DType, shape: ShapeT
    # ) -> None: ...
    # @overload
    # def __init__[ShapeT: _Shape](
    #     self: LinearOperator[np.complex128, ShapeT], /, dtype: onp.AnyComplex128DType, shape: ShapeT
    # ) -> None: ...
    # @overload
    # def __init__[ShapeT: _Shape](self: LinearOperator[Any, ShapeT], /, dtype: type | str | None, shape: ShapeT) -> None: ...

    @override
    def __getstate__(self, /) -> dict[str, Any]: ...
    def __setstate__(self, state: dict[str, Any], /) -> None: ...

    #
    @property
    def H[LinearOperatorT: LinearOperator[Any, Any]](self: _CanAdjoint[LinearOperatorT], /) -> LinearOperatorT: ...
    @final
    def adjoint[LinearOperatorT: LinearOperator[Any, Any]](self: _CanAdjoint[LinearOperatorT], /) -> LinearOperatorT: ...
    def _adjoint(self, /) -> LinearOperator[_SCT_co, _ShapeT_co]: ...

    #
    @property
    def T(self, /) -> _TransposedLinearOperator[_SCT_co, _ShapeT_co]: ...
    def transpose(self, /) -> _TransposedLinearOperator[_SCT_co, _ShapeT_co]: ...

    #
    @overload  # float matrix
    def matvec(self, /, x: onp.Matrix[_Real]) -> onp.Matrix[_SCT_co]: ...
    @overload  # complex matrix
    def matvec(self, /, x: onp.Matrix[np.complex128]) -> onp.Matrix[np.complex128]: ...
    @overload  # float array
    def matvec(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload  # complex array
    def matvec(self, /, x: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload  # unknown array
    def matvec(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    rmatvec = matvec

    #
    @overload
    def matmat(self, /, X: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def matmat(self, /, X: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def matmat(self, /, X: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    rmatmat = matmat

    #
    @overload
    def dot[SCT: _Scalar, ShapeT: _Shape](
        self, /, x: LinearOperator[SCT, ShapeT]
    ) -> _ProductLinearOperator[_SCT_co, SCT, ShapeT]: ...
    @overload
    def dot(self, /, x: onp.ToFloat) -> _ScaledLinearOperator[_SCT_co, _ShapeT_co]: ...
    @overload
    def dot(self, /, x: onp.ToJustComplex128) -> _ScaledLinearOperator[np.complex128, _ShapeT_co]: ...
    @overload
    def dot(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def dot(self, /, x: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def dot(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    __mul__ = dot

    # keep in sync with `dot`
    @overload
    def rdot[SCT: _Scalar, ShapeT: _Shape](
        self, /, x: LinearOperator[SCT, ShapeT]
    ) -> _ProductLinearOperator[_SCT_co, SCT, ShapeT]: ...
    @overload
    def rdot(self, /, x: onp.ToFloat) -> _ScaledLinearOperator[_SCT_co, _ShapeT_co]: ...
    @overload
    def rdot(self, /, x: onp.ToJustComplex128) -> _ScaledLinearOperator[np.complex128, _ShapeT_co]: ...
    @overload
    def rdot(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def rdot(self, /, x: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def rdot(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    __rmul__ = rdot

    #
    @overload
    def __matmul__[SCT: _Scalar, ShapeT: _Shape](
        self, /, x: LinearOperator[SCT, ShapeT]
    ) -> _ProductLinearOperator[_SCT_co, SCT, ShapeT]: ...
    @overload
    def __matmul__(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def __matmul__(self, /, x: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def __matmul__(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    __call__ = __matmul__

    # TODO(@jorenham): shape-typing
    @overload
    def __rmatmul__[SCT: _Scalar, ShapeT: _Shape](
        self, /, x: LinearOperator[SCT, ShapeT]
    ) -> _ProductLinearOperator[_SCT_co, SCT, ShapeT]: ...
    @overload
    def __rmatmul__(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def __rmatmul__(self, /, x: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
    @overload
    def __rmatmul__(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...

    #
    @overload
    def __truediv__(self, other: onp.ToFloat, /) -> _ScaledLinearOperator[_SCT_co, _ShapeT_co]: ...
    @overload
    def __truediv__(self, other: onp.ToJustComplex128, /) -> _ScaledLinearOperator[np.complex128, _ShapeT_co]: ...
    @overload
    def __truediv__(self, other: onp.ToComplex, /) -> _ScaledLinearOperator[Any, _ShapeT_co]: ...

    #
    def __neg__(self, /) -> _ScaledLinearOperator[_SCT_co, _ShapeT_co]: ...
    def __add__[SCT: _Scalar](self, x: LinearOperator[SCT], /) -> _SumLinearOperator[_SCT_co, SCT, _ShapeT_co]: ...
    __sub__ = __add__
    def __pow__(self, p: onp.ToInt, /) -> _PowerLinearOperator[_SCT_co, _ShapeT_co]: ...

@final
class _CustomLinearOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    args: tuple[()]

    #
    @overload  # no dtype
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.int8 | Any, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: None = None,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype known (positional)
    def __init__(
        self,
        /,
        shape: _ShapeT_co,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.ToDType[_SCT_co],
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype known (keyword)
    def __init__(
        self,
        /,
        shape: _ShapeT_co,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.ToDType[_SCT_co],
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like int_ (positional)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.int_, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like int_ (keyword)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.int_, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like float64 (positional)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.float64, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like float64 (keyword)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.float64, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like complex128 (positional)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.complex128, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # dtype-like complex128 (keyword)
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[np.complex128, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # unknown dtype
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[Any, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: type | str | None = None,
        rmatmat: _FunMatMat | None = None,
        xp: None = None,
    ) -> None: ...
    @overload  # xp given
    def __init__[ShapeT: _Shape](
        self: _CustomLinearOperator[Any, ShapeT],
        /,
        shape: ShapeT,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: Incomplete | None = None,
        rmatmat: _FunMatMat | None = None,
        *,
        xp: ModuleType,
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

@type_check_only
class _UnaryLinearOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    A: LinearOperator[_SCT_co, _ShapeT_co]
    args: tuple[LinearOperator[_SCT_co, _ShapeT_co]]

    #
    def __new__(cls, A: LinearOperator[_SCT_co, _ShapeT_co], xp: ModuleType | None = None) -> Self: ...
    def __init__(self, /, A: LinearOperator[_SCT_co, _ShapeT_co], xp: ModuleType | None = None) -> None: ...

    #
    @override
    def _adjoint(self, /) -> _AdjointLinearOperator[_SCT_co, _ShapeT_co]: ...

@final
class _AdjointLinearOperator(_UnaryLinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]): ...

@final
class _TransposedLinearOperator(_UnaryLinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]): ...

@final
class _SumLinearOperator(LinearOperator[_SCT1_co | _SCT2_co, _ShapeT_co], Generic[_SCT1_co, _SCT2_co, _ShapeT_co]):
    args: tuple[LinearOperator[_SCT1_co, _ShapeT_co], LinearOperator[_SCT2_co, _ShapeT_co]]

    def __new__(
        cls, A: LinearOperator[_SCT1_co, _ShapeT_co], B: LinearOperator[_SCT2_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> Self: ...
    def __init__(
        self, /, A: LinearOperator[_SCT1_co, _ShapeT_co], B: LinearOperator[_SCT2_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

@final
class _ProductLinearOperator(LinearOperator[_SCT1_co | _SCT2_co, _ShapeT_co], Generic[_SCT1_co, _SCT2_co, _ShapeT_co]):
    args: tuple[LinearOperator[_SCT1_co, _ShapeT_co], LinearOperator[_SCT2_co, _ShapeT_co]]

    #
    def __new__(
        cls, A: LinearOperator[_SCT1_co, _ShapeT_co], B: LinearOperator[_SCT2_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> Self: ...
    def __init__(
        self, /, A: LinearOperator[_SCT1_co, _ShapeT_co], B: LinearOperator[_SCT2_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

# mypy reports a false positive `overload-overlap` error with numpy<2.5
# mypy: disable-error-code=overload-overlap

@final
class _ScaledLinearOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    args: tuple[LinearOperator[_SCT_co, _ShapeT_co], _SCT_co | complex]

    #
    @overload
    def __new__(cls, A: LinearOperator[_SCT_co, _ShapeT_co], alpha: _SCT_co | complex, xp: ModuleType | None = None) -> Self: ...
    @overload
    def __new__[ShapeT: _Shape](
        cls, A: LinearOperator[npc.floating, ShapeT], alpha: onp.ToFloat64, xp: ModuleType | None = None
    ) -> _ScaledLinearOperator[np.float64, ShapeT]: ...
    @overload
    def __new__[ShapeT: _Shape](
        cls, A: LinearOperator[npc.complexfloating, ShapeT], alpha: onp.ToComplex128, xp: ModuleType | None = None
    ) -> _ScaledLinearOperator[np.complex128, ShapeT]: ...

    #
    @overload
    def __init__(
        self, /, A: LinearOperator[_SCT_co, _ShapeT_co], alpha: _SCT_co | complex, xp: ModuleType | None = None
    ) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](
        self: _ScaledLinearOperator[np.float64, ShapeT],
        /,
        A: LinearOperator[npc.floating, ShapeT],
        alpha: onp.ToFloat64,
        xp: ModuleType | None = None,
    ) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](
        self: _ScaledLinearOperator[np.complex128, ShapeT],
        /,
        A: LinearOperator[Any, ShapeT],
        alpha: onp.ToComplex128,
        xp: ModuleType | None = None,
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

@final
class _PowerLinearOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    args: tuple[LinearOperator[_SCT_co, _ShapeT_co], SupportsIndex]

    @override
    def __new__(cls, A: LinearOperator[_SCT_co, _ShapeT_co], p: SupportsIndex, xp: ModuleType | None = None) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, A: LinearOperator[_SCT_co, _ShapeT_co], p: SupportsIndex, xp: ModuleType | None = None) -> None: ...  # pyrefly:ignore[bad-override]

    #
    @override
    def _adjoint(self, /) -> Self: ...

class MatrixLinearOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    A: _spbase[_SCT_co, _ShapeT_co] | onp.ArrayND[_SCT_co, _ShapeT_co]
    args: tuple[_spbase[_SCT_co, _ShapeT_co] | onp.ArrayND[_SCT_co, _ShapeT_co]]

    def __new__(
        cls, /, A: _spbase[_SCT_co, _ShapeT_co] | onp.ArrayND[_SCT_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> Self: ...
    def __init__(
        self, /, A: _spbase[_SCT_co, _ShapeT_co] | onp.ArrayND[_SCT_co, _ShapeT_co], xp: ModuleType | None = None
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> _AdjointMatrixOperator[_SCT_co, _ShapeT_co]: ...

@final
class _AdjointMatrixOperator(MatrixLinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    # pyrefly: ignore [bad-override]
    args: tuple[LinearOperator[_SCT_co, _ShapeT_co]]  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleVariableOverride]

    #
    @override
    def __new__(cls, A: LinearOperator[_SCT_co, _ShapeT_co], xp: ModuleType | None = None) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, A: LinearOperator[_SCT_co, _ShapeT_co], xp: ModuleType | None = None) -> None: ...  # pyrefly:ignore[bad-override]

    #
    @override
    def _adjoint(self, /) -> MatrixLinearOperator[_SCT_co, _ShapeT_co]: ...  # type: ignore[override] # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]

class IdentityOperator(LinearOperator[_SCT_co, _ShapeT_co], Generic[_SCT_co, _ShapeT_co]):
    @overload
    def __new__(cls, shape: _ShapeT_co, dtype: onp.ToDType[_SCT_co], xp: None = None) -> Self: ...
    @overload
    def __new__[ShapeT: _Shape](
        cls, shape: ShapeT, dtype: onp.AnyFloat64DType | None = None, xp: None = None
    ) -> IdentityOperator[np.float64, ShapeT]: ...
    @overload
    def __new__[ShapeT: _Shape](
        cls, shape: ShapeT, dtype: onp.AnyComplex128DType, xp: None = None
    ) -> IdentityOperator[np.complex128, ShapeT]: ...
    @overload
    def __new__[ShapeT: _Shape](cls, shape: ShapeT, dtype: str, xp: None = None) -> IdentityOperator[Any, ShapeT]: ...
    @overload
    def __new__[ShapeT: _Shape](
        cls, shape: ShapeT, dtype: Incomplete | None = None, *, xp: ModuleType
    ) -> IdentityOperator[Any, ShapeT]: ...

    #
    @overload
    def __init__(self, /, shape: _ShapeT_co, dtype: onp.ToDType[_SCT_co], xp: None = None) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](
        self: IdentityOperator[np.float64, ShapeT], /, shape: ShapeT, dtype: onp.AnyFloat64DType | None = None, xp: None = None
    ) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](
        self: IdentityOperator[np.complex128, ShapeT], /, shape: ShapeT, dtype: onp.AnyComplex128DType, xp: None = None
    ) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](self: IdentityOperator[Any, ShapeT], /, shape: ShapeT, dtype: str, xp: None = None) -> None: ...
    @overload
    def __init__[ShapeT: _Shape](
        self: IdentityOperator[Any, ShapeT], /, shape: ShapeT, dtype: Incomplete | None = None, *, xp: ModuleType
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

#
@overload
def aslinearoperator[ScalarT: _Scalar, ShapeT: _Shape](
    A: nptc.CanArray[ShapeT, np.dtype[ScalarT]],
) -> MatrixLinearOperator[ScalarT, ShapeT]: ...
@overload
def aslinearoperator[ScalarT: _Scalar, ShapeT: _Shape](A: _spbase[ScalarT, ShapeT]) -> MatrixLinearOperator[ScalarT, ShapeT]: ...
@overload
def aslinearoperator[ScalarT: npc.inexact, ShapeT: _Shape](
    A: _HasShapeAndDTypeAndMatVec[ScalarT, ShapeT],
) -> MatrixLinearOperator[ScalarT, ShapeT]: ...
@overload
def aslinearoperator[ScalarT: npc.inexact, ShapeT: _Shape](
    A: _HasShapeAndMatVec[ScalarT, ShapeT],
) -> MatrixLinearOperator[ScalarT, ShapeT]: ...
