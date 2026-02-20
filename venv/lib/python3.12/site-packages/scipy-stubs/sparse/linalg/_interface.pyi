# mypy: disable-error-code="override"
import types
from collections.abc import Callable, Iterable
from typing import Any, ClassVar, Final, Generic, Protocol, Self, TypeAlias, final, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

__all__ = ["LinearOperator", "aslinearoperator"]

_SCT = TypeVar("_SCT", bound=npc.number | np.bool_)
_SCT_co = TypeVar("_SCT_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_SCT1_co = TypeVar("_SCT1_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_SCT2_co = TypeVar("_SCT2_co", bound=npc.number | np.bool_, default=_SCT1_co, covariant=True)
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_FunMatVecT_co = TypeVar("_FunMatVecT_co", bound=_FunMatVec, default=_FunMatVec, covariant=True)

_LinearOperatorT = TypeVar("_LinearOperatorT", bound=LinearOperator[Any])
_LinearOperatorT_co = TypeVar("_LinearOperatorT_co", bound=LinearOperator[Any], covariant=True)

_ToShape: TypeAlias = Iterable[op.CanIndex]
_Real: TypeAlias = np.bool_ | npc.integer | npc.floating
_FunMatVec: TypeAlias = Callable[[onp.ArrayND[Any]], onp.ToComplex1D | onp.ToComplex2D]
_FunMatMat: TypeAlias = Callable[[onp.Array2D[Any]], onp.ToComplex2D]

@type_check_only
class _CanAdjoint(Protocol[_LinearOperatorT_co]):
    def _adjoint(self, /) -> _LinearOperatorT_co: ...

@type_check_only
class _HasShapeAndMatVec(Protocol[_SCT_co]):
    shape: tuple[int, int]
    @overload
    def matvec(self, /, x: onp.CanArray1D[np.float64]) -> onp.CanArray1D[_SCT_co]: ...
    @overload
    def matvec(self, /, x: onp.CanArray2D[np.float64]) -> onp.CanArray2D[_SCT_co]: ...
    @overload
    def matvec(self, /, x: onp.CanArray1D[np.complex128]) -> onp.ToComplex1D: ...
    @overload
    def matvec(self, /, x: onp.CanArray2D[np.complex128]) -> onp.ToComplex2D: ...

@type_check_only
class _HasShapeAndDTypeAndMatVec(Protocol[_SCT_co]):
    shape: tuple[int, int]
    @property
    def dtype(self, /) -> np.dtype[_SCT_co]: ...
    @overload
    def matvec(self, /, x: onp.CanArray1D[np.float64] | onp.CanArray1D[np.complex128]) -> onp.ToComplex1D: ...
    @overload
    def matvec(self, /, x: onp.CanArray2D[np.float64] | onp.CanArray2D[np.complex128]) -> onp.ToComplex2D: ...

###

class LinearOperator(Generic[_SCT_co]):
    __array_ufunc__: ClassVar[None] = None
    ndim: ClassVar[int] = 2

    shape: Final[tuple[int, int]]
    dtype: np.dtype[_SCT_co]

    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    # keep in sync with `_CustomLinearOperator.__init__`
    @overload  # no dtype
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: None = None,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.int8 | Any]: ...
    @overload  # dtype known (positional)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.ToDType[_SCT],
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[_SCT]: ...
    @overload  # dtype known (keyword)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.ToDType[_SCT],
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[_SCT]: ...
    @overload  # dtype-like int_ (positional)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.int_]: ...
    @overload  # dtype-like int_ (keyword)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.int_]: ...
    @overload  # dtype-like float64 (positional)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.float64]: ...
    @overload  # dtype-like float64 (positional)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.float64]: ...
    @overload  # dtype-like complex128 (positional)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.complex128]: ...
    @overload  # dtype-like complex128 (keyword)
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[np.complex128]: ...
    @overload  # unknown dtype
    def __new__(
        cls,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: type | str | None = None,
        rmatmat: _FunMatMat | None = None,
    ) -> _CustomLinearOperator[Any]: ...

    # NOTE: the `__init__` method cannot be annotated, because it will cause mypy to `__new__`:
    # https://github.com/python/mypy/issues/17251

    # ruff: noqa: ERA001

    # @overload
    # def __init__(self, /, dtype: onp.ToDType[_SCT_co], shape: _ToShape) -> None: ...
    # @overload
    # def __init__(self: LinearOperator[np.int_], /, dtype: onp.AnyIntDType, shape: _ToShape) -> None: ...
    # @overload
    # def __init__(self: LinearOperator[np.float64], /, dtype: onp.AnyFloat64DType, shape: _ToShape) -> None: ...
    # @overload
    # def __init__(self: LinearOperator[np.complex128], /, dtype: onp.AnyComplex128DType, shape: _ToShape) -> None: ...
    # @overload
    # def __init__(self: LinearOperator[Any], /, dtype: type | str | None, shape: _ToShape) -> None: ...

    #
    @property
    def H(self: _CanAdjoint[_LinearOperatorT], /) -> _LinearOperatorT: ...
    @final
    def adjoint(self: _CanAdjoint[_LinearOperatorT], /) -> _LinearOperatorT: ...
    def _adjoint(self, /) -> LinearOperator[_SCT_co]: ...

    #
    @property
    def T(self, /) -> _TransposedLinearOperator[_SCT_co]: ...
    def transpose(self, /) -> _TransposedLinearOperator[_SCT_co]: ...

    #
    @overload  # float array 1d
    def matvec(self, /, x: onp.ToFloatStrict1D) -> onp.Array1D[_SCT_co]: ...
    @overload  # float matrix
    def matvec(self, /, x: onp.Matrix[_Real]) -> onp.Matrix[_SCT_co]: ...
    @overload  # complex matrix
    def matvec(self, /, x: onp.Matrix[np.complex128]) -> onp.Matrix[np.complex128]: ...
    @overload  # float array 2d
    def matvec(self, /, x: onp.ToFloatStrict2D) -> onp.Array2D[_SCT_co]: ...
    @overload  # complex array 1d
    def matvec(self, /, x: onp.ToJustComplex128Strict1D) -> onp.Array1D[np.complex128]: ...
    @overload  # complex array 2d
    def matvec(self, /, x: onp.ToJustComplex128Strict2D) -> onp.Array2D[np.complex128]: ...
    @overload  # float array
    def matvec(self, /, x: onp.ToFloat2D) -> onp.ArrayND[_SCT_co]: ...
    @overload  # complex array
    def matvec(self, /, x: onp.ToJustComplex128_2D) -> onp.ArrayND[np.complex128]: ...
    @overload  # unknown array
    def matvec(self, /, x: onp.ToComplex2D) -> onp.ArrayND[Any]: ...
    rmatvec = matvec

    #
    @overload
    def matmat(self, /, X: onp.ToFloat2D) -> onp.Array2D[_SCT_co]: ...
    @overload
    def matmat(self, /, X: onp.ToJustComplex128_2D) -> onp.Array2D[np.complex128]: ...
    @overload
    def matmat(self, /, X: onp.ToComplex2D) -> onp.Array2D[Any]: ...
    rmatmat = matmat

    #
    @overload
    def dot(self, /, x: LinearOperator[_SCT]) -> _ProductLinearOperator[_SCT_co, _SCT]: ...
    @overload
    def dot(self, /, x: onp.ToFloat) -> _ScaledLinearOperator[_SCT_co]: ...
    @overload
    def dot(self, /, x: onp.ToJustComplex128) -> _ScaledLinearOperator[np.complex128]: ...
    @overload
    def dot(self, /, x: onp.ToFloatStrict1D) -> onp.Array1D[_SCT_co]: ...
    @overload
    def dot(self, /, x: onp.ToJustComplex128Strict1D) -> onp.Array1D[np.complex128]: ...
    @overload
    def dot(self, /, x: onp.ToFloatStrict2D) -> onp.Array2D[_SCT_co]: ...
    @overload
    def dot(self, /, x: onp.ToJustComplex128Strict2D) -> onp.Array2D[np.complex128]: ...
    @overload
    def dot(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def dot(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    __mul__ = dot
    __rmul__ = dot
    __call__ = dot

    #
    @overload
    def __matmul__(self, /, x: LinearOperator[_SCT]) -> _ProductLinearOperator[_SCT_co, _SCT]: ...
    @overload
    def __matmul__(self, /, x: onp.ToFloatStrict1D) -> onp.Array1D[_SCT_co]: ...
    @overload
    def __matmul__(self, /, x: onp.ToJustComplex128Strict1D) -> onp.Array1D[np.complex128]: ...
    @overload
    def __matmul__(self, /, x: onp.ToFloatStrict2D) -> onp.Array2D[_SCT_co]: ...
    @overload
    def __matmul__(self, /, x: onp.ToJustComplex128Strict2D) -> onp.Array2D[np.complex128]: ...
    @overload
    def __matmul__(self, /, x: onp.ToFloatND) -> onp.ArrayND[_SCT_co]: ...
    @overload
    def __matmul__(self, /, x: onp.ToComplexND) -> onp.ArrayND[Any]: ...
    __rmatmul__ = __matmul__

    #
    @overload
    def __truediv__(self, other: onp.ToFloat, /) -> _ScaledLinearOperator[_SCT_co]: ...
    @overload
    def __truediv__(self, other: onp.ToJustComplex128, /) -> _ScaledLinearOperator[np.complex128]: ...
    @overload
    def __truediv__(self, other: onp.ToComplex, /) -> _ScaledLinearOperator[Any]: ...

    #
    def __neg__(self, /) -> _ScaledLinearOperator[_SCT_co]: ...
    def __add__(self, x: LinearOperator[_SCT], /) -> _SumLinearOperator[_SCT_co, _SCT]: ...
    __sub__ = __add__
    def __pow__(self, p: onp.ToInt, /) -> _PowerLinearOperator[_SCT_co]: ...

@final
class _CustomLinearOperator(LinearOperator[_SCT_co], Generic[_SCT_co, _FunMatVecT_co]):
    args: tuple[()]

    #
    @overload  # no dtype
    def __init__(
        self: _CustomLinearOperator[np.int8 | Any],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: None = None,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype known (positional)
    def __init__(
        self,
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.ToDType[_SCT_co],
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype known (keyword)
    def __init__(
        self,
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.ToDType[_SCT_co],
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like int_ (positional)
    def __init__(
        self: _CustomLinearOperator[np.int_],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like int_ (keyword)
    def __init__(
        self: _CustomLinearOperator[np.int_],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyIntDType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like float64 (positional)
    def __init__(
        self: _CustomLinearOperator[np.float64],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like float64 (keyword)
    def __init__(
        self: _CustomLinearOperator[np.float64],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like complex128 (positional)
    def __init__(
        self: _CustomLinearOperator[np.complex128],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None,
        matmat: _FunMatMat | None,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # dtype-like complex128 (keyword)
    def __init__(
        self: _CustomLinearOperator[np.complex128],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...
    @overload  # unknown dtype
    def __init__(
        self: _CustomLinearOperator[Any],
        /,
        shape: _ToShape,
        matvec: _FunMatVec,
        rmatvec: _FunMatVec | None = None,
        matmat: _FunMatMat | None = None,
        dtype: type | str | None = None,
        rmatmat: _FunMatMat | None = None,
    ) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

@type_check_only
class _UnaryLinearOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    A: LinearOperator[_SCT_co]
    args: tuple[LinearOperator[_SCT_co]]

    #
    def __new__(cls, A: LinearOperator[_SCT_co]) -> Self: ...
    def __init__(self, /, A: LinearOperator[_SCT_co]) -> None: ...
    @override
    def _adjoint(self, /) -> _AdjointLinearOperator[_SCT_co]: ...

@final
class _AdjointLinearOperator(_UnaryLinearOperator[_SCT_co], Generic[_SCT_co]): ...

@final
class _TransposedLinearOperator(_UnaryLinearOperator[_SCT_co], Generic[_SCT_co]): ...

@final
class _SumLinearOperator(LinearOperator[_SCT1_co | _SCT2_co], Generic[_SCT1_co, _SCT2_co]):
    args: tuple[LinearOperator[_SCT1_co], LinearOperator[_SCT2_co]]

    def __new__(cls, A: LinearOperator[_SCT1_co], B: LinearOperator[_SCT2_co]) -> Self: ...
    def __init__(self, /, A: LinearOperator[_SCT1_co], B: LinearOperator[_SCT2_co]) -> None: ...
    @override
    def _adjoint(self, /) -> Self: ...

@final
class _ProductLinearOperator(LinearOperator[_SCT1_co | _SCT2_co], Generic[_SCT1_co, _SCT2_co]):
    args: tuple[LinearOperator[_SCT1_co], LinearOperator[_SCT2_co]]

    #
    def __new__(cls, A: LinearOperator[_SCT1_co], B: LinearOperator[_SCT2_co]) -> Self: ...
    def __init__(self, /, A: LinearOperator[_SCT1_co], B: LinearOperator[_SCT2_co]) -> None: ...
    @override
    def _adjoint(self, /) -> Self: ...

@final
class _ScaledLinearOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    args: tuple[LinearOperator[_SCT_co], _SCT_co | complex]

    #
    @overload
    def __new__(cls, A: LinearOperator[_SCT_co], alpha: _SCT_co | complex) -> Self: ...  # type: ignore[overload-overlap]
    @overload
    def __new__(cls, A: LinearOperator[npc.floating], alpha: onp.ToFloat64) -> _ScaledLinearOperator[np.float64]: ...
    @overload
    def __new__(cls, A: LinearOperator[npc.complexfloating], alpha: onp.ToComplex128) -> _ScaledLinearOperator[np.complex128]: ...

    #
    @overload
    def __init__(self, /, A: LinearOperator[_SCT_co], alpha: _SCT_co | complex) -> None: ...
    @overload
    def __init__(self: _ScaledLinearOperator[np.float64], /, A: LinearOperator[npc.floating], alpha: onp.ToFloat64) -> None: ...
    @overload
    def __init__(self: _ScaledLinearOperator[np.complex128], /, A: LinearOperator, alpha: onp.ToComplex128) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

@final
class _PowerLinearOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    args: tuple[LinearOperator[_SCT_co], op.CanIndex]

    def __new__(cls, A: LinearOperator[_SCT_co], p: op.CanIndex) -> Self: ...
    def __init__(self, /, A: LinearOperator[_SCT_co], p: op.CanIndex) -> None: ...
    @override
    def _adjoint(self, /) -> Self: ...

class MatrixLinearOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    A: _spbase | onp.Array2D[_SCT_co]
    args: tuple[_spbase | onp.Array2D[_SCT_co]]

    def __new__(cls, A: _spbase | onp.ArrayND[_SCT_co]) -> Self: ...
    def __init__(self, /, A: _spbase | onp.ArrayND[_SCT_co]) -> None: ...
    @override
    def _adjoint(self, /) -> _AdjointMatrixOperator[_SCT_co]: ...

@final
class _AdjointMatrixOperator(MatrixLinearOperator[_SCT_co], Generic[_SCT_co]):
    args: tuple[MatrixLinearOperator[_SCT_co]]  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleVariableOverride]

    @property
    @override
    def dtype(self, /) -> np.dtype[_SCT_co]: ...  # pyright: ignore[reportIncompatibleVariableOverride]  # pyrefly: ignore[bad-override]

    #
    def __new__(cls, adjoint_array: LinearOperator[_SCT_co]) -> Self: ...
    def __init__(self, /, adjoint_array: LinearOperator[_SCT_co]) -> None: ...
    @override
    def _adjoint(self, /) -> MatrixLinearOperator[_SCT_co]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override]

class IdentityOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    @overload
    def __new__(cls, shape: _ToShape, dtype: onp.ToDType[_SCT_co]) -> Self: ...
    @overload
    def __new__(cls, shape: _ToShape, dtype: onp.AnyFloat64DType | None = None) -> IdentityOperator[np.float64]: ...
    @overload
    def __new__(cls, shape: _ToShape, dtype: onp.AnyComplex128DType) -> IdentityOperator[np.complex128]: ...
    @overload
    def __new__(cls, shape: _ToShape, dtype: str) -> IdentityOperator[Any]: ...

    #
    @overload
    def __init__(self, /, shape: _ToShape, dtype: onp.ToDType[_SCT_co]) -> None: ...
    @overload
    def __init__(self: IdentityOperator[np.float64], /, shape: _ToShape, dtype: onp.AnyFloat64DType | None = None) -> None: ...
    @overload
    def __init__(self: IdentityOperator[np.complex128], /, shape: _ToShape, dtype: onp.AnyComplex128DType) -> None: ...
    @overload
    def __init__(self: IdentityOperator[Any], /, shape: _ToShape, dtype: str) -> None: ...

    #
    @override
    def _adjoint(self, /) -> Self: ...

#
@overload
def aslinearoperator(A: onp.CanArrayND[_InexactT]) -> MatrixLinearOperator[_InexactT]: ...
@overload
def aslinearoperator(A: _spbase[_InexactT]) -> MatrixLinearOperator[_InexactT]: ...
@overload
def aslinearoperator(
    A: onp.ArrayND[np.bool_ | npc.integer | np.float64] | _spbase[np.bool_ | npc.integer | np.float64],
) -> MatrixLinearOperator[np.float64]: ...
@overload
def aslinearoperator(A: _HasShapeAndDTypeAndMatVec[_InexactT]) -> MatrixLinearOperator[_InexactT]: ...
@overload
def aslinearoperator(A: _HasShapeAndMatVec[_InexactT]) -> MatrixLinearOperator[_InexactT]: ...
