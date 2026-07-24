from collections.abc import Sequence
from typing import Any, Generic, Literal, Self, SupportsIndex, overload, override, type_check_only
from typing_extensions import TypeIs, TypeVar

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._coo import coo_array, coo_matrix
from ._csr import csr_matrix
from ._index import IndexMixin
from ._matrix import spmatrix
from ._typing import _ToShape2D

__all__ = ["isspmatrix_lil", "lil_array", "lil_matrix"]

#

type _ToMatrixPy[_T] = list[_T] | list[list[_T]]  # intentionally invariant
type _ToMatrix[_ScalarT: npc.number | np.bool] = (
    _spbase[_ScalarT] | onp.CanArrayND[_ScalarT] | Sequence[onp.CanArrayND[_ScalarT]] | _ToMatrixPy[_ScalarT]
)
type _ToAnyLIL = _ToShape2D | _ToMatrix[npc.number | np.bool]

_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool, default=Any, covariant=True)

###

class _lil_base(_spbase[_ScalarT_co, tuple[int, int]], IndexMixin[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    dtype: np.dtype[_ScalarT_co]
    data: onp.Array1D[np.object_]
    rows: onp.Array1D[np.object_]

    @property
    @override
    def format(self, /) -> Literal["lil"]: ...
    @property
    @override
    def ndim(self, /) -> Literal[2]: ...
    @property
    @override
    def shape(self, /) -> tuple[int, int]: ...

    #
    def __init__(
        self,
        /,
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def __iadd__(self, other: onp.ToFalse | _spbase | onp.ArrayND[npc.number | np.bool], /) -> Self: ...  # type:ignore[override]
    @override
    def __isub__(self, other: onp.ToFalse | _spbase | onp.ArrayND[npc.number | np.bool], /) -> Self: ...  # type:ignore[override]
    @override
    def __imul__(self, other: onp.ToComplex, /) -> Self: ...  # type: ignore[override]
    @override
    def __itruediv__(self, other: onp.ToComplex, /) -> Self: ...  # type: ignore[override]

    #
    @override
    def tolil(self, /, copy: bool = False) -> Self: ...  # type: ignore[override]
    @override
    def resize(self, /, *shape: int) -> None: ...  # pyright: ignore[reportIncompatibleMethodOverride] # pyrefly: ignore[bad-override] # ty: ignore[invalid-method-override] # zuban: ignore[override]

    #
    @override
    @overload
    def count_nonzero(self, /, axis: None = None) -> np.intp: ...
    @overload
    def count_nonzero(self, /, axis: SupportsIndex) -> onp.Array1D[np.intp]: ...

    #
    def getrowview(self, /, i: int) -> Self: ...
    def getrow(self, /, i: onp.ToJustInt) -> lil_array[_ScalarT_co] | csr_matrix[_ScalarT_co]: ...

class lil_array(_lil_base[_ScalarT_co], sparray[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    # NOTE: These two methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_array[_ScalarT_co, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_array[_ScalarT, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_as_any__(self, /) -> lil_array[Any]: ...

    # NOTE: keep the in sync with `lil_matrix.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT_co] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like | None
    def __init__(
        self: lil_array[np.float64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (positional)
    def __init__(
        self: lil_array[np.bool],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (keyword)
    def __init__(
        self: lil_array[np.bool],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (positional)
    def __init__(
        self: lil_array[np.int_],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (keyword)
    def __init__(
        self: lil_array[np.int_],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like (positional)
    def __init__(
        self: lil_array[np.float64],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like (keyword)
    def __init__(
        self: lil_array[np.float64],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: lil_array[np.complex128],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (keyword)
    def __init__(
        self: lil_array[np.complex128],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: lil_array[np.bool],
        /,
        arg1: _ToMatrixPy[bool],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.int, dtype: int-like | None
    def __init__(
        self: lil_array[np.int_],
        /,
        arg1: _ToMatrixPy[int],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: lil_array[np.float64],
        /,
        arg1: _ToMatrixPy[float],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: lil_array[np.complex128],
        /,
        arg1: _ToMatrixPy[complex],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <unknown>
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def getrow(self, /, i: onp.ToJustInt) -> lil_array[_ScalarT_co]: ...

class lil_matrix(_lil_base[_ScalarT_co], spmatrix[_ScalarT_co], Generic[_ScalarT_co]):
    # NOTE: These two methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_matrix[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> coo_matrix[_ScalarT]: ...
    @override
    @type_check_only
    def __assoc_as_any__(self, /) -> lil_matrix[Any]: ...

    # NOTE: keep the in sync with `lil_array.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT_co] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like | None
    def __init__(
        self: lil_matrix[np.float64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (positional)
    def __init__(
        self: lil_matrix[np.bool],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (keyword)
    def __init__(
        self: lil_matrix[np.bool],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (positional)
    def __init__(
        self: lil_matrix[np.int_],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (keyword)
    def __init__(
        self: lil_matrix[np.int_],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like (positional)
    def __init__(
        self: lil_matrix[np.float64],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like (keyword)
    def __init__(
        self: lil_matrix[np.float64],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: lil_matrix[np.complex128],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (keyword)
    def __init__(
        self: lil_matrix[np.complex128],
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: lil_matrix[np.bool],
        /,
        arg1: _ToMatrixPy[bool],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.int, dtype: int-like | None
    def __init__(
        self: lil_matrix[np.int_],
        /,
        arg1: _ToMatrixPy[int],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: lil_matrix[np.float64],
        /,
        arg1: _ToMatrixPy[float],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: lil_matrix[np.complex128],
        /,
        arg1: _ToMatrixPy[complex],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: <unknown>
    def __init__(
        self,
        /,
        arg1: _ToAnyLIL,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def getrow(self, /, i: onp.ToJustInt) -> csr_matrix[_ScalarT_co]: ...

    # NOTE: using `@override` together with `@overload` causes stubtest to crash...
    @override
    @overload
    def getnnz(self, /, axis: None = None) -> int: ...
    @overload
    def getnnz(self, /, axis: SupportsIndex) -> onp.Array1D[np.int32]: ...

#
def isspmatrix_lil(x: object) -> TypeIs[lil_matrix]: ...
