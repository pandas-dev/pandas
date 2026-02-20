from collections.abc import Sequence
from typing import Any, ClassVar, Generic, Literal, TypeAlias, overload, type_check_only
from typing_extensions import TypeAliasType, TypeIs, TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import sparray
from ._compressed import _cs_matrix
from ._csr import _csr_base, csr_array, csr_matrix
from ._matrix import spmatrix
from ._typing import _Sparse2D, _ToShape2D

__all__ = ["csc_array", "csc_matrix", "isspmatrix_csc"]

_T = TypeVar("_T")

_Scalar: TypeAlias = npc.number | np.bool_
_ScalarT = TypeVar("_ScalarT", bound=_Scalar)
_ScalarT_co = TypeVar("_ScalarT_co", bound=_Scalar, default=Any, covariant=True)

_Seq2D: TypeAlias = Sequence[Sequence[_T]]

_ToIndices: TypeAlias = onp.CanArrayND[npc.integer] | Sequence[int]

_RawCSC = TypeAliasType(
    # `(data, (row_ind, col_ind))` or `(data, indices, indptr)`
    "_RawCSC",
    tuple[_T, tuple[_ToIndices, _ToIndices]] | tuple[_T, _ToIndices, _ToIndices],
    type_params=(_T,),
)
_ToCSC = TypeAliasType(
    "_ToCSC",
    (
        _Sparse2D[_ScalarT]
        | onp.CanArrayND[_ScalarT]
        | _RawCSC[onp.CanArrayND[_ScalarT] | Sequence[_ScalarT]]
        | Sequence[onp.CanArrayND[_ScalarT] | Sequence[_ScalarT]]
    ),
    type_params=(_ScalarT,),
)
_ToAnyCSC = TypeAliasType(
    "_ToAnyCSC", _ToShape2D | _Sparse2D[_Scalar] | onp.ToArray2D[complex, _Scalar] | _RawCSC[onp.ToComplex1D]
)
_ToBoolCSC: TypeAlias = _Seq2D[bool] | _RawCSC[Sequence[bool]]
_ToIntCSC: TypeAlias = _Seq2D[op.JustInt] | _RawCSC[Sequence[op.JustInt]]
_ToFloatCSC: TypeAlias = _Seq2D[op.JustFloat] | _RawCSC[Sequence[op.JustFloat]] | _ToShape2D
_ToComplexCSC: TypeAlias = _Seq2D[op.JustComplex] | _RawCSC[Sequence[op.JustComplex]]

###

class _csc_base(_cs_matrix[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    _format: ClassVar = "csc"

    @property
    @override
    def format(self, /) -> Literal["csc"]: ...
    @property
    @override
    def ndim(self, /) -> Literal[2]: ...
    @property
    @override
    def shape(self, /) -> tuple[int, int]: ...

    #
    @override
    def transpose(  # type: ignore[override]
        self, /, axes: tuple[Literal[1, -1], Literal[0]] | None = None, copy: bool = False
    ) -> _csr_base[_ScalarT_co, tuple[int, int]]: ...

    #
    @override
    @overload
    def count_nonzero(self, /, axis: None = None) -> np.intp: ...
    @overload
    def count_nonzero(self, /, axis: op.CanIndex) -> onp.Array1D[np.intp]: ...

class csc_array(_csc_base[_ScalarT_co], sparray[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> csc_array[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> csc_array[_ScalarT]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> csc_array[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> csc_array[np.float64]: ...

    # NOTE: keep in sync with `csc_matrix.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToCSC[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT_co] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like bool, dtype: bool-like | None
    def __init__(
        self: csc_array[np.bool_],
        /,
        arg1: _ToBoolCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~int, dtype: int-like | None
    def __init__(
        self: csc_array[np.int_],
        /,
        arg1: _ToIntCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~float, dtype: float64-like | None
    def __init__(
        self: csc_array[np.float64],
        /,
        arg1: _ToFloatCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~complex, dtype: complex128-like | None
    def __init__(
        self: csc_array[np.complex128],
        /,
        arg1: _ToComplexCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: bool-like (positional)
    def __init__(
        self: csc_array[np.bool_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: bool-like (keyword)
    def __init__(
        self: csc_array[np.bool_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: int-like (positional)
    def __init__(
        self: csc_array[np.int_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: int-like (keyword)
    def __init__(
        self: csc_array[np.int_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: float-like (positional)
    def __init__(
        self: csc_array[np.float64],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: float-like (keyword)
    def __init__(
        self: csc_array[np.float64],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: complex128-like (positional)
    def __init__(
        self: csc_array[np.complex128],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: complex128-like (keyword)
    def __init__(
        self: csc_array[np.complex128],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def transpose(  # type: ignore[override]
        self, /, axes: tuple[Literal[1, -1], Literal[0]] | None = None, copy: bool = False
    ) -> csr_array[_ScalarT_co, tuple[int, int]]: ...

class csc_matrix(_csc_base[_ScalarT_co], spmatrix[_ScalarT_co], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> csc_matrix[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> csc_matrix[_ScalarT]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> csc_matrix[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> csc_matrix[np.float64]: ...

    # NOTE: keep in sync with `csc_array.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToCSC[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: onp.ToDType[_ScalarT_co] | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like bool, dtype: bool-like | None
    def __init__(
        self: csc_matrix[np.bool_],
        /,
        arg1: _ToBoolCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyBoolDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~int, dtype: int-like | None
    def __init__(
        self: csc_matrix[np.int_],
        /,
        arg1: _ToIntCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~float, dtype: float64-like | None
    def __init__(
        self: csc_matrix[np.float64],
        /,
        arg1: _ToFloatCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d array-like ~complex, dtype: complex128-like | None
    def __init__(
        self: csc_matrix[np.complex128],
        /,
        arg1: _ToComplexCSC,
        shape: _ToShape2D | None = None,
        dtype: onp.AnyComplex128DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: bool-like (positional)
    def __init__(
        self: csc_matrix[np.bool_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: bool-like (keyword)
    def __init__(
        self: csc_matrix[np.bool_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: int-like (positional)
    def __init__(
        self: csc_matrix[np.int_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: int-like (keyword)
    def __init__(
        self: csc_matrix[np.int_],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: float-like (positional)
    def __init__(
        self: csc_matrix[np.float64],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: float-like (keyword)
    def __init__(
        self: csc_matrix[np.float64],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyFloat64DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: complex128-like (positional)
    def __init__(
        self: csc_matrix[np.complex128],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # dtype: complex128-like (keyword)
    def __init__(
        self: csc_matrix[np.complex128],
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (positional)
    def __init__(
        self,
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-D, dtype: <known> (keyword)
    def __init__(
        self,
        /,
        arg1: _ToAnyCSC,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...

    #
    @override
    def transpose(  # type: ignore[override]
        self, /, axes: tuple[Literal[1, -1], Literal[0]] | None = None, copy: bool = False
    ) -> csr_matrix[_ScalarT_co]: ...

    #
    @override
    @overload
    def getnnz(self, /, axis: None = None) -> int: ...
    @overload
    def getnnz(self, /, axis: op.CanIndex) -> onp.Array1D[np.int32]: ...

def isspmatrix_csc(x: object) -> TypeIs[csc_matrix[Any]]: ...
