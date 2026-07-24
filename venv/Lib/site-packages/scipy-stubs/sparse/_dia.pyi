from collections.abc import Sequence
from typing import Any, ClassVar, Generic, Literal, overload, override, type_check_only
from typing_extensions import TypeIs, TypeVar

import numpy as np
import numpy.typing as npt
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._coo import coo_array, coo_matrix
from ._data import _data_matrix
from ._matrix import spmatrix
from ._typing import _ToShape2D

__all__ = ["dia_array", "dia_matrix", "isspmatrix_dia"]

type _ToMatrixPy[T] = Sequence[T] | Sequence[Sequence[T]]
type _ToMatrix[ScalarT: npc.number | np.bool] = (
    _spbase[ScalarT] | onp.CanArrayND[ScalarT] | Sequence[onp.CanArrayND[ScalarT]] | _ToMatrixPy[ScalarT]
)
type _ToData[ScalarT: npc.number | np.bool] = tuple[onp.ArrayND[ScalarT], onp.ArrayND[npc.integer]]

_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool, default=Any, covariant=True)

###

class _dia_base(_data_matrix[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    _format: ClassVar = "dia"

    data: onp.Array2D[_ScalarT_co]
    offsets: onp.Array1D[np.int32]

    @property
    @override
    def format(self, /) -> Literal["dia"]: ...
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

class dia_array(_dia_base[_ScalarT_co], sparray[_ScalarT_co, tuple[int, int]], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_array[_ScalarT_co, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__[ScalarT: npc.number | np.bool](self, sctype: ScalarT, /) -> coo_array[ScalarT, tuple[int, int]]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> dia_array[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> dia_array[np.float64]: ...
    @override
    @type_check_only
    def __assoc_as_any__(self, /) -> dia_array[Any]: ...

    # NOTE: keep in sync with `dia_matrix.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co] | _ToData[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: float64-like | None
    def __init__(
        self: dia_array[np.float64],
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
        self: dia_array[np.bool],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: bool-like (keyword)
    def __init__(
        self: dia_array[np.bool],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyBoolDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (positional)
    def __init__(
        self: dia_array[np.int64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: int-like (keyword)
    def __init__(
        self: dia_array[np.int64],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyIntDType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (positional)
    def __init__(
        self: dia_array[np.complex128],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: complex128-like (keyword)
    def __init__(
        self: dia_array[np.complex128],
        /,
        arg1: _ToShape2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.AnyComplex128DType,
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: dia_array[np.bool],
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
        self: dia_array[np.int_],
        /,
        arg1: _ToMatrixPy[op.JustInt],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: dia_array[np.float64],
        /,
        arg1: _ToMatrixPy[op.JustFloat],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: dia_array[np.complex128],
        /,
        arg1: _ToMatrixPy[op.JustComplex],
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
        arg1: onp.ToComplex2D,
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
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        *,
        dtype: onp.ToDType[_ScalarT_co],
        copy: bool = False,
        maxprint: int | None = None,
    ) -> None: ...

class dia_matrix(_dia_base[_ScalarT_co], spmatrix[_ScalarT_co], Generic[_ScalarT_co]):
    # NOTE: These four methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @override
    @type_check_only
    def __assoc_stacked__(self, /) -> coo_matrix[_ScalarT_co]: ...
    @override
    @type_check_only
    def __assoc_stacked_as__[ScalarT: npc.number | np.bool](self, sctype: ScalarT, /) -> coo_matrix[ScalarT]: ...
    @override
    @type_check_only
    def __assoc_as_float32__(self, /) -> dia_matrix[np.float32]: ...
    @override
    @type_check_only
    def __assoc_as_float64__(self, /) -> dia_matrix[np.float64]: ...
    @override
    @type_check_only
    def __assoc_as_any__(self, /) -> dia_matrix[Any]: ...

    # NOTE: keep in sync with `dia_array.__init__`
    @overload  # matrix-like (known dtype), dtype: None
    def __init__(
        self,
        /,
        arg1: _ToMatrix[_ScalarT_co] | _ToData[_ScalarT_co],
        shape: _ToShape2D | None = None,
        dtype: None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # 2-d shape-like, dtype: None
    def __init__(
        self: dia_matrix[np.float64],
        /,
        arg1: _ToShape2D,
        shape: None = None,
        dtype: None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.bool, dtype: bool-like | None
    def __init__(
        self: dia_matrix[np.bool],
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
        self: dia_matrix[np.int_],
        /,
        arg1: _ToMatrixPy[op.JustInt],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyIntDType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.float, dtype: float64-like | None
    def __init__(
        self: dia_matrix[np.float64],
        /,
        arg1: _ToMatrixPy[op.JustFloat],
        shape: _ToShape2D | None = None,
        dtype: onp.AnyFloat64DType | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...
    @overload  # matrix-like builtins.complex, dtype: complex128-like | None
    def __init__(
        self: dia_matrix[np.complex128],
        /,
        arg1: _ToMatrixPy[op.JustComplex],
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
        arg1: onp.ToComplex2D,
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
        arg1: onp.ToComplex2D,
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
        arg1: onp.ToComplex2D,
        shape: _ToShape2D | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

def isspmatrix_dia(x: object) -> TypeIs[dia_matrix]: ...
